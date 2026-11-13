"""Bounded whole-locus fitting and cross-fitting for evidence backends."""
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from loguru import logger
from contextlib import nullcontext

from ....utils.exceptions import IPyradError
from .sfs_prior import folded_saf, fit_folded_sfs, sfs_posterior_moments


class CrossFittedSFS:
    """Every evaluation locus uses a prior trained outside its locus fold.

    Folds use sorted locus rank modulo fold count. Fitting loci are evenly
    spaced within each fold, with a fixed total budget. This is locus holdout,
    not individual holdout or proof of biological population membership.
    """
    def __init__(self, loci, groups, *, fit_loci=1024, folds=5, max_iter=10000):
        self.loci = np.unique(np.asarray(loci, dtype=np.int64))
        if folds < 2 or fit_loci < folds or max_iter < 1 or len(self.loci) < folds:
            raise IPyradError('SFS requires >=2 folds, at least one fitting locus per fold, positive iterations, and enough retained loci.')
        self.folds = folds
        self.max_iter = max_iter
        self.groups = groups
        if any(len(ix) > 128 for ix in groups.values()):
            raise IPyradError('SFS currently supports at most 128 individuals per group.')
        chosen = []
        budget = min(fit_loci, len(self.loci))
        for fold in range(folds):
            pool = self.loci[fold::folds]
            size = min(len(pool), budget//folds + (fold < budget % folds))
            chosen.extend(pool[np.linspace(0, len(pool)-1, size, dtype=int)])
        self.fit_loci = np.sort(chosen)
        self.parts = {g: [] for g in groups}
        self.labels = []
        self.bytes = 0
        self.models = {}
        self.rows = []
        self.diagnostics = {}

    def fold_ids(self, loci):
        loci = np.asarray(loci)
        ranks = np.searchsorted(self.loci, loci)
        if np.any(ranks >= len(self.loci)) or np.any(self.loci[ranks] != loci):
            raise ValueError('Unknown locus in SFS cross-fitting.')
        return ranks % self.folds

    def add(self, likelihoods, valid, loci):
        """Add only selected fitting loci after the backend's usual filters."""
        keep = np.isin(loci, self.fit_loci)
        if not keep.any():
            return
        for group, indices in self.groups.items():
            f = folded_saf(likelihoods[indices][:, keep], valid[indices][:, keep])
            self.bytes += f.nbytes
            if self.bytes > 256*1024**2:
                raise IPyradError('SFS fitting exceeds 256 MiB; reduce --sfs-fit-loci.')
            self.parts[group].append(f)
        self.labels.append(self.fold_ids(np.asarray(loci)[keep]))

    def fit(self, cores=1):
        if not self.labels:
            raise IPyradError('No eligible SFS fitting sites.')
        labels = np.concatenate(self.labels)
        matrices = {g: np.concatenate(parts) for g, parts in self.parts.items()}
        self.parts.clear()
        self.labels.clear()
        def task(key):
            group, fold = key
            try:
                weights, diagnostics = fit_folded_sfs(matrices[group][labels != fold], max_iter=self.max_iter)
            except ValueError as exc:
                raise IPyradError(f'SFS group {group}, held-out fold {fold}: {exc}') from exc
            return key, weights, diagnostics
        keys = [(g, f) for g in self.groups for f in range(self.folds)]
        # Each independent group/fold fit uses one BLAS thread. Bound outer work.
        try:
            from threadpoolctl import threadpool_limits
            limits = threadpool_limits(limits=1)
            workers = max(1, min(int(cores or 1), len(keys), 4))
        except ImportError:
            # Optional analysis dependency; minimal installations remain usable.
            limits, workers = nullcontext(), 1
        with limits, ThreadPoolExecutor(max_workers=workers) as pool:
            for (g, fold), w, d in pool.map(task, keys):
                self.models[g, fold] = w
                self.diagnostics[f'{g}:{fold}'] = d
                for k, weight in enumerate(w):
                    self.rows.append(dict(population=g, held_out_fold=fold,
                        n_individuals=len(self.groups[g]), minor_allele_count=k,
                        weight=float(weight), **d))
                logger.info('SFS {} fold {}: {} informative sites, {} iterations', g, fold, d['informative_sites'], d['iterations'])

    def moments(self, group, likelihoods, valid, loci):
        labels = self.fold_ids(loci)
        result = None
        for fold in np.unique(labels):
            take = np.flatnonzero(labels == fold)
            part = sfs_posterior_moments(likelihoods[:, take], valid[:, take], self.models[group, int(fold)])
            if result is None:
                result = {k: np.empty((len(valid), len(labels), 3) if k == 'GP' else
                    (len(valid), len(labels)) if k == 'dosage' else (len(labels),), dtype=v.dtype)
                    for k, v in part.items()}
            for key, value in part.items():
                if key in ('GP', 'dosage'):
                    result[key][:, take] = value
                else:
                    result[key][take] = value
        return result

    def table(self):
        return pd.DataFrame(self.rows)

    def metadata(self):
        return dict(model='folded finite-sample SFS; exchangeable chromosome assignments conditional on count',
            fitting='regularized EM, Dirichlet alpha=1.5 per folded count bin; flat SAF rows omitted',
            holdout='whole-locus cross-fitting; sorted locus rank modulo folds; not individual holdout',
            folds=self.folds, fitting_loci=self.fit_loci.tolist(), max_iterations=self.max_iter,
            fitting_locus_selection='evenly spaced within each fold', diagnostics=self.diagnostics,
            missing_policy='constant likelihood in fixed-cohort SAF; hypergeometric prior projection for observed-cell summaries',
            ascertainment='spectrum conditional on input scope, allele model and eligibility filters; SNP spectrum is not an all-base spectrum',
            uncertainty='fitted-prior uncertainty is not included in per-site posterior variance; small groups and low depth can be weakly identified')
