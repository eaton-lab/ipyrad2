#!/usr/bin/env python

"""Sklearn-backed sNMF-style clustering on SNP HDF5 inputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import warnings

from loguru import logger
import numpy as np
import pandas as pd

from ...utils.exceptions import IPyradError
from ..extracters.snps_extracter import _MISSING_GENO
from .common import (
    build_imputed_sample_data_summary,
    count_linkage_blocks,
    ensure_output_paths,
    get_numerical_input,
    impute_genotype_matrix,
    log_snp_imputation_summary,
    log_snp_view_summary,
    marker_ids_from_view,
    normalize_impute_method,
    parse_k_range,
    require_hdf5_input,
    require_sklearn,
    run_snps_extracter_for_method,
    summarize_prepared_snp_view,
    write_assignments,
    write_marker_cluster_matrix,
    write_membership,
    write_sample_data_summary,
    write_stats_file,
)


_DEFAULT_ALPHA_W = 1e-4
_DEFAULT_ALPHA_H = "same"
_DEFAULT_L1_RATIO = 1.0
_DEFAULT_N_INIT = 10
_DEFAULT_CV_REPLICATES = 5
_DEFAULT_CV_HOLDOUT = 0.1
_NMF_TOL = 1e-3
_MAX_ITER = 3000
_EPSILON = 1e-12


@dataclass
class SNMFFit:
    """Selected sNMF fit for one value of K."""

    membership: np.ndarray
    genotype_frequencies: np.ndarray
    allele_frequencies: np.ndarray
    reconstruction_err: float
    n_iter: int
    hit_max_iter: bool


@dataclass
class SNMFCrossEntropyScore:
    """Cross-entropy score summary for one value of K."""

    mean_cross_entropy: float
    sd_cross_entropy: float
    capped_fit_count: int
    total_fit_count: int


def _normalize_alpha_h(alpha_h: float | str) -> float | str:
    """Normalize alpha_H, allowing floats or sklearn's `'same'` sentinel."""
    if alpha_h == "same":
        return "same"
    if isinstance(alpha_h, str):
        text = alpha_h.strip().lower()
        if text == "same":
            return "same"
        try:
            alpha_h = float(text)
        except ValueError as exc:
            raise IPyradError("sNMF alpha_H must be a non-negative float or 'same'.") from exc
    if float(alpha_h) < 0:
        raise IPyradError("sNMF alpha_H must be non-negative.")
    return float(alpha_h)


def _validate_snmf_configuration(
    *,
    k: int | None,
    k_range: str | None,
    nsamples: int | None,
    alpha_w: float,
    alpha_h: float | str,
    l1_ratio: float,
    n_init: int,
    cv_replicates: int,
    cv_holdout: float,
) -> None:
    """Validate CLI and runtime sNMF configuration."""
    if (k is None) == (k_range is None):
        raise IPyradError("Specify exactly one of -k or --k-range.")
    if alpha_w < 0:
        raise IPyradError("sNMF alpha_W must be non-negative.")
    _normalize_alpha_h(alpha_h)
    if not 0.0 <= l1_ratio <= 1.0:
        raise IPyradError("sNMF l1_ratio must be between 0 and 1.")
    if n_init < 1:
        raise IPyradError("sNMF n_init must be at least 1.")
    if cv_replicates < 1:
        raise IPyradError("sNMF cv_replicates must be at least 1.")
    if not 0.0 < cv_holdout < 1.0:
        raise IPyradError("sNMF cv_holdout must be between 0 and 1.")
    if nsamples is None:
        return
    if k_range is not None:
        lower, upper = parse_k_range(k_range)
        if lower < 2:
            raise IPyradError("K ranges for sNMF must start at 2 or greater.")
        if upper >= nsamples:
            raise IPyradError("Maximum K must be smaller than the number of retained samples.")
    else:
        if k is None or k < 2:
            raise IPyradError("K for sNMF must be 2 or greater.")
        if k >= nsamples:
            raise IPyradError("K must be smaller than the number of retained samples.")


def _encode_genotypes_disjunctive(matrix: np.ndarray) -> np.ndarray:
    """Expand a diploid genotype dosage matrix into three indicator columns per SNP."""
    if matrix.ndim != 2:
        raise IPyradError("sNMF requires a 2D genotype matrix after imputation.")
    if np.any((matrix < 0) | (matrix > 2)):
        raise IPyradError("sNMF genotype encoding expects diploid genotypes in {0, 1, 2}.")
    encoded = np.eye(3, dtype=np.float64)[matrix.astype(np.int64, copy=False)]
    return encoded.reshape(matrix.shape[0], matrix.shape[1] * 3)


def _normalize_membership(weights: np.ndarray) -> np.ndarray:
    """Normalize NMF weights into ancestry coefficients that sum to one per sample."""
    row_sums = weights.sum(axis=1, keepdims=True)
    membership = np.divide(
        weights,
        row_sums,
        out=np.full_like(weights, 1.0 / max(weights.shape[1], 1), dtype=np.float64),
        where=row_sums != 0,
    )
    return membership


def _normalize_genotype_frequencies(components: np.ndarray) -> np.ndarray:
    """Normalize NMF components into valid genotype-frequency blocks per SNP."""
    if components.shape[1] % 3 != 0:
        raise IPyradError("sNMF components do not align to 3-state genotype blocks.")
    blocks = components.reshape(components.shape[0], components.shape[1] // 3, 3)
    sums = blocks.sum(axis=2, keepdims=True)
    genotype_frequencies = np.divide(
        blocks,
        sums,
        out=np.full_like(blocks, 1.0 / 3.0, dtype=np.float64),
        where=sums != 0,
    )
    return genotype_frequencies


def _derive_allele_frequencies(genotype_frequencies: np.ndarray) -> np.ndarray:
    """Convert genotype probabilities to derived-allele frequencies."""
    return genotype_frequencies[:, :, 1] / 2.0 + genotype_frequencies[:, :, 2]


def _init_schedule(seed: int | None, n_init: int) -> list[tuple[str, int]]:
    """Return one deterministic NNDSVD start plus optional random restarts."""
    base = int(seed) if seed is not None else int(np.random.default_rng().integers(2**31))
    rng = np.random.default_rng(base)
    schedule = [("nndsvda", base)]
    for _ in range(max(n_init - 1, 0)):
        schedule.append(("random", int(rng.integers(2**31))))
    return schedule


def _fit_snmf_once(
    encoded: np.ndarray,
    *,
    k: int,
    init: str,
    seed: int,
    alpha_w: float,
    alpha_h: float | str,
    l1_ratio: float,
) -> SNMFFit:
    """Fit one sparse NMF model and return normalized ancestry outputs."""
    decomposition, _KMeans, _LDA, _TSNE = require_sklearn()
    from sklearn.exceptions import ConvergenceWarning

    model = decomposition.NMF(
        n_components=k,
        init=init,
        solver="cd",
        random_state=seed,
        tol=_NMF_TOL,
        max_iter=_MAX_ITER,
        alpha_W=alpha_w,
        alpha_H=alpha_h,
        l1_ratio=l1_ratio,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        weights = model.fit_transform(encoded)
    membership = _normalize_membership(weights.astype(np.float64, copy=False))
    genotype_frequencies = _normalize_genotype_frequencies(
        model.components_.astype(np.float64, copy=False)
    )
    return SNMFFit(
        membership=membership,
        genotype_frequencies=genotype_frequencies,
        allele_frequencies=_derive_allele_frequencies(genotype_frequencies),
        reconstruction_err=float(model.reconstruction_err_),
        n_iter=int(model.n_iter_),
        hit_max_iter=(
            int(model.n_iter_) >= _MAX_ITER
            or any(issubclass(w.category, ConvergenceWarning) for w in caught)
        ),
    )


def _fit_snmf(
    matrix: np.ndarray,
    *,
    k: int,
    seed: int | None,
    alpha_w: float,
    alpha_h: float | str,
    l1_ratio: float,
    n_init: int,
) -> SNMFFit:
    """Fit sNMF across one or more initializations and keep the best objective."""
    encoded = _encode_genotypes_disjunctive(matrix)
    best_fit = None
    capped_init_count = 0
    for init, fit_seed in _init_schedule(seed, n_init):
        fit = _fit_snmf_once(
            encoded,
            k=k,
            init=init,
            seed=fit_seed,
            alpha_w=alpha_w,
            alpha_h=alpha_h,
            l1_ratio=l1_ratio,
        )
        logger.debug(
            "sNMF K={} init={} seed={} reconstruction_err={} n_iter={}",
            k,
            init,
            fit_seed,
            fit.reconstruction_err,
            fit.n_iter,
        )
        capped_init_count += int(fit.hit_max_iter)
        if best_fit is None or fit.reconstruction_err < best_fit.reconstruction_err:
            best_fit = fit
    if best_fit is None:
        raise IPyradError("sNMF failed to produce any NMF fit.")
    logger.debug(
        "sNMF K={} capped initializations={}/{} selected_fit_capped={}",
        k,
        capped_init_count,
        n_init,
        best_fit.hit_max_iter,
    )
    return best_fit


def _log_k_convergence_warning(
    *,
    k: int,
    capped_fit_count: int,
    total_fit_count: int,
    selected_fit_hit_max_iter: bool,
) -> None:
    """Emit at most one convergence warning for one evaluated K."""
    if capped_fit_count == 0:
        return
    suffix = " (selected fit also capped)." if selected_fit_hit_max_iter else "."
    logger.warning(
        "sNMF K={} convergence warning: {}/{} fits reached max_iter={}{}",
        k,
        capped_fit_count,
        total_fit_count,
        _MAX_ITER,
        suffix,
    )


def _predict_genotype_probabilities(
    membership: np.ndarray,
    genotype_frequencies: np.ndarray,
    heldout: np.ndarray,
) -> np.ndarray:
    """Predict genotype-state probabilities for held-out sample/SNP pairs."""
    sample_indices = heldout[:, 0]
    snp_indices = heldout[:, 1]
    q = membership[sample_indices, :]
    g = np.moveaxis(genotype_frequencies[:, snp_indices, :], 0, 1)
    predicted = np.einsum("mk,mkg->mg", q, g)
    sums = predicted.sum(axis=1, keepdims=True)
    return np.divide(
        predicted,
        sums,
        out=np.full_like(predicted, 1.0 / 3.0, dtype=np.float64),
        where=sums != 0,
    )


def _mean_cross_entropy(truth: np.ndarray, predicted: np.ndarray) -> float:
    """Return mean genotype cross-entropy against predicted genotype probabilities."""
    probs = predicted[np.arange(truth.size), truth.astype(np.int64, copy=False)]
    probs = np.clip(probs, _EPSILON, 1.0)
    return float(-np.mean(np.log(probs)))


def _cross_entropy_for_k(
    raw_matrix: np.ndarray,
    extracter,
    *,
    k: int,
    impute_method: str | None,
    seed: int | None,
    alpha_w: float,
    alpha_h: float | str,
    l1_ratio: float,
    n_init: int,
    cv_replicates: int,
    holdout_fraction: float,
) -> SNMFCrossEntropyScore:
    """Score one K by masked-genotype cross-entropy."""
    rng = np.random.default_rng(seed)
    observed = np.argwhere(raw_matrix != _MISSING_GENO)
    if observed.size == 0:
        raise IPyradError("sNMF cross-entropy requires at least one observed genotype.")

    scores = []
    capped_fit_count = 0
    for rep in range(cv_replicates):
        nholdout = max(1, int(np.ceil(observed.shape[0] * holdout_fraction)))
        chosen = rng.choice(observed.shape[0], size=nholdout, replace=False)
        heldout = observed[chosen]

        training = raw_matrix.copy()
        training[heldout[:, 0], heldout[:, 1]] = _MISSING_GENO
        imputed = impute_genotype_matrix(
            training,
            extracter,
            impute_method=impute_method,
            random_seed=int(rng.integers(2**31)),
        )
        fit = _fit_snmf(
            imputed,
            k=k,
            seed=int(rng.integers(2**31)),
            alpha_w=alpha_w,
            alpha_h=alpha_h,
            l1_ratio=l1_ratio,
            n_init=n_init,
        )
        capped_fit_count += int(fit.hit_max_iter)
        if fit.hit_max_iter:
            logger.debug(
                "sNMF K={} CV replicate {} reached max_iter={}",
                k,
                rep,
                _MAX_ITER,
            )
        truth = raw_matrix[heldout[:, 0], heldout[:, 1]]
        predicted = _predict_genotype_probabilities(
            fit.membership,
            fit.genotype_frequencies,
            heldout,
        )
        scores.append(_mean_cross_entropy(truth, predicted))

    return SNMFCrossEntropyScore(
        mean_cross_entropy=float(np.mean(scores)),
        sd_cross_entropy=float(np.std(scores)),
        capped_fit_count=capped_fit_count,
        total_fit_count=cv_replicates,
    )


def run_snmf_method(
    *,
    data: Path | str,
    name: str,
    outdir: Path | str,
    k: int | None,
    k_range: str | None,
    min_sample_coverage: float,
    max_sample_missing: float,
    min_minor_allele_frequency: float,
    imap,
    minmap,
    exclude,
    include_reference: bool,
    impute_method: str | None,
    subsample: bool,
    random_seed: int | None,
    cores: int,
    force: bool,
    alpha_w: float = _DEFAULT_ALPHA_W,
    alpha_h: float | str = _DEFAULT_ALPHA_H,
    l1_ratio: float = _DEFAULT_L1_RATIO,
    n_init: int = _DEFAULT_N_INIT,
    cv_replicates: int = _DEFAULT_CV_REPLICATES,
    cv_holdout: float = _DEFAULT_CV_HOLDOUT,
    log_level: str = "INFO",
    min_genotype_depth: int = 0,
    min_site_qual: float = 0.0,
) -> None:
    """CLI entrypoint for sparse NMF-style clustering with cross-entropy K scoring."""
    require_hdf5_input(data, "snmf")
    normalized_impute = normalize_impute_method(impute_method)
    normalized_alpha_h = _normalize_alpha_h(alpha_h)
    _validate_snmf_configuration(
        k=k,
        k_range=k_range,
        nsamples=None,
        alpha_w=alpha_w,
        alpha_h=normalized_alpha_h,
        l1_ratio=l1_ratio,
        n_init=n_init,
        cv_replicates=cv_replicates,
        cv_holdout=cv_holdout,
    )

    outdir = Path(outdir).expanduser().absolute()
    paths = {
        "membership": outdir / f"{name}.membership.tsv",
        "allele_frequencies": outdir / f"{name}.allele_frequencies.tsv",
        "assignments": outdir / f"{name}.assignments.tsv",
        "k_scan": outdir / f"{name}.k_scan.tsv",
        "sample_data_summary": outdir / f"{name}.sample_data_summary.tsv",
        "stats": outdir / f"{name}.stats.txt",
    }
    ensure_output_paths(paths.values(), force=force)
    outdir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "sNMF regularization: alpha_W={} alpha_H={} l1_ratio={} n_init={}",
        alpha_w,
        normalized_alpha_h,
        l1_ratio,
        n_init,
    )
    logger.info(
        "sNMF K scoring: cv_replicates={} cv_holdout={}",
        cv_replicates,
        cv_holdout,
    )

    extracter = run_snps_extracter_for_method(
        data=data,
        min_sample_coverage=min_sample_coverage,
        max_sample_missing=max_sample_missing,
        min_minor_allele_frequency=min_minor_allele_frequency,
        imap=imap,
        minmap=minmap,
        min_genotype_depth=min_genotype_depth,
        min_site_qual=min_site_qual,
        exclude=exclude,
        include_reference=include_reference,
        cores=cores,
        log_level=log_level,
    )
    _validate_snmf_configuration(
        k=k,
        k_range=k_range,
        nsamples=len(extracter.snames),
        alpha_w=alpha_w,
        alpha_h=normalized_alpha_h,
        l1_ratio=l1_ratio,
        n_init=n_init,
        cv_replicates=cv_replicates,
        cv_holdout=cv_holdout,
    )

    prepared = get_numerical_input(
        extracter,
        subsample=subsample,
        random_seed=random_seed,
        impute_method=normalized_impute,
        log_level=log_level,
    )
    if prepared.matrix.shape[1] == 0:
        raise IPyradError("sNMF requires at least one SNP after filtering.")
    log_snp_imputation_summary("snmf", prepared.imputation, subsample=subsample)
    log_snp_view_summary(
        "snmf",
        summarize_prepared_snp_view(extracter, prepared.view, subsample=subsample),
        view_label="prepared",
    )
    raw_matrix = prepared.view.genos.copy()

    if k_range is not None:
        lower, upper = parse_k_range(k_range)
        k_values = list(range(lower, upper + 1))
    else:
        k_values = [int(k)]

    run_rows = []
    selected = None
    seed_rng = np.random.default_rng(random_seed)
    for current_k in k_values:
        fit_seed = int(seed_rng.integers(2**31))
        cv_seed = int(seed_rng.integers(2**31))
        fit = _fit_snmf(
            prepared.matrix,
            k=current_k,
            seed=fit_seed,
            alpha_w=alpha_w,
            alpha_h=normalized_alpha_h,
            l1_ratio=l1_ratio,
            n_init=n_init,
        )
        cv_score = _cross_entropy_for_k(
            raw_matrix,
            extracter,
            k=current_k,
            impute_method=normalized_impute,
            seed=cv_seed,
            alpha_w=alpha_w,
            alpha_h=normalized_alpha_h,
            l1_ratio=l1_ratio,
            n_init=n_init,
            cv_replicates=cv_replicates,
            holdout_fraction=cv_holdout,
        )
        capped_fit_count = int(fit.hit_max_iter) + cv_score.capped_fit_count
        total_fit_count = 1 + cv_score.total_fit_count
        _log_k_convergence_warning(
            k=current_k,
            capped_fit_count=capped_fit_count,
            total_fit_count=total_fit_count,
            selected_fit_hit_max_iter=fit.hit_max_iter,
        )
        row = {
            "k": current_k,
            "mean_cross_entropy": cv_score.mean_cross_entropy,
            "sd_cross_entropy": cv_score.sd_cross_entropy,
            "best_reconstruction_err": fit.reconstruction_err,
            "best_n_iter": fit.n_iter,
            "selected": False,
        }
        run_rows.append(row)
        logger.info(
            "sNMF K={} mean_cross_entropy={} sd_cross_entropy={} reconstruction_err={}",
            current_k,
            cv_score.mean_cross_entropy,
            cv_score.sd_cross_entropy,
            fit.reconstruction_err,
        )
        if selected is None or cv_score.mean_cross_entropy < selected[0]["mean_cross_entropy"]:
            selected = (row, fit, capped_fit_count, total_fit_count)

    if selected is None:
        raise IPyradError("sNMF failed to evaluate any K values.")
    selected_row, selected_fit, selected_capped_fit_count, selected_total_fit_count = selected
    selected_row["selected"] = True
    selected_k = int(selected_row["k"])
    logger.info(
        "selected sNMF K={} by minimum mean cross-entropy {}",
        selected_k,
        selected_row["mean_cross_entropy"],
    )

    write_membership(paths["membership"], extracter.snames, selected_fit.membership)
    write_assignments(paths["assignments"], extracter.snames, selected_fit.membership)
    write_marker_cluster_matrix(
        paths["allele_frequencies"],
        marker_ids_from_view(prepared.view),
        selected_fit.allele_frequencies,
    )
    pd.DataFrame.from_records(run_rows).to_csv(paths["k_scan"], sep="\t", index=False)
    sample_summary = build_imputed_sample_data_summary(
        samples=extracter.snames,
        matrix=prepared.view.genos,
        impute_method=normalized_impute,
    )
    write_sample_data_summary(paths["sample_data_summary"], sample_summary)
    write_stats_file(
        paths["stats"],
        tool="snmf",
        extracter=extracter,
        subsample=subsample,
        random_seed=random_seed,
        impute_method=impute_method,
        summary={
            "k_selected": selected_k,
            "k_range": k_range if k_range is not None else "NA",
            "selected_cross_entropy": float(selected_row["mean_cross_entropy"]),
            "selected_cross_entropy_sd": float(selected_row["sd_cross_entropy"]),
            "alpha_W": alpha_w,
            "alpha_H": normalized_alpha_h,
            "l1_ratio": l1_ratio,
            "n_init": n_init,
            "cv_replicates": cv_replicates,
            "cv_holdout": cv_holdout,
            "nmf_tol": _NMF_TOL,
            "nmf_max_iter": _MAX_ITER,
            "best_reconstruction_err": selected_fit.reconstruction_err,
            "best_n_iter": selected_fit.n_iter,
            "selected_fit_hit_max_iter": selected_fit.hit_max_iter,
            "capped_fits_for_selected_k": selected_capped_fit_count,
            "total_fits_for_selected_k": selected_total_fit_count,
            "linked_post_filter_snps": int(extracter.stats["post_filter_snps"]),
            "linked_post_filter_snp_containing_linkage_blocks": int(
                extracter.stats["post_filter_snp_containing_linkage_blocks"]
            ),
            "exported_snps": int(
                count_linkage_blocks(prepared.view)
                if subsample
                else prepared.view.snpsmap.shape[0]
            ),
            "exported_snp_containing_linkage_blocks": count_linkage_blocks(prepared.view),
            "samples_retained": len(extracter.snames),
        },
    )
    logger.debug(
        "sNMF output files: {}",
        {key: str(path) for key, path in paths.items()},
    )
    logger.info("wrote sNMF outputs to {}", outdir)


__all__ = ["run_snmf_method"]
