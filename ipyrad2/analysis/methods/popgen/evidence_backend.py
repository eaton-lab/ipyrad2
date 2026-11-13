"""Evidence-aware summaries of a fixed biallelic SNP catalogue."""
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

from ....utils.exceptions import IPyradError
from ...evidence import GenotypeEvidenceReader
from ...extracters.snps_extracter import SNPsExtracter
from ..common import build_sample_data_summary
from .likelihoods import likelihoods_from_pl, posterior_moments
from .models import PopgenResult
from .sfs_fitting import CrossFittedSFS


def _check_individuals(samples, mapping):
    if mapping is None:
        return dict(zip(samples, samples))
    if not isinstance(mapping, dict):
        pairs = []
        for line in Path(mapping).read_text().splitlines():
            if line.strip() and not line.lstrip().startswith('#'):
                fields = line.split()
                if len(fields) != 2:
                    raise IPyradError('Individual map requires sample and individual columns.')
                pairs.append(fields)
        if len({p[0] for p in pairs}) != len(pairs):
            raise IPyradError('Individual map contains duplicate sample names.')
        mapping = dict(pairs)
    if any(s not in mapping for s in samples):
        raise IPyradError('Individual map must identify every selected sample.')
    result = {s: str(mapping[s]) for s in samples}
    if len(set(result.values())) != len(result):
        raise IPyradError('Likelihood popgen requires one row per individual; combine replicate assays before this analysis.')
    return result


def run_evidence_popgen(*, data, requested_stats, min_sample_coverage,
                       max_sample_missing, min_minor_allele_frequency, imap,
                       minmap, exclude, include_reference, subsample_unlinked,
                       random_seed, cores, log_level, min_genotype_depth=0,
                       min_site_qual=0., evidence_prior='uniform', individual_map=None,
                       sfs_fit_loci=1024, sfs_folds=5, sfs_max_iter=10000):
    if set(requested_stats) != {'heterozygosity'}:
        raise IPyradError('Likelihood mode currently supports heterozygosity (including SNP-catalogue gene diversity), not pi/Fst/SFS or inbreeding statistics.')
    if evidence_prior not in ('uniform', 'population', 'sfs'):
        raise IPyradError('Unknown evidence prior.')
    if evidence_prior in ('population', 'sfs') and not imap:
        raise IPyradError('Population evidence prior requires an explicit --imap; groups are biological assumptions.')
    if include_reference:
        raise IPyradError('Likelihood popgen cannot use a synthetic reference sample.')
    if min_minor_allele_frequency:
        raise IPyradError('MAF filtering is not supported in likelihood mode; use a fixed SNP catalogue.')
    if not 0 <= max_sample_missing <= 1 or min_sample_coverage < 1 or min_genotype_depth < 0 or min_site_qual < 0:
        raise IPyradError('Invalid evidence coverage, missingness, depth, or QUAL filter.')
    # Reuse sample/group parsing only, never the hard-GT filtering path.
    tool = SNPsExtracter(data=Path(data), min_sample_coverage=min_sample_coverage,
                        max_sample_missing=max_sample_missing, min_minor_allele_frequency=0,
                        imap=imap, minmap=minmap, exclude=exclude, include_reference=False,
                        cores=cores, min_genotype_depth=0, min_site_qual=min_site_qual)
    samples = list(tool.snames)
    if not samples:
        raise IPyradError('No samples selected for likelihood popgen.')
    identities = _check_individuals(samples, individual_map)
    # SNPsExtracter expands sample patterns and enforces unique group membership.
    with GenotypeEvidenceReader(data) as reader:
        if any(s not in reader.names for s in samples):
            raise IPyradError('Likelihood popgen requires observed samples with native evidence.')
        nv = len(reader.snpsmap)
        candidates = []
        available = np.zeros(len(samples), dtype=np.int64)
        for start in range(0, nv, 512):
            sites = np.arange(start, min(nv, start+512))
            d = reader.read(sites, samples)
            keep = d['allele_count'] == 2
            if min_site_qual:
                qual = reader.io['site_qual'][sites]
                keep &= np.isfinite(qual) & (qual >= min_site_qual)
            good = d['usable_PL'] & (d['DP'] >= max(1, min_genotype_depth))
            available += good[:, keep].sum(axis=1)
            candidates.extend(sites[keep])
        candidates = np.asarray(candidates, dtype=np.int64)
        missing = 1-available/len(candidates) if len(candidates) else np.ones(len(samples))
        dropped = [s for s, rate in zip(samples, missing) if rate > max_sample_missing]
        samples = [s for s in samples if s not in dropped]
        if not samples:
            raise IPyradError('No samples remain after evidence missingness filtering.')
        populations = {g: [s for s in members if s in samples] for g, members in tool.imap.items()}
        populations = {g: members for g, members in populations.items() if members}
        indices = {g: np.array([samples.index(s) for s in members]) for g, members in populations.items()}
        selected = []
        for start in range(0, len(candidates), 512):
            sites = candidates[start:start+512]
            d = reader.read(sites, samples)
            good = d['usable_PL'] & (d['DP'] >= max(1, min_genotype_depth))
            keep = good.sum(axis=0) >= min_sample_coverage
            for group, ix in indices.items():
                keep &= good[ix].sum(axis=0) >= tool.minmap[group]
            selected.extend(sites[keep])
        selected = np.asarray(selected, dtype=np.int64)
        linked_count = len(selected)
        if subsample_unlinked and len(selected):
            rng = np.random.default_rng(random_seed)
            loci = reader.snpsmap[selected, 0]
            # First draw per randomly ordered site; sorted output preserves coordinates.
            shuffled = rng.permutation(len(selected))
            _, first = np.unique(loci[shuffled], return_index=True)
            selected = np.sort(selected[shuffled[first]])
        sfs = None
        if evidence_prior == 'sfs':
            locus_ids = reader.snpsmap[selected, 0]
            sfs = CrossFittedSFS(locus_ids, indices, fit_loci=sfs_fit_loci,
                                folds=sfs_folds, max_iter=sfs_max_iter)
            fitting = selected[np.isin(locus_ids, sfs.fit_loci)]
            for start in range(0, len(fitting), 512):
                sites = fitting[start:start+512]
                d = reader.read(sites, samples)
                good = d['usable_PL'] & (d['DP'] >= max(1, min_genotype_depth))
                likelihood, good = likelihoods_from_pl(d['PL'][..., :3], good)
                sfs.add(likelihood, good, reader.snpsmap[sites, 0])
            sfs.fit(cores)
        ns = len(samples)
        count = np.zeros(ns, dtype=np.int64)
        het = np.zeros(ns)
        dosage = np.zeros(ns)
        hard_count = np.zeros(ns, dtype=np.int64)
        hard_het = np.zeros(ns, dtype=np.int64)
        aggregates = {g: dict(sites=0, cells=0, het=0., diversity=0., frequency=0., frequency_variance=0.) for g in populations}
        for start in range(0, len(selected), 128):
            sites = selected[start:start+128]
            d = reader.read(sites, samples)
            good = d['usable_PL'] & (d['DP'] >= max(1, min_genotype_depth))
            likelihood, good = likelihoods_from_pl(d['PL'][..., :3], good)
            gt = d['GT']
            hard = good & np.all(gt != 255, axis=-1)
            hard_count += hard.sum(axis=1)
            hard_het += (hard & (gt[..., 0] != gt[..., 1])).sum(axis=1)
            for group, ix in indices.items():
                result = (sfs.moments(group, likelihood[ix], good[ix], reader.snpsmap[sites, 0])
                          if sfs is not None else posterior_moments(likelihood[ix], good[ix], evidence_prior))
                gp = result['GP']
                count[ix] += good[ix].sum(axis=1)
                het[ix] += np.nansum(gp[..., 1], axis=1)
                dosage[ix] += np.nansum(result['dosage'], axis=1)
                eligible = result['n_individuals'] >= max(1, tool.minmap[group])
                a = aggregates[group]
                a['sites'] += int(eligible.sum())
                a['cells'] += int(result['n_individuals'][eligible].sum())
                with np.errstate(invalid='ignore', divide='ignore'):
                    per_site_het = np.nansum(gp[..., 1], axis=0)/result['n_individuals']
                for key, value in [('het', per_site_het), ('diversity', result['gene_diversity']),
                                   ('frequency', result['allele_frequency']), ('frequency_variance', result['allele_frequency_variance'])]:
                    a[key] += float(np.nansum(value[eligible]))
    with np.errstate(invalid='ignore', divide='ignore'):
        hmean = het/count
        dsmean = dosage/count
        rates = 1-count/len(selected) if len(selected) else np.ones(ns)
    membership = {s: group for group, members in populations.items() for s in members}
    sample_stats = pd.DataFrame(dict(sample=samples, population=[membership[s] for s in samples],
        catalogue_sites=len(selected), sites_with_likelihood=count,
        posterior_expected_heterozygotes=het, posterior_heterozygosity=hmean,
        mean_posterior_alt_dosage=dsmean, canonical_called_sites_on_same_cells=hard_count,
        canonical_heterozygotes_on_same_cells=hard_het))
    rows = []
    for g, a in aggregates.items():
        denom = a['sites'] or np.nan
        rows.append(dict(population=g, n_samples=len(populations[g]), catalogue_sites_used=a['sites'],
            sample_cells_with_likelihood=a['cells'], posterior_heterozygosity=a['het']/denom,
            posterior_snp_gene_diversity=a['diversity']/denom,
            mean_posterior_alt_frequency=a['frequency']/denom,
            mean_posterior_frequency_variance=a['frequency_variance']/denom))
    summary = dict(input_backend='snp_likelihoods', evidence_prior=evidence_prior,
        scope='fixed biallelic SNP catalogue; not per-base diversity',
        model=('independent uniform genotype prior' if evidence_prior == 'uniform' else 'HWE conditional on frequency; Beta(1,1) frequency integrated jointly'),
        missing_policy='no posterior for zero-read, hard-masked, unmatched, or incomplete PL cells',
        genotype_filter='does not require a nonmissing hard GT',
        diversity_formula='E[2*K*(M-K)/(M*(M-1)) | reads], M=2*n usable individuals; includes joint count variance',
        uncertainty_note='frequency variance is per-site model uncertainty, not a locus-block confidence interval',
        samples_selected_initial=tool.initial_snames, samples_dropped_by_max_missing=dropped,
        samples_final=samples, individual_ids=identities, imap=populations, minmap=tool.minmap,
        original_catalogue_sites=nv, biallelic_qual_candidates=len(candidates),
        linked_post_filter_snps=linked_count, exported_snps=len(selected),
        subsample_unlinked=subsample_unlinked, random_seed=random_seed)
    if sfs is not None:
        summary['model'] = 'cross-fitted folded finite-sample SFS'
        summary['sfs_prior'] = sfs.metadata()
    logger.info('Computed evidence-aware SNP-catalogue summaries for {} populations', len(populations))
    return PopgenResult(sample_data_summary=build_sample_data_summary(samples=samples,
        missing_fraction=pd.Series(rates, index=samples), post_imputation_missing_fraction=pd.Series(rates, index=samples),
        imputation_algorithm='not-imputed'), sample_stats=sample_stats,
        population_stats=pd.DataFrame(rows), summary=summary, sfs=sfs.table() if sfs is not None else None)
