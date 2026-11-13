"""Per-callable-base diversity under an explicit biallelic locus model.

A reference/unseen native pileup record supplies uncertainty about one latent
alternative allele. At an observed biallelic record we condition on its two
nucleotides and use their three native PL entries. We do not reconstruct a
four-nucleotide likelihood vector from a symbolic allele.
"""
from collections import Counter

import numpy as np
import pandas as pd
from loguru import logger

from ....utils.exceptions import IPyradError
from ...rad_evidence import RadSiteEvidenceReader
from ...extracters.snps_extracter import SNPsExtracter
from ..common import build_sample_data_summary
from .evidence_backend import _check_individuals
from .likelihoods import likelihoods_from_pl, posterior_moments
from .models import PopgenResult
from .sfs_fitting import CrossFittedSFS


CHUNK_BASES = 2048


def biallelic_locus_chunk(reader, sites, samples, min_depth=1):
    """Select native PL triads, retaining one explicit exclusion per base.

    Returns one column per physical base (never per record). `kind` is 1 for
    reference/unseen, 2 for an observed nucleotide pair, and 0 if excluded.
    Complete native PL and the reader's original hard masks remain mandatory.
    """
    d = reader.read(sites, samples)
    ns, nv = len(samples), len(sites)
    pl = np.full((ns, nv, 3), -1, dtype=np.int32)
    good = np.zeros((ns, nv), dtype=bool)
    kind = np.zeros(nv, dtype=np.uint8)
    reasons = Counter()
    records = {}
    for j, site in enumerate(d['record_site']):
        records.setdefault(int(site), []).append(j)
    for col, site in enumerate(sites):
        js = records.get(int(site), [])
        if int(d['reference'][col]) not in (65, 67, 71, 84) or d['site_status'][col] == 2:
            reasons['non_acgt_reference'] += 1
            continue
        if not js:
            reasons['no_record'] += 1
            continue
        if len(js) != 1:
            reasons['multiple_records'] += 1
            continue
        j = js[0]
        alleles = d['alleles'][j]
        concrete = [i for i, a in enumerate(alleles) if a in ('A', 'C', 'G', 'T')]
        symbols = [i for i, a in enumerate(alleles) if a == '<*>']
        if (not alleles or alleles[0] != chr(int(d['reference'][col]))
                or len(set(alleles)) != len(alleles)
                or len(concrete) + len(symbols) != len(alleles) or len(symbols) > 1):
            reasons['indel_or_unsupported_alleles'] += 1
            continue
        if len(concrete) > 2:
            reasons['multiallelic_nucleotides'] += 1
            continue
        if len(concrete) == 1 and len(symbols) == 1:
            a, b = concrete[0], symbols[0]
            kind[col] = 1
        elif len(concrete) == 2:
            a, b = concrete
            kind[col] = 2
        else:
            # REF-only PL has no likelihood for a heterozygous alternative.
            reasons['no_alternative_likelihood'] += 1
            continue
        a, b = sorted((a, b))
        cols = [a*(a+1)//2+a, b*(b+1)//2+a, b*(b+1)//2+b]
        pl[:, col] = d['PL'][:, j, cols]
        good[:, col] = d['usable_PL'][:, j] & (d['DP'][:, j] >= max(1, min_depth))
    return pl, good, kind, reasons


def run_locus_evidence_popgen(*, data, requested_stats, min_sample_coverage,
                             max_sample_missing, min_minor_allele_frequency,
                             imap, minmap, exclude, include_reference,
                             subsample_unlinked, random_seed, cores, log_level,
                             min_genotype_depth=0, min_site_qual=0.,
                             evidence_prior='uniform', individual_map=None,
                             sfs_fit_loci=1024, sfs_folds=5, sfs_max_iter=10000):
    """Average posterior chromosome-pair differences over callable RAD bases."""
    if not set(requested_stats) <= {'pi', 'heterozygosity'}:
        raise IPyradError('Locus likelihoods currently support pi and heterozygosity only.')
    if evidence_prior not in ('uniform', 'population', 'sfs'):
        raise IPyradError('Unknown evidence prior.')
    if evidence_prior in ('population', 'sfs') and not imap:
        raise IPyradError('Population evidence prior requires an explicit --imap.')
    if include_reference:
        raise IPyradError('Locus likelihoods cannot use the synthetic reference sample.')
    if min_minor_allele_frequency or min_site_qual or subsample_unlinked or random_seed is not None:
        raise IPyradError('Locus likelihoods require all callable bases: MAF, variant QUAL, SNP subsampling and seed options are unsupported.')
    if not 0 <= max_sample_missing <= 1 or min_sample_coverage < 1 or min_genotype_depth < 0:
        raise IPyradError('Invalid locus evidence coverage, missingness or depth filter.')
    # Sample names, pattern expansion and population parsing only; no SNP data
    # or hard-call filters determine the all-base site domain.
    selection = SNPsExtracter(data=data, min_sample_coverage=min_sample_coverage,
        max_sample_missing=max_sample_missing, min_minor_allele_frequency=0,
        imap=imap, minmap=minmap, exclude=exclude, include_reference=False,
        cores=cores, min_genotype_depth=0, min_site_qual=0)
    initial = list(selection.snames)
    identities = _check_individuals(initial, individual_map)
    with RadSiteEvidenceReader(data) as reader:
        if not initial or any(s not in reader.names for s in initial):
            raise IPyradError('Locus likelihoods require observed samples with native evidence.')
        samples = initial
        dropped = []
        if max_sample_missing < 1:
            available = np.zeros(len(initial), dtype=np.int64)
            acgt = 0
            for start in range(0, reader.nsites, CHUNK_BASES):
                sites = np.arange(start, min(start+CHUNK_BASES, reader.nsites))
                _, good, _, _ = biallelic_locus_chunk(reader, sites, initial, min_genotype_depth)
                available += good.sum(axis=1)
                acgt += int(np.isin(reader.group['reference'][start:start+len(sites)], [65,67,71,84]).sum())
            missing = 1-available/acgt if acgt else np.ones(len(initial))
            dropped = [s for s, value in zip(initial, missing) if value > max_sample_missing]
            samples = [s for s in initial if s not in dropped]
        if not samples:
            raise IPyradError('No samples remain after locus-evidence missingness filtering.')
        populations = {g: [s for s in members if s in samples] for g, members in selection.imap.items()}
        populations = {g: members for g, members in populations.items() if members}
        ix = {g: np.array([samples.index(s) for s in members]) for g,members in populations.items()}
        sfs = None
        if evidence_prior == 'sfs':
            intervals = reader.group['intervals'][:]
            offsets = np.r_[0, np.cumsum(intervals[:, 2]-intervals[:, 1])]
            locus_ids = reader.group['site_map'][offsets[:-1], 2]
            sfs = CrossFittedSFS(locus_ids, ix, fit_loci=sfs_fit_loci,
                                folds=sfs_folds, max_iter=sfs_max_iter)
            for pos in np.flatnonzero(np.isin(locus_ids, sfs.fit_loci)):
                for start in range(int(offsets[pos]), int(offsets[pos+1]), CHUNK_BASES):
                    sites = np.arange(start, min(start+CHUNK_BASES, offsets[pos+1]))
                    pl, good, _, _ = biallelic_locus_chunk(reader, sites, samples, min_genotype_depth)
                    keep = good.sum(axis=0) >= min_sample_coverage
                    for group, indices in ix.items():
                        keep &= good[indices].sum(axis=0) >= selection.minmap[group]
                    if keep.any():
                        likelihood, good = likelihoods_from_pl(pl[:, keep], good[:, keep])
                        sfs.add(likelihood, good, np.full(keep.sum(), locus_ids[pos]))
            sfs.fit(cores)
        counts = np.zeros(len(samples), dtype=np.int64)
        het = np.zeros(len(samples))
        agg = {g: dict(callable_bases=0, callable_sample_cells=0,
                      callable_unseen_alt_bases=0, callable_observed_alt_bases=0,
                      posterior_pairwise_differences_sum=0., heterozygosity_sum=0.) for g in populations}
        exclusions = Counter()
        acgt = model_bases = selected_bases = 0
        for start in range(0, reader.nsites, CHUNK_BASES):
            sites = np.arange(start, min(start+CHUNK_BASES, reader.nsites))
            pl, good, kind, reasons = biallelic_locus_chunk(reader, sites, samples, min_genotype_depth)
            exclusions.update(reasons)
            acgt += int(np.isin(reader.group['reference'][start:start+len(sites)], [65,67,71,84]).sum())
            model_bases += int((kind > 0).sum())
            keep = good.sum(axis=0) >= min_sample_coverage
            for group, indices in ix.items():
                keep &= good[indices].sum(axis=0) >= selection.minmap[group]
            selected_bases += int(keep.sum())
            if not keep.any():
                continue
            likelihood, good = likelihoods_from_pl(pl[:, keep], good[:, keep])
            for group, indices in ix.items():
                moments = (sfs.moments(group, likelihood[indices], good[indices],
                           reader.group['site_map'][sites, 2][keep]) if sfs is not None
                           else posterior_moments(likelihood[indices], good[indices], evidence_prior))
                gp = moments['GP']
                counts[indices] += good[indices].sum(axis=1)
                het[indices] += np.nansum(gp[..., 1], axis=1)
                n = moments['n_individuals']
                valid = n >= max(1, selection.minmap[group])
                a = agg[group]
                a['callable_bases'] += int(valid.sum())
                a['callable_sample_cells'] += int(n[valid].sum())
                a['callable_unseen_alt_bases'] += int((valid & (kind[keep] == 1)).sum())
                a['callable_observed_alt_bases'] += int((valid & (kind[keep] == 2)).sum())
                a['posterior_pairwise_differences_sum'] += float(np.nansum(moments['gene_diversity'][valid]))
                a['heterozygosity_sum'] += float(np.sum(np.nansum(gp[...,1], axis=0)[valid]/n[valid]))
            if start % (CHUNK_BASES*100) == 0:
                logger.info('Locus likelihoods: {} / {} bases', start+len(sites), reader.nsites)
        total = reader.nsites
    with np.errstate(divide='ignore', invalid='ignore'):
        h = het/counts
    membership = {s:g for g,members in populations.items() for s in members}
    sample_stats = pd.DataFrame(dict(sample=samples, population=[membership[s] for s in samples],
        retained_locus_bases=total, acgt_locus_bases=acgt, bases_with_likelihood=counts,
        posterior_expected_heterozygous_bases=het, posterior_heterozygosity_per_base=h))
    rows = []
    for g,a in agg.items():
        denom = a['callable_bases'] or np.nan
        row = dict(population=g, n_samples=len(populations[g]), retained_locus_bases=total,
                   acgt_locus_bases=acgt, **a)
        row['posterior_pi'] = a['posterior_pairwise_differences_sum']/denom
        row['posterior_heterozygosity_per_base'] = row.pop('heterozygosity_sum')/denom
        rows.append(row)
    rates = 1-counts/acgt if acgt else np.ones(len(samples))
    summary = dict(input_backend='locus_likelihoods', evidence_scope='loci', evidence_prior=evidence_prior,
        scope='per callable base within final retained locus intervals; not whole-genome diversity',
        allele_model='biallelic: REF versus one latent unseen alternative at REF/<*> sites; conditional on the observed nucleotide pair at biallelic sites',
        symbolic_policy='No symbolic allele is assigned a nucleotide; no four-base likelihood reconstruction. Unseen third-allele states at observed biallelic sites are outside this conditional model.',
        prior_model=('independent uniform genotype prior' if evidence_prior=='uniform' else 'HWE conditional on frequency; Beta(1,1) frequency integrated jointly'),
        pi_formula='mean over callable bases of E[2*K*(M-K)/(M*(M-1)) | reads], M=2*n usable individuals; includes joint count covariance',
        denominator='Each eligible physical base contributes once per population with >=max(1,minmap) usable individuals, after global coverage/minmap filters. Missing, masked and excluded bases contribute neither numerator nor denominator.',
        singleton_note='With one usable diploid individual, pi equals posterior heterozygosity of its two chromosomes.',
        uncertainty_note='Model-based posterior mean, not an unbiasedness guarantee or confidence interval. Low-depth estimates can depend strongly on the prior.',
        missing_policy='positive DP, complete native PL, PASS status and no hard mask; no missing-cell posterior imputation',
        samples_selected_initial=initial, samples_dropped_by_max_missing=dropped, samples_final=samples,
        individual_ids=identities, imap=populations, minmap=selection.minmap,
        retained_locus_bases=total, acgt_locus_bases=acgt, eligible_model_bases=model_bases,
        globally_callable_bases=selected_bases, structural_exclusions=dict(exclusions),
        sample_missingness_denominator='all retained ACGT bases before global/population coverage filters')
    if sfs is not None:
        summary['prior_model'] = 'cross-fitted folded finite-sample SFS'
        summary['sfs_prior'] = sfs.metadata()
    return PopgenResult(sample_data_summary=build_sample_data_summary(samples=samples,
        missing_fraction=pd.Series(rates,index=samples),post_imputation_missing_fraction=pd.Series(rates,index=samples),
        imputation_algorithm='not-imputed'),sample_stats=sample_stats,
        population_stats=pd.DataFrame(rows),summary=summary, sfs=sfs.table() if sfs is not None else None)
