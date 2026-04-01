#!/usr/bin/env python

"""SNP-HDF5 backend for genome-wide population-genetic statistics."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from loguru import logger

from ....utils.exceptions import IPyradError
from ...extractors.snp_extractor import SNPExtractor
from ..common import build_sample_data_summary
from ..common import calculate_sample_missing_fraction
from .common import build_pairwise_stats_dataframe
from .common import build_population_stats_dataframe
from .common import build_sfs_dataframe
from .common import build_global_stats_dataframe
from .common import build_sample_stats_dataframe
from .estimators import fis_from_heterozygosity
from .models import PopgenResult
from .estimators import SNP_STATS
from .estimators import hudson_fst_components_block
from .estimators import summarize_genotype_block


def run_snp_popgen(
    *,
    data: Path,
    requested_stats: list[str],
    min_sample_coverage: float,
    max_sample_missing: float,
    min_minor_allele_frequency: float,
    imap,
    minmap,
    exclude,
    include_reference: bool,
    subsample_unlinked: bool,
    random_seed: int | None,
    cores: int,
    log_level: str,
) -> PopgenResult:
    """Compute SNP-backed population-genetic summaries from filtered genotypes."""
    unsupported = sorted(set(requested_stats).difference(SNP_STATS))
    if unsupported:
        raise IPyradError(
            "The following popgen statistics require sequence HDF5 in this phase: "
            + ", ".join(unsupported)
        )

    tool = SNPExtractor(
        data=Path(data),
        min_sample_coverage=min_sample_coverage,
        max_sample_missing=max_sample_missing,
        min_minor_allele_frequency=min_minor_allele_frequency,
        imap=imap,
        minmap=minmap,
        exclude=exclude,
        include_reference=include_reference,
        cores=cores,
    )
    tool.run(log_level=log_level)
    view = tool.get_view(
        subsample=subsample_unlinked,
        random_seed=random_seed,
        log_level=log_level,
    )

    global_stats_df = None
    if "fit" in requested_stats:
        global_summary = summarize_genotype_block(
            view.genos,
            include_minor_allele_count=False,
        )
        global_valid = (
            (global_summary.called_samples >= int(min_sample_coverage))
            & (global_summary.chromosome_count >= 2)
        )
        global_sites_used = int(np.sum(global_valid))
        global_observed = (
            float(np.nansum(global_summary.observed_heterozygosity[global_valid])) / global_sites_used
            if global_sites_used
            else np.nan
        )
        global_expected = (
            float(np.nansum(global_summary.expected_heterozygosity[global_valid])) / global_sites_used
            if global_sites_used
            else np.nan
        )
        global_stats_df = build_global_stats_dataframe(
            sites_used_heterozygosity=global_sites_used,
            observed_heterozygosity=global_observed,
            expected_heterozygosity_total=global_expected,
        )

    pop_names = list(tool.imap)
    sample_index = {name: idx for idx, name in enumerate(tool.snames)}
    pop_indices = {
        pop: np.array([sample_index[name] for name in names], dtype=np.int64)
        for pop, names in tool.imap.items()
    }
    pop_stats = {
        pop: {
            "population": pop,
            "n_samples": len(names),
            "sites_used_heterozygosity": 0,
            "observed_heterozygosity_sum": 0.0,
            "expected_heterozygosity_sum": 0.0,
        }
        for pop, names in tool.imap.items()
    }
    pairwise_keys = [
        (pop1, pop2)
        for idx, pop1 in enumerate(pop_names)
        for pop2 in pop_names[idx + 1 :]
    ]
    pairwise_stats = {
        key: {
            "population1": key[0],
            "population2": key[1],
            "sites_used": 0,
            "fst_num_sum": 0.0,
            "fst_den_sum": 0.0,
        }
        for key in pairwise_keys
    }
    max_minor_allele_count = (
        max((2 * len(names) for names in tool.imap.values()), default=0)
    )
    sfs_counts = np.zeros((len(pop_names), max_minor_allele_count + 1), dtype=np.int64)
    pop_summaries = {}

    for pop_idx, pop in enumerate(pop_names):
        summary = summarize_genotype_block(
            view.genos[pop_indices[pop], :],
            include_minor_allele_count="sfs" in requested_stats,
        )
        pop_summaries[pop] = summary
        valid = (
            (summary.called_samples >= int(tool.minmap[pop]))
            & (summary.chromosome_count >= 2)
        )
        if any(stat in requested_stats for stat in ("heterozygosity", "fis")):
            pop_stats[pop]["sites_used_heterozygosity"] = int(np.sum(valid))
            pop_stats[pop]["observed_heterozygosity_sum"] = float(
                np.nansum(summary.observed_heterozygosity[valid])
            )
            pop_stats[pop]["expected_heterozygosity_sum"] = float(
                np.nansum(summary.expected_heterozygosity[valid])
            )
        if "sfs" in requested_stats and summary.minor_allele_count is not None:
            sfs_valid = valid & summary.biallelic
            if np.any(sfs_valid):
                counts = np.bincount(
                    summary.minor_allele_count[sfs_valid],
                    minlength=sfs_counts.shape[1],
                )
                sfs_counts[pop_idx, : counts.size] += counts

    for pop1, pop2 in pairwise_keys:
        summary1 = pop_summaries[pop1]
        summary2 = pop_summaries[pop2]
        valid = (
            (summary1.called_samples >= int(tool.minmap[pop1]))
            & (summary2.called_samples >= int(tool.minmap[pop2]))
            & (summary1.chromosome_count > 0)
            & (summary2.chromosome_count > 0)
        )
        if not np.any(valid):
            continue
        pairwise_stats[(pop1, pop2)]["sites_used"] = int(np.sum(valid))
        fst_num, fst_den = hudson_fst_components_block(
            summary1.allele_counts[valid],
            summary2.allele_counts[valid],
            summary1.pi[valid],
            summary2.pi[valid],
        )
        pairwise_stats[(pop1, pop2)]["fst_num_sum"] = float(np.nansum(fst_num))
        pairwise_stats[(pop1, pop2)]["fst_den_sum"] = float(np.nansum(fst_den))

    population_rows = []
    for pop in pop_names:
        row = {
            "population": pop,
            "n_samples": pop_stats[pop]["n_samples"],
        }
        if any(stat in requested_stats for stat in ("heterozygosity", "fis")):
            sites_used = pop_stats[pop]["sites_used_heterozygosity"]
            observed = (
                pop_stats[pop]["observed_heterozygosity_sum"] / sites_used if sites_used else np.nan
            )
            expected = (
                pop_stats[pop]["expected_heterozygosity_sum"] / sites_used if sites_used else np.nan
            )
            row["sites_used_heterozygosity"] = sites_used
            row["observed_heterozygosity"] = observed
            row["expected_heterozygosity"] = expected
            if "fis" in requested_stats:
                row["fis"] = fis_from_heterozygosity(observed, expected)
        population_rows.append(row)
    population_df = build_population_stats_dataframe(
        population_rows,
        requested_stats=requested_stats,
        include_window_metadata=False,
    )

    pairwise_rows = []
    if "fst" in requested_stats and pairwise_keys:
        for pop1, pop2 in pairwise_keys:
            denom = pairwise_stats[(pop1, pop2)]["fst_den_sum"]
            pairwise_rows.append(
                {
                    "population1": pop1,
                    "population2": pop2,
                    "sites_used": pairwise_stats[(pop1, pop2)]["sites_used"],
                    "fst": (
                        pairwise_stats[(pop1, pop2)]["fst_num_sum"] / denom
                        if denom > 0
                        else np.nan
                    ),
                }
            )
    pairwise_df = build_pairwise_stats_dataframe(
        pairwise_rows,
        requested_stats=requested_stats,
        include_window_metadata=False,
    )

    sfs_rows = []
    if "sfs" in requested_stats:
        for pop_idx, pop in enumerate(pop_names):
            for mac, count in enumerate(sfs_counts[pop_idx]):
                if count <= 0:
                    continue
                sfs_rows.append(
                    {
                        "population": pop,
                        "minor_allele_count": mac,
                        "site_count": int(count),
                    }
                )
    sfs_df = build_sfs_dataframe(sfs_rows)

    missing_fraction = calculate_sample_missing_fraction(view.genos, tool.snames)
    sample_data_summary = build_sample_data_summary(
        samples=tool.snames,
        missing_fraction=missing_fraction,
        post_imputation_missing_fraction=missing_fraction,
        imputation_algorithm="not-imputed",
    )
    missing_mask = view.genos == 255
    sites_missing = np.sum(missing_mask, axis=1).astype(np.int64, copy=False)
    sites_total = int(view.genos.shape[1])
    sites_called = np.full(len(tool.snames), sites_total, dtype=np.int64) - sites_missing
    heterozygous_sites = np.sum(view.genos == 1, axis=1).astype(np.int64, copy=False)
    sample_stats = build_sample_stats_dataframe(
        samples=tool.snames,
        imap=tool.imap,
        sites_total=sites_total,
        sites_called=sites_called,
        sites_missing=sites_missing,
        heterozygous_sites=heterozygous_sites,
    )
    summary = {
        "input_backend": "snp",
        "requested_stats": list(requested_stats),
        "samples_selected_initial": tool.initial_snames,
        "samples_dropped_by_max_missing": tool.dropped_samples_by_missing,
        "samples_final": tool.snames,
        "imap": tool.imap,
        "minmap": tool.minmap,
        "include_reference": tool.include_reference,
        "max_sample_missing": tool.max_sample_missing,
        "subsample_unlinked": subsample_unlinked,
        "random_seed": random_seed,
        "min_minor_allele_frequency": min_minor_allele_frequency,
        "linked_post_filter_snps": int(tool.stats["post_filter_snps"]),
        "exported_snps": int(view.genos.shape[1]),
    }
    logger.info("computed SNP-backed popgen statistics for {} population(s)", len(tool.imap))
    return PopgenResult(
        sample_data_summary=sample_data_summary,
        sample_stats=sample_stats,
        global_stats=global_stats_df,
        population_stats=population_df,
        pairwise_stats=pairwise_df,
        sfs=sfs_df,
        summary=summary,
    )
