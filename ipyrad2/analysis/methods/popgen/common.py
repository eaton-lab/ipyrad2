#!/usr/bin/env python

"""Shared table builders for popgen outputs."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ....utils.exceptions import IPyradError
from .estimators import fit_from_heterozygosity


def build_sample_stats_dataframe(
    *,
    samples: list[str],
    imap: dict[str, list[str]],
    sites_total: int | np.ndarray | list[int],
    sites_called: np.ndarray | list[int],
    sites_missing: np.ndarray | list[int],
    heterozygous_sites: np.ndarray | list[int],
) -> pd.DataFrame:
    """Build the shared per-sample population-genetic summary table."""
    nsamples = len(samples)
    sample_to_population = {}
    for population, members in imap.items():
        for sample in members:
            if sample in sample_to_population:
                raise IPyradError(
                    f"Sample {sample!r} appears multiple times in popgen IMAP groups."
                )
            sample_to_population[sample] = population

    populations = [sample_to_population.get(sample, "all") for sample in samples]
    if any(sample not in sample_to_population for sample in samples):
        missing = [sample for sample in samples if sample not in sample_to_population]
        raise IPyradError(
            "Per-sample popgen summaries require every retained sample to belong to one "
            f"population. Missing assignments: {', '.join(missing)}"
        )

    called = np.asarray(sites_called, dtype=np.int64)
    missing = np.asarray(sites_missing, dtype=np.int64)
    heterozygous = np.asarray(heterozygous_sites, dtype=np.int64)
    if np.isscalar(sites_total):
        total = np.full(nsamples, int(sites_total), dtype=np.int64)
    else:
        total = np.asarray(sites_total, dtype=np.int64)

    for label, values in (
        ("sites_total", total),
        ("sites_called", called),
        ("sites_missing", missing),
        ("heterozygous_sites", heterozygous),
    ):
        if values.ndim != 1 or values.shape[0] != nsamples:
            raise IPyradError(
                f"Per-sample popgen {label} must have one value per retained sample."
            )

    observed = np.full(nsamples, np.nan, dtype=float)
    valid = called > 0
    observed[valid] = heterozygous[valid] / called[valid]
    called_fraction = np.zeros(nsamples, dtype=float)
    missing_fraction = np.zeros(nsamples, dtype=float)
    valid_total = total > 0
    called_fraction[valid_total] = called[valid_total] / total[valid_total]
    missing_fraction[valid_total] = missing[valid_total] / total[valid_total]
    homozygous = called - heterozygous
    return pd.DataFrame(
        {
            "sample": samples,
            "population": populations,
            "sites_total": total,
            "sites_called": called,
            "called_fraction": called_fraction,
            "sites_missing": missing,
            "missing_fraction": missing_fraction,
            "homozygous_sites": homozygous,
            "heterozygous_sites": heterozygous,
            "observed_heterozygosity": observed,
        }
    )


def build_global_stats_dataframe(
    *,
    sites_used_heterozygosity: int,
    observed_heterozygosity: float,
    expected_heterozygosity_total: float,
) -> pd.DataFrame:
    """Build the shared one-row global population-genetic summary table."""
    sites_used = int(sites_used_heterozygosity)
    if sites_used < 0:
        raise IPyradError("Global popgen sites_used_heterozygosity cannot be negative.")
    return pd.DataFrame(
        {
            "sites_used_heterozygosity": [sites_used],
            "observed_heterozygosity": [observed_heterozygosity],
            "expected_heterozygosity_total": [expected_heterozygosity_total],
            "fit": [
                fit_from_heterozygosity(
                    observed_heterozygosity,
                    expected_heterozygosity_total,
                )
            ],
        }
    )


def build_population_stats_dataframe(
    rows: list[dict[str, object]],
    *,
    requested_stats: list[str],
    include_window_metadata: bool,
) -> pd.DataFrame | None:
    """Return one ordered population-stats table from row dictionaries."""
    if not rows:
        return None
    columns: list[str] = []
    if include_window_metadata:
        columns.extend(
            [
                "window_id",
                "window_mode",
                "scaffold",
                "start",
                "end",
                "first_locus",
                "last_locus",
                "nloci",
                "sites_total",
            ]
        )
    columns.extend(["population", "n_samples"])
    if "pi" in requested_stats:
        columns.extend(["sites_used_pi", "pi"])
    if "theta_w" in requested_stats or "tajima_d" in requested_stats:
        columns.extend(["sites_used_theta", "segregating_sites"])
        if "theta_w" in requested_stats:
            columns.append("theta_w")
        if "tajima_d" in requested_stats:
            columns.append("tajima_d")
    if any(stat in requested_stats for stat in ("heterozygosity", "fis")):
        columns.extend(
            [
                "sites_used_heterozygosity",
                "observed_heterozygosity",
                "expected_heterozygosity",
            ]
        )
        if "fis" in requested_stats:
            columns.append("fis")
    return pd.DataFrame(rows, columns=columns)


def build_pairwise_stats_dataframe(
    rows: list[dict[str, object]],
    *,
    requested_stats: list[str],
    include_window_metadata: bool,
) -> pd.DataFrame | None:
    """Return one ordered pairwise-stats table from row dictionaries."""
    if not rows:
        return None
    columns: list[str] = []
    if include_window_metadata:
        columns.extend(
            [
                "window_id",
                "window_mode",
                "scaffold",
                "start",
                "end",
                "first_locus",
                "last_locus",
                "nloci",
                "sites_total",
            ]
        )
    columns.extend(["population1", "population2", "sites_used"])
    if "dxy" in requested_stats:
        columns.append("dxy")
    if "fst" in requested_stats:
        columns.append("fst")
    return pd.DataFrame(rows, columns=columns)


def build_sfs_dataframe(rows: list[dict[str, object]]) -> pd.DataFrame | None:
    """Return one ordered folded-SFS table."""
    if not rows:
        return None
    return pd.DataFrame(rows, columns=["population", "minor_allele_count", "site_count"])
