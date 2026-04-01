#!/usr/bin/env python

"""Estimator helpers for population-genetic statistics."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .models import ORDERED_STATS
from .models import SEQUENCE_ONLY_STATS
from .models import SEQUENCE_STATS
from .models import SNP_STATS
from .models import normalize_stats

_SEQ_ALLELE_COUNTS = np.zeros((256, 4), dtype=np.int16)
_SEQ_CALLED = np.zeros(256, dtype=bool)
_SEQ_HET = np.zeros(256, dtype=bool)

for base, counts, hetero in (
    ("A", (2, 0, 0, 0), False),
    ("C", (0, 2, 0, 0), False),
    ("G", (0, 0, 2, 0), False),
    ("T", (0, 0, 0, 2), False),
    ("R", (1, 0, 1, 0), True),
    ("Y", (0, 1, 0, 1), True),
    ("S", (0, 1, 1, 0), True),
    ("W", (1, 0, 0, 1), True),
    ("K", (0, 0, 1, 1), True),
    ("M", (1, 1, 0, 0), True),
):
    code = ord(base)
    _SEQ_ALLELE_COUNTS[code] = np.array(counts, dtype=np.int16)
    _SEQ_CALLED[code] = True
    _SEQ_HET[code] = hetero

_MISSING_SEQ_CODES = {ord("N"), ord("-")}
_MISSING_GENO = 255


@dataclass
class SiteSummary:
    """Derived per-site allele-count summary for one population."""

    called_samples: int
    chromosome_count: int
    allele_counts: np.ndarray
    observed_heterozygosity: float
    expected_heterozygosity: float
    pi: float
    segregating: bool
    biallelic: bool
    minor_allele_count: int | None


@dataclass
class BlockSummary:
    """Vectorized per-site allele-count summary for one alignment/genotype block."""

    called_samples: np.ndarray
    chromosome_count: np.ndarray
    allele_counts: np.ndarray
    observed_heterozygosity: np.ndarray
    expected_heterozygosity: np.ndarray
    pi: np.ndarray
    segregating: np.ndarray
    biallelic: np.ndarray
    minor_allele_count: np.ndarray | None

def gene_diversity_block(counts: np.ndarray) -> np.ndarray:
    """Return unbiased gene diversity for one `(nsites, nalleles)` count matrix."""
    if counts.ndim != 2:
        raise ValueError("gene_diversity_block expects a 2D count matrix")
    totals = np.sum(counts, axis=1)
    result = np.full(totals.shape, np.nan, dtype=np.float64)
    valid = totals > 1
    if not np.any(valid):
        return result

    valid_counts = counts[valid].astype(np.float64, copy=False)
    valid_totals = totals[valid].astype(np.float64, copy=False)
    freqs = valid_counts / valid_totals[:, None]
    base = 1.0 - np.sum(freqs * freqs, axis=1)
    result[valid] = (valid_totals / (valid_totals - 1.0)) * base
    return result


def summarize_sequence_site(values: np.ndarray) -> SiteSummary:
    """Return allele counts and diversity summary from one sequence column."""
    called_mask = _SEQ_CALLED[values]
    called_samples = int(np.sum(called_mask))
    if called_samples == 0:
        return SiteSummary(
            called_samples=0,
            chromosome_count=0,
            allele_counts=np.zeros(4, dtype=np.int64),
            observed_heterozygosity=float("nan"),
            expected_heterozygosity=float("nan"),
            pi=float("nan"),
            segregating=False,
            biallelic=False,
            minor_allele_count=None,
        )
    called = values[called_mask]
    counts = _SEQ_ALLELE_COUNTS[called].sum(axis=0).astype(np.int64, copy=False)
    hetero_count = int(np.sum(_SEQ_HET[called]))
    chroms = int(np.sum(counts))
    nonzero = counts[counts > 0]
    return SiteSummary(
        called_samples=called_samples,
        chromosome_count=chroms,
        allele_counts=counts,
        observed_heterozygosity=hetero_count / called_samples,
        expected_heterozygosity=gene_diversity(counts),
        pi=gene_diversity(counts),
        segregating=bool(nonzero.size > 1),
        biallelic=bool(nonzero.size == 2),
        minor_allele_count=int(np.min(nonzero)) if nonzero.size == 2 else None,
    )


def summarize_sequence_block(
    values: np.ndarray,
    *,
    include_minor_allele_count: bool,
) -> BlockSummary:
    """Return vectorized site summaries from one `(samples, sites)` sequence block."""
    if values.ndim != 2:
        raise ValueError("summarize_sequence_block expects a 2D block")
    if values.shape[1] == 0:
        empty = np.zeros(0, dtype=np.int64)
        empty_bool = np.zeros(0, dtype=bool)
        empty_float = np.zeros(0, dtype=np.float64)
        return BlockSummary(
            called_samples=empty.copy(),
            chromosome_count=empty.copy(),
            allele_counts=np.zeros((0, 4), dtype=np.int64),
            observed_heterozygosity=empty_float.copy(),
            expected_heterozygosity=empty_float.copy(),
            pi=empty_float.copy(),
            segregating=empty_bool.copy(),
            biallelic=empty_bool.copy(),
            minor_allele_count=empty.copy() if include_minor_allele_count else None,
        )

    called_samples = np.sum(_SEQ_CALLED[values], axis=0).astype(np.int64, copy=False)
    allele_counts = np.sum(_SEQ_ALLELE_COUNTS[values], axis=0).astype(np.int64, copy=False)
    hetero_count = np.sum(_SEQ_HET[values], axis=0).astype(np.int64, copy=False)
    chromosome_count = np.sum(allele_counts, axis=1).astype(np.int64, copy=False)
    pi = gene_diversity_block(allele_counts)
    expected_heterozygosity = pi.copy()

    observed_heterozygosity = np.full(called_samples.shape, np.nan, dtype=np.float64)
    called_mask = called_samples > 0
    observed_heterozygosity[called_mask] = (
        hetero_count[called_mask] / called_samples[called_mask]
    )

    positive = allele_counts > 0
    nonzero_count = np.sum(positive, axis=1)
    segregating = nonzero_count > 1
    biallelic = nonzero_count == 2

    minor_allele_count = None
    if include_minor_allele_count:
        masked = np.where(positive, allele_counts, np.iinfo(np.int64).max)
        minor_allele_count = np.min(masked, axis=1).astype(np.int64, copy=False)
        minor_allele_count[~biallelic] = 0

    return BlockSummary(
        called_samples=called_samples,
        chromosome_count=chromosome_count,
        allele_counts=allele_counts,
        observed_heterozygosity=observed_heterozygosity,
        expected_heterozygosity=expected_heterozygosity,
        pi=pi,
        segregating=segregating,
        biallelic=biallelic,
        minor_allele_count=minor_allele_count,
    )


def summarize_genotype_site(values: np.ndarray) -> SiteSummary:
    """Return allele counts and diversity summary from one diploid genotype column."""
    called_mask = values != _MISSING_GENO
    called = values[called_mask]
    called_samples = int(called.size)
    if called_samples == 0:
        return SiteSummary(
            called_samples=0,
            chromosome_count=0,
            allele_counts=np.zeros(2, dtype=np.int64),
            observed_heterozygosity=float("nan"),
            expected_heterozygosity=float("nan"),
            pi=float("nan"),
            segregating=False,
            biallelic=False,
            minor_allele_count=None,
        )
    hom_ref = int(np.sum(called == 0))
    het = int(np.sum(called == 1))
    hom_alt = int(np.sum(called == 2))
    counts = np.array([2 * hom_ref + het, 2 * hom_alt + het], dtype=np.int64)
    nonzero = counts[counts > 0]
    return SiteSummary(
        called_samples=called_samples,
        chromosome_count=int(np.sum(counts)),
        allele_counts=counts,
        observed_heterozygosity=het / called_samples,
        expected_heterozygosity=gene_diversity(counts),
        pi=gene_diversity(counts),
        segregating=bool(hom_ref and (het or hom_alt) or hom_alt and (het or hom_ref)),
        biallelic=bool(nonzero.size == 2),
        minor_allele_count=int(np.min(nonzero)) if nonzero.size == 2 else None,
    )


def summarize_genotype_block(
    values: np.ndarray,
    *,
    include_minor_allele_count: bool,
) -> BlockSummary:
    """Return vectorized site summaries from one `(samples, sites)` genotype block."""
    if values.ndim != 2:
        raise ValueError("summarize_genotype_block expects a 2D block")
    if values.shape[1] == 0:
        empty = np.zeros(0, dtype=np.int64)
        empty_bool = np.zeros(0, dtype=bool)
        empty_float = np.zeros(0, dtype=np.float64)
        return BlockSummary(
            called_samples=empty.copy(),
            chromosome_count=empty.copy(),
            allele_counts=np.zeros((0, 2), dtype=np.int64),
            observed_heterozygosity=empty_float.copy(),
            expected_heterozygosity=empty_float.copy(),
            pi=empty_float.copy(),
            segregating=empty_bool.copy(),
            biallelic=empty_bool.copy(),
            minor_allele_count=empty.copy() if include_minor_allele_count else None,
        )

    called_samples = np.sum(values != _MISSING_GENO, axis=0).astype(np.int64, copy=False)
    hom_ref = np.sum(values == 0, axis=0).astype(np.int64, copy=False)
    het = np.sum(values == 1, axis=0).astype(np.int64, copy=False)
    hom_alt = np.sum(values == 2, axis=0).astype(np.int64, copy=False)
    allele_counts = np.column_stack((2 * hom_ref + het, 2 * hom_alt + het)).astype(
        np.int64,
        copy=False,
    )
    chromosome_count = np.sum(allele_counts, axis=1).astype(np.int64, copy=False)
    pi = gene_diversity_block(allele_counts)
    expected_heterozygosity = pi.copy()

    observed_heterozygosity = np.full(called_samples.shape, np.nan, dtype=np.float64)
    called_mask = called_samples > 0
    observed_heterozygosity[called_mask] = het[called_mask] / called_samples[called_mask]

    positive = allele_counts > 0
    nonzero_count = np.sum(positive, axis=1)
    segregating = nonzero_count > 1
    biallelic = nonzero_count == 2

    minor_allele_count = None
    if include_minor_allele_count:
        minor_allele_count = np.min(allele_counts, axis=1).astype(np.int64, copy=False)
        minor_allele_count[~biallelic] = 0

    return BlockSummary(
        called_samples=called_samples,
        chromosome_count=chromosome_count,
        allele_counts=allele_counts,
        observed_heterozygosity=observed_heterozygosity,
        expected_heterozygosity=expected_heterozygosity,
        pi=pi,
        segregating=segregating,
        biallelic=biallelic,
        minor_allele_count=minor_allele_count,
    )


def gene_diversity(counts: np.ndarray) -> float:
    """Return unbiased gene diversity from allele counts."""
    total = int(np.sum(counts))
    if total <= 1:
        return float("nan")
    freqs = counts / total
    base = 1.0 - float(np.sum(freqs * freqs))
    return float((total / (total - 1)) * base)


def _inbreeding_from_heterozygosity(
    observed_heterozygosity: float,
    expected_heterozygosity: float,
) -> float:
    """Return one inbreeding coefficient as `1 - Ho/He`."""
    if np.isnan(observed_heterozygosity) or np.isnan(expected_heterozygosity):
        return float("nan")
    if expected_heterozygosity <= 0:
        return float("nan")
    return float(1.0 - (observed_heterozygosity / expected_heterozygosity))


def fis_from_heterozygosity(
    observed_heterozygosity: float,
    expected_heterozygosity: float,
) -> float:
    """Return within-population Fis as `1 - Ho/He`."""
    return _inbreeding_from_heterozygosity(
        observed_heterozygosity,
        expected_heterozygosity,
    )


def fit_from_heterozygosity(
    observed_heterozygosity: float,
    expected_heterozygosity_total: float,
) -> float:
    """Return global Fit as `1 - Ho/Ht`."""
    return _inbreeding_from_heterozygosity(
        observed_heterozygosity,
        expected_heterozygosity_total,
    )


def sitewise_dxy(counts1: np.ndarray, counts2: np.ndarray) -> float:
    """Return sitewise dxy between two populations."""
    total1 = int(np.sum(counts1))
    total2 = int(np.sum(counts2))
    if total1 == 0 or total2 == 0:
        return float("nan")
    freqs1 = counts1 / total1
    freqs2 = counts2 / total2
    return float(1.0 - np.sum(freqs1 * freqs2))


def sitewise_dxy_block(counts1: np.ndarray, counts2: np.ndarray) -> np.ndarray:
    """Return vectorized sitewise dxy for two `(nsites, nalleles)` count matrices."""
    if counts1.shape != counts2.shape:
        raise ValueError("sitewise_dxy_block requires equal-shaped count matrices")
    total1 = np.sum(counts1, axis=1)
    total2 = np.sum(counts2, axis=1)
    result = np.full(total1.shape, np.nan, dtype=np.float64)
    valid = (total1 > 0) & (total2 > 0)
    if not np.any(valid):
        return result
    freqs1 = counts1[valid].astype(np.float64, copy=False) / total1[valid, None]
    freqs2 = counts2[valid].astype(np.float64, copy=False) / total2[valid, None]
    result[valid] = 1.0 - np.sum(freqs1 * freqs2, axis=1)
    return result


def hudson_fst_components(counts1: np.ndarray, counts2: np.ndarray) -> tuple[float, float]:
    """Return numerator and denominator components for Hudson-style Fst."""
    dxy = sitewise_dxy(counts1, counts2)
    if math.isnan(dxy):
        return float("nan"), float("nan")
    pi1 = gene_diversity(counts1)
    pi2 = gene_diversity(counts2)
    if math.isnan(pi1):
        pi1 = 0.0
    if math.isnan(pi2):
        pi2 = 0.0
    return float(dxy - ((pi1 + pi2) / 2.0)), float(dxy)


def hudson_fst_components_block(
    counts1: np.ndarray,
    counts2: np.ndarray,
    pi1: np.ndarray,
    pi2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return vectorized Hudson Fst numerator and denominator arrays."""
    dxy = sitewise_dxy_block(counts1, counts2)
    numerator = dxy - ((np.nan_to_num(pi1, nan=0.0) + np.nan_to_num(pi2, nan=0.0)) / 2.0)
    numerator[np.isnan(dxy)] = np.nan
    return numerator, dxy


def wattersons_theta(segregating_sites: int, chromosome_count: int, sites_used: int) -> float:
    """Return per-site Watterson's theta for fixed chromosome count."""
    if segregating_sites <= 0 or chromosome_count <= 1 or sites_used <= 0:
        return float("nan")
    a1 = np.sum(1.0 / np.arange(1, chromosome_count))
    return float((segregating_sites / a1) / sites_used)


def tajimas_d(pi_total: float, segregating_sites: int, chromosome_count: int) -> float:
    """Return Tajima's D from total pairwise differences and segregating sites."""
    if chromosome_count <= 1 or segregating_sites <= 0:
        return float("nan")
    n = float(chromosome_count)
    a1 = np.sum(1.0 / np.arange(1, chromosome_count))
    a2 = np.sum(1.0 / (np.arange(1, chromosome_count) ** 2))
    b1 = (n + 1.0) / (3.0 * (n - 1.0))
    b2 = (2.0 * (n * n + n + 3.0)) / (9.0 * n * (n - 1.0))
    c1 = b1 - (1.0 / a1)
    c2 = b2 - ((n + 2.0) / (a1 * n)) + (a2 / (a1 * a1))
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    theta_total = segregating_sites / a1
    var = (e1 * segregating_sites) + (e2 * segregating_sites * (segregating_sites - 1))
    if var <= 0:
        return float("nan")
    return float((pi_total - theta_total) / math.sqrt(var))
