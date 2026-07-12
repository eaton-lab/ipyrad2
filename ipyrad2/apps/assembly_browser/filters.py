"""Interactive SNP/locus filtering for the assembly browser."""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil
from typing import Iterable

import numpy as np
import pandas as pd

from .data import AssemblyStore
from .data import SNPSMAP_COLUMNS


MISSING_GENO = 255
MISSING_SNP = 78


@dataclass(frozen=True)
class FilterParams:
    """User-adjustable filters for browser summaries."""

    samples: tuple[str, ...] = ()
    max_sample_missing: float = 1.0
    max_site_missing: float = 1.0
    min_sample_coverage: float = 0.0
    min_minor_allele_frequency: float = 0.0
    min_genotype_depth: int = 0
    min_site_qual: float = 0.0
    chunk_size: int = 50_000


@dataclass(frozen=True)
class FilterResult:
    """Filtered summary payload for charts and tables."""

    params: FilterParams
    input_samples: list[str]
    retained_samples: list[str]
    dropped_samples: list[str]
    sample_summary: pd.DataFrame
    site_summary: pd.DataFrame
    locus_summary: pd.DataFrame
    filter_counts: pd.DataFrame
    totals: dict[str, int | float]


def _normalize_fraction(value: float) -> float:
    return min(1.0, max(0.0, float(value)))


def _coverage_threshold(min_sample_coverage: float, nsamples: int) -> int:
    """Convert an absolute or fractional sample coverage threshold to sample count."""
    if min_sample_coverage <= 1:
        return int(ceil(_normalize_fraction(min_sample_coverage) * nsamples))
    return int(min_sample_coverage)


def _selected_samples(store: AssemblyStore, requested: Iterable[str]) -> tuple[list[str], list[int]]:
    names = store.genotype_sample_names()
    if requested:
        selected = [name for name in requested if name in names]
    else:
        selected = list(names)
    if not selected:
        raise ValueError("No selected samples are present in the HDF5 genotype matrix.")
    index = store.sample_name_to_index()
    return selected, [index[name] for name in selected]


def _chunk_site_metrics(
    chunk: dict[str, np.ndarray],
    *,
    min_genotype_depth: int,
) -> dict[str, np.ndarray | int]:
    """Compute per-site metrics for one HDF5 chunk."""
    genos = chunk["genos"]
    sample_dp = chunk["sample_dp"]
    site_qual = chunk["site_qual"]
    snpsmap = chunk["snpsmap"]

    genos = genos.copy()
    depth_masked = 0
    if min_genotype_depth > 0 and sample_dp is not None:
        called = np.any(genos != MISSING_GENO, axis=2)
        low_depth = called & (sample_dp < min_genotype_depth)
        if np.any(low_depth):
            genos[low_depth] = MISSING_GENO
            depth_masked = int(np.count_nonzero(low_depth))

    called_genotype = np.any(genos != MISSING_GENO, axis=2)
    called_samples = called_genotype.sum(axis=0).astype(np.int32)
    missing_fraction = 1.0 - (called_samples / max(1, genos.shape[0]))
    missing_genotypes_by_sample = (~called_genotype).sum(axis=1).astype(np.int64)

    bad_allele = np.any((genos == 2) | (genos == 3), axis=(0, 2))
    genomask = np.ma.array(genos, mask=(genos == MISSING_GENO))
    invariant = np.all(genomask.sum(axis=2).mean(axis=0).round().astype(int).data == genomask.sum(axis=2).data, axis=0)

    ref_alleles = (genomask == 0).sum(axis=2).sum(axis=0).data.astype(float)
    alt_alleles = (genomask == 1).sum(axis=2).sum(axis=0).data.astype(float)
    denom = ref_alleles + alt_alleles
    with np.errstate(divide="ignore", invalid="ignore"):
        alt_freq = np.divide(alt_alleles, denom, out=np.zeros_like(alt_alleles), where=denom > 0)
    maf = np.minimum(alt_freq, 1.0 - alt_freq)

    return {
        "loc": snpsmap[:, 0].astype(np.int64),
        "loc_idx": snpsmap[:, 1].astype(np.int64),
        "loc_pos": snpsmap[:, 2].astype(np.int64),
        "scaff": snpsmap[:, 3].astype(np.int64),
        "pos": snpsmap[:, 4].astype(np.int64),
        "called_samples": called_samples,
        "missing_fraction": missing_fraction.astype(np.float32),
        "maf": maf.astype(np.float32),
        "site_qual": (
            site_qual.astype(np.float32)
            if site_qual is not None
            else np.full(suffix_len(snpsmap), np.nan, dtype=np.float32)
        ),
        "bad_allele": bad_allele.astype(bool),
        "invariant": invariant.astype(bool),
        "missing_by_sample": missing_genotypes_by_sample,
        "total_genotypes_by_sample": np.full(genos.shape[0], genos.shape[1], dtype=np.int64),
        "depth_masked": depth_masked,
    }


def suffix_len(array: np.ndarray) -> int:
    """Return first-axis length with a name that reads clearly at call sites."""
    return int(array.shape[0])


def _collect_site_metrics(
    store: AssemblyStore,
    sample_indices: list[int],
    params: FilterParams,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, int]:
    chunks = []
    missing_by_sample = np.zeros(len(sample_indices), dtype=np.int64)
    total_by_sample = np.zeros(len(sample_indices), dtype=np.int64)
    depth_masked = 0

    for chunk in store.site_chunks(sample_indices, chunk_size=params.chunk_size):
        metrics = _chunk_site_metrics(
            chunk,
            min_genotype_depth=params.min_genotype_depth,
        )
        depth_masked += int(metrics["depth_masked"])
        missing_by_sample += metrics["missing_by_sample"]
        total_by_sample += metrics["total_genotypes_by_sample"]
        chunks.append(
            pd.DataFrame(
                {
                    "loc": metrics["loc"],
                    "loc_idx": metrics["loc_idx"],
                    "loc_pos": metrics["loc_pos"],
                    "scaff": metrics["scaff"],
                    "pos": metrics["pos"],
                    "called_samples": metrics["called_samples"],
                    "missing_fraction": metrics["missing_fraction"],
                    "maf": metrics["maf"],
                    "site_qual": metrics["site_qual"],
                    "bad_allele": metrics["bad_allele"],
                    "invariant": metrics["invariant"],
                }
            )
        )

    if not chunks:
        columns = SNPSMAP_COLUMNS + [
            "called_samples",
            "missing_fraction",
            "maf",
            "site_qual",
            "bad_allele",
            "invariant",
        ]
        return pd.DataFrame(columns=columns), missing_by_sample, total_by_sample, depth_masked
    return pd.concat(chunks, ignore_index=True), missing_by_sample, total_by_sample, depth_masked


def _sample_summary(names: list[str], missing: np.ndarray, totals: np.ndarray) -> pd.DataFrame:
    with np.errstate(divide="ignore", invalid="ignore"):
        missing_fraction = np.divide(
            missing,
            totals,
            out=np.zeros(missing.shape, dtype=float),
            where=totals > 0,
        )
    return pd.DataFrame(
        {
            "sample": names,
            "missing_genotypes": missing.astype(int),
            "total_genotypes": totals.astype(int),
            "missing_fraction": missing_fraction,
        }
    )


def _apply_site_filters(site_df: pd.DataFrame, params: FilterParams, nsamples: int) -> pd.DataFrame:
    min_called = _coverage_threshold(params.min_sample_coverage, nsamples)
    pass_mask = (
        (~site_df["bad_allele"])
        & (~site_df["invariant"])
        & (site_df["called_samples"] >= min_called)
        & (site_df["missing_fraction"] <= _normalize_fraction(params.max_site_missing))
        & (site_df["maf"] >= max(0.0, float(params.min_minor_allele_frequency)))
    )
    if "site_qual" in site_df and not site_df["site_qual"].isna().all():
        pass_mask &= site_df["site_qual"].fillna(-np.inf) >= float(params.min_site_qual)

    out = site_df.copy()
    out["pass_filter"] = pass_mask.to_numpy(dtype=bool)
    return out


def _locus_summary(site_df: pd.DataFrame) -> pd.DataFrame:
    if site_df.empty:
        return pd.DataFrame(
            columns=[
                "loc",
                "scaff",
                "snps",
                "retained_snps",
                "mean_missing_fraction",
                "mean_maf",
                "pass_filter",
            ]
        )
    grouped = site_df.groupby("loc", sort=True)
    summary = grouped.agg(
        scaff=("scaff", "first"),
        snps=("loc", "size"),
        retained_snps=("pass_filter", "sum"),
        mean_missing_fraction=("missing_fraction", "mean"),
        mean_maf=("maf", "mean"),
    ).reset_index()
    summary["pass_filter"] = summary["retained_snps"] > 0
    return summary


def _filter_counts(site_df: pd.DataFrame, params: FilterParams, nsamples: int) -> pd.DataFrame:
    min_called = _coverage_threshold(params.min_sample_coverage, nsamples)
    reasons = {
        "non_biallelic_or_complex": site_df["bad_allele"],
        "invariant": site_df["invariant"],
        "below_min_sample_coverage": site_df["called_samples"] < min_called,
        "above_max_site_missing": site_df["missing_fraction"] > _normalize_fraction(params.max_site_missing),
        "below_maf": site_df["maf"] < max(0.0, float(params.min_minor_allele_frequency)),
    }
    if "site_qual" in site_df and not site_df["site_qual"].isna().all():
        reasons["below_min_site_qual"] = site_df["site_qual"].fillna(-np.inf) < float(params.min_site_qual)
    rows = [{"filter": name, "sites": int(np.count_nonzero(mask))} for name, mask in reasons.items()]
    rows.append({"filter": "retained", "sites": int(np.count_nonzero(site_df["pass_filter"]))})
    return pd.DataFrame(rows)


def apply_filters(store: AssemblyStore, params: FilterParams) -> FilterResult:
    """Apply browser filters and return compact summaries."""
    input_names, input_indices = _selected_samples(store, params.samples)
    site_df, missing, totals, depth_masked = _collect_site_metrics(store, input_indices, params)
    sample_df = _sample_summary(input_names, missing, totals)

    max_sample_missing = _normalize_fraction(params.max_sample_missing)
    keep_sample = sample_df["missing_fraction"] <= max_sample_missing
    if not keep_sample.all() and keep_sample.any():
        retained_names = sample_df.loc[keep_sample, "sample"].tolist()
        retained_names, retained_indices = _selected_samples(store, retained_names)
        site_df, missing, totals, depth_masked = _collect_site_metrics(store, retained_indices, params)
        sample_df = _sample_summary(retained_names, missing, totals)
    else:
        retained_names = input_names

    dropped_samples = [name for name in input_names if name not in set(retained_names)]
    site_df = _apply_site_filters(site_df, params, len(retained_names))
    locus_df = _locus_summary(site_df)
    counts_df = _filter_counts(site_df, params, len(retained_names))

    totals_dict = {
        "input_samples": len(input_names),
        "retained_samples": len(retained_names),
        "dropped_samples": len(dropped_samples),
        "input_snps": int(site_df.shape[0]),
        "retained_snps": int(site_df["pass_filter"].sum()) if not site_df.empty else 0,
        "input_loci": int(site_df["loc"].nunique()) if not site_df.empty else 0,
        "retained_loci": int(locus_df["pass_filter"].sum()) if not locus_df.empty else 0,
        "depth_masked_genotypes": int(depth_masked),
    }

    return FilterResult(
        params=params,
        input_samples=input_names,
        retained_samples=retained_names,
        dropped_samples=dropped_samples,
        sample_summary=sample_df,
        site_summary=site_df,
        locus_summary=locus_df,
        filter_counts=counts_df,
        totals=totals_dict,
    )
