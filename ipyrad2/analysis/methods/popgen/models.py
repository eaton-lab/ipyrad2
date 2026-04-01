#!/usr/bin/env python

"""Typed internal models for popgen requests, stats, and results."""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class StatSpec:
    """One popgen statistic definition."""

    name: str
    backends: frozenset[str]
    output_families: tuple[str, ...]
    manifest_formula: str | None = None
    supports_windowing: bool = True


STAT_SPECS: tuple[StatSpec, ...] = (
    StatSpec("pi", frozenset({"sequence"}), ("population",)),
    StatSpec("dxy", frozenset({"sequence"}), ("pairwise",)),
    StatSpec("fst", frozenset({"sequence", "snp"}), ("pairwise",)),
    StatSpec("tajima_d", frozenset({"sequence"}), ("population",)),
    StatSpec("theta_w", frozenset({"sequence"}), ("population",)),
    StatSpec("heterozygosity", frozenset({"sequence", "snp"}), ("population",)),
    StatSpec(
        "fis",
        frozenset({"sequence", "snp"}),
        ("population",),
        manifest_formula="fis = 1 - Ho/He",
    ),
    StatSpec(
        "fit",
        frozenset({"sequence", "snp"}),
        ("global",),
        manifest_formula="fit = 1 - Ho/Ht_total",
    ),
    StatSpec(
        "sfs",
        frozenset({"sequence", "snp"}),
        ("sfs",),
        supports_windowing=False,
    ),
)
STAT_SPEC_BY_NAME = {spec.name: spec for spec in STAT_SPECS}
ORDERED_STATS = tuple(spec.name for spec in STAT_SPECS)
SEQUENCE_STATS = {spec.name for spec in STAT_SPECS if "sequence" in spec.backends}
SNP_STATS = {spec.name for spec in STAT_SPECS if "snp" in spec.backends}
SEQUENCE_ONLY_STATS = SEQUENCE_STATS.difference(SNP_STATS)


def normalize_stats(stats) -> list[str]:
    """Return normalized statistic names in canonical order."""
    if stats is None:
        return []
    if isinstance(stats, str):
        tokens = [token.strip().lower() for token in stats.split(",")]
    else:
        tokens = [str(token).strip().lower() for token in stats]
    tokens = [token for token in tokens if token]
    if not tokens:
        return []
    if tokens == ["all"]:
        return ["all"]
    unknown = sorted(set(tokens).difference(ORDERED_STATS))
    if unknown:
        raise ValueError(f"unknown popgen statistics: {', '.join(unknown)}")
    return [name for name in ORDERED_STATS if name in tokens]


@dataclass(frozen=True)
class WindowingConfig:
    """Normalized windowing options for one popgen run."""

    mode: str | None = None
    window_size: int | None = None
    step_size: int | None = None
    loci_per_window: int | None = None
    locus_step: int | None = None

    def as_manifest_dict(self) -> dict[str, Any]:
        """Return ordered manifest-ready windowing fields."""
        return {
            "window_mode": self.mode or "none",
            "window_size": self.window_size,
            "step_size": self.step_size,
            "loci_per_window": self.loci_per_window,
            "locus_step": self.locus_step,
        }


@dataclass(frozen=True)
class PopgenRequest:
    """Resolved execution request for one popgen run."""

    data: Path
    name: str
    outdir: Path
    requested_stats: tuple[str, ...]
    backend: str
    has_sequence: bool
    has_snp: bool
    min_sample_coverage: float
    max_sample_missing: float
    min_minor_allele_frequency: float
    imap: Any
    minmap: Any
    exclude: tuple[str, ...]
    include_reference: bool
    subsample_unlinked: bool
    random_seed: int | None
    cores: int
    force: bool
    log_level: str
    windowing: WindowingConfig = field(default_factory=WindowingConfig)

    @property
    def requested_stat_list(self) -> list[str]:
        """Return requested stats as a list for manifest serialization."""
        return list(self.requested_stats)

    def requested_stat_formulas(self) -> dict[str, str]:
        """Return manifest formula notes for requested stats."""
        formulas = {
            f"{name}_formula": spec.manifest_formula
            for name, spec in STAT_SPEC_BY_NAME.items()
            if name in self.requested_stats and spec.manifest_formula is not None
        }
        return formulas


@dataclass(frozen=True)
class OutputTableSpec:
    """One tabular popgen output mapping."""

    key: str
    suffix: str


OUTPUT_TABLE_SPECS: tuple[OutputTableSpec, ...] = (
    OutputTableSpec("sample_stats", "sample_stats.tsv"),
    OutputTableSpec("global_stats", "global_stats.tsv"),
    OutputTableSpec("population_stats", "population_stats.tsv"),
    OutputTableSpec("pairwise_stats", "pairwise_stats.tsv"),
    OutputTableSpec("sfs", "sfs.tsv"),
    OutputTableSpec("window_population_stats", "window_population_stats.tsv"),
    OutputTableSpec("window_pairwise_stats", "window_pairwise_stats.tsv"),
)


@dataclass
class PopgenResult:
    """Typed result object returned by one popgen backend."""

    sample_data_summary: pd.DataFrame
    sample_stats: pd.DataFrame
    summary: dict[str, Any]
    global_stats: pd.DataFrame | None = None
    population_stats: pd.DataFrame | None = None
    pairwise_stats: pd.DataFrame | None = None
    sfs: pd.DataFrame | None = None
    window_population_stats: pd.DataFrame | None = None
    window_pairwise_stats: pd.DataFrame | None = None

    def get_output_table(self, key: str) -> pd.DataFrame | None:
        """Return one named tabular output."""
        return getattr(self, key)

    def iter_output_tables(self) -> list[tuple[str, pd.DataFrame]]:
        """Return ordered non-empty tabular outputs written as files."""
        tables: list[tuple[str, pd.DataFrame]] = [("sample_stats", self.sample_stats)]
        for spec in OUTPUT_TABLE_SPECS[1:]:
            table = self.get_output_table(spec.key)
            if table is not None and not table.empty:
                tables.append((spec.key, table))
        return tables
