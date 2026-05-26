#!/usr/bin/env python

"""Typed internal models for baba requests and outputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class QuartetDefinition:
    """One resolved quartet test."""

    source: str
    p1: str
    p2: str
    p3: str
    p4: str

    def as_dict(self) -> dict[str, str]:
        """Return one manifest/table-ready row."""
        return {
            "source": self.source,
            "p1": self.p1,
            "p2": self.p2,
            "p3": self.p3,
            "p4": self.p4,
        }


@dataclass(frozen=True)
class BabaRequest:
    """Resolved execution request for one baba run."""

    data: Path
    name: str
    outdir: Path
    tests: Path | None
    tree: Path | None
    imap: Path | None
    minmap: Path | None
    min_sample_coverage: int
    min_genotype_depth: int
    min_site_qual: float
    exclude: tuple[str, ...]
    include_reference: bool
    resampling: str
    bootstrap_replicates: int
    jackknife_block_bp: int
    jackknife_block_loci: int
    seed: int | None
    f_branch: bool
    f_branch_p_threshold: float
    write_block_table: bool
    clustering_stats: bool
    cores: int
    force: bool
    log_level: str
    logged_command: str | None = None


@dataclass
class BabaResult:
    """Typed result object returned by one baba run."""

    quartets: pd.DataFrame
    rooted: pd.DataFrame
    resolved_tests: pd.DataFrame
    manifest: str
    summary_json: str
    f_branch: pd.DataFrame | None = None
    f_branch_matrix: pd.DataFrame | None = None
    f_branch_z: pd.DataFrame | None = None
    f_branch_p: pd.DataFrame | None = None
    blocks: pd.DataFrame | None = None
    tree_text: str | None = None
