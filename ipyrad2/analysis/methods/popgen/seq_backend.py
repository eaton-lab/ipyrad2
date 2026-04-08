#!/usr/bin/env python

"""Sequence-HDF5 backend for genome-wide population-genetic statistics."""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
from loguru import logger

from ....utils.exceptions import IPyradError
from ....utils.parallel import run_with_pool
from ...extracters.sequence_common import build_sequence_imap_minmap
from ...extracters.sequence_common import normalize_sequence_population_inputs
from ...extracters.sequence_common import resolve_sequence_sample_subset
from ...extracters.sequence_common import sync_sequence_imap_after_sample_drop
from ..common import build_sample_data_summary
from .common import build_pairwise_stats_dataframe
from .common import build_population_stats_dataframe
from .common import build_sfs_dataframe
from .common import build_global_stats_dataframe
from .common import build_sample_stats_dataframe
from .estimators import fis_from_heterozygosity
from .estimators import hudson_fst_components_block
from .models import PopgenResult
from .estimators import summarize_sequence_block
from .estimators import tajimas_d
from .estimators import wattersons_theta


SEQUENCE_CHUNK_SITES = 5000
_HET_LOOKUP = np.zeros(256, dtype=bool)
_MISSING_LOOKUP = np.zeros(256, dtype=bool)
for _base in "RYSWKM":
    _HET_LOOKUP[ord(_base)] = True
for _base in "N-":
    _MISSING_LOOKUP[ord(_base)] = True


@dataclass(frozen=True)
class WindowSpec:
    """One planned popgen window over sequence-HDF5 spans."""

    window_id: int
    window_mode: str
    scaffold: str
    start: int | None
    end: int | None
    first_locus: int | None
    last_locus: int | None
    nloci: int
    sites_total: int
    spans: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class WindowChunkView:
    """One window projected into local chunk coordinates."""

    window_id: int
    local_spans: tuple[tuple[int, int], ...]


@dataclass
class SequenceAccumulator:
    """Typed additive accumulator for one set of sequence summaries."""

    missing_counts: np.ndarray
    heterozygous_counts: np.ndarray
    total_sites: int
    sites_used_pi: np.ndarray
    pi_sum: np.ndarray
    sites_used_theta: np.ndarray
    segregating_sites: np.ndarray
    theta_pi_total: np.ndarray
    sites_used_heterozygosity: np.ndarray
    observed_heterozygosity_sum: np.ndarray
    expected_heterozygosity_sum: np.ndarray
    global_sites_used_heterozygosity: int
    global_observed_heterozygosity_sum: float
    global_expected_heterozygosity_sum: float
    pair_sites_used: np.ndarray
    pair_dxy_sum: np.ndarray
    pair_fst_num_sum: np.ndarray
    pair_fst_den_sum: np.ndarray
    sfs_counts: np.ndarray

    @classmethod
    def zeros(
        cls,
        *,
        npops: int,
        npairs: int,
        max_minor_allele_count: int,
        nsamples: int,
    ) -> "SequenceAccumulator":
        """Return a zero-valued accumulator with stable array shapes."""
        return cls(
            missing_counts=np.zeros(nsamples, dtype=np.int64),
            heterozygous_counts=np.zeros(nsamples, dtype=np.int64),
            total_sites=0,
            sites_used_pi=np.zeros(npops, dtype=np.int64),
            pi_sum=np.zeros(npops, dtype=np.float64),
            sites_used_theta=np.zeros(npops, dtype=np.int64),
            segregating_sites=np.zeros(npops, dtype=np.int64),
            theta_pi_total=np.zeros(npops, dtype=np.float64),
            sites_used_heterozygosity=np.zeros(npops, dtype=np.int64),
            observed_heterozygosity_sum=np.zeros(npops, dtype=np.float64),
            expected_heterozygosity_sum=np.zeros(npops, dtype=np.float64),
            global_sites_used_heterozygosity=0,
            global_observed_heterozygosity_sum=0.0,
            global_expected_heterozygosity_sum=0.0,
            pair_sites_used=np.zeros(npairs, dtype=np.int64),
            pair_dxy_sum=np.zeros(npairs, dtype=np.float64),
            pair_fst_num_sum=np.zeros(npairs, dtype=np.float64),
            pair_fst_den_sum=np.zeros(npairs, dtype=np.float64),
            sfs_counts=np.zeros((npops, max_minor_allele_count + 1), dtype=np.int64),
        )

    def merge(self, other: "SequenceAccumulator") -> None:
        """Add another accumulator into this one."""
        self.missing_counts += other.missing_counts
        self.heterozygous_counts += other.heterozygous_counts
        self.total_sites += int(other.total_sites)
        self.sites_used_pi += other.sites_used_pi
        self.pi_sum += other.pi_sum
        self.sites_used_theta += other.sites_used_theta
        self.segregating_sites += other.segregating_sites
        self.theta_pi_total += other.theta_pi_total
        self.sites_used_heterozygosity += other.sites_used_heterozygosity
        self.observed_heterozygosity_sum += other.observed_heterozygosity_sum
        self.expected_heterozygosity_sum += other.expected_heterozygosity_sum
        self.global_sites_used_heterozygosity += int(other.global_sites_used_heterozygosity)
        self.global_observed_heterozygosity_sum += float(other.global_observed_heterozygosity_sum)
        self.global_expected_heterozygosity_sum += float(other.global_expected_heterozygosity_sum)
        self.pair_sites_used += other.pair_sites_used
        self.pair_dxy_sum += other.pair_dxy_sum
        self.pair_fst_num_sum += other.pair_fst_num_sum
        self.pair_fst_den_sum += other.pair_fst_den_sum
        self.sfs_counts += other.sfs_counts


@dataclass
class SequenceChunkSummary:
    """One chunk's genome-wide and window-local summary products."""

    genome: SequenceAccumulator
    windows: dict[int, SequenceAccumulator] = field(default_factory=dict)


def _plan_sequence_chunk_spans(
    data: Path,
    *,
    target_sites: int = SEQUENCE_CHUNK_SITES,
) -> list[tuple[tuple[int, int], ...]]:
    """Group phymap spans into moderate-size chunks for I/O and parallel work."""
    chunks: list[tuple[tuple[int, int], ...]] = []
    current: list[tuple[int, int]] = []
    current_sites = 0

    with h5py.File(data, "r") as io5:
        for row in io5["phymap"]:
            start = int(row[1])
            end = int(row[2])
            width = end - start
            if current and current_sites + width > target_sites:
                chunks.append(tuple(current))
                current = []
                current_sites = 0
            if current and current[-1][1] == start:
                current[-1] = (current[-1][0], end)
            else:
                current.append((start, end))
            current_sites += width

    if current:
        chunks.append(tuple(current))
    return chunks


def _load_sequence_chunk_from_phy(
    phy: h5py.Dataset,
    sidxs: list[int],
    spans: tuple[tuple[int, int], ...],
) -> np.ndarray:
    """Load one `(samples, sites)` sequence chunk from one or more contiguous spans."""
    total_sites = sum(end - start for start, end in spans)
    block = np.empty((len(sidxs), total_sites), dtype=np.uint8)
    offset = 0
    for start, end in spans:
        width = end - start
        block[:, offset : offset + width] = phy[sidxs, start:end]
        offset += width
    block[block == 45] = 78
    return block


def _count_sequence_missing(block: np.ndarray) -> np.ndarray:
    """Return per-sample missing-call counts for one loaded sequence block."""
    return np.sum(_MISSING_LOOKUP[block], axis=1).astype(np.int64, copy=False)


def _slice_block_by_local_spans(
    block: np.ndarray,
    spans: tuple[tuple[int, int], ...],
) -> np.ndarray:
    """Return one block subset from local `(start, end)` spans."""
    if not spans:
        return block[:, 0:0]
    if len(spans) == 1:
        start, end = spans[0]
        return block[:, start:end]
    return np.concatenate(
        [block[:, start:end] for start, end in spans],
        axis=1,
    )


def _summarize_sequence_missing_chunk(
    *,
    data: Path,
    sidxs: list[int],
    spans: tuple[tuple[int, int], ...],
) -> tuple[np.ndarray, int]:
    """Return missing-call totals for one chunk."""
    with h5py.File(data, "r") as io5:
        block = _load_sequence_chunk_from_phy(io5["phy"], sidxs, spans)
    return _count_sequence_missing(block), int(block.shape[1])


def _calculate_sample_missing(
    data: Path,
    sidxs: list[int],
    snames: list[str],
    chunk_spans: list[tuple[tuple[int, int], ...]],
    *,
    cores: int,
    log_level: str,
) -> tuple[pd.Series, int]:
    """Return per-sample missingness across the full selected sequence matrix."""
    missing = np.zeros(len(sidxs), dtype=np.int64)
    total_sites = 0

    if cores > 1 and len(chunk_spans) > 1:
        jobs = {
            idx: (
                _summarize_sequence_missing_chunk,
                {"data": data, "sidxs": sidxs, "spans": spans},
            )
            for idx, spans in enumerate(chunk_spans)
        }
        results = run_with_pool(
            jobs,
            log_level,
            cores,
            msg="Calculating popgen sample missingness",
        )
        for chunk_missing, chunk_sites in results.values():
            missing += chunk_missing
            total_sites += chunk_sites
    else:
        with h5py.File(data, "r") as io5:
            phy = io5["phy"]
            for spans in chunk_spans:
                block = _load_sequence_chunk_from_phy(phy, sidxs, spans)
                missing += _count_sequence_missing(block)
                total_sites += block.shape[1]

    values = np.zeros(len(sidxs), dtype=float) if total_sites == 0 else missing / total_sites
    return pd.Series(values, index=snames, dtype=float), total_sites


def _init_sequence_accumulators(
    npops: int,
    npairs: int,
    max_minor_allele_count: int,
    nsamples: int,
) -> SequenceAccumulator:
    """Create additive accumulators for one or more sequence chunks."""
    return SequenceAccumulator.zeros(
        npops=npops,
        npairs=npairs,
        max_minor_allele_count=max_minor_allele_count,
        nsamples=nsamples,
    )


def _summarize_sequence_block_into(
    accum: SequenceAccumulator,
    block: np.ndarray,
    *,
    pop_indices: list[np.ndarray],
    pop_sizes: np.ndarray,
    pop_minmap: np.ndarray,
    pair_indices: list[tuple[int, int]],
    track_sample_stats: bool,
    need_pi: bool,
    need_theta: bool,
    need_heterozygosity: bool,
    need_global_fit: bool,
    global_min_sample_coverage: int,
    need_pairwise: bool,
    need_sfs: bool,
) -> None:
    """Accumulate one block of sequence sites into additive numeric summaries."""
    accum.total_sites += int(block.shape[1])
    if track_sample_stats:
        accum.missing_counts += _count_sequence_missing(block)
        accum.heterozygous_counts += np.sum(_HET_LOOKUP[block], axis=1).astype(
            np.int64,
            copy=False,
        )

    if need_global_fit:
        global_summary = summarize_sequence_block(
            block,
            include_minor_allele_count=False,
        )
        global_valid = (
            (global_summary.called_samples >= global_min_sample_coverage)
            & (global_summary.chromosome_count >= 2)
        )
        accum.global_sites_used_heterozygosity += int(np.sum(global_valid))
        accum.global_observed_heterozygosity_sum += float(
            np.nansum(global_summary.observed_heterozygosity[global_valid])
        )
        accum.global_expected_heterozygosity_sum += float(
            np.nansum(global_summary.expected_heterozygosity[global_valid])
        )

    pop_summaries = []
    for pop_idx, indices in enumerate(pop_indices):
        summary = summarize_sequence_block(
            block[indices, :],
            include_minor_allele_count=need_sfs,
        )
        pop_summaries.append(summary)

        valid = (
            (summary.called_samples >= pop_minmap[pop_idx])
            & (summary.chromosome_count >= 2)
        )

        if need_pi:
            accum.sites_used_pi[pop_idx] += int(np.sum(valid))
            accum.pi_sum[pop_idx] += float(np.nansum(summary.pi[valid]))

        if need_heterozygosity:
            accum.sites_used_heterozygosity[pop_idx] += int(np.sum(valid))
            accum.observed_heterozygosity_sum[pop_idx] += float(
                np.nansum(summary.observed_heterozygosity[valid])
            )
            accum.expected_heterozygosity_sum[pop_idx] += float(
                np.nansum(summary.expected_heterozygosity[valid])
            )

        if need_sfs and summary.minor_allele_count is not None:
            sfs_valid = valid & summary.biallelic
            if np.any(sfs_valid):
                counts = np.bincount(
                    summary.minor_allele_count[sfs_valid],
                    minlength=accum.sfs_counts.shape[1],
                )
                accum.sfs_counts[pop_idx, : counts.size] += counts

        if need_theta:
            full = (
                (summary.called_samples == pop_sizes[pop_idx])
                & (summary.chromosome_count >= 2)
            )
            accum.sites_used_theta[pop_idx] += int(np.sum(full))
            accum.theta_pi_total[pop_idx] += float(np.nansum(summary.pi[full]))
            accum.segregating_sites[pop_idx] += int(np.sum(summary.segregating[full]))

    if not need_pairwise:
        return

    for pair_idx, (idx1, idx2) in enumerate(pair_indices):
        summary1 = pop_summaries[idx1]
        summary2 = pop_summaries[idx2]
        valid = (
            (summary1.called_samples >= pop_minmap[idx1])
            & (summary2.called_samples >= pop_minmap[idx2])
            & (summary1.chromosome_count > 0)
            & (summary2.chromosome_count > 0)
        )
        if not np.any(valid):
            continue

        accum.pair_sites_used[pair_idx] += int(np.sum(valid))
        fst_num, fst_den = hudson_fst_components_block(
            summary1.allele_counts[valid],
            summary2.allele_counts[valid],
            summary1.pi[valid],
            summary2.pi[valid],
        )
        accum.pair_dxy_sum[pair_idx] += float(np.nansum(fst_den))
        accum.pair_fst_num_sum[pair_idx] += float(np.nansum(fst_num))
        accum.pair_fst_den_sum[pair_idx] += float(np.nansum(fst_den))


def _summarize_sequence_chunk(
    *,
    data: Path,
    sidxs: list[int],
    spans: tuple[tuple[int, int], ...],
    pop_indices: list[np.ndarray],
    pop_sizes: np.ndarray,
    pop_minmap: np.ndarray,
    pair_indices: list[tuple[int, int]],
    max_minor_allele_count: int,
    track_sample_stats: bool,
    need_pi: bool,
    need_theta: bool,
    need_heterozygosity: bool,
    need_global_fit: bool,
    global_min_sample_coverage: int,
    need_pairwise: bool,
    need_sfs: bool,
    window_views: tuple[WindowChunkView, ...] = (),
) -> SequenceChunkSummary:
    """Load and summarize one sequence chunk for serial or pooled reduction."""
    with h5py.File(data, "r") as io5:
        block = _load_sequence_chunk_from_phy(io5["phy"], sidxs, spans)

    accum = _init_sequence_accumulators(
        len(pop_indices),
        len(pair_indices),
        max_minor_allele_count,
        len(sidxs),
    )
    _summarize_sequence_block_into(
        accum,
        block,
        pop_indices=pop_indices,
        pop_sizes=pop_sizes,
        pop_minmap=pop_minmap,
        pair_indices=pair_indices,
        track_sample_stats=track_sample_stats,
        need_pi=need_pi,
        need_theta=need_theta,
        need_heterozygosity=need_heterozygosity,
        need_global_fit=need_global_fit,
        global_min_sample_coverage=global_min_sample_coverage,
        need_pairwise=need_pairwise,
        need_sfs=need_sfs,
    )
    window_accums: dict[int, SequenceAccumulator] = {}
    if window_views:
        for view in window_views:
            subblock = _slice_block_by_local_spans(block, view.local_spans)
            window_accum = _init_sequence_accumulators(
                len(pop_indices),
                len(pair_indices),
                0,
                0,
            )
            _summarize_sequence_block_into(
                window_accum,
                subblock,
                pop_indices=pop_indices,
                pop_sizes=pop_sizes,
                pop_minmap=pop_minmap,
                pair_indices=pair_indices,
                track_sample_stats=False,
                need_pi=need_pi,
                need_theta=need_theta,
                need_heterozygosity=need_heterozygosity,
                need_global_fit=False,
                global_min_sample_coverage=1,
                need_pairwise=need_pairwise,
                need_sfs=False,
            )
            window_accums[view.window_id] = window_accum
    return SequenceChunkSummary(genome=accum, windows=window_accums)


def _merge_sequence_accumulators(
    target: SequenceAccumulator,
    source: SequenceAccumulator,
) -> None:
    """Merge one chunk accumulator into the running totals."""
    target.merge(source)


def _collect_sequence_accumulators(
    *,
    data: Path,
    sidxs: list[int],
    chunk_spans: list[tuple[tuple[int, int], ...]],
    pop_indices: list[np.ndarray],
    pop_sizes: np.ndarray,
    pop_minmap: np.ndarray,
    pair_indices: list[tuple[int, int]],
    track_sample_stats: bool,
    need_pi: bool,
    need_theta: bool,
    need_heterozygosity: bool,
    need_global_fit: bool,
    global_min_sample_coverage: int,
    need_pairwise: bool,
    need_sfs: bool,
    chunk_window_views: dict[int, tuple[WindowChunkView, ...]] | None,
    cores: int,
    log_level: str,
) -> tuple[SequenceAccumulator, dict[int, SequenceAccumulator]]:
    """Return merged numeric accumulators for the selected sequence samples."""
    max_minor_allele_count = int(2 * np.max(pop_sizes)) if pop_sizes.size else 0
    merged = _init_sequence_accumulators(
        len(pop_indices),
        len(pair_indices),
        max_minor_allele_count,
        len(sidxs),
    )
    window_totals: dict[int, SequenceAccumulator] = {}

    if cores > 1 and len(chunk_spans) > 1:
        jobs = {
            idx: (
                _summarize_sequence_chunk,
                {
                    "data": data,
                    "sidxs": sidxs,
                    "spans": spans,
                    "pop_indices": pop_indices,
                    "pop_sizes": pop_sizes,
                    "pop_minmap": pop_minmap,
                    "pair_indices": pair_indices,
                    "max_minor_allele_count": max_minor_allele_count,
                    "track_sample_stats": track_sample_stats,
                    "need_pi": need_pi,
                    "need_theta": need_theta,
                    "need_heterozygosity": need_heterozygosity,
                    "need_global_fit": need_global_fit,
                    "global_min_sample_coverage": global_min_sample_coverage,
                    "need_pairwise": need_pairwise,
                    "need_sfs": need_sfs,
                    "window_views": tuple() if chunk_window_views is None else chunk_window_views.get(idx, tuple()),
                },
            )
            for idx, spans in enumerate(chunk_spans)
        }
        results = run_with_pool(
            jobs,
            log_level,
            cores,
            msg="Computing popgen summaries",
        )
        for chunk_result in results.values():
            _merge_sequence_accumulators(merged, chunk_result.genome)
            for window_id, window_accum in chunk_result.windows.items():
                target = window_totals.get(window_id)
                if target is None:
                    target = _init_sequence_accumulators(
                        len(pop_indices),
                        len(pair_indices),
                        0,
                        0,
                    )
                    window_totals[window_id] = target
                target.merge(window_accum)
        return merged, window_totals

    with h5py.File(data, "r") as io5:
        phy = io5["phy"]
        for idx, spans in enumerate(chunk_spans):
            block = _load_sequence_chunk_from_phy(phy, sidxs, spans)
            _summarize_sequence_block_into(
                merged,
                block,
                pop_indices=pop_indices,
                pop_sizes=pop_sizes,
                pop_minmap=pop_minmap,
                pair_indices=pair_indices,
                track_sample_stats=track_sample_stats,
                need_pi=need_pi,
                need_theta=need_theta,
                need_heterozygosity=need_heterozygosity,
                need_global_fit=need_global_fit,
                global_min_sample_coverage=global_min_sample_coverage,
                need_pairwise=need_pairwise,
                need_sfs=need_sfs,
            )
            for view in tuple() if chunk_window_views is None else chunk_window_views.get(idx, tuple()):
                subblock = _slice_block_by_local_spans(block, view.local_spans)
                target = window_totals.get(view.window_id)
                if target is None:
                    target = _init_sequence_accumulators(
                        len(pop_indices),
                        len(pair_indices),
                        0,
                        0,
                    )
                    window_totals[view.window_id] = target
                _summarize_sequence_block_into(
                    target,
                    subblock,
                    pop_indices=pop_indices,
                    pop_sizes=pop_sizes,
                    pop_minmap=pop_minmap,
                    pair_indices=pair_indices,
                    track_sample_stats=False,
                    need_pi=need_pi,
                    need_theta=need_theta,
                    need_heterozygosity=need_heterozygosity,
                    need_global_fit=False,
                    global_min_sample_coverage=1,
                    need_pairwise=need_pairwise,
                    need_sfs=False,
                )
    return merged, window_totals


def _decode_hdf5_strings(values: np.ndarray) -> list[str]:
    """Return decoded strings from one HDF5 string array."""
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def _merge_spans(spans: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    """Merge adjacent spans for efficient chunk loading."""
    if not spans:
        return tuple()
    merged = [spans[0]]
    for start, end in spans[1:]:
        prev_start, prev_end = merged[-1]
        if prev_end == start:
            merged[-1] = (prev_start, end)
        else:
            merged.append((start, end))
    return tuple(merged)


def _load_sequence_window_context(
    data: Path,
) -> tuple[np.ndarray, list[str], list[int]]:
    """Load phymap rows and scaffold metadata for window planning."""
    with h5py.File(data, "r") as io5:
        phymap = np.asarray(io5["phymap"])
        scaffold_names = _decode_hdf5_strings(io5.attrs["scaffold_names"])
        scaffold_lengths = [int(value) for value in io5.attrs["scaffold_lengths"]]
    return phymap, scaffold_names, scaffold_lengths


def _plan_genomic_windows(
    phymap: np.ndarray,
    scaffold_names: list[str],
    scaffold_lengths: list[int],
    *,
    window_size: int,
    step_size: int,
) -> list[WindowSpec]:
    """Plan scaffold-coordinate windows over phymap-delimited sequence loci."""
    specs: list[WindowSpec] = []
    row_numbers = np.arange(phymap.shape[0], dtype=np.int64)
    window_id = 1

    for scaff_idx, (scaffold, scaffold_length) in enumerate(
        zip(scaffold_names, scaffold_lengths)
    ):
        scaffold_mask = phymap[:, 0] == scaff_idx
        scaffold_rows = phymap[scaffold_mask]
        scaffold_row_numbers = row_numbers[scaffold_mask]
        locus_start0 = (
            scaffold_rows[:, 3].astype(np.int64, copy=False) - 1
            if scaffold_rows.size
            else np.zeros(0, dtype=np.int64)
        )
        locus_end0 = (
            scaffold_rows[:, 4].astype(np.int64, copy=False)
            if scaffold_rows.size
            else np.zeros(0, dtype=np.int64)
        )

        for start0 in range(0, scaffold_length, step_size):
            end0 = min(start0 + window_size, scaffold_length)
            spans: list[tuple[int, int]] = []
            first_locus = None
            last_locus = None
            nloci = 0

            if scaffold_rows.size:
                overlap_mask = (locus_end0 > start0) & (locus_start0 < end0)
                selected_rows = scaffold_rows[overlap_mask]
                selected_row_numbers = scaffold_row_numbers[overlap_mask]
                nloci = int(selected_rows.shape[0])
                if nloci:
                    first_locus = int(selected_row_numbers[0]) + 1
                    last_locus = int(selected_row_numbers[-1]) + 1
                for row in selected_rows:
                    locus_start = int(row[3]) - 1
                    locus_end = int(row[4])
                    overlap_start = max(locus_start, start0)
                    overlap_end = min(locus_end, end0)
                    if overlap_start >= overlap_end:
                        continue
                    phy_start = int(row[1]) + (overlap_start - locus_start)
                    phy_end = phy_start + (overlap_end - overlap_start)
                    spans.append((phy_start, phy_end))

            merged_spans = _merge_spans(spans)
            sites_total = sum(stop - start for start, stop in merged_spans)
            specs.append(
                WindowSpec(
                    window_id=window_id,
                    window_mode="genomic",
                    scaffold=scaffold,
                    start=start0 + 1,
                    end=end0,
                    first_locus=first_locus,
                    last_locus=last_locus,
                    nloci=nloci,
                    sites_total=sites_total,
                    spans=merged_spans,
                )
            )
            window_id += 1

    return specs


def _plan_locus_windows(
    phymap: np.ndarray,
    scaffold_names: list[str],
    *,
    loci_per_window: int,
    locus_step: int,
) -> list[WindowSpec]:
    """Plan anonymous RAD windows over consecutive phymap loci."""
    specs: list[WindowSpec] = []
    window_id = 1

    for start_idx in range(0, phymap.shape[0], locus_step):
        rows = phymap[start_idx : start_idx + loci_per_window]
        if rows.shape[0] == 0:
            break
        merged_spans = _merge_spans(
            [(int(row[1]), int(row[2])) for row in rows]
        )
        scaffold_ids = rows[:, 0].astype(np.int64, copy=False)
        same_scaffold = bool(np.all(scaffold_ids == scaffold_ids[0]))
        scaffold = (
            scaffold_names[int(scaffold_ids[0])]
            if same_scaffold
            else "multiple"
        )
        start = int(rows[0][3]) if same_scaffold else None
        end = int(rows[-1][4]) if same_scaffold else None
        specs.append(
            WindowSpec(
                window_id=window_id,
                window_mode="locus",
                scaffold=scaffold,
                start=start,
                end=end,
                first_locus=start_idx + 1,
                last_locus=start_idx + int(rows.shape[0]),
                nloci=int(rows.shape[0]),
                sites_total=sum(stop - begin for begin, stop in merged_spans),
                spans=merged_spans,
            )
        )
        window_id += 1

    return specs


def _build_population_rows_from_accum(
    accum: SequenceAccumulator,
    *,
    pop_names: list[str],
    pop_sizes: np.ndarray,
    requested_stats: list[str],
    metadata: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Build ordered population rows from additive accumulators."""
    metadata = {} if metadata is None else metadata
    need_pi = "pi" in requested_stats
    need_theta_w = "theta_w" in requested_stats
    need_tajima_d = "tajima_d" in requested_stats
    need_heterozygosity = any(
        stat in requested_stats for stat in ("heterozygosity", "fis")
    )
    need_fis = "fis" in requested_stats
    rows: list[dict[str, Any]] = []

    for pop_idx, pop in enumerate(pop_names):
        row = dict(metadata)
        row["population"] = pop
        row["n_samples"] = int(pop_sizes[pop_idx])
        if need_pi:
            sites_used = int(accum.sites_used_pi[pop_idx])
            row["sites_used_pi"] = sites_used
            row["pi"] = accum.pi_sum[pop_idx] / sites_used if sites_used else np.nan
        if need_theta_w or need_tajima_d:
            sites_used = int(accum.sites_used_theta[pop_idx])
            segregating_sites = int(accum.segregating_sites[pop_idx])
            row["sites_used_theta"] = sites_used
            row["segregating_sites"] = segregating_sites
            if need_theta_w:
                row["theta_w"] = wattersons_theta(
                    segregating_sites,
                    2 * int(pop_sizes[pop_idx]),
                    sites_used,
                )
            if need_tajima_d:
                row["tajima_d"] = tajimas_d(
                    float(accum.theta_pi_total[pop_idx]),
                    segregating_sites,
                    2 * int(pop_sizes[pop_idx]),
                )
        if need_heterozygosity:
            sites_used = int(accum.sites_used_heterozygosity[pop_idx])
            observed = (
                accum.observed_heterozygosity_sum[pop_idx] / sites_used
                if sites_used
                else np.nan
            )
            expected = (
                accum.expected_heterozygosity_sum[pop_idx] / sites_used
                if sites_used
                else np.nan
            )
            row["sites_used_heterozygosity"] = sites_used
            row["observed_heterozygosity"] = observed
            row["expected_heterozygosity"] = expected
            if need_fis:
                row["fis"] = fis_from_heterozygosity(observed, expected)
        rows.append(row)

    return rows


def _build_pairwise_rows_from_accum(
    accum: SequenceAccumulator,
    *,
    pop_names: list[str],
    pair_indices: list[tuple[int, int]],
    requested_stats: list[str],
    metadata: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Build ordered pairwise rows from additive accumulators."""
    metadata = {} if metadata is None else metadata
    rows: list[dict[str, Any]] = []
    need_dxy = "dxy" in requested_stats
    need_fst = "fst" in requested_stats

    for pair_idx, (idx1, idx2) in enumerate(pair_indices):
        row = dict(metadata)
        row["population1"] = pop_names[idx1]
        row["population2"] = pop_names[idx2]
        row["sites_used"] = int(accum.pair_sites_used[pair_idx])
        if need_dxy:
            row["dxy"] = (
                accum.pair_dxy_sum[pair_idx] / row["sites_used"]
                if row["sites_used"]
                else np.nan
            )
        if need_fst:
            denom = accum.pair_fst_den_sum[pair_idx]
            row["fst"] = (
                accum.pair_fst_num_sum[pair_idx] / denom
                if denom > 0
                else np.nan
            )
        rows.append(row)

    return rows


def _build_sfs_rows_from_accum(
    accum: SequenceAccumulator,
    *,
    pop_names: list[str],
) -> list[dict[str, Any]]:
    """Build ordered folded-SFS rows from additive accumulators."""
    rows: list[dict[str, Any]] = []
    for pop_idx, pop in enumerate(pop_names):
        for mac, count in enumerate(accum.sfs_counts[pop_idx]):
            if count <= 0:
                continue
            rows.append(
                {
                    "population": pop,
                    "minor_allele_count": mac,
                    "site_count": int(count),
                }
            )
    return rows


def _window_metadata(spec: WindowSpec) -> dict[str, Any]:
    """Return ordered metadata fields for one window."""
    return {
        "window_id": spec.window_id,
        "window_mode": spec.window_mode,
        "scaffold": spec.scaffold,
        "start": spec.start,
        "end": spec.end,
        "first_locus": spec.first_locus,
        "last_locus": spec.last_locus,
        "nloci": spec.nloci,
        "sites_total": spec.sites_total,
    }


def _window_has_retained_sites(
    accum: SequenceAccumulator,
    *,
    requested_stats: list[str],
    pair_indices: list[tuple[int, int]],
) -> bool:
    """Return True if one window retained sites for any requested scalar stat."""
    if "pi" in requested_stats and np.any(accum.sites_used_pi > 0):
        return True
    if any(stat in requested_stats for stat in ("theta_w", "tajima_d")) and np.any(
        accum.sites_used_theta > 0
    ):
        return True
    if any(stat in requested_stats for stat in ("heterozygosity", "fis")) and np.any(
        accum.sites_used_heterozygosity > 0
    ):
        return True
    if pair_indices and any(stat in requested_stats for stat in ("dxy", "fst")) and np.any(
        accum.pair_sites_used > 0
    ):
        return True
    return False


def _localize_window_to_chunk(
    *,
    chunk_spans: tuple[tuple[int, int], ...],
    window_spans: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, int], ...]:
    """Project one window's phy-coordinate spans into local chunk offsets."""
    local_spans: list[tuple[int, int]] = []
    offset = 0
    for chunk_start, chunk_end in chunk_spans:
        for window_start, window_end in window_spans:
            overlap_start = max(chunk_start, window_start)
            overlap_end = min(chunk_end, window_end)
            if overlap_start >= overlap_end:
                continue
            local_spans.append(
                (
                    offset + (overlap_start - chunk_start),
                    offset + (overlap_end - chunk_start),
                )
            )
        offset += chunk_end - chunk_start
    return tuple(local_spans)


def _plan_sequence_windows(
    *,
    data: Path,
    window_size: int | None,
    step_size: int | None,
    loci_per_window: int | None,
    locus_step: int | None,
) -> tuple[list[WindowSpec], dict[str, Any]]:
    """Return planned windows and manifest summary metadata."""
    mode = "genomic" if window_size is not None else "locus"
    phymap, scaffold_names, scaffold_lengths = _load_sequence_window_context(data)
    if mode == "genomic":
        specs = _plan_genomic_windows(
            phymap,
            scaffold_names,
            scaffold_lengths,
            window_size=window_size,
            step_size=step_size,
        )
    else:
        specs = _plan_locus_windows(
            phymap,
            scaffold_names,
            loci_per_window=loci_per_window,
            locus_step=locus_step,
        )

    summary = {
        "window_mode": mode,
        "window_size": window_size,
        "step_size": step_size,
        "loci_per_window": loci_per_window,
        "locus_step": locus_step,
        "windows_planned": len(specs),
        "windows_written": 0,
        "windows_skipped": len(specs),
    }
    return specs, summary


def _plan_chunk_window_views(
    chunk_spans: list[tuple[tuple[int, int], ...]],
    specs: list[WindowSpec],
) -> dict[int, tuple[WindowChunkView, ...]]:
    """Return per-chunk local window projections for one set of window specs."""
    mapping: dict[int, tuple[WindowChunkView, ...]] = {}
    planned_specs = [spec for spec in specs if spec.spans]
    for chunk_idx, spans in enumerate(chunk_spans):
        views: list[WindowChunkView] = []
        for spec in planned_specs:
            local_spans = _localize_window_to_chunk(
                chunk_spans=spans,
                window_spans=spec.spans,
            )
            if local_spans:
                views.append(
                    WindowChunkView(
                        window_id=spec.window_id,
                        local_spans=local_spans,
                    )
                )
        if views:
            mapping[chunk_idx] = tuple(views)
    return mapping

def _build_window_output_tables(
    *,
    specs: list[WindowSpec],
    window_totals: dict[int, SequenceAccumulator],
    pop_names: list[str],
    pop_sizes: np.ndarray,
    pair_indices: list[tuple[int, int]],
    requested_stats: list[str],
    summary: dict[str, Any],
) -> tuple[pd.DataFrame | None, pd.DataFrame | None, dict[str, Any]]:
    """Return ordered window tables from one-pass window accumulators."""
    need_pairwise = any(stat in requested_stats for stat in ("dxy", "fst")) and bool(pair_indices)
    population_rows: list[dict[str, Any]] = []
    pairwise_rows: list[dict[str, Any]] = []
    windows_written = 0
    windows_skipped = 0
    for spec in specs:
        accum = window_totals.get(spec.window_id)
        if accum is None or not _window_has_retained_sites(
            accum,
            requested_stats=requested_stats,
            pair_indices=pair_indices,
        ):
            windows_skipped += 1
            continue
        metadata = _window_metadata(spec)
        population_rows.extend(
            _build_population_rows_from_accum(
                accum,
                pop_names=pop_names,
                pop_sizes=pop_sizes,
                requested_stats=requested_stats,
                metadata=metadata,
            )
        )
        if need_pairwise:
            pairwise_rows.extend(
                _build_pairwise_rows_from_accum(
                    accum,
                    pop_names=pop_names,
                    pair_indices=pair_indices,
                    requested_stats=requested_stats,
                    metadata=metadata,
                )
            )
        windows_written += 1

    summary["windows_written"] = windows_written
    summary["windows_skipped"] = windows_skipped
    return (
        build_population_stats_dataframe(
            population_rows,
            requested_stats=requested_stats,
            include_window_metadata=True,
        ),
        build_pairwise_stats_dataframe(
            pairwise_rows,
            requested_stats=requested_stats,
            include_window_metadata=True,
        ),
        summary,
    )


def run_sequence_popgen(
    *,
    data: Path,
    requested_stats: list[str],
    min_sample_coverage: float,
    max_sample_missing: float,
    imap,
    minmap,
    exclude,
    include_reference: bool,
    cores: int,
    log_level: str,
    window_size: int | None,
    step_size: int | None,
    loci_per_window: int | None,
    locus_step: int | None,
) -> dict[str, Any]:
    """Compute genome-wide sequence-backed population-genetic summaries."""
    raw_imap, raw_minmap = normalize_sequence_population_inputs(imap, minmap)
    user_imap = bool(raw_imap)
    chunk_spans = _plan_sequence_chunk_spans(data)

    initial_snames, initial_sidxs, _exclude_set = resolve_sequence_sample_subset(
        data,
        exclude=exclude,
        include_reference=include_reference,
        imap=raw_imap,
    )
    final_imap, final_minmap = build_sequence_imap_minmap(
        initial_snames,
        min_sample_coverage=min_sample_coverage,
        imap=raw_imap,
        minmap=raw_minmap,
    )

    if max_sample_missing < 1.0:
        initial_missing, total_sites = _calculate_sample_missing(
            data,
            initial_sidxs,
            initial_snames,
            chunk_spans,
            cores=cores,
            log_level=log_level,
        )
        keep_names = [
            name for name in initial_snames if float(initial_missing.loc[name]) <= max_sample_missing
        ]
    else:
        initial_missing = None
        total_sites = 0
        keep_names = list(initial_snames)

    if not keep_names:
        raise IPyradError("No samples passed max_sample_missing filter.")

    dropped = [name for name in initial_snames if name not in keep_names]
    final_snames = list(keep_names)
    sname_to_index = {name: idx for idx, name in enumerate(initial_snames)}
    final_sidxs = [initial_sidxs[sname_to_index[name]] for name in final_snames]
    final_imap, final_minmap = sync_sequence_imap_after_sample_drop(
        final_snames,
        user_imap=user_imap,
        imap=final_imap,
        minmap=final_minmap,
        min_sample_coverage=min_sample_coverage,
    )

    pop_names = list(final_imap)
    pop_name_to_index = {name: idx for idx, name in enumerate(final_snames)}
    pop_indices = [
        np.array([pop_name_to_index[name] for name in final_imap[pop]], dtype=np.int64)
        for pop in pop_names
    ]
    pop_sizes = np.array([len(final_imap[pop]) for pop in pop_names], dtype=np.int64)
    pop_minmap = np.array([int(final_minmap[pop]) for pop in pop_names], dtype=np.int64)
    pair_indices = [
        (idx1, idx2)
        for idx1 in range(len(pop_names))
        for idx2 in range(idx1 + 1, len(pop_names))
    ]

    need_pi = "pi" in requested_stats
    need_theta = any(stat in requested_stats for stat in ("theta_w", "tajima_d"))
    need_heterozygosity = any(
        stat in requested_stats for stat in ("heterozygosity", "fis")
    )
    need_global_fit = "fit" in requested_stats
    need_pairwise = any(stat in requested_stats for stat in ("dxy", "fst")) and bool(pair_indices)
    need_sfs = "sfs" in requested_stats
    scalar_stats_requested = any(stat != "sfs" for stat in requested_stats)
    window_specs: list[WindowSpec] = []
    window_summary: dict[str, Any] = {
        "window_mode": "none",
        "window_size": None,
        "step_size": None,
        "loci_per_window": None,
        "locus_step": None,
        "windows_planned": 0,
        "windows_written": 0,
        "windows_skipped": 0,
    }
    chunk_window_views: dict[int, tuple[WindowChunkView, ...]] | None = None
    if window_size is not None or loci_per_window is not None:
        window_specs, window_summary = _plan_sequence_windows(
            data=data,
            window_size=window_size,
            step_size=step_size,
            loci_per_window=loci_per_window,
            locus_step=locus_step,
        )
        if scalar_stats_requested and window_specs:
            chunk_window_views = _plan_chunk_window_views(chunk_spans, window_specs)

    accum, window_totals = _collect_sequence_accumulators(
        data=data,
        sidxs=final_sidxs,
        chunk_spans=chunk_spans,
        pop_indices=pop_indices,
        pop_sizes=pop_sizes,
        pop_minmap=pop_minmap,
        pair_indices=pair_indices,
        track_sample_stats=True,
        need_pi=need_pi,
        need_theta=need_theta,
        need_heterozygosity=need_heterozygosity,
        need_global_fit=need_global_fit,
        global_min_sample_coverage=int(min_sample_coverage),
        need_pairwise=need_pairwise,
        need_sfs=need_sfs,
        chunk_window_views=chunk_window_views,
        cores=cores,
        log_level=log_level,
    )

    if initial_missing is None:
        total_sites = int(accum.total_sites)
        values = np.zeros(len(final_snames), dtype=float)
        if total_sites:
            values = accum.missing_counts / total_sites
        sample_missing = pd.Series(values, index=final_snames, dtype=float)
    else:
        sample_missing = initial_missing.loc[final_snames]
    sample_data_summary = build_sample_data_summary(
        samples=final_snames,
        missing_fraction=sample_missing,
        post_imputation_missing_fraction=sample_missing,
        imputation_algorithm="not-imputed",
    )
    sample_sites_total = int(accum.total_sites)
    sample_sites_missing = accum.missing_counts.astype(np.int64, copy=False)
    sample_sites_called = (
        np.full(len(final_snames), sample_sites_total, dtype=np.int64)
        - sample_sites_missing
    )
    sample_stats = build_sample_stats_dataframe(
        samples=final_snames,
        imap=final_imap,
        sites_total=sample_sites_total,
        sites_called=sample_sites_called,
        sites_missing=sample_sites_missing,
        heterozygous_sites=accum.heterozygous_counts,
    )
    global_stats_df = None
    if need_global_fit:
        global_sites_used = int(accum.global_sites_used_heterozygosity)
        global_observed = (
            float(accum.global_observed_heterozygosity_sum) / global_sites_used
            if global_sites_used
            else np.nan
        )
        global_expected = (
            float(accum.global_expected_heterozygosity_sum) / global_sites_used
            if global_sites_used
            else np.nan
        )
        global_stats_df = build_global_stats_dataframe(
            sites_used_heterozygosity=global_sites_used,
            observed_heterozygosity=global_observed,
            expected_heterozygosity_total=global_expected,
        )

    population_df = build_population_stats_dataframe(
        _build_population_rows_from_accum(
            accum,
            pop_names=pop_names,
            pop_sizes=pop_sizes,
            requested_stats=requested_stats,
        ),
        requested_stats=requested_stats,
        include_window_metadata=False,
    )
    pairwise_df = build_pairwise_stats_dataframe(
        _build_pairwise_rows_from_accum(
            accum,
            pop_names=pop_names,
            pair_indices=pair_indices,
            requested_stats=requested_stats,
        ) if need_pairwise else [],
        requested_stats=requested_stats,
        include_window_metadata=False,
    )
    sfs_df = build_sfs_dataframe(
        _build_sfs_rows_from_accum(accum, pop_names=pop_names)
    ) if need_sfs else None

    window_population_df = None
    window_pairwise_df = None
    if window_specs and scalar_stats_requested:
        (
            window_population_df,
            window_pairwise_df,
            window_summary,
        ) = _build_window_output_tables(
            specs=window_specs,
            window_totals=window_totals,
            pop_names=pop_names,
            pop_sizes=pop_sizes,
            pair_indices=pair_indices,
            requested_stats=requested_stats,
            summary=window_summary,
        )

    summary = {
        "input_backend": "sequence",
        "requested_stats": list(requested_stats),
        "samples_selected_initial": initial_snames,
        "samples_dropped_by_max_missing": dropped,
        "samples_final": final_snames,
        "imap": final_imap,
        "minmap": final_minmap,
        "include_reference": include_reference,
        "max_sample_missing": max_sample_missing,
        "sites_considered_for_missingness": total_sites,
        **window_summary,
    }
    if window_summary["window_mode"] != "none" and "sfs" in requested_stats:
        summary["window_sfs_note"] = "Windowed SFS is not written in this phase."
    logger.info("computed sequence-backed popgen statistics for {} population(s)", len(final_imap))
    return PopgenResult(
        sample_data_summary=sample_data_summary,
        sample_stats=sample_stats,
        global_stats=global_stats_df,
        population_stats=population_df,
        pairwise_stats=pairwise_df,
        sfs=sfs_df,
        window_population_stats=window_population_df,
        window_pairwise_stats=window_pairwise_df,
        summary=summary,
    )
