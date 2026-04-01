#!/usr/bin/env python

"""Stats helpers for ipyrad2 map final BAM reporting."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math
from pathlib import Path
import re
import sys
from typing import Iterable, Iterator, Sequence

import pandas as pd

from ..utils.parallel import stream_pipeline_lines


BIN = Path(sys.prefix) / "bin"
BIN_SAMTOOLS = str(BIN / "samtools")

MAPQ_REPORT_THRESHOLD = 20
SOFT_CLIP_REPORT_THRESHOLD = 25
NM_REPORT_THRESHOLD = 50
ABS_TLEN_REPORT_THRESHOLD = 2000


@dataclass(frozen=True)
class MappingJobResult:
    """Summary returned from one mapping worker job."""

    sname: str
    bam_path: Path
    is_paired: bool
    nreads_processed: int
    nreads_filtered_before_bam_by_unmapped_or_nonprimary: int
    nreads_written_before_duplicate_removal: int
    duplicate_stats: dict[str, int]
    nreads_filtered_before_bam_by_mate_unmapped_or_cross_scaffold: int = 0


@dataclass(frozen=True)
class _SamRecord:
    """Minimal parsed SAM fields required for map reporting stats."""

    qname: str
    is_read1: bool
    is_read2: bool
    rname: str
    mapq: int
    soft_clip: int
    nm: int
    abs_tlen: int


def _normalize_stat_key(label: str) -> str:
    """Normalize one markdup stats label to snake_case."""
    return re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")


def parse_markdup_report(path: Path) -> dict[str, int]:
    """Parse one samtools markdup stats report into a normalized dict."""
    stats: dict[str, int] = {}
    with path.open("r") as infile:
        for raw_line in infile:
            line = raw_line.strip()
            if not line:
                continue
            if ":" in line:
                label, value = line.rsplit(":", 1)
            else:
                label, value = line.rsplit(maxsplit=1)
            try:
                stats[_normalize_stat_key(label)] = int(value)
            except ValueError:
                continue
    return stats


def _duplicate_stat(stats: dict[str, int], *candidates: str) -> int:
    """Return the first duplicate-removal stat found among candidate keys."""
    for key in candidates:
        if key in stats:
            return int(stats[key])
    return 0


def _soft_clip_from_cigar(cigar: str) -> int:
    """Return the number of soft-clipped bases in one CIGAR string."""
    if cigar == "*" or not cigar:
        return 0
    total = 0
    digits: list[str] = []
    for char in cigar:
        if char.isdigit():
            digits.append(char)
            continue
        length = int("".join(digits)) if digits else 0
        if char == "S":
            total += length
        digits.clear()
    return total


def _nm_from_fields(fields: Sequence[str]) -> int:
    """Extract the NM tag from SAM optional fields."""
    for field in fields[11:]:
        if field.startswith("NM:i:"):
            return int(field[5:])
    return 0


def _parse_sam_record(line: str) -> _SamRecord:
    """Parse the minimal SAM fields needed for mapper reporting stats."""
    fields = line.split("\t")
    flag = int(fields[1])
    return _SamRecord(
        qname=fields[0],
        is_read1=bool(flag & 0x40),
        is_read2=bool(flag & 0x80),
        rname=fields[2],
        mapq=int(fields[4]),
        soft_clip=_soft_clip_from_cigar(fields[5]),
        nm=_nm_from_fields(fields),
        abs_tlen=abs(int(fields[8])),
    )


def _iter_templates_from_qname_group(
    records: Sequence[_SamRecord],
) -> Iterator[tuple[_SamRecord, ...]]:
    """Yield pair templates and singleton templates for one qname group."""
    read1s = []
    read2s = []
    other = []
    for record in records:
        if record.is_read1 and not record.is_read2:
            read1s.append(record)
        elif record.is_read2 and not record.is_read1:
            read2s.append(record)
        else:
            other.append(record)

    npairs = min(len(read1s), len(read2s))
    for idx in range(npairs):
        yield (read1s[idx], read2s[idx])
    for record in read1s[npairs:]:
        yield (record,)
    for record in read2s[npairs:]:
        yield (record,)
    for record in other:
        yield (record,)


def _iter_templates_from_collated_lines(
    lines: Iterable[str],
) -> Iterator[tuple[_SamRecord, ...]]:
    """Yield templates from a collated SAM stream."""
    current_qname: str | None = None
    current_records: list[_SamRecord] = []
    for line in lines:
        record = _parse_sam_record(line)
        if current_qname is None:
            current_qname = record.qname
        if record.qname != current_qname:
            yield from _iter_templates_from_qname_group(current_records)
            current_qname = record.qname
            current_records = [record]
            continue
        current_records.append(record)
    if current_records:
        yield from _iter_templates_from_qname_group(current_records)


def _safe_ratio(numer: int, denom: int) -> float:
    """Return numer/denom or NaN when the denominator is zero."""
    if denom == 0:
        return float("nan")
    return float(numer / denom)


def _median_from_histogram(hist: Counter[int], total: int) -> float:
    """Return the median represented by a sparse integer histogram."""
    left_rank = (total - 1) // 2
    right_rank = total // 2
    cumulative = 0
    left_value = None
    right_value = None
    for value in sorted(hist):
        cumulative += hist[value]
        if left_value is None and cumulative > left_rank:
            left_value = value
        if right_value is None and cumulative > right_rank:
            right_value = value
            break
    assert left_value is not None and right_value is not None
    return float((left_value + right_value) / 2)


def _histogram_summary(hist: Counter[int]) -> tuple[float, float, float]:
    """Return mean, median, and stdev for a sparse histogram."""
    total = sum(hist.values())
    if total == 0:
        return float("nan"), float("nan"), float("nan")
    mean = sum(value * count for value, count in hist.items()) / total
    median = _median_from_histogram(hist, total)
    variance = sum(((value - mean) ** 2) * count for value, count in hist.items()) / total
    return float(mean), float(median), float(math.sqrt(variance))


def collect_single_end_bam_stats(job_result: MappingJobResult) -> dict:
    """Collect read-level reporting stats from one final single-end BAM."""
    cmd = [
        BIN_SAMTOOLS, "view",
        str(job_result.bam_path),
    ]

    mapq_hist: Counter[int] = Counter()
    soft_hist: Counter[int] = Counter()
    nm_hist: Counter[int] = Counter()
    nreads_in_final_bam = 0
    reads_below_mapq_20 = 0
    reads_above_soft_clip_25 = 0
    reads_at_or_above_nm_50 = 0
    reads_passing_reporting_thresholds = 0

    for line in stream_pipeline_lines([cmd]):
        record = _parse_sam_record(line)
        nreads_in_final_bam += 1
        mapq_hist[record.mapq] += 1
        soft_hist[record.soft_clip] += 1
        nm_hist[record.nm] += 1

        mapq_fail = record.mapq < MAPQ_REPORT_THRESHOLD
        soft_fail = record.soft_clip > SOFT_CLIP_REPORT_THRESHOLD
        nm_fail = record.nm >= NM_REPORT_THRESHOLD
        reads_below_mapq_20 += int(mapq_fail)
        reads_above_soft_clip_25 += int(soft_fail)
        reads_at_or_above_nm_50 += int(nm_fail)
        if not any((mapq_fail, soft_fail, nm_fail)):
            reads_passing_reporting_thresholds += 1

    read_mapq_mean, read_mapq_median, read_mapq_stdev = _histogram_summary(mapq_hist)
    read_soft_clip_mean, read_soft_clip_median, read_soft_clip_stdev = _histogram_summary(soft_hist)
    read_nm_mean, read_nm_median, read_nm_stdev = _histogram_summary(nm_hist)
    return {
        "nreads_processed": int(job_result.nreads_processed),
        "nreads_filtered_before_bam_by_unmapped_or_nonprimary": int(
            job_result.nreads_filtered_before_bam_by_unmapped_or_nonprimary
        ),
        "nreads_filtered_before_bam_by_mate_unmapped_or_cross_scaffold": int(
            job_result.nreads_filtered_before_bam_by_mate_unmapped_or_cross_scaffold
        ),
        "nreads_written_before_duplicate_removal": int(
            job_result.nreads_written_before_duplicate_removal
        ),
        "nreads_in_final_bam": int(nreads_in_final_bam),
        "proportion_retained_in_final_bam": _safe_ratio(
            nreads_in_final_bam,
            job_result.nreads_processed,
        ),
        "reads_below_mapq_20": int(reads_below_mapq_20),
        "reads_above_soft_clip_25": int(reads_above_soft_clip_25),
        "reads_at_or_above_nm_50": int(reads_at_or_above_nm_50),
        "reads_passing_reporting_thresholds": int(reads_passing_reporting_thresholds),
        "proportion_reads_passing_reporting_thresholds": _safe_ratio(
            reads_passing_reporting_thresholds,
            nreads_in_final_bam,
        ),
        "read_mapq_mean": float(read_mapq_mean),
        "read_mapq_median": float(read_mapq_median),
        "read_mapq_stdev": float(read_mapq_stdev),
        "read_soft_clip_mean": float(read_soft_clip_mean),
        "read_soft_clip_median": float(read_soft_clip_median),
        "read_soft_clip_stdev": float(read_soft_clip_stdev),
        "read_nm_mean": float(read_nm_mean),
        "read_nm_median": float(read_nm_median),
        "read_nm_stdev": float(read_nm_stdev),
    }


def collect_paired_bam_stats(job_result: MappingJobResult) -> dict:
    """Collect exact pair-level reporting stats from one final paired-end BAM."""
    tmp_prefix = job_result.bam_path.parent / "tmpdir" / f"{job_result.sname}.stats_collate"
    cmd1 = [
        BIN_SAMTOOLS, "collate",
        "-O",
        "-u",
        "-T", str(tmp_prefix),
        str(job_result.bam_path),
    ]
    cmd2 = [
        BIN_SAMTOOLS, "view",
        "-",
    ]

    min_mapq_hist: Counter[int] = Counter()
    max_soft_hist: Counter[int] = Counter()
    max_nm_hist: Counter[int] = Counter()
    abs_tlen_hist: Counter[int] = Counter()

    npairs_evaluated_in_final_bam = 0
    npairs_with_both_mates_in_final_bam = 0
    nsingletons_in_final_bam = 0
    pairs_with_a_read_below_mapq_20 = 0
    pairs_with_a_read_above_soft_clip_25 = 0
    pairs_with_a_read_at_or_above_nm_50 = 0
    pairs_failing_same_reference_pairing = 0
    pairs_above_abs_tlen_2000 = 0
    pairs_passing_reporting_thresholds = 0

    try:
        for template in _iter_templates_from_collated_lines(stream_pipeline_lines([cmd1, cmd2])):
            npairs_evaluated_in_final_bam += 1
            mapq_values = [record.mapq for record in template]
            soft_values = [record.soft_clip for record in template]
            nm_values = [record.nm for record in template]
            min_mapq = min(mapq_values)
            max_soft = max(soft_values)
            max_nm = max(nm_values)
            min_mapq_hist[min_mapq] += 1
            max_soft_hist[max_soft] += 1
            max_nm_hist[max_nm] += 1

            has_both_mates = len(template) == 2
            if has_both_mates:
                npairs_with_both_mates_in_final_bam += 1
            else:
                nsingletons_in_final_bam += 1

            mapq_fail = min_mapq < MAPQ_REPORT_THRESHOLD
            soft_fail = max_soft > SOFT_CLIP_REPORT_THRESHOLD
            nm_fail = max_nm >= NM_REPORT_THRESHOLD
            same_reference_fail = (not has_both_mates) or (template[0].rname != template[1].rname)
            tlen_fail = False
            if has_both_mates and not same_reference_fail:
                pair_abs_tlen = max(record.abs_tlen for record in template)
                abs_tlen_hist[pair_abs_tlen] += 1
                tlen_fail = pair_abs_tlen > ABS_TLEN_REPORT_THRESHOLD

            pairs_with_a_read_below_mapq_20 += int(mapq_fail)
            pairs_with_a_read_above_soft_clip_25 += int(soft_fail)
            pairs_with_a_read_at_or_above_nm_50 += int(nm_fail)
            pairs_failing_same_reference_pairing += int(same_reference_fail)
            pairs_above_abs_tlen_2000 += int(tlen_fail)
            if not any((mapq_fail, soft_fail, nm_fail, same_reference_fail, tlen_fail)):
                pairs_passing_reporting_thresholds += 1
    finally:
        for path in tmp_prefix.parent.glob(tmp_prefix.name + "*"):
            try:
                path.unlink()
            except FileNotFoundError:
                pass

    pair_min_mapq_mean, pair_min_mapq_median, pair_min_mapq_stdev = _histogram_summary(min_mapq_hist)
    pair_max_soft_clip_mean, pair_max_soft_clip_median, pair_max_soft_clip_stdev = _histogram_summary(max_soft_hist)
    pair_max_nm_mean, pair_max_nm_median, pair_max_nm_stdev = _histogram_summary(max_nm_hist)
    pair_abs_tlen_mean, pair_abs_tlen_median, pair_abs_tlen_stdev = _histogram_summary(abs_tlen_hist)
    return {
        "npairs_processed": int(job_result.nreads_processed // 2),
        "nreads_processed": int(job_result.nreads_processed),
        "nreads_filtered_before_bam_by_unmapped_or_nonprimary": int(
            job_result.nreads_filtered_before_bam_by_unmapped_or_nonprimary
        ),
        "nreads_filtered_before_bam_by_mate_unmapped_or_cross_scaffold": int(
            job_result.nreads_filtered_before_bam_by_mate_unmapped_or_cross_scaffold
        ),
        "npairs_written_before_duplicate_removal": int(
            job_result.nreads_written_before_duplicate_removal // 2
        ),
        "nreads_written_before_duplicate_removal": int(
            job_result.nreads_written_before_duplicate_removal
        ),
        "npairs_evaluated_in_final_bam": int(npairs_evaluated_in_final_bam),
        "npairs_with_both_mates_in_final_bam": int(npairs_with_both_mates_in_final_bam),
        "nsingletons_in_final_bam": int(nsingletons_in_final_bam),
        "duplicate_records_removed_total": _duplicate_stat(
            job_result.duplicate_stats,
            "duplicate_total",
        ),
        "duplicate_pairs_removed": _duplicate_stat(
            job_result.duplicate_stats,
            "duplicate_pair",
            "duplicate_pairs",
        ),
        "duplicate_singletons_removed": _duplicate_stat(
            job_result.duplicate_stats,
            "duplicate_single",
            "duplicate_singleton",
            "duplicate_singletons",
        ),
        "pairs_with_a_read_below_mapq_20": int(pairs_with_a_read_below_mapq_20),
        "pairs_with_a_read_above_soft_clip_25": int(pairs_with_a_read_above_soft_clip_25),
        "pairs_with_a_read_at_or_above_nm_50": int(pairs_with_a_read_at_or_above_nm_50),
        "pairs_failing_same_reference_pairing": int(pairs_failing_same_reference_pairing),
        "pairs_above_abs_tlen_2000": int(pairs_above_abs_tlen_2000),
        "pairs_passing_reporting_thresholds": int(pairs_passing_reporting_thresholds),
        "proportion_pairs_with_both_mates_in_final_bam": _safe_ratio(
            npairs_with_both_mates_in_final_bam,
            npairs_evaluated_in_final_bam,
        ),
        "proportion_pairs_passing_reporting_thresholds": _safe_ratio(
            pairs_passing_reporting_thresholds,
            npairs_evaluated_in_final_bam,
        ),
        "pair_min_mapq_mean": float(pair_min_mapq_mean),
        "pair_min_mapq_median": float(pair_min_mapq_median),
        "pair_min_mapq_stdev": float(pair_min_mapq_stdev),
        "pair_max_soft_clip_mean": float(pair_max_soft_clip_mean),
        "pair_max_soft_clip_median": float(pair_max_soft_clip_median),
        "pair_max_soft_clip_stdev": float(pair_max_soft_clip_stdev),
        "pair_max_nm_mean": float(pair_max_nm_mean),
        "pair_max_nm_median": float(pair_max_nm_median),
        "pair_max_nm_stdev": float(pair_max_nm_stdev),
        "pair_abs_tlen_mean": float(pair_abs_tlen_mean),
        "pair_abs_tlen_median": float(pair_abs_tlen_median),
        "pair_abs_tlen_stdev": float(pair_abs_tlen_stdev),
    }


def _report_header(is_paired: bool) -> str:
    """Return a short header block describing map reporting semantics."""
    mode_text = (
        "pair-level thresholds for paired-end datasets"
        if is_paired
        else "read-level thresholds for single-end datasets"
    )
    paired_note = (
        "# Paired-end final BAMs keep only mapped mates on the same scaffold.\n"
        if is_paired
        else ""
    )
    return (
        "# ipyrad2 map stats\n"
        "# Final BAMs are coordinate sorted and indexed.\n"
        "# The thresholds below are reporting only and are not applied during mapping.\n"
        f"{paired_note}"
        f"# Reporting mode: {mode_text}\n"
        f"# MAPQ threshold: {MAPQ_REPORT_THRESHOLD}\n"
        f"# Soft-clipped bases threshold: {SOFT_CLIP_REPORT_THRESHOLD}\n"
        f"# NM threshold: {NM_REPORT_THRESHOLD}\n"
        f"# Absolute TLEN threshold: {ABS_TLEN_REPORT_THRESHOLD}\n\n"
    )


def render_map_stats_report(stats: dict[str, dict], is_paired: bool) -> str:
    """Render the final mapper stats report."""
    df = pd.DataFrame({name: stats[name] for name in sorted(stats)}).T
    float_columns = {
        column
        for column in df.columns
        if column.startswith("proportion_")
        or column.endswith("_mean")
        or column.endswith("_median")
        or column.endswith("_stdev")
    }
    formatted = df.copy()
    for column in formatted.columns:
        if column in float_columns:
            formatted[column] = formatted[column].map(lambda x: f"{float(x):.3f}")
        else:
            formatted[column] = formatted[column].map(lambda x: f"{int(x)}")
    return _report_header(is_paired) + formatted.to_string() + "\n"
