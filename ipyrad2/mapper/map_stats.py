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
        "input_reads": int(job_result.nreads_processed),
        "reads_removed_unmapped_or_nonprimary": int(
            job_result.nreads_filtered_before_bam_by_unmapped_or_nonprimary
        ),
        "reads_in_final_bam": int(nreads_in_final_bam),
        "fraction_input_reads_retained_in_final_bam": _safe_ratio(
            nreads_in_final_bam,
            job_result.nreads_processed,
        ),
        "reads_failing_min_mapq_20": int(reads_below_mapq_20),
        "reads_failing_max_softclip_25": int(reads_above_soft_clip_25),
        "reads_failing_max_nm_50": int(reads_at_or_above_nm_50),
        "reads_passing_all_preview_filters": int(reads_passing_reporting_thresholds),
        "fraction_reads_passing_all_preview_filters": _safe_ratio(
            reads_passing_reporting_thresholds,
            nreads_in_final_bam,
        ),
        "mapq_mean": float(read_mapq_mean),
        "mapq_median": float(read_mapq_median),
        "mapq_stdev": float(read_mapq_stdev),
        "softclip_mean": float(read_soft_clip_mean),
        "softclip_median": float(read_soft_clip_median),
        "softclip_stdev": float(read_soft_clip_stdev),
        "nm_mean": float(read_nm_mean),
        "nm_median": float(read_nm_median),
        "nm_stdev": float(read_nm_stdev),
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
        "input_templates": int(job_result.nreads_processed // 2),
        "reads_removed_unmapped_or_nonprimary": int(
            job_result.nreads_filtered_before_bam_by_unmapped_or_nonprimary
        ),
        "reads_removed_same_scaffold_pairing": int(
            job_result.nreads_filtered_before_bam_by_mate_unmapped_or_cross_scaffold
        ),
        "templates_in_final_bam": int(npairs_evaluated_in_final_bam),
        "fraction_input_templates_retained_in_final_bam": _safe_ratio(
            npairs_evaluated_in_final_bam,
            job_result.nreads_processed // 2,
        ),
        "duplicate_records_removed": _duplicate_stat(
            job_result.duplicate_stats,
            "duplicate_total",
        ),
        "templates_failing_min_mapq_20": int(pairs_with_a_read_below_mapq_20),
        "templates_failing_max_softclip_25": int(pairs_with_a_read_above_soft_clip_25),
        "templates_failing_max_nm_50": int(pairs_with_a_read_at_or_above_nm_50),
        "templates_failing_max_abs_tlen_2000": int(pairs_above_abs_tlen_2000),
        "templates_passing_all_preview_filters": int(pairs_passing_reporting_thresholds),
        "fraction_templates_passing_all_preview_filters": _safe_ratio(
            pairs_passing_reporting_thresholds,
            npairs_evaluated_in_final_bam,
        ),
        "min_mapq_mean": float(pair_min_mapq_mean),
        "min_mapq_median": float(pair_min_mapq_median),
        "min_mapq_stdev": float(pair_min_mapq_stdev),
        "max_softclip_mean": float(pair_max_soft_clip_mean),
        "max_softclip_median": float(pair_max_soft_clip_median),
        "max_softclip_stdev": float(pair_max_soft_clip_stdev),
        "max_nm_mean": float(pair_max_nm_mean),
        "max_nm_median": float(pair_max_nm_median),
        "max_nm_stdev": float(pair_max_nm_stdev),
        "abs_tlen_mean": float(pair_abs_tlen_mean),
        "abs_tlen_median": float(pair_abs_tlen_median),
        "abs_tlen_stdev": float(pair_abs_tlen_stdev),
    }


def _report_header(is_paired: bool) -> str:
    """Return a short header block describing map reporting semantics."""
    paired_note = (
        "# Paired-end final BAMs keep only mapped mates on the same scaffold.\n"
        if is_paired
        else ""
    )
    return (
        "# ipyrad2 map stats\n"
        "# Final BAMs are coordinate sorted and indexed.\n"
        f"{paired_note}"
        "\n"
    )


def _format_frame(df: pd.DataFrame) -> str:
    """Format one stats table for display."""
    float_columns = {
        column
        for column in df.columns
        if column.startswith("fraction_")
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
    return formatted.to_string() + "\n"


def _section_frame(
    stats: dict[str, dict],
    columns: list[str],
    *,
    float_defaults: set[str] | None = None,
) -> pd.DataFrame:
    """Build a frame for one report section with stable columns and defaults."""
    float_defaults = float_defaults or set()
    rows = {}
    for sname in sorted(stats):
        sample_stats = stats[sname]
        rows[sname] = {
            column: sample_stats.get(column, float("nan") if column in float_defaults else 0)
            for column in columns
        }
    return pd.DataFrame.from_dict(rows, orient="index", columns=columns)


def build_map_stats_payload(
    stats: dict[str, dict],
    is_paired: bool,
    logged_command: str | None = None,
) -> dict[str, object]:
    """Build the structured run-level mapper stats payload."""
    if is_paired:
        applied_columns = [
            "input_templates",
            "reads_removed_unmapped_or_nonprimary",
            "reads_removed_same_scaffold_pairing",
            "duplicate_records_removed",
            "templates_in_final_bam",
            "fraction_input_templates_retained_in_final_bam",
        ]
        preview_effect_columns = [
            "templates_failing_min_mapq_20",
            "templates_failing_max_softclip_25",
            "templates_failing_max_nm_50",
            "templates_failing_max_abs_tlen_2000",
            "templates_passing_all_preview_filters",
            "fraction_templates_passing_all_preview_filters",
        ]
        preview_summary_columns = [
            "min_mapq_mean",
            "min_mapq_median",
            "min_mapq_stdev",
            "max_softclip_mean",
            "max_softclip_median",
            "max_softclip_stdev",
            "max_nm_mean",
            "max_nm_median",
            "max_nm_stdev",
            "abs_tlen_mean",
            "abs_tlen_median",
            "abs_tlen_stdev",
        ]
        preview_flags = "-qm/--min-map-q, -ms/--max-softclip, -me/--max-nm, -mt/--max-tlen"
        preview_mode_note = "# Preview mode: pair-level thresholds evaluated on final BAM templates.\n"
    else:
        applied_columns = [
            "input_reads",
            "reads_removed_unmapped_or_nonprimary",
            "reads_in_final_bam",
            "fraction_input_reads_retained_in_final_bam",
        ]
        preview_effect_columns = [
            "reads_failing_min_mapq_20",
            "reads_failing_max_softclip_25",
            "reads_failing_max_nm_50",
            "reads_passing_all_preview_filters",
            "fraction_reads_passing_all_preview_filters",
        ]
        preview_summary_columns = [
            "mapq_mean",
            "mapq_median",
            "mapq_stdev",
            "softclip_mean",
            "softclip_median",
            "softclip_stdev",
            "nm_mean",
            "nm_median",
            "nm_stdev",
        ]
        preview_flags = "-qm/--min-map-q, -ms/--max-softclip, -me/--max-nm"
        preview_mode_note = "# Preview mode: read-level thresholds evaluated on final BAM reads.\n"

    applied_frame = _section_frame(stats, applied_columns, float_defaults={
        column for column in applied_columns if column.startswith("fraction_")
    })
    preview_effect_frame = _section_frame(stats, preview_effect_columns, float_defaults={
        column for column in preview_effect_columns if column.startswith("fraction_")
    })
    preview_summary_frame = _section_frame(stats, preview_summary_columns, float_defaults=set(preview_summary_columns))

    payload: dict[str, object] = {
        "is_paired": is_paired,
        "applied_mapping_summary": (
            applied_frame.rename_axis("sample").reset_index().to_dict(orient="records")
        ),
        "assemble_read_filter_preview": {
            "description": "These preview thresholds were not applied during mapping.",
            "flags": preview_flags,
            "mode_note": preview_mode_note.strip(),
            "mapq_threshold": MAPQ_REPORT_THRESHOLD,
            "soft_clipped_bases_threshold": SOFT_CLIP_REPORT_THRESHOLD,
            "nm_threshold": NM_REPORT_THRESHOLD,
            "absolute_tlen_threshold": ABS_TLEN_REPORT_THRESHOLD if is_paired else None,
            "filter_effects": (
                preview_effect_frame.rename_axis("sample").reset_index().to_dict(orient="records")
            ),
            "metric_summaries": (
                preview_summary_frame.rename_axis("sample").reset_index().to_dict(orient="records")
            ),
        },
    }
    if logged_command:
        payload["command"] = logged_command
    return payload


def _frame_from_payload_records(records: list[dict]) -> pd.DataFrame:
    """Build one report frame from payload records with a stable sample index."""
    if not records:
        return pd.DataFrame(index=pd.Index([], name="sample"))
    return pd.DataFrame.from_records(records).set_index("sample")


def render_map_stats_payload_report(payload: dict[str, object]) -> str:
    """Render the final mapper stats report from a structured payload."""
    is_paired = bool(payload["is_paired"])
    applied_frame = _frame_from_payload_records(payload["applied_mapping_summary"])
    preview = payload["assemble_read_filter_preview"]
    preview_effect_frame = _frame_from_payload_records(preview["filter_effects"])
    preview_summary_frame = _frame_from_payload_records(preview["metric_summaries"])

    command_header = ""
    if payload.get("command"):
        command_header = f"CMD: {payload['command']}\n\n"

    return (
        command_header
        + _report_header(is_paired)
        + "## Applied mapping summary\n"
        + "# These counts describe filters already applied during ipyrad2 map.\n\n"
        + _format_frame(applied_frame)
        + "\n## Assemble read-filter preview (not applied during mapping)\n"
        + f"# {preview['description']}\n"
        + f"# Use them to guide ipyrad2 assemble read filters: {preview['flags']}.\n"
        + f"{preview['mode_note']}\n"
        + f"# MAPQ threshold: {preview['mapq_threshold']}\n"
        + f"# Soft-clipped bases threshold: {preview['soft_clipped_bases_threshold']}\n"
        + f"# NM threshold: {preview['nm_threshold']}\n"
        + (
            f"# Absolute TLEN threshold: {preview['absolute_tlen_threshold']}\n"
            if preview["absolute_tlen_threshold"] is not None
            else ""
        )
        + "\n### Preview filter effects\n"
        + _format_frame(preview_effect_frame)
        + "\n### Preview metric summaries\n"
        + _format_frame(preview_summary_frame)
        + "\n"
    )


def render_map_stats_report(
    stats: dict[str, dict],
    is_paired: bool,
    logged_command: str | None = None,
) -> str:
    """Render the final mapper stats report."""
    payload = build_map_stats_payload(
        stats,
        is_paired,
        logged_command=logged_command,
    )
    return render_map_stats_payload_report(payload)
