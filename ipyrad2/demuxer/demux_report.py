#!/usr/bin/env python

"""Reporting helpers for demux logging and stats output."""

from __future__ import annotations

from collections import Counter
from math import ceil
from pathlib import Path
from typing import Dict, List, Tuple

from loguru import logger
import pandas as pd

from ipyrad2.utils.kmers import InferredJunctionSet


DEMUX_STATS_PREFIX = "ipyrad_demux_stats_"
PRESERVED_WARNING_PREVIEW = 3
SUSPECTED_BARCODE_REPORT_LIMIT = 25
SUSPECTED_BARCODE_MIN_COUNT = 50
SUSPECTED_BARCODE_MIN_FRACTION = 0.0001


def format_motif_set(junction: InferredJunctionSet | None) -> str:
    """Return a detailed human-readable description of a motif set."""
    if junction is None:
        return "<none>"
    motifs = ", ".join(junction.motifs) if junction.motifs else "<none>"
    if junction.position_mode == "barcode_boundary" and junction.boundary_supports:
        return f"[{motifs}] from barcode boundaries {junction.position_summary}"
    return f"[{motifs}] at offset {junction.offset}"


def format_logged_motif_set(junction: InferredJunctionSet | None) -> str:
    """Return a concise motif-set description for demux INFO logging."""
    if junction is None:
        return "<none>"
    motifs = ", ".join(junction.motifs) if junction.motifs else "<none>"
    if junction.position_mode == "barcode_boundary" and junction.boundary_supports:
        return f"[{motifs}] inferred from barcode boundaries"
    return f"[{motifs}] at offset {junction.offset}"


def warn_multi_motif_inference(
    read_label: str,
    junction: InferredJunctionSet,
    max_reads_kmer: int,
) -> None:
    """Warn when multiple motifs were auto-detected on one read end."""
    if len(junction.motifs) <= 1:
        return
    logger.warning(
        "{} cutsite motif inference found multiple motifs {}. "
        "This can reflect low-quality data or too few sampled reads; consider increasing "
        "--max-reads-kmer from {}. It can also be valid multi-enzyme data such as 3RAD. "
        "In that case, reads will be demultiplexed using all detected motifs. "
        "You can enter cutsite motifs manually to suppress this warning.",
        read_label,
        junction.motifs,
        max_reads_kmer,
    )


def _append_motif_report_rows(
    rows: List[List[object]],
    *,
    read_end: str,
    role: str,
    source: str,
    decision: str,
    inferred: InferredJunctionSet | None,
) -> None:
    """Append motif report rows for one selected or detected motif set."""
    if inferred is None or not inferred.motifs:
        return
    for motif, count, fraction in zip(
        inferred.motifs,
        inferred.motif_counts,
        inferred.motif_support_fractions,
    ):
        position = (
            inferred.position_summary
            if inferred.position_mode == "barcode_boundary"
            else f"offset {inferred.offset}"
        )
        rows.append([read_end, role, source, decision, position, motif, count, fraction])


def format_preserved_file_preview(paths: List[Path]) -> str:
    """Return a short deterministic preview of preserved outdir files."""
    preview = [path.name for path in paths[:PRESERVED_WARNING_PREVIEW]]
    rendered = ", ".join(preview)
    if len(paths) > PRESERVED_WARNING_PREVIEW:
        rendered += f", ... and {len(paths) - PRESERVED_WARNING_PREVIEW} more"
    return rendered


def next_demux_stats_path(outdir: Path) -> Path:
    """Return the next available numbered demux stats path in outdir."""
    idx = 0
    while True:
        outfile = outdir / f"{DEMUX_STATS_PREFIX}{idx}.txt"
        if outfile.exists():
            idx += 1
            continue
        return outfile


def aggregate_file_stat_counter(
    file_stats: Dict[str, Tuple],
    stat_idx: int,
) -> Counter:
    """Aggregate one counter-like column across all raw file stats."""
    aggregate = Counter()
    for stats in file_stats.values():
        aggregate.update(stats[stat_idx])
    return aggregate


def aggregate_suspected_barcode_stats(
    file_stats: Dict[str, Tuple],
) -> Dict[Tuple[str, bytes], Tuple[int, int]]:
    """Aggregate bounded suspected-barcode summaries across raw files."""
    aggregate: Dict[Tuple[str, bytes], Tuple[int, int]] = {}
    for stats in file_stats.values():
        if len(stats) < 5:
            continue
        for key, (estimate, error) in stats[4].items():
            current_estimate, current_error = aggregate.get(key, (0, 0))
            aggregate[key] = (
                current_estimate + int(estimate),
                current_error + int(error),
            )
    return aggregate


def total_raw_reads_from_file_stats(file_stats: Dict[str, Tuple]) -> int:
    """Return total raw reads represented by demux per-file stats."""
    total = 0
    for stats in file_stats.values():
        total += sum(stats[0].values())
        total += sum(stats[1].values())
        total += sum(stats[3].values())
    return total


def format_observed_barcode(barcode: bytes) -> str | Tuple[str, ...]:
    """Return a human-readable barcode observation for stats output."""
    if b"_" in barcode:
        return tuple(part.decode() for part in barcode.split(b"_"))
    return barcode.decode()


def _expected_barcodes_by_read_end(
    names_to_barcodes: Dict[str, Tuple[str, str]],
) -> Dict[str, Tuple[str, ...]]:
    """Return user-entered expected barcode strings grouped by read end."""
    expected = {"R1": set(), "R2": set(), "i7": set()}
    for barcode1, barcode2 in names_to_barcodes.values():
        if barcode1:
            expected["R1"].add(barcode1)
            expected["i7"].add(barcode1)
        if barcode2:
            expected["R2"].add(barcode2)
    return {
        read_end: tuple(sorted(values))
        for read_end, values in expected.items()
    }


def _hamming_distance(left: str, right: str) -> int:
    """Return Hamming distance for equal-length strings."""
    return sum(a != b for a, b in zip(left, right))


def _nearest_expected_barcode(
    observed: str,
    expected: Tuple[str, ...],
) -> Tuple[str, int | None, str]:
    """Return nearest expected barcode, mismatch count, and relationship label."""
    leading_deletion = [
        barcode
        for barcode in expected
        if len(barcode) == len(observed) + 1 and barcode[1:] == observed
    ]
    if leading_deletion:
        return sorted(leading_deletion)[0], 0, "leading_base_deletion"

    same_length = [barcode for barcode in expected if len(barcode) == len(observed)]
    if not same_length:
        return "", None, ""
    nearest = min(same_length, key=lambda barcode: (_hamming_distance(observed, barcode), barcode))
    return nearest, _hamming_distance(observed, nearest), "same_length"


def _suspected_barcode_report_rows(
    *,
    file_stats: Dict[str, Tuple],
    names_to_barcodes: Dict[str, Tuple[str, str]],
) -> List[List[object]]:
    """Return report rows for frequent suspected missing barcodes."""
    raw_reads = total_raw_reads_from_file_stats(file_stats)
    if raw_reads <= 0:
        return []

    threshold = max(
        SUSPECTED_BARCODE_MIN_COUNT,
        ceil(raw_reads * SUSPECTED_BARCODE_MIN_FRACTION),
    )
    expected_by_read_end = _expected_barcodes_by_read_end(names_to_barcodes)
    rows = []
    for (read_end, barcode), (estimate, error) in aggregate_suspected_barcode_stats(file_stats).items():
        min_count = max(0, estimate - error)
        if min_count < threshold:
            continue
        observed = barcode.decode(errors="replace")
        nearest, distance, relationship = _nearest_expected_barcode(
            observed,
            expected_by_read_end.get(read_end, ()),
        )
        rows.append(
            [
                read_end,
                observed,
                estimate,
                min_count,
                error,
                estimate / raw_reads,
                nearest,
                "" if distance is None else distance,
                relationship,
            ]
        )

    return sorted(
        rows,
        key=lambda row: (-int(row[3]), row[0], row[1]),
    )[:SUSPECTED_BARCODE_REPORT_LIMIT]


def write_demux_stats(
    *,
    outdir: Path,
    file_stats: Dict[str, Tuple],
    sample_stats: Dict[str, int],
    names_to_barcodes: Dict[str, Tuple[str, str]],
    barcodes_to_names: Dict[bytes, str],
    i7: bool,
    re1_source: str | None,
    re1_inference: InferredJunctionSet | None,
    re1_detected_inference: InferredJunctionSet | None,
    re1_motif_decision: str | None,
    re2_source: str | None,
    re2_inference: InferredJunctionSet | None,
    re2_detected_inference: InferredJunctionSet | None,
    re2_motif_decision: str | None,
    barcode_boundary_collisions: List[Dict[str, str]],
) -> Path:
    """Write the numbered demux stats report and return its path."""
    stats_file = next_demux_stats_path(outdir)
    with stats_file.open("w", encoding="utf-8") as outfile:
        outfile.write("# Raw file statistics\n######################\n")
        file_df = pd.DataFrame(
            index=sorted(file_stats),
            columns=["total_reads", "cut_found", "bar_matched", "bar_ambiguous"],
        )
        for key in sorted(file_stats):
            stats = file_stats[key]
            not_cut = sum(stats[0].values())
            matched = sum(stats[1].values())
            ambiguous = sum(stats[3].values())
            total = not_cut + matched + ambiguous
            file_df.loc[key, :] = total, matched + ambiguous, matched, ambiguous
        outfile.write(file_df.to_string() + "\n\n")

        outfile.write("# Sample demux statistics\n######################\n")
        sample_df = pd.DataFrame(
            index=sorted(sample_stats),
            columns=["reads_raw"],
            data=[sample_stats[i] for i in sorted(sample_stats)],
        )
        outfile.write(sample_df.to_string() + "\n\n")

        if not i7:
            outfile.write("# Restriction motif inference\n######################\n")
            motif_rows = []
            for read_end, source, selected, detected, decision in (
                (
                    "R1",
                    re1_source,
                    re1_inference,
                    re1_detected_inference,
                    re1_motif_decision,
                ),
                (
                    "R2",
                    re2_source,
                    re2_inference,
                    re2_detected_inference,
                    re2_motif_decision,
                ),
            ):
                decision = decision or ""
                if source == "manual":
                    _append_motif_report_rows(
                        motif_rows,
                        read_end=read_end,
                        role="detected",
                        source="auto",
                        decision=decision,
                        inferred=detected,
                    )
                _append_motif_report_rows(
                    motif_rows,
                    read_end=read_end,
                    role="selected",
                    source=source or "",
                    decision=decision,
                    inferred=selected,
                )
            if motif_rows:
                motif_df = pd.DataFrame(
                    motif_rows,
                    columns=[
                        "read_end",
                        "role",
                        "source",
                        "decision",
                        "position",
                        "motif",
                        "support",
                        "support_fraction",
                    ],
                )
                outfile.write(
                    motif_df.to_string(
                        index=False,
                        float_format=lambda value: f"{value:.6f}",
                    ) + "\n\n"
                )
            else:
                outfile.write("none\n\n")

            outfile.write("# Barcode boundary collisions\n######################\n")
            if barcode_boundary_collisions:
                collision_df = pd.DataFrame(barcode_boundary_collisions)
                outfile.write(collision_df.to_string(index=False) + "\n\n")
            else:
                outfile.write("none\n\n")

            outfile.write("# Barcode boundary ambiguity statistics\n######################\n")
            outfile.write(
                "Reads listed in this section matched multiple sample/barcode-boundary "
                "candidates and were not assigned or written to any sample output file.\n"
            )
            ambiguous_obs = aggregate_file_stat_counter(file_stats, 3)
            if ambiguous_obs:
                ambiguity_df = pd.DataFrame(
                    {
                        "observed_boundary_candidates": [
                            key.decode(errors="replace") for key in ambiguous_obs
                        ],
                        "N_records": [ambiguous_obs[key] for key in ambiguous_obs],
                    }
                )
                outfile.write(ambiguity_df.to_string(index=False) + "\n\n")
            else:
                outfile.write("none\n\n")

        outfile.write("# Suspected missing barcode statistics\n######################\n")
        outfile.write(
            "This section reports frequent unassigned barcode-like observations with "
            "bounded-memory estimated counts. min_records is a conservative lower bound.\n"
        )
        suspected_rows = _suspected_barcode_report_rows(
            file_stats=file_stats,
            names_to_barcodes=names_to_barcodes,
        )
        if suspected_rows:
            suspected_df = pd.DataFrame(
                suspected_rows,
                columns=[
                    "read_end",
                    "observed_barcode",
                    "estimated_records",
                    "min_records",
                    "max_error",
                    "fraction_raw_reads_est",
                    "nearest_expected_barcode",
                    "nearest_expected_mismatches",
                    "nearest_expected_relationship",
                ],
            )
            outfile.write(
                suspected_df.to_string(
                    index=False,
                    float_format=lambda value: f"{value:.6f}",
                ) + "\n\n"
            )
        else:
            outfile.write("none\n\n")

        outfile.write("# Barcode detection statistics\n######################\n")
        data = []
        bar_obs = aggregate_file_stat_counter(file_stats, 1)
        sorted_bar_obs = sorted(bar_obs, key=lambda x: bar_obs[x], reverse=True)

        for name in sorted(names_to_barcodes):
            truebar = names_to_barcodes[name]
            for foundbar in sorted_bar_obs:
                if name != barcodes_to_names[foundbar]:
                    continue
                count = bar_obs[foundbar]
                if count:
                    data.append([name, truebar, format_observed_barcode(foundbar), count])

        bad_bars = aggregate_file_stat_counter(file_stats, 0)
        bad_bar_obs = sorted(bad_bars, key=lambda x: bad_bars[x], reverse=True)
        for badbar in bad_bar_obs:
            data.append(["no_match", "", format_observed_barcode(badbar), bad_bars[badbar]])

        barcodes_df = pd.DataFrame(
            index=[i[0] for i in data],
            columns=["true_bar", "observed_bar", "N_records"],
            data=[i[1:] for i in data],
        )
        outfile.write(barcodes_df.to_string() + "\n")

    logger.info(f"demultiplexing statistics written to {stats_file}")
    return stats_file
