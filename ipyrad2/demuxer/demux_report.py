#!/usr/bin/env python

"""Reporting helpers for demux logging and stats output."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

from loguru import logger
import pandas as pd

from ipyrad2.utils.kmers import InferredJunctionSet


DEMUX_STATS_PREFIX = "ipyrad_demux_stats_"
PRESERVED_WARNING_PREVIEW = 3


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
    file_stats: Dict[str, Tuple[Dict[str, int], Dict[bytes, int], Dict[str, int], Dict[bytes, int]]],
    stat_idx: int,
) -> Counter:
    """Aggregate one counter-like column across all raw file stats."""
    aggregate = Counter()
    for stats in file_stats.values():
        aggregate.update(stats[stat_idx])
    return aggregate


def format_observed_barcode(barcode: bytes) -> str | Tuple[str, ...]:
    """Return a human-readable barcode observation for stats output."""
    if b"_" in barcode:
        return tuple(part.decode() for part in barcode.split(b"_"))
    return barcode.decode()


def write_demux_stats(
    *,
    outdir: Path,
    file_stats: Dict[str, Tuple[Dict[str, int], Dict[bytes, int], Dict[str, int], Dict[bytes, int]]],
    sample_stats: Dict[str, int],
    names_to_barcodes: Dict[str, Tuple[str, str]],
    barcodes_to_names: Dict[bytes, str],
    i7: bool,
    re1_source: str | None,
    re1_inference: InferredJunctionSet | None,
    re2_source: str | None,
    re2_inference: InferredJunctionSet | None,
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
            for read_end, source, inferred in (
                ("R1", re1_source, re1_inference),
                ("R2", re2_source, re2_inference),
            ):
                if inferred is None or not inferred.motifs:
                    continue
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
                    motif_rows.append([read_end, source or "", position, motif, count, fraction])
            if motif_rows:
                motif_df = pd.DataFrame(
                    motif_rows,
                    columns=[
                        "read_end",
                        "source",
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
