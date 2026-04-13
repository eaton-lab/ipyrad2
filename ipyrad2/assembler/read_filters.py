#!/usr/bin/env python

"""Helpers to prepare assemble-time filtered BAMs."""

from __future__ import annotations

from pathlib import Path
import sys

from loguru import logger

from ..utils.exceptions import IPyradError
from ..utils.parallel import run_pipeline
from ..utils.parallel import stream_pipeline_lines


BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")

_PRIMARY_MAPPED_EXCLUDE_FLAGS = 0x904
_LAYOUT_PROBE_READ_LIMIT = 1000


def _iter_bam_view_lines(
    bam_file: Path,
    *,
    exclude_flags: int | None = None,
):
    """Yield `samtools view` output lines for one BAM under optional exclusions."""
    cmd = [BIN_SAM, "view"]
    if exclude_flags:
        cmd.extend(["-F", hex(int(exclude_flags))])
    cmd.append(str(bam_file))
    return stream_pipeline_lines([cmd])


def _parse_sam_flag(line: str, bam_file: Path) -> int:
    """Return the FLAG column from one SAM alignment line."""
    fields = line.split("\t", 2)
    if len(fields) < 2:
        raise IPyradError(f"Could not parse BAM record while probing layout: {bam_file}")
    try:
        return int(fields[1])
    except ValueError as exc:
        raise IPyradError(f"Could not parse BAM FLAG while probing layout: {bam_file}") from exc


def _sample_primary_mapped_layout(bam_file: Path) -> str | None:
    """Classify layout from the first sampled primary mapped reads, if any."""
    saw_paired = False
    saw_single = False
    sampled = 0
    with _iter_bam_view_lines(
        bam_file,
        exclude_flags=_PRIMARY_MAPPED_EXCLUDE_FLAGS,
    ) as lines:
        for line in lines:
            flag = _parse_sam_flag(line, bam_file)
            if flag & 0x1:
                saw_paired = True
            else:
                saw_single = True
            sampled += 1
            if saw_paired and saw_single:
                raise IPyradError(
                    "BAM contains mixed single-end and paired-end primary mapped reads; "
                    "mixed layouts are supported across samples, not within one BAM: "
                    f"{bam_file}"
                )
            if sampled >= _LAYOUT_PROBE_READ_LIMIT:
                break

    if sampled == 0:
        return None
    return "paired" if saw_paired else "single"


def _sample_any_paired_record(bam_file: Path) -> bool:
    """Return True when any sampled alignment record advertises paired layout."""
    sampled = 0
    with _iter_bam_view_lines(bam_file) as lines:
        for line in lines:
            if _parse_sam_flag(line, bam_file) & 0x1:
                return True
            sampled += 1
            if sampled >= _LAYOUT_PROBE_READ_LIMIT:
                break
    return False


def classify_bam_layout(bam_file: Path) -> str:
    """Return `paired` or `single` for one BAM from a sampled read probe."""
    primary_layout = _sample_primary_mapped_layout(bam_file)
    if primary_layout is not None:
        return primary_layout

    # Preserve legacy behavior for pathological/empty inputs where no primary
    # mapped reads are present by falling back to the broader paired-read probe.
    return "paired" if _sample_any_paired_record(bam_file) else "single"


def bam_appears_paired(bam_file: Path) -> bool:
    """Return True when the BAM appears to contain paired reads."""
    return classify_bam_layout(bam_file) == "paired"


def build_mapped_read_filter_expr(
    *,
    is_paired: bool,
    max_tlen: int | None,
    max_softclip: int | None,
    max_nm: int | None,
) -> str | None:
    """Return one samtools view -e expression for assemble-time BAM filtering."""
    clauses: list[str] = []

    if is_paired:
        pair_terms = [
            "((flag&4)==0)",
            "((flag&8)==0)",
            '(rnext=="=" || rnext==rname)',
        ]
        if max_tlen is not None:
            pair_terms.append(f"(tlen>={-max_tlen} && tlen<={max_tlen})")
        clauses.append("(" + " && ".join(pair_terms) + ")")

    if max_softclip is not None:
        clauses.append(f"(sclen <= {max_softclip})")

    if max_nm is not None:
        clauses.append(f"([NM] <= {max_nm})")

    if not clauses:
        return None
    return " && ".join(clauses)


def get_analysis_bam_path(tmpdir: Path, sname: str) -> Path:
    """Return the temp filtered BAM path for one sample."""
    return tmpdir / "analysis_bams" / f"{sname}.analysis.filtered.bam"


def prepare_filtered_analysis_bam(
    *,
    sname: str,
    bam_file: Path,
    is_paired: bool,
    tmpdir: Path,
    min_map_q: int,
    max_tlen: int | None,
    max_softclip: int | None,
    max_nm: int | None,
    threads: int,
) -> Path:
    """Write and index one assemble-time filtered BAM, then return its path."""
    out_bam = get_analysis_bam_path(tmpdir, sname)
    out_bam.parent.mkdir(parents=True, exist_ok=True)

    if not is_paired and max_tlen is not None:
        logger.debug(
            "ignoring paired-read assemble filters for single-end BAM: {}",
            bam_file.name,
        )

    cmd = [
        BIN_SAM, "view",
        "-b",
        "-@", str(max(1, threads)),
        "-q", str(min_map_q),
    ]
    expr = build_mapped_read_filter_expr(
        is_paired=is_paired,
        max_tlen=max_tlen,
        max_softclip=max_softclip,
        max_nm=max_nm,
    )
    if expr is not None:
        cmd.extend(["-e", expr])
    cmd.extend([
        "-o", str(out_bam),
        str(bam_file),
    ])
    run_pipeline([cmd])

    cmd = [
        BIN_SAM, "index",
        "-c",
        "-@", str(max(1, threads)),
        str(out_bam),
    ]
    run_pipeline([cmd])
    return out_bam
