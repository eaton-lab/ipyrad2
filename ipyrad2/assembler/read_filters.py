#!/usr/bin/env python

"""Helpers to prepare assemble-time filtered BAMs."""

from __future__ import annotations

from pathlib import Path
import sys

from loguru import logger

from ..utils.exceptions import IPyradError
from ..utils.parallel import run_pipeline


BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")

_PRIMARY_MAPPED_EXCLUDE_FLAGS = 0x904
_PRIMARY_MAPPED_SINGLE_EXCLUDE_FLAGS = 0x905


def _count_bam_records(
    bam_file: Path,
    *,
    require_flags: int | None = None,
    exclude_flags: int | None = None,
) -> int:
    """Return one `samtools view -c` count for a BAM under the requested flag masks."""
    cmd = [BIN_SAM, "view", "-c"]
    if exclude_flags:
        cmd.extend(["-F", hex(int(exclude_flags))])
    if require_flags:
        cmd.extend(["-f", hex(int(require_flags))])
    cmd.append(str(bam_file))
    _, out, _ = run_pipeline([cmd])
    if isinstance(out, bytes):
        out = out.decode()
    text = out.strip() if isinstance(out, str) else str(out).strip()
    return int(text or "0")


def classify_bam_layout(bam_file: Path) -> str:
    """Return `paired` or `single` for one BAM, rejecting hybrid primary layouts."""
    paired_primary = _count_bam_records(
        bam_file,
        require_flags=0x1,
        exclude_flags=_PRIMARY_MAPPED_EXCLUDE_FLAGS,
    )
    single_primary = _count_bam_records(
        bam_file,
        exclude_flags=_PRIMARY_MAPPED_SINGLE_EXCLUDE_FLAGS,
    )
    if paired_primary and single_primary:
        raise IPyradError(
            "BAM contains mixed single-end and paired-end primary mapped reads; "
            "mixed layouts are supported across samples, not within one BAM: "
            f"{bam_file}"
        )
    if paired_primary:
        return "paired"
    if single_primary:
        return "single"

    # Preserve legacy behavior for pathological/empty inputs where no primary
    # mapped reads are present by falling back to the broader paired-read probe.
    any_paired = _count_bam_records(bam_file, require_flags=0x1)
    return "paired" if any_paired else "single"


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
