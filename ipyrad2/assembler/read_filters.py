#!/usr/bin/env python

"""Helpers to prepare assemble-time filtered BAMs."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
import json
import re
import sys

from loguru import logger

from ..utils.exceptions import IPyradError
from ..utils.parallel import run_pipeline
from ..utils.parallel import stream_pipeline_lines


BIN = Path(sys.prefix) / "bin"
BIN_SAM = str(BIN / "samtools")

_PRIMARY_MAPPED_EXCLUDE_FLAGS = 0x904
_LAYOUT_PROBE_READ_LIMIT = 1000


@dataclass(frozen=True)
class FilteredAnalysisBamResult:
    """Filtered assemble BAM plus per-sample read-count diagnostics."""

    bam_path: Path
    reads_before_filtering: int
    reads_after_filtering: int
    filter_diagnostics: dict = field(default_factory=dict)


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
    min_aligned_len: int | None,
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

    if min_aligned_len is not None:
        clauses.append(f"((qlen - sclen) >= {min_aligned_len})")

    if not clauses:
        return None
    return " && ".join(clauses)


def get_analysis_bam_path(tmpdir: Path, sname: str) -> Path:
    """Return the temp pre-paralog filtered BAM path for one sample."""
    return tmpdir / "analysis_bams" / f"{sname}.analysis.filtered.bam"


def get_calling_bam_path(tmpdir: Path, sname: str) -> Path:
    """Return the temp post-paralog calling BAM path used only for joint calling."""
    return tmpdir / "calling_bams" / f"{sname}.variant.filtered.bam"


def get_paralog_bam_path(tmpdir: Path, sname: str) -> Path:
    """Return the temp loci-restricted BAM path used for paralog scoring."""
    return tmpdir / "paralog_bams" / f"{sname}.paralog.filtered.bam"


def _count_bam_records(bam_file: Path, threads: int) -> int:
    """Return the number of alignment records in one BAM."""
    cmd = [
        BIN_SAM,
        "view",
        "-c",
        "-@",
        str(max(1, threads)),
        str(bam_file),
    ]
    _, out, _ = run_pipeline([cmd])
    return int(out.decode().strip() or "0")


def summarize_read_filter_diagnostics(
    lines,
    *,
    min_map_q: int,
    is_paired: bool,
    max_tlen: int | None,
    max_softclip: int | None,
    max_nm: int | None,
    min_aligned_len: int | None,
) -> dict:
    """Count diagnostic failures in one SAM stream, without selecting reads.

    Independent counts may overlap. First-failure counts use the declared order
    and sum to the total rejected records. Query length includes insertions and
    soft clips, matching samtools qlen; aligned query length is qlen - sclen.
    The existing samtools command remains authoritative for filtering.
    """
    settings = {
        "min_map_q": min_map_q,
        "is_paired": is_paired,
        "max_tlen": max_tlen,
        "max_softclip": max_softclip,
        "max_nm": max_nm,
        "min_aligned_len": min_aligned_len,
    }
    order = ["mapping_quality"]
    if is_paired:
        order.append("pair_geometry")
        if max_tlen is not None:
            order.append("template_length")
    for key, value in (
        ("softclip", max_softclip),
        ("mismatch_count", max_nm),
        ("aligned_query_length", min_aligned_len),
    ):
        if value is not None:
            order.append(key)
    independent = dict.fromkeys(order, 0)
    first = dict.fromkeys(order, 0)
    total = accepted = missing_nm = 0
    flag_counts = Counter()
    for line in lines:
        if line.startswith("@") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        flag, mq, tlen = int(f[1]), int(f[4]), int(f[8])
        cigar = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5])]
        qlen = sum(n for n, op in cigar if op in "MIS=X")
        sclen = sum(n for n, op in cigar if op == "S")
        nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), None)
        missing_nm += nm is None
        failures = {
            "mapping_quality": mq < min_map_q,
            "pair_geometry": bool(flag & 0xC) or f[6] not in ("=", f[2]),
            "template_length": max_tlen is not None and abs(tlen) > max_tlen,
            "softclip": max_softclip is not None and sclen > max_softclip,
            "mismatch_count": max_nm is not None and (nm is None or nm > max_nm),
            "aligned_query_length": min_aligned_len is not None
            and qlen - sclen < min_aligned_len,
        }
        failed = [key for key in order if failures[key]]
        for key in failed:
            independent[key] += 1
        total += 1
        if failed:
            first[failed[0]] += 1
        else:
            accepted += 1
        for name, bit in (
            ("unmapped", 4),
            ("secondary", 0x100),
            ("qc_fail", 0x200),
            ("duplicate", 0x400),
            ("supplementary", 0x800),
        ):
            flag_counts[name] += bool(flag & bit)
    return {
        "unit": "alignment_records",
        "settings": settings,
        "filter_order": order,
        "input_records": total,
        "predicted_retained_records": accepted,
        "rejected_records": total - accepted,
        "independent_failures": independent,
        "first_failure": first,
        "records_without_nm": missing_nm,
        "input_flag_counts": dict(flag_counts),
        "semantics": "Independent failures overlap; first_failure is exclusive in filter_order. "
        "Flags are descriptive, not additional filters. Aligned query length includes insertions. "
        "Missing NM fails an enabled max_nm. Counts precede overlap quality masking and paralog/locus filtering.",
    }


def collect_read_filter_diagnostics(
    bam_file: Path, *, threads: int = 1, **settings
) -> dict:
    """Replace the input count scan with a single bounded-memory SAM scan."""
    cmd = [BIN_SAM, "view", "-@", str(max(1, threads)), str(bam_file)]
    with stream_pipeline_lines([cmd]) as lines:
        return summarize_read_filter_diagnostics(lines, **settings)


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
    min_aligned_len: int | None,
    threads: int,
) -> FilteredAnalysisBamResult:
    """Write and index one pre-paralog assemble BAM plus its filter counts."""
    out_bam = get_analysis_bam_path(tmpdir, sname)
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    diagnostics = collect_read_filter_diagnostics(
        bam_file, threads=threads, min_map_q=min_map_q, is_paired=is_paired, max_tlen=max_tlen,
        max_softclip=max_softclip, max_nm=max_nm, min_aligned_len=min_aligned_len,
    )
    reads_before_filtering = diagnostics["input_records"]

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
        min_aligned_len=min_aligned_len,
    )
    if expr is not None:
        cmd.extend(["-e", expr])
    if is_paired:
        # HTSlib does not correct overlaps for non-PROPER_PAIR reads even with
        # --count-orphans. Correct only those pairs, retaining original flags.
        cmd.extend(["-u", str(bam_file)])
        run_pipeline([
            cmd,
            [BIN_SAM, "collate", "-u", "-O", "-T", str(out_bam.with_suffix(".collate")), "-"],
            [BIN_SAM, "view", "-h", "-"],
            [sys.executable, "-m", "ipyrad2.assembler.pair_overlaps",
             "--stats", str(out_bam.with_suffix(".overlaps.json"))],
            [BIN_SAM, "sort", "-@", str(max(1, threads)), "-o", str(out_bam), "-"],
        ])
    else:
        cmd.extend(["-o", str(out_bam), str(bam_file)])
        run_pipeline([cmd])

    cmd = [
        BIN_SAM, "index",
        "-c",
        "-@", str(max(1, threads)),
        str(out_bam),
    ]
    run_pipeline([cmd])
    reads_after_filtering = _count_bam_records(out_bam, threads)
    diagnostics["observed_retained_records"] = reads_after_filtering
    diagnostics["counts_reconcile"] = (
        diagnostics["predicted_retained_records"] == reads_after_filtering
    )
    if not diagnostics["counts_reconcile"]:
        logger.warning("read-filter diagnostic counts differ from samtools output for {}", sname)
    out_bam.with_suffix(".filters.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    return FilteredAnalysisBamResult(
        bam_path=out_bam,
        reads_before_filtering=reads_before_filtering,
        reads_after_filtering=reads_after_filtering,
        filter_diagnostics=diagnostics,
    )


def prepare_variant_call_bam(
    *,
    sname: str,
    bam_file: Path,
    keep_bed: Path,
    tmpdir: Path,
    threads: int,
) -> Path:
    """Write and index one post-paralog calling BAM used only for joint calling."""
    if not keep_bed.exists():
        raise IPyradError(f"Retained per-sample loci BED not found for {sname}: {keep_bed}")

    out_bam = get_calling_bam_path(tmpdir, sname)
    out_bam.parent.mkdir(parents=True, exist_ok=True)

    if keep_bed.stat().st_size == 0:
        # Preserve the sample in the joint call with a valid header-only BAM
        # when no loci survive per-sample paralog filtering.
        cmd = [
            BIN_SAM,
            "view",
            "-b",
            "-H",
            "-@",
            str(max(1, threads)),
            "-o",
            str(out_bam),
            str(bam_file),
        ]
    else:
        cmd = [
            BIN_SAM,
            "view",
            "-b",
            "-h",
            "-@",
            str(max(1, threads)),
            "-L",
            str(keep_bed),
            "-o",
            str(out_bam),
            str(bam_file),
        ]
    run_pipeline([cmd])

    cmd = [
        BIN_SAM,
        "index",
        "-c",
        "-@",
        str(max(1, threads)),
        str(out_bam),
    ]
    run_pipeline([cmd])
    return out_bam


def prepare_paralog_bam(
    *,
    sname: str,
    bam_file: Path,
    regions_bed: Path,
    tmpdir: Path,
    threads: int,
) -> Path:
    """Write and index one loci-restricted BAM used only for paralog scoring."""
    if not regions_bed.exists():
        raise IPyradError(
            f"Shared loci BED not found for paralog BAM prep on {sname}: {regions_bed}"
        )

    out_bam = get_paralog_bam_path(tmpdir, sname)
    out_bam.parent.mkdir(parents=True, exist_ok=True)

    if regions_bed.stat().st_size == 0:
        cmd = [
            BIN_SAM,
            "view",
            "-b",
            "-H",
            "-@",
            str(max(1, threads)),
            "-o",
            str(out_bam),
            str(bam_file),
        ]
    else:
        cmd = [
            BIN_SAM,
            "view",
            "-b",
            "-h",
            "-@",
            str(max(1, threads)),
            "-L",
            str(regions_bed),
            "-o",
            str(out_bam),
            str(bam_file),
        ]
    run_pipeline([cmd])

    cmd = [
        BIN_SAM,
        "index",
        "-c",
        "-@",
        str(max(1, threads)),
        str(out_bam),
    ]
    run_pipeline([cmd])
    return out_bam
