#!/usr/bin/env python

"""Orchestrate the active end-to-end `ipyrad2 assemble` workflow."""

from typing import List
import shutil
from pathlib import Path
from loguru import logger
import pandas as pd
from .beds import (
    get_name_from_bam,
    get_names_from_bams,
    get_reference_sort_order,
    get_coverage_bed_graphs,
    get_across_sample_loci_bed,
    get_sample_depth_stats_in_final_loci,
    sort_bed_by_reference_order,
    write_callable_regions_bed,
)
from .loci import (
    write_sam_faidx,
    get_reference_in_loci_beds,
    make_lowdepth_mask,
    make_paralog_mask,
    merge_final_vcf_mask_beds,
    merge_sample_mask_beds,
    get_consensus,
    build_locus_fasta_database,
    write_final_outputs,
    write_assemble_stats_report,
)
from .paralogs import (
    aggregate_across_samples,
    get_sample_paralog_tables,
    write_per_sample_final_good,
)
from .read_filters import BIN_SAM, classify_bam_layout, prepare_filtered_analysis_bam
from .variants import (
    get_chunked_loci_beds,
    get_group_called_variants_in_vcf_chunks,
    get_concat_chunk_vcfs,
    get_filtered_vcf,
    apply_wgs_het_allele_balance_mask,
    write_variant_postfilter_stats,
    load_variant_postfilter_stats,
    summarize_variant_support_by_sample_type,
    get_vcf_with_indels_resolved,
    compact_resolved_vcf_to_final_loci_contigs,
    load_variant_resolution_stats,
    write_vcf,
)
from .write_snps import write_snps_hdf5
from ..utils.parallel import run_pipeline, run_with_pool, run_with_pool_iter
from ..utils.exceptions import IPyradError
from ..utils.pops import expand_imap_patterns, parse_imap, parse_pops_file
from ..utils.profiling import profile_stage


def existing_results_force_or_raise(outdir, tmpdir, name, force):
    """Apply assemble overwrite policy for the current output prefix."""
    if (outdir / f"{name}.loci.gz").exists() or tmpdir.exists():
        if not force:
            raise IPyradError(
                f"outfiles with prefix {name} already exist in {outdir}. Use --force to overwrite."
            )
        else:
            # collect relevant files and rm
            logger.debug(f"removing previous ipyrad assemble files from {outdir}")
            if tmpdir.exists():
                shutil.rmtree(tmpdir)
            rfiles = [
                outdir / f"{name}.loci.txt",
                outdir / f"{name}.loci.gz",
                outdir / f"{name}.bed",
                outdir / f"{name}.vcf.gz",
                outdir / f"{name}.vcf.gz.csi",
                outdir / f"{name}.hdf5",
                outdir / f"{name}.stats.txt",
                outdir / f"{name}.stats_counts.tsv",
                outdir / f"{name}.stats_sample_cov.txt",
                outdir / f"{name}.stats_locus_coverage.txt",
            ]
            for r in rfiles:
                if r.exists():
                    r.unlink()


def _log_mapped_read_filter_settings(
    *,
    min_map_q: int,
    max_tlen: int | None,
    max_softclip: int | None,
    max_nm: int | None,
    min_aligned_len: int | None,
) -> None:
    """Log the assemble-time mapped-read filter settings."""
    parts = [f"MAPQ>={min_map_q}", "same scaffold pairs only"]
    if max_tlen is not None:
        parts.append(f"abs(TLEN)<={max_tlen}")
    if max_softclip is not None:
        parts.append(f"softclip<={max_softclip}")
    if max_nm is not None:
        parts.append(f"NM<={max_nm}")
    if min_aligned_len is not None:
        parts.append(f"aligned_len>={min_aligned_len}")
    logger.info("filtering mapped reads before assembly: {}", ", ".join(parts))


def _count_nonempty_lines(path: Path) -> int:
    """Return the number of non-empty lines in a text file."""
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def _get_mixed_paralog_summary_path(phase_dir: Path) -> Path:
    """Return the mixed-mode per-locus paralog QC summary path."""
    return phase_dir / "paralogs.mixed_summary.tsv"


def _get_mixed_paralog_counts_path(phase_dir: Path) -> Path:
    """Return the mixed-mode aggregate paralog count summary path."""
    return phase_dir / "paralogs.mixed.counts.tsv"


def _write_mixed_paralog_counts(phase_dir: Path, counts: dict[str, int]) -> Path:
    """Persist compact mixed-mode paralog counts for the final stats report."""
    path = _get_mixed_paralog_counts_path(phase_dir)
    lines = [f"{key}\t{int(value)}" for key, value in counts.items()]
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return path


def _load_mixed_paralog_counts(tmpdir: Path) -> dict[str, int]:
    """Load mixed-mode aggregate paralog counts when present."""
    defaults = {
        "loci_fail_paralog_rad": 0,
        "loci_fail_paralog_wgs": 0,
        "loci_fail_paralog_both": 0,
        "loci_pass_paralog_rad_fail_paralog_wgs": 0,
    }
    path = _get_mixed_paralog_counts_path(tmpdir / "phase")
    if not path.exists():
        return defaults
    stats = dict(defaults)
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if not raw_line.strip():
            continue
        key, value = raw_line.split("\t", 1)
        if key in stats:
            stats[key] = int(value)
    return stats


def _write_mixed_paralog_summary(
    *,
    phase_dir: Path,
    rad_metrics: pd.DataFrame,
    wgs_metrics: pd.DataFrame,
) -> dict[str, int]:
    """Write one mixed-mode per-locus QC table and return aggregate counts."""
    rad_cols = [
        "rid",
        "chrom",
        "start",
        "end",
        "n_data",
        "n_good",
        "n_fail",
        "fail_frac_among_data",
        "good_frac_among_data",
        "drop_global",
        "keep_global",
    ]
    wgs_cols = list(rad_cols)

    rad = rad_metrics[[col for col in rad_cols if col in rad_metrics.columns]].copy()
    wgs = wgs_metrics[[col for col in wgs_cols if col in wgs_metrics.columns]].copy()

    if "rid" not in rad.columns or "rid" not in wgs.columns:
        raise IPyradError(
            "Mixed RAD/WGS paralog summaries require per-locus rid values."
        )

    rad = rad.rename(
        columns={
            "n_data": "rad_n_data",
            "n_good": "rad_n_good",
            "n_fail": "rad_n_fail",
            "fail_frac_among_data": "rad_fail_frac_among_data",
            "good_frac_among_data": "rad_good_frac_among_data",
            "drop_global": "rad_drop_global",
            "keep_global": "rad_keep_global",
        }
    )
    wgs = wgs.rename(
        columns={
            "chrom": "wgs_chrom",
            "start": "wgs_start",
            "end": "wgs_end",
            "n_data": "wgs_n_data",
            "n_good": "wgs_n_good",
            "n_fail": "wgs_n_fail",
            "fail_frac_among_data": "wgs_fail_frac_among_data",
            "good_frac_among_data": "wgs_good_frac_among_data",
            "drop_global": "wgs_drop_global",
            "keep_global": "wgs_keep_global",
        }
    )

    merged = rad.merge(wgs, on="rid", how="outer")
    if "chrom" not in merged.columns and "wgs_chrom" in merged.columns:
        merged["chrom"] = merged["wgs_chrom"]
    if "start" not in merged.columns and "wgs_start" in merged.columns:
        merged["start"] = merged["wgs_start"]
    if "end" not in merged.columns and "wgs_end" in merged.columns:
        merged["end"] = merged["wgs_end"]

    fill_int = [
        "rad_n_data",
        "rad_n_good",
        "rad_n_fail",
        "wgs_n_data",
        "wgs_n_good",
        "wgs_n_fail",
    ]
    fill_float = [
        "rad_fail_frac_among_data",
        "rad_good_frac_among_data",
        "wgs_fail_frac_among_data",
        "wgs_good_frac_among_data",
    ]
    fill_bool = [
        "rad_drop_global",
        "rad_keep_global",
        "wgs_drop_global",
        "wgs_keep_global",
    ]
    for col in fill_int:
        if col in merged.columns:
            merged[col] = merged[col].fillna(0).astype("int64")
    for col in fill_float:
        if col in merged.columns:
            merged[col] = merged[col].fillna(0.0).astype(float)
    for col in fill_bool:
        if col in merged.columns:
            merged[col] = merged[col].fillna(False).astype(bool)

    merged["controlling_keep_global"] = merged.get("rad_keep_global", False)
    merged["controlling_drop_global"] = merged.get("rad_drop_global", False)
    merged = merged.sort_values(["chrom", "start", "end", "rid"], kind="stable")
    summary_path = _get_mixed_paralog_summary_path(phase_dir)
    merged.to_csv(summary_path, sep="\t", index=False)

    counts = {
        "loci_fail_paralog_rad": int(merged["rad_drop_global"].sum()),
        "loci_fail_paralog_wgs": int(merged["wgs_drop_global"].sum()),
        "loci_fail_paralog_both": int(
            (merged["rad_drop_global"] & merged["wgs_drop_global"]).sum()
        ),
        "loci_pass_paralog_rad_fail_paralog_wgs": int(
            (merged["rad_keep_global"] & merged["wgs_drop_global"]).sum()
        ),
    }
    _write_mixed_paralog_counts(phase_dir, counts)
    logger.info("mixed RAD/WGS paralog QC summary written to {}", summary_path)
    return counts


def _load_reference_scaffold_order(tmpdir: Path) -> dict[str, int]:
    """Return reference scaffold order from the assemble REF_info.txt file."""
    ref_info = tmpdir / "REF_info.txt"
    order: dict[str, int] = {}
    with ref_info.open("r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if not parts[0]:
                continue
            order[parts[0]] = idx
    if not order:
        raise IPyradError(f"Reference scaffold order file is empty: {ref_info}")
    return order


def _load_reference_scaffold_lengths(tmpdir: Path) -> dict[str, int]:
    """Return reference scaffold lengths from the assemble REF_info.txt file."""
    ref_info = tmpdir / "REF_info.txt"
    lengths: dict[str, int] = {}
    with ref_info.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2 or not parts[0]:
                continue
            try:
                lengths[parts[0]] = int(parts[1])
            except ValueError as exc:
                raise IPyradError(
                    f"Reference scaffold length file is malformed: {ref_info}"
                ) from exc
    if not lengths:
        raise IPyradError(f"Reference scaffold length file is empty: {ref_info}")
    return lengths


def _load_reference_scaffold_records(tmpdir: Path) -> list[tuple[str, int]]:
    """Return ordered `(scaffold, length)` records from the assemble REF_info.txt file."""
    ref_info = tmpdir / "REF_info.txt"
    records: list[tuple[str, int]] = []
    with ref_info.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2 or not parts[0]:
                continue
            try:
                records.append((parts[0], int(parts[1])))
            except ValueError as exc:
                raise IPyradError(
                    f"Reference scaffold length file is malformed: {ref_info}"
                ) from exc
    if not records:
        raise IPyradError(f"Reference scaffold length file is empty: {ref_info}")
    return records


def _get_bam_header_reference_records(
    bam_file: Path,
) -> list[tuple[str, int]]:
    """Return ordered `(contig, length)` records from BAM `@SQ` header lines."""
    cmd = [BIN_SAM, "view", "-H", str(bam_file)]
    _, out, _ = run_pipeline([cmd])
    text = out.decode() if isinstance(out, bytes) else str(out)
    records: list[tuple[str, int]] = []
    for line_no, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line or not line.startswith("@SQ"):
            continue
        contig = None
        contig_len = None
        for field in line.split("\t")[1:]:
            if field.startswith("SN:"):
                contig = field[3:]
            elif field.startswith("LN:"):
                contig_len = field[3:]
        if not contig or contig_len is None:
            raise IPyradError(
                f"BAM header @SQ line is malformed for {bam_file} on line {line_no}: {raw_line}"
            )
        try:
            records.append((contig, int(contig_len)))
        except ValueError as exc:
            raise IPyradError(
                f"BAM header @SQ line is malformed for {bam_file} on line {line_no}: {raw_line}"
            ) from exc
    if not records:
        raise IPyradError(f"BAM header contains no @SQ records: {bam_file}")
    return records


def _describe_bam_reference_mismatch(
    bam_records: list[tuple[str, int]],
    reference_records: list[tuple[str, int]],
) -> str | None:
    """Return one compact description of the first BAM-header/reference mismatch."""
    if len(bam_records) != len(reference_records):
        return (
            f"BAM header has {len(bam_records)} contigs, reference has "
            f"{len(reference_records)}"
        )
    for bam_record, reference_record in zip(bam_records, reference_records):
        bam_name, bam_len = bam_record
        ref_name, ref_len = reference_record
        if bam_name != ref_name:
            return f"first differing @SQ contig is BAM {bam_name}, reference {ref_name}"
        if bam_len != ref_len:
            return (
                f"first differing @SQ length is {bam_name} "
                f"(BAM {bam_len}, reference {ref_len})"
            )
    return None


def _validate_bam_header_records_match_reference(
    sample_bam_records: dict[str, list[tuple[str, int]]],
    tmpdir: Path,
    reference: Path,
) -> None:
    """Reject BAM headers whose reference dictionary differs from `-r`."""
    reference_records = _load_reference_scaffold_records(tmpdir)
    affected_samples: list[str] = []
    for sname, bam_records in sample_bam_records.items():
        mismatch = _describe_bam_reference_mismatch(bam_records, reference_records)
        if mismatch is not None:
            affected_samples.append(f"{sname}: {mismatch}")

    if not affected_samples:
        return

    sample_summaries = " | ".join(affected_samples[:5])
    if len(affected_samples) > 5:
        sample_summaries += f" | ... ({len(affected_samples) - 5} more samples)"
    raise IPyradError(
        "BAM headers do not match the current reference passed "
        f"to -r ({reference}). These BAMs were mapped against a different reference "
        "dictionary. If you reused the same reference path during mapping, stale "
        "bwa-mem2 sidecar index files (.ann/.amb/.pac/.0123/.bwt.2bit.64) are a "
        "likely cause. Remap these BAMs against that exact reference before "
        f"running assemble. Affected samples: {sample_summaries}"
    )


def _validate_analysis_bams_match_reference(
    sample_bams: dict[str, Path],
    tmpdir: Path,
    reference: Path,
) -> None:
    """Reject BAMs whose header reference dictionary differs from `-r`."""
    sample_bam_records = {
        sname: _get_bam_header_reference_records(bam_file)
        for sname, bam_file in sample_bams.items()
    }
    _validate_bam_header_records_match_reference(
        sample_bam_records,
        tmpdir,
        reference,
    )


def _normalize_user_loci_bed(loci_bed: Path, tmpdir: Path) -> tuple[Path, int]:
    """Validate and normalize a user-provided loci BED into the assemble tmpdir."""
    if not loci_bed.exists():
        raise IPyradError(f"--loci-bed file not found: {loci_bed}")
    if not loci_bed.is_file():
        raise IPyradError(f"--loci-bed must point to a file: {loci_bed}")

    scaffold_order = _load_reference_scaffold_order(tmpdir)
    scaffold_lengths = _load_reference_scaffold_lengths(tmpdir)
    records: list[tuple[str, int, int, int]] = []
    ignored_extra_cols = False

    with loci_bed.open("r", encoding="utf-8") as handle:
        for lineno, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if (
                not line
                or line.startswith("#")
                or line.startswith("track ")
                or line.startswith("browser ")
            ):
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) < 3:
                raise IPyradError(
                    f"--loci-bed line {lineno} must have at least 3 tab-delimited columns."
                )
            chrom = parts[0]
            if chrom not in scaffold_order:
                raise IPyradError(
                    f"--loci-bed contains scaffold not present in reference: {chrom}"
                )
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError as exc:
                raise IPyradError(
                    f"--loci-bed line {lineno} has non-integer start/end coordinates."
                ) from exc
            if start < 0:
                raise IPyradError(f"--loci-bed line {lineno} has start < 0.")
            if end <= start:
                raise IPyradError(f"--loci-bed line {lineno} must satisfy end > start.")
            if end > scaffold_lengths[chrom]:
                raise IPyradError(
                    f"--loci-bed line {lineno} exceeds reference length for {chrom}: {end} > {scaffold_lengths[chrom]}"
                )
            if len(parts) > 3:
                ignored_extra_cols = True
            records.append((chrom, start, end, lineno))

    if not records:
        raise IPyradError(f"--loci-bed contains no loci: {loci_bed}")

    records.sort(key=lambda item: (scaffold_order[item[0]], item[1], item[2], item[3]))

    by_chrom_last_end: dict[str, int] = {}
    for chrom, start, end, _lineno in records:
        prev_end = by_chrom_last_end.get(chrom)
        if prev_end is not None and start < prev_end:
            raise IPyradError(f"--loci-bed contains overlapping intervals on {chrom}.")
        by_chrom_last_end[chrom] = end

    out_bed = tmpdir / "beds" / "loci.raw.bed"
    with out_bed.open("w", encoding="utf-8") as out:
        for chrom, start, end, _lineno in records:
            out.write(f"{chrom}\t{start}\t{end}\n")

    if ignored_extra_cols:
        logger.debug("ignored extra columns beyond BED3 while normalizing --loci-bed")
    return out_bed, len(records)


def _normalize_bam_rename_file(
    rename_bams: Path, bam_paths: list[Path]
) -> dict[str, str]:
    """Parse and validate explicit BAM-basename sample-name overrides."""
    rename_bams = rename_bams.expanduser().absolute()
    if not rename_bams.exists():
        raise IPyradError(f"--rename-bams file not found: {rename_bams}")
    if not rename_bams.is_file():
        raise IPyradError(f"--rename-bams must point to a file: {rename_bams}")

    basenames = [path.name for path in bam_paths]
    basename_counts: dict[str, int] = {}
    for name in basenames:
        basename_counts[name] = basename_counts.get(name, 0) + 1
    duplicate_inputs = sorted(
        name for name, count in basename_counts.items() if count > 1
    )
    if duplicate_inputs:
        raise IPyradError(
            "--rename-bams cannot be used when input BAM basenames are duplicated: "
            + ", ".join(duplicate_inputs)
        )

    rename_map: dict[str, str] = {}
    with rename_bams.open("r", encoding="utf-8") as handle:
        for lineno, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                raise IPyradError(
                    f"--rename-bams line {lineno} must contain exactly 2 columns."
                )
            bam_name, sample_name = parts
            if bam_name in rename_map:
                raise IPyradError(
                    f"--rename-bams assigns BAM basename multiple times: {bam_name}"
                )
            rename_map[bam_name] = sample_name

    if not rename_map:
        raise IPyradError(f"--rename-bams contains no rename mappings: {rename_bams}")

    extra = sorted(set(rename_map).difference(basenames))
    if extra:
        raise IPyradError(
            "--rename-bams contains BAM basenames not present in this assemble run: "
            + ", ".join(extra)
        )
    return rename_map


def _collect_named_bams(
    bam_paths: list[Path],
    rename_map: dict[str, str],
) -> tuple[dict[str, Path], list[tuple[str, str]]]:
    """Resolve final sample names for BAM inputs from header names plus overrides."""
    bam_dict: dict[str, Path] = {}
    renamed: list[tuple[str, str]] = []
    unresolved = [bam_file for bam_file in bam_paths if bam_file.name not in rename_map]
    if len(unresolved) == 1:
        unresolved_names = {unresolved[0]: get_name_from_bam(unresolved[0])}
    else:
        unresolved_names = get_names_from_bams(unresolved)
    for bam_file in bam_paths:
        sname = rename_map.get(bam_file.name)
        if sname is None:
            sname = unresolved_names[bam_file]
        else:
            renamed.append((sname, bam_file.name))
        if sname in bam_dict:
            raise IPyradError(f"Multiple input files of sample name {sname}")
        bam_dict[sname] = bam_file
    return bam_dict, renamed


def _probe_bam_metadata(bam_file: Path) -> dict[str, object]:
    """Return sampled layout and ordered `@SQ` records for one BAM."""
    return {
        "layout": classify_bam_layout(bam_file),
        "header_records": _get_bam_header_reference_records(bam_file),
    }


def _collect_bam_metadata(
    bam_dict: dict[str, Path],
    log_level: str,
    max_workers: int,
) -> dict[str, dict[str, object]]:
    """Collect startup BAM metadata in parallel without progress logging."""
    jobs_iter = (
        (sname, (_probe_bam_metadata, {"bam_file": bam_file}))
        for sname, bam_file in bam_dict.items()
    )
    return {
        sname: result
        for sname, result in run_with_pool_iter(
            jobs_iter,
            log_level,
            max_workers=max_workers,
            max_inflight=max_workers,
        )
    }


def _normalize_populations_file(
    populations: Path,
    tmpdir: Path,
    sample_names: list[str],
) -> tuple[Path, dict[str, list[str]], dict[str, int] | None]:
    """Validate and normalize grouped-calling sample assignments for assemble."""
    populations = populations.expanduser().absolute()
    if not populations.exists():
        raise IPyradError(f"--populations file not found: {populations}")
    if not populations.is_file():
        raise IPyradError(f"--populations must point to a file: {populations}")

    parsed_minmap: dict[str, int] | None = None
    try:
        imap, parsed_minmap = parse_pops_file(populations)
    except IPyradError as pops_exc:
        try:
            imap = parse_imap(populations)
        except IPyradError as imap_exc:
            raise pops_exc from imap_exc

    if not imap:
        raise IPyradError(
            f"--populations contains no sample assignments: {populations}"
        )
    imap, _unmatched = expand_imap_patterns(
        imap,
        sample_names,
        mapping_name="--populations",
        available_name="this assemble run",
    )

    sample_to_group: dict[str, str] = {}
    for group, names in imap.items():
        for name in names:
            sample_to_group[name] = group

    assembled = set(sample_names)
    assigned = set(sample_to_group)
    missing = sorted(assembled.difference(assigned))
    if missing:
        raise IPyradError(
            "--populations is missing assembled sample(s): " + ", ".join(missing)
        )
    out_path = tmpdir / "populations.normalized.tsv"
    with out_path.open("w", encoding="utf-8") as out:
        for sample in sample_names:
            out.write(f"{sample}\t{sample_to_group[sample]}\n")

    normalized_imap = {
        group: [sample for sample in sample_names if sample_to_group[sample] == group]
        for group in imap
    }
    return out_path, normalized_imap, parsed_minmap


def _prepare_analysis_bams(
    *,
    bam_dict: dict[str, Path],
    sample_layouts: dict[str, bool],
    tmpdir: Path,
    min_map_q: int,
    max_tlen: int | None,
    max_softclip: int | None,
    max_nm: int | None,
    min_aligned_len: int | None,
    threads: int,
    workers: int,
    log_level: str,
) -> dict[str, Path]:
    """Write temp filtered BAMs for all assemble inputs and return their paths."""
    _log_mapped_read_filter_settings(
        min_map_q=min_map_q,
        max_tlen=max_tlen,
        max_softclip=max_softclip,
        max_nm=max_nm,
        min_aligned_len=min_aligned_len,
    )
    jobs = {}
    for sname, bam_file in bam_dict.items():
        jobs[sname] = (
            prepare_filtered_analysis_bam,
            {
                # These BAMs are the single filtered inputs reused by all later
                # assemble stages, so the same mapped-read thresholds are applied
                # once here instead of being reimplemented downstream.
                "sname": sname,
                "bam_file": bam_file,
                "is_paired": sample_layouts[sname],
                "tmpdir": tmpdir,
                "min_map_q": min_map_q,
                "max_tlen": max_tlen,
                "max_softclip": max_softclip,
                "max_nm": max_nm,
                "min_aligned_len": min_aligned_len,
                "threads": threads,
            },
        )
    return run_with_pool(jobs, log_level, workers, msg="Filtering mapped reads")


def _run_paralog_stage(
    *,
    sample_bams: dict[str, Path],
    regions_bed: Path,
    reference: Path,
    bed_dir: Path,
    phase_dir: Path,
    min_map_q: int,
    min_base_q: int,
    softclip_len_threshold: int,
    softclip_frac_max: float,
    depth_z_max: float,
    third_frac_cut: float,
    min_3allele_sites: int,
    maf_threshold: float,
    max_sites_above_maf: int,
    paralog_fail_frac_max: float,
    workers: int,
    log_level: str,
    rad_sample_names: list[str] | None = None,
    wgs_sample_names: list[str] | None = None,
) -> Path:
    """Run the active per-sample and across-sample paralog filtering stage."""
    rad_sample_names = sorted(rad_sample_names or [])
    wgs_sample_names = sorted(wgs_sample_names or [])
    mixed_mode = bool(rad_sample_names and wgs_sample_names)
    wgs_sample_name_set = set(wgs_sample_names)

    logger.info("filtering paralogs within samples")
    callable_regions_bed = write_callable_regions_bed(
        regions_bed,
        reference,
        phase_dir / "loci.callable.paralog.bed",
    )
    logger.info(
        "excluding non-ACGT reference positions from within-sample paralog variant calling"
    )
    if mixed_mode and (softclip_len_threshold is not None or softclip_frac_max is not None):
        logger.info(
            "mixed RAD/WGS assembly detected; skipping softclip-based paralog failure for WGS samples"
        )

    # Score every sample against the shared RAD-defined loci BED. RAD samples
    # still define the loci, but WGS samples are also evaluated here so their
    # sample-specific masks and the shared drop decision use the same evidence.
    kwargs = dict(
        regions_bed=regions_bed,
        callable_regions_bed=callable_regions_bed,
        reference_fasta=reference,
        tmpdir=phase_dir,
        min_map_q=min_map_q,
        min_base_q=min_base_q,
        indel_pad_bp=10,
        min_allele_depth=2,
        max_abs_dp_z_max=depth_z_max,
        third_frac_cut=third_frac_cut,
        min_3allele_sites=min_3allele_sites,
        maf_threshold=maf_threshold,
        max_sites_above_maf=max_sites_above_maf,
        softclip_len_threshold=softclip_len_threshold,
        softclip_frac_max=softclip_frac_max,
    )
    jobs = {}
    for sname, bam_file in sample_bams.items():
        sample_kwargs = kwargs
        if mixed_mode and sname in wgs_sample_name_set:
            sample_kwargs = sample_kwargs | {
                "softclip_len_threshold": None,
                "softclip_frac_max": None,
            }
        ikwargs = sample_kwargs | dict(bam=bam_file, prefix=sname)
        jobs[sname] = (get_sample_paralog_tables, ikwargs)
    run_with_pool(jobs, log_level, workers, msg="Filtering paralogs within samples")

    logger.info("aggregating paralog filtering across samples")

    if mixed_mode:
        logger.info(
            "mixed RAD/WGS assembly detected; using RAD samples to control shared paralog locus retention and WGS samples for QC only"
        )
        rad_prefix = phase_dir / "paralogs.rad"
        rad_metrics = aggregate_across_samples(
            regions_bed=regions_bed,
            sample_prefixes=rad_sample_names,
            in_dir=phase_dir,
            out_prefix=rad_prefix,
            fail_frac_max=paralog_fail_frac_max,
            min_data_samples=1,
        )
        wgs_prefix = phase_dir / "paralogs.wgs"
        wgs_metrics = aggregate_across_samples(
            regions_bed=regions_bed,
            sample_prefixes=wgs_sample_names,
            in_dir=phase_dir,
            out_prefix=wgs_prefix,
            fail_frac_max=paralog_fail_frac_max,
            min_data_samples=1,
        )
        mixed_counts = _write_mixed_paralog_summary(
            phase_dir=phase_dir,
            rad_metrics=rad_metrics,
            wgs_metrics=wgs_metrics,
        )
        shared_good_bed = Path(f"{rad_prefix}.shared_good.final.bed")
        n_keep = int(rad_metrics["keep_global"].sum())
        n_total_loci = int(rad_metrics.shape[0])
        logger.info(
            "mixed paralog summary: RAD-fail={} WGS-fail={} both-fail={} RAD-pass/WGS-fail={}",
            mixed_counts["loci_fail_paralog_rad"],
            mixed_counts["loci_fail_paralog_wgs"],
            mixed_counts["loci_fail_paralog_both"],
            mixed_counts["loci_pass_paralog_rad_fail_paralog_wgs"],
        )
    else:
        # Then combine the per-sample calls while counting failures only among
        # samples that actually have read data in a locus. This avoids punishing
        # no-SNP / low-information loci that still had real coverage.
        shared_prefix = phase_dir / "paralogs"
        metrics = aggregate_across_samples(
            regions_bed=regions_bed,
            sample_prefixes=sorted(sample_bams),
            in_dir=phase_dir,
            out_prefix=shared_prefix,
            fail_frac_max=paralog_fail_frac_max,
            min_data_samples=1,
        )
        shared_good_bed = Path(f"{shared_prefix}.shared_good.final.bed")
        n_keep = int(metrics["keep_global"].sum())
        n_total_loci = int(metrics.shape[0])

    # Keep a named copy of the shared paralog-filtered BED for debugging while
    # also promoting it to the canonical loci.bed used by downstream stages.
    final_loci_bed = bed_dir / "loci.paralog_filtered.bed"
    ref_info = bed_dir.parent / "REF_info.txt"
    sort_bed_by_reference_order(shared_good_bed, final_loci_bed, ref_info)
    sort_bed_by_reference_order(shared_good_bed, bed_dir / "loci.bed", ref_info)

    # Also materialize sample-specific final BEDs that already respect the
    # global shared filter, so later sample-level masking can reuse them.
    write_per_sample_final_good(
        sample_prefixes=sorted(sample_bams),
        in_dir=phase_dir,
        shared_good_bed=shared_good_bed,
        out_dir=bed_dir,
    )

    logger.info("paralog filtering retained {}/{} shared loci", n_keep, n_total_loci)
    if (not final_loci_bed.exists()) or final_loci_bed.stat().st_size == 0:
        raise IPyradError("No loci passed paralog filtering.")
    return final_loci_bed


def _run_variant_stage(
    *,
    tmpdir: Path,
    reference: Path,
    bam_dict: dict[str, Path],
    group_samples_file: Path | None,
    min_map_q: int,
    min_base_q: int,
    min_sample_depth: int,
    min_geno_q: int,
    min_site_q: int,
    cores: int,
    threads: int,
    log_level: str,
    wgs_samples: list[str] | None = None,
) -> Path:
    """Call, filter, and resolve joint variants inside the canonical loci BED."""
    vcf_dir = tmpdir / "vcfs"
    vcf_dir.mkdir(parents=True, exist_ok=True)
    callable_loci_bed = write_callable_regions_bed(
        tmpdir / "beds" / "loci.bed",
        reference,
        tmpdir / "beds" / "loci.callable.variant.bed",
    )
    if _count_nonempty_lines(callable_loci_bed) == 0:
        raise IPyradError(
            "No callable A/C/G/T reference positions remain in the final loci BED after excluding non-ACGT bases."
        )

    # Keep the number of concurrent joint-calling jobs aligned with the normal
    # assemble worker budget. Each mpileup/call job is memory-heavy because it
    # scans all BAMs in one loci chunk, so oversubscribing here drives the
    # highest peak RAM in typical assemble runs.
    chunk_threads = max(1, min(2, threads))
    variant_workers = max(1, cores // max(1, threads))

    # Still create more chunks than inflight jobs so long loci sets distribute
    # better across workers without forcing more concurrent mpileup processes.
    chunk_count = max(8, 2 * variant_workers)
    chunk_beds = get_chunked_loci_beds(
        tmpdir, chunk_count, source_bed=callable_loci_bed
    )

    jobs = {}
    for chunk in chunk_beds:
        jobs[str(chunk)] = (
            get_group_called_variants_in_vcf_chunks,
            dict(
                tmpdir=tmpdir,
                reference=reference,
                bam_files=list(bam_dict.values()),
                group_samples_file=group_samples_file,
                min_base_q=min_base_q,
                min_map_q=min_map_q,
                locus_chunk=chunk,
                threads=chunk_threads,
            ),
        )
    run_with_pool(jobs, log_level, variant_workers, msg="Calling variants")

    # Collapse chunk-level calls back to one project VCF, then apply the shared
    # genotype/site filters before resolving SNP/indel conflicts.
    get_concat_chunk_vcfs(tmpdir, threads)
    logger.info("filtering variants")
    get_filtered_vcf(tmpdir, min_sample_depth, min_geno_q, min_site_q, threads)
    if wgs_samples:
        logger.info("masking WGS heterozygous genotypes by allele balance")
        ab_stats = apply_wgs_het_allele_balance_mask(
            tmpdir / "vcfs" / "loci.filtered.vcf.gz",
            wgs_samples,
            low=0.20,
            high=0.80,
        )
        write_variant_postfilter_stats(tmpdir, **ab_stats)
        logger.info(
            "masked {} / {} WGS heterozygous genotypes outside allele-balance range [0.20, 0.80]",
            ab_stats["wgs_het_genotypes_masked_by_allele_balance"],
            ab_stats["wgs_het_genotypes_examined_for_allele_balance"],
        )
    logger.info("resolving indels and snps")
    return get_vcf_with_indels_resolved(tmpdir, reference, threads)


def _build_sample_masks(
    *,
    tmpdir: Path,
    all_dict: dict[str, Path],
    min_sample_depth: int,
    workers: int,
    log_level: str,
) -> dict[str, Path]:
    """Build low-depth and sample-specific masks, then merge them per sample."""
    # Low-depth masks apply to every analysis BAM because both RAD and WGS
    # samples contribute consensus sequence and final genotype outputs.
    jobs = {}
    for sname in all_dict:
        jobs[sname] = (
            make_lowdepth_mask,
            dict(sname=sname, min_sample_depth=min_sample_depth, tmpdir=tmpdir),
        )
    run_with_pool(jobs, log_level, workers, msg="Building low-depth masks")

    # Sample-specific paralog masks now apply to every sample that was scored in
    # the shared RAD-defined loci, including WGS samples.
    paralog_masks: dict[str, Path] = {}
    if all_dict:
        jobs = {}
        for sname in all_dict:
            jobs[sname] = (make_paralog_mask, dict(sname=sname, tmpdir=tmpdir))
        paralog_masks = run_with_pool(
            jobs,
            log_level,
            workers,
            msg="Building sample-specific paralog masks",
        )

    # Merge low-depth, paralog, and any overlapping-indel-cluster exclusions
    # into the final mask BED that the consensus step uses to write Ns for
    # sample-specific filtered positions.
    jobs = {}
    for sname in all_dict:
        jobs[sname] = (merge_sample_mask_beds, dict(sname=sname, tmpdir=tmpdir))
    run_with_pool(jobs, log_level, workers, msg="Merging sample masks")
    return paralog_masks


def _write_consensus_and_outputs(
    *,
    name: str,
    outdir: Path,
    tmpdir: Path,
    snames: List[str],
    reference: Path,
    masks: List[str] | None,
    sample_masks: dict[str, Path] | None,
    shared_loci_after_delimiting: int,
    shared_loci_after_paralog_filtering: int,
    min_locus_sample_coverage: int,
    min_locus_trim_sample_coverage: int,
    min_locus_length: int,
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    max_sample_hetero_frequency: float,
    consensus_workers: int,
    final_vcf_mask_workers: int,
    workers: int,
    threads: int,
    log_level: str,
    cores: int | None = None,
    rad_samples: list[str] | None = None,
    wgs_samples: list[str] | None = None,
) -> None:
    """Write consensus sequences, final locus outputs, and the SNP database."""
    # Slice the reference to the final shared loci so every downstream output is
    # anchored to the same canonical set of assembled windows.
    logger.info("extracting reference sequence in locus beds")
    write_sam_faidx(tmpdir)
    reference_fasta = get_reference_in_loci_beds(tmpdir, reference)

    # Consensus calling applies the merged per-sample mask beds on top of the
    # resolved project VCF to create sample FASTAs for the final database.
    logger.info("extracting consensus sequences")
    jobs = {}
    for sname in snames:
        jobs[sname] = (
            get_consensus,
            dict(
                sname=sname,
                reference_fasta=reference_fasta,
                tmpdir=tmpdir,
                keep_insertions=False,
            ),
        )
    logger.debug(
        "consensus stage using up to {} workers across {} samples",
        consensus_workers,
        len(snames),
    )
    with profile_stage("consensus extraction"):
        run_with_pool(
            jobs, log_level, consensus_workers, msg="Extracting consensus sequences"
        )

    # Build one FASTA database spanning all consensus sequences. The final .loci
    # and HDF5 writers both consume this database to stay coordinate-consistent.
    logger.info("assembling loci")
    build_locus_fasta_database(name, snames, reference, tmpdir, masks)

    logger.info("writing outfiles (.loci, .hdf5, .bed, .stats.txt)")
    loci_summary = write_final_outputs(
        snames=snames,
        name=name,
        outdir=outdir,
        tmpdir=tmpdir,
        reference=reference,
        min_locus_sample_coverage=min_locus_sample_coverage,
        min_locus_trim_sample_coverage=min_locus_trim_sample_coverage,
        min_locus_length=min_locus_length,
        max_locus_hetero_frequency=max_locus_hetero_frequency,
        max_locus_variant_frequency=max_locus_variant_frequency,
        max_sample_hetero_frequency=max_sample_hetero_frequency,
        cores=max(1, int(cores or workers)),
        log_level=log_level,
    )
    if loci_summary["nloci_after_filtering"] == 0:
        raise IPyradError("No loci passed final trimming/filtering.")
    final_loci_bed = outdir / f"{name}.bed"

    # The final VCF is filtered to the trimmed/retained outdir BED, then the
    # SNP dataset is appended to the same output HDF5 for downstream analyses.
    logger.info("writing variants file (.vcf.gz)")
    with profile_stage("final VCF writing"):
        compact_resolved_vcf_to_final_loci_contigs(tmpdir, reference, final_loci_bed)
        final_vcf_masks = {}
        if snames:
            jobs = {
                sname: (merge_final_vcf_mask_beds, dict(sname=sname, tmpdir=tmpdir))
                for sname in snames
            }
            final_vcf_masks = run_with_pool(
                jobs,
                log_level,
                final_vcf_mask_workers,
                msg="Merging final VCF mask BEDs",
            )
        final_vcf = write_vcf(
            name,
            outdir,
            tmpdir,
            threads,
            sample_masks=final_vcf_masks,
            cores=cores,
            log_level=log_level,
        )

    mixed_run_summary: dict[str, int] | None = None
    rad_samples = sorted(rad_samples or [])
    wgs_samples = sorted(wgs_samples or [])
    if rad_samples and wgs_samples:
        mixed_run_summary = {
            "rad_samples": len(rad_samples),
            "wgs_samples": len(wgs_samples),
        }
        mixed_run_summary.update(_load_mixed_paralog_counts(tmpdir))
        mixed_run_summary.update(load_variant_postfilter_stats(tmpdir))
        support_stats = summarize_variant_support_by_sample_type(
            final_vcf, rad_samples, wgs_samples
        )
        mixed_run_summary.update(support_stats)
        logger.info(
            "mixed RAD/WGS summary: sites RAD-only={} WGS-only={} both={} WGS het masks={}",
            mixed_run_summary["sites_supported_rad_only"],
            mixed_run_summary["sites_supported_wgs_only"],
            mixed_run_summary["sites_supported_both"],
            mixed_run_summary["wgs_het_genotypes_masked_by_allele_balance"],
        )
        logger.info(
            "mixed RAD/WGS paralog summary: RAD-fail={} WGS-fail={} both-fail={} RAD-pass/WGS-fail={}",
            mixed_run_summary["loci_fail_paralog_rad"],
            mixed_run_summary["loci_fail_paralog_wgs"],
            mixed_run_summary["loci_fail_paralog_both"],
            mixed_run_summary["loci_pass_paralog_rad_fail_paralog_wgs"],
        )
        if mixed_run_summary.get("sites_supported_neither", 0):
            logger.debug(
                "mixed RAD/WGS final VCF retained {} site(s) with no ALT support after masking",
                mixed_run_summary["sites_supported_neither"],
            )

    logger.info("writing snps database (.hdf5)")
    with profile_stage("SNP HDF5 writing"):
        nsnps_written = write_snps_hdf5(name, outdir, snames, reference)

    logger.info("computing final sample depth summary")
    jobs = {}
    loci_bed = final_loci_bed
    for sname in snames:
        jobs[sname] = (
            get_sample_depth_stats_in_final_loci,
            dict(sname=sname, loci_bed=loci_bed, tmpdir=tmpdir),
        )
    sample_depth_stats = run_with_pool(
        jobs,
        log_level,
        workers,
        msg="Summarizing final sample depths",
    )

    # Write the final human-readable assemble report after all final outputs
    # exist, so its counts reflect the exact BED, VCF, and HDF5 products.
    write_assemble_stats_report(
        name=name,
        outdir=outdir,
        snames=snames,
        shared_loci_after_delimiting=shared_loci_after_delimiting,
        shared_loci_after_paralog_filtering=shared_loci_after_paralog_filtering,
        loci_summary=loci_summary,
        sample_depth_stats=sample_depth_stats,
        nsnps_written=nsnps_written,
        overlap_stats=load_variant_resolution_stats(tmpdir),
        mixed_run_summary=mixed_run_summary,
    )


def run_assembler(
    rad_bams: List[Path],
    wgs_bams: List[Path] | None,
    reference: Path,
    outdir: Path,
    name: str,
    loci_bed: Path | None,
    min_map_q: int,
    max_tlen: int | None,
    max_softclip: int | None,
    max_nm: int | None,
    min_site_q: int,
    min_geno_q: int,
    min_base_q: int,
    min_sample_depth: int,  # sample must have depth cov or site is masked.
    min_locus_sample_coverage: int,  # locus must have data for N samples (used in locus delim)
    min_locus_trim_sample_coverage: int,  # trim r/l to region with at least N samples data (default 4)
    min_locus_length: int,
    min_locus_merge_distance: int,  # merge loci within this distance
    max_locus_hetero_frequency: float,
    max_locus_variant_frequency: float,
    max_sample_hetero_frequency: float,
    softclip_len_threshold: int,
    softclip_frac_max: float,
    depth_z_max: float,
    third_frac_cut: float,
    min_3allele_sites: int,
    maf_threshold: float,
    max_sites_above_maf: int,
    paralog_fail_frac_max: float,
    populations: Path | None,
    rename_bams: Path | None,
    masks: List[str] | None,
    cores: int,
    threads: int,
    force: bool,
    log_level: str,
    min_aligned_len: int | None = None,
):
    # Normalize the top-level input/output paths first so later stages can
    # treat everything as concrete local files.
    loci_bed = loci_bed.expanduser().absolute() if loci_bed else None
    reference = reference.expanduser().absolute()
    outdir = outdir.expanduser().absolute()
    tmpdir = outdir / f"{name}_tmpdir"

    # Assemble runs multi-threaded worker jobs, so derive the parallel job
    # count up front and validate the simple numeric filters before any work.
    workers = max(1, cores // threads)

    if max_tlen is not None and max_tlen < 0:
        raise IPyradError("max_tlen must be >= 0 when provided.")
    if max_softclip is not None and max_softclip < 0:
        raise IPyradError("max_softclip must be >= 0 when provided.")
    if max_nm is not None and max_nm < 0:
        raise IPyradError("max_nm must be >= 0 when provided.")
    if min_aligned_len is not None and min_aligned_len < 0:
        raise IPyradError("min_aligned_len must be >= 0 when provided.")
    if not 0 <= max_sample_hetero_frequency <= 1:
        raise IPyradError("max_sample_hetero_frequency must be between 0 and 1.")
    if min_3allele_sites < 0:
        raise IPyradError("min_3allele_sites must be >= 0.")
    if max_sites_above_maf < 0:
        raise IPyradError("max_sites_above_maf must be >= 0.")

    # Apply overwrite policy before creating any temp outputs for this run.
    existing_results_force_or_raise(outdir, tmpdir, name, force)

    # Create the working directory layout used by the active assemble path.
    outdir.mkdir(exist_ok=True, parents=True)
    tmpdir.mkdir(exist_ok=True)
    bed_dir = tmpdir / "beds"
    bed_dir.mkdir(exist_ok=True)
    phase_dir = tmpdir / "phase"
    phase_dir.mkdir(exist_ok=True)

    expanded_rad_bams = (
        [bam_file.expanduser().absolute() for bam_file in rad_bams] if rad_bams else []
    )
    expanded_wgs_bams = (
        [bam_file.expanduser().absolute() for bam_file in wgs_bams] if wgs_bams else []
    )

    # Validate that some BAM inputs exist before parsing optional rename maps.
    # This keeps the main missing-input errors stable instead of surfacing
    # rename-file validation noise first.
    bam_dict: dict[str, Path] = {}
    wgs_dict: dict[str, Path] = {}
    if not expanded_rad_bams and loci_bed is None:
        raise IPyradError(
            "No RAD bam files found. These are required unless --loci-bed is provided."
        )
    if not expanded_rad_bams and not expanded_wgs_bams:
        raise IPyradError(
            "No input BAM files found. Provide --rad-bams and/or --wgs-bams."
        )

    rename_map: dict[str, str] = {}
    if rename_bams is not None:
        rename_map = _normalize_bam_rename_file(
            rename_bams,
            expanded_rad_bams + expanded_wgs_bams,
        )

    # Load the RAD BAMs that define loci when no explicit loci BED is supplied.
    logger.info("loading sample meta-data")
    bam_dict, rad_renamed = _collect_named_bams(expanded_rad_bams, rename_map)
    if bam_dict:
        logger.info(f"loaded {len(bam_dict)} RAD samples")

    # Optional WGS BAMs are still normalized and filtered here so later
    # assemble milestones can reuse the same prepared analysis BAMs.
    wgs_dict, wgs_renamed = _collect_named_bams(expanded_wgs_bams, rename_map)
    if wgs_dict:
        logger.info(f"loaded {len(wgs_dict)} WGS samples")

    renamed_pairs = sorted(rad_renamed + wgs_renamed)
    if renamed_pairs:
        shown = min(10, len(renamed_pairs))
        logger.info(
            "renamed {} BAM sample name(s) from --rename-bams", len(renamed_pairs)
        )
        logger.info(
            "showing first {}/{} BAM rename mappings", shown, len(renamed_pairs)
        )
        max_len = max(len(sample_name) for sample_name, _ in renamed_pairs[:shown])
        for sample_name, bam_name in renamed_pairs[:shown]:
            logger.info("{} <- {}", sample_name.ljust(max_len), bam_name)

    duplicate_names = sorted(set(bam_dict) & set(wgs_dict))
    if duplicate_names:
        joined = ", ".join(duplicate_names)
        raise IPyradError(
            f"RAD and WGS inputs resolve to duplicate sample names: {joined}"
        )

    # Prepare one combined sample map, then replace the original BAMs with the
    # filtered analysis BAMs that downstream assemble stages should consume.
    all_dict = wgs_dict | bam_dict
    snames = sorted(all_dict)
    startup_workers = max(1, min(cores, len(snames)))
    consensus_workers = max(1, min(cores, len(snames)))
    final_vcf_mask_workers = max(1, min(cores, len(snames)))
    all_dict = {i: all_dict[i] for i in snames}
    group_samples_file = None
    if populations is not None:
        (
            group_samples_file,
            population_imap,
            parsed_minmap,
        ) = _normalize_populations_file(populations, tmpdir, snames)
        logger.info("grouped variant calling enabled from --populations")
        logger.info(
            "loaded {} population group(s): {}",
            len(population_imap),
            ", ".join(
                f"{group}={len(names)}" for group, names in population_imap.items()
            ),
        )
        if parsed_minmap is not None:
            logger.info(
                "ignoring per-population minmap thresholds in --populations; assemble currently uses this file for grouped calling only"
            )
        logger.debug("normalized grouped-calling samples file: {}", group_samples_file)

    bam_metadata = _collect_bam_metadata(all_dict, log_level, startup_workers)
    sample_layouts = {
        sname: (bam_metadata[sname]["layout"] == "paired")
        for sname in snames
    }
    n_paired = sum(sample_layouts.values())
    n_single = len(sample_layouts) - n_paired
    logger.info(
        "assemble BAM layout summary: {} paired-end, {} single-end",
        n_paired,
        n_single,
    )
    if n_paired and n_single:
        logger.info(
            "mixed single-end and paired-end BAMs detected across samples; paired-end samples use fragment-span coverage and single-end samples use read-span coverage during locus delimiting"
        )

    logger.debug("fetching reference scaffold order")
    get_reference_sort_order(reference, tmpdir)
    logger.info("validating input BAMs against the current reference")
    _validate_bam_header_records_match_reference(
        {
            sname: bam_metadata[sname]["header_records"]
            for sname in snames
        },
        tmpdir,
        reference,
    )

    all_dict = _prepare_analysis_bams(
        bam_dict=all_dict,
        sample_layouts=sample_layouts,
        tmpdir=tmpdir,
        min_map_q=min_map_q,
        max_tlen=max_tlen,
        max_softclip=max_softclip,
        max_nm=max_nm,
        min_aligned_len=min_aligned_len,
        threads=threads,
        workers=workers,
        log_level=log_level,
    )
    all_dict = {sname: all_dict[sname] for sname in snames}
    bam_dict = {sname: all_dict[sname] for sname in sorted(bam_dict)}
    wgs_dict = {sname: all_dict[sname] for sname in sorted(wgs_dict)}

    # Record runtime settings and initialize reference metadata before locus
    # delimiting starts.
    logger.info(
        f"using up to {cores} cores (up to {workers} multi-threaded jobs using {threads} threads)"
    )
    normalized_loci_bed = None
    if loci_bed is not None:
        normalized_loci_bed, input_locus_count = _normalize_user_loci_bed(
            loci_bed, tmpdir
        )
        logger.info(
            "using provided loci BED instead of RAD-based shared-locus delimiting"
        )
        logger.info("loaded {} loci from {}", input_locus_count, normalized_loci_bed)
        logger.debug(
            "ignoring RAD-delimiting options because --loci-bed was provided: min_locus_sample_coverage={}, min_locus_length={}, min_locus_merge_distance={}",
            min_locus_sample_coverage,
            min_locus_length,
            min_locus_merge_distance,
        )

    # [1] Locus delimiting:
    # Build per-sample coverage intervals for all analysis BAMs so both RAD and
    # WGS samples can later generate low-depth masks, but define the shared loci
    # BED only from the RAD subset.
    logger.info("delimiting sample coverage beds")
    kwargs = dict(
        reference=reference,
        tmpdir=tmpdir,
        min_map_q=min_map_q,
        min_sample_depth=min_sample_depth,
        min_merge_distance=min_locus_merge_distance,
        threads=threads,
    )
    jobs = {}
    for sname, bam_file in all_dict.items():
        ikwargs = kwargs | dict(
            sname=sname,
            bam_file=bam_file,
            is_paired=sample_layouts[sname],
        )
        jobs[sname] = (get_coverage_bed_graphs, ikwargs)
    run_with_pool(jobs, log_level, workers, msg="Delimiting sample coverage beds")

    if normalized_loci_bed is None:
        logger.info("delimiting shared coverage beds (loci)")
        args = (
            list(bam_dict),
            min_locus_sample_coverage,
            min_locus_merge_distance,
            min_locus_length,
            ".fragments.merged.bed",
            tmpdir,
        )
        loci_bed = get_across_sample_loci_bed(*args)

        # Preserve the pre-paralog shared loci so the raw RAD-defined windows stay
        # available for debugging after the filtered BED becomes canonical.
        shutil.copy2(loci_bed, bed_dir / "loci.raw.bed")
        shared_loci_after_delimiting = _count_nonempty_lines(bed_dir / "loci.raw.bed")
    else:
        loci_bed = normalized_loci_bed
        shared_loci_after_delimiting = _count_nonempty_lines(loci_bed)

    # [2] Paralog filtering:
    # Score every sample against the shared RAD-defined loci BED, reduce those
    # calls across samples, and promote the filtered shared BED to the
    # canonical loci.bed consumed by the rest of the assemble workflow.
    final_loci_bed = _run_paralog_stage(
        sample_bams=all_dict,
        regions_bed=loci_bed,
        reference=reference,
        bed_dir=bed_dir,
        phase_dir=phase_dir,
        min_map_q=min_map_q,
        min_base_q=min_base_q,
        softclip_len_threshold=softclip_len_threshold,
        softclip_frac_max=softclip_frac_max,
        depth_z_max=depth_z_max,
        third_frac_cut=third_frac_cut,
        min_3allele_sites=min_3allele_sites,
        maf_threshold=maf_threshold,
        max_sites_above_maf=max_sites_above_maf,
        paralog_fail_frac_max=paralog_fail_frac_max,
        workers=workers,
        log_level=log_level,
        rad_sample_names=sorted(bam_dict),
        wgs_sample_names=sorted(wgs_dict),
    )
    logger.info("paralog-filtered shared BED written to {}", final_loci_bed)
    shared_loci_after_paralog_filtering = _count_nonempty_lines(final_loci_bed)

    # [3] Variant calling:
    # Jointly call variants across every filtered analysis BAM inside the
    # canonical shared loci, then apply project-wide genotype/site filtering.
    with profile_stage("variant calling"):
        _run_variant_stage(
            tmpdir=tmpdir,
            reference=reference,
            bam_dict=all_dict,
            group_samples_file=group_samples_file,
            min_map_q=min_map_q,
            min_base_q=min_base_q,
            min_sample_depth=min_sample_depth,
            min_geno_q=min_geno_q,
            min_site_q=min_site_q,
            cores=cores,
            threads=threads,
            log_level=log_level,
            wgs_samples=sorted(wgs_dict) if (bam_dict and wgs_dict) else None,
        )

    # [4] Sample masks:
    # Build low-depth masks for all samples, add sample-specific RAD paralog
    # masks where available, and build the merged consensus BED masks that
    # sample-specific FASTA generation consumes directly.
    with profile_stage("sample mask building"):
        paralog_masks = _build_sample_masks(
            tmpdir=tmpdir,
            all_dict=all_dict,
            min_sample_depth=min_sample_depth,
            workers=workers,
            log_level=log_level,
        )

    # [5] Final outputs:
    # Write consensus FASTAs, build the locus database, and then materialize
    # the final loci/bed/vcf/hdf5 outputs from the fully filtered locus set.
    _write_consensus_and_outputs(
        name=name,
        outdir=outdir,
        tmpdir=tmpdir,
        snames=snames,
        reference=reference,
        masks=masks,
        sample_masks=paralog_masks,
        shared_loci_after_delimiting=shared_loci_after_delimiting,
        shared_loci_after_paralog_filtering=shared_loci_after_paralog_filtering,
        min_locus_sample_coverage=min_locus_sample_coverage,
        min_locus_trim_sample_coverage=min_locus_trim_sample_coverage,
        min_locus_length=min_locus_length,
        max_locus_hetero_frequency=max_locus_hetero_frequency,
        max_locus_variant_frequency=max_locus_variant_frequency,
        max_sample_hetero_frequency=max_sample_hetero_frequency,
        consensus_workers=consensus_workers,
        final_vcf_mask_workers=final_vcf_mask_workers,
        workers=workers,
        threads=threads,
        log_level=log_level,
        cores=cores,
        rad_samples=sorted(bam_dict),
        wgs_samples=sorted(wgs_dict),
    )
    logger.info("assemble outputs written to {}", outdir)
    return
