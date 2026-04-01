#!/usr/bin/env python

"""Map reads to a reference and write coordinate-sorted BAMs."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Dict, Tuple

from loguru import logger

from ..utils.exceptions import IPyradError
from ..utils.names import expand_path
from ..utils.parallel import run_pipeline
from ..utils.parallel import run_with_pool
from .map_samples_prep import prepare_map_samples
from .map_stats import MappingJobResult
from .map_stats import collect_paired_bam_stats
from .map_stats import collect_single_end_bam_stats
from .map_stats import parse_markdup_report
from .map_stats import render_map_stats_report


BIN = Path(sys.prefix) / "bin"
BIN_BWA = str(BIN / "bwa-mem2")
BIN_SAMTOOLS = str(BIN / "samtools")


def _require(condition: bool, message: str) -> None:
    """Raise a user-facing mapper error when a precondition is not met."""
    if not condition:
        raise IPyradError(message)


def _output_bam_path(sname: str, outdir: Path) -> Path:
    """Return the final BAM output path for one sample."""
    return outdir / f"{sname}.filtered.bam"


def _output_index_path(bam_path: Path) -> Path:
    """Return the CSI index path for one BAM path."""
    return bam_path.with_suffix(bam_path.suffix + ".csi")


def _remove_output_artifacts(bam_path: Path) -> None:
    """Remove one BAM and its CSI index when they exist."""
    for path in (bam_path, _output_index_path(bam_path)):
        try:
            path.unlink()
        except FileNotFoundError:
            pass


def _cleanup_prefix(prefix: Path) -> None:
    """Remove temporary files created with one samtools prefix."""
    for path in prefix.parent.glob(prefix.name + "*"):
        try:
            path.unlink()
        except FileNotFoundError:
            pass


def _load_save_counts(path: Path) -> dict[str, int]:
    """Load one samtools --save-counts JSON file."""
    with path.open("r", encoding="utf-8") as infile:
        data = json.load(infile)
    return {key: int(value) for key, value in data.items()}


def _finalize_indexed_bam(tmp_bam: Path, final_bam: Path) -> None:
    """Rename a temporary indexed BAM into its final output path."""
    tmp_index = _output_index_path(tmp_bam)
    final_index = _output_index_path(final_bam)
    os.replace(tmp_bam, final_bam)
    os.replace(tmp_index, final_index)


def _check_mapper_dependencies() -> None:
    """Ensure required mapper tool binaries are available."""
    for label, tool_path in (("bwa-mem2", BIN_BWA), ("samtools", BIN_SAMTOOLS)):
        path = Path(tool_path)
        _require(path.exists(), f"Required mapper binary not found: {tool_path}")
        _require(os.access(path, os.X_OK), f"Required mapper binary is not executable: {tool_path}")


def _validate_runtime_settings(
    cores: int,
    threads: int,
    is_paired: bool,
    mark_dups_by_coords: bool,
    mark_dups_by_umis: bool,
) -> None:
    """Validate mapper runtime settings before launching jobs."""
    _require(cores > 0, "cores must be a positive integer.")
    _require(threads > 0, "threads must be a positive integer.")
    _require(threads <= cores, "threads cannot exceed cores.")
    if mark_dups_by_coords and mark_dups_by_umis:
        raise IPyradError("you cannot select both mark_dups_by_coords and mark_dups_by_umis.")
    if (mark_dups_by_coords or mark_dups_by_umis) and not is_paired:
        raise IPyradError("Data do not appear to be paired. Cannot remove duplicates for SE data.")


def _check_duplicate_mode_warnings(
    mark_dups_by_coords: bool,
    mark_dups_by_umis: bool,
) -> None:
    """Emit duplicate-removal mode warnings."""
    if mark_dups_by_coords:
        logger.warning(
            "removing PCR duplicates by coordinates. Be sure this run includes only WGS samples, not RAD"
        )
    if mark_dups_by_umis:
        logger.warning(
            "removing PCR duplicates by i5 UMIs. Be sure you ran `ipyrad2 trim` with `-U` to store i5 tags for these samples"
        )


def _index_ref_with_bwa(reference: Path) -> None:
    """Index the reference with bwa-mem2 unless it is already indexed."""
    _require(reference.exists(), f"reference path {reference} does not exist.")

    suffixes = [".pac", ".ann", ".amb", ".0123", ".bwt.2bit.64"]
    index_paths = [reference.with_suffix(reference.suffix + suffix) for suffix in suffixes]
    if all(path.exists() for path in index_paths):
        logger.debug(f"reference is already bwa indexed: {reference}")
        return

    if not os.access(reference.parent, os.W_OK | os.X_OK):
        raise IPyradError(
            "cannot index reference because you do not have write access to its directory."
        )

    logger.info(f"indexing reference: {reference.name}")
    cmd = [
        BIN_BWA, "index",
        str(reference),
    ]
    run_pipeline([cmd])


def _bwa_mem_cmd(
    sname: str,
    fastqs: Tuple[Path, Path | None],
    reference: Path,
    bwa_threads: int,
) -> list[str]:
    """Build the bwa-mem2 alignment command for one sample."""
    cmd = [
        BIN_BWA, "mem",
        "-Y",
        "-T", "20",
        "-R", f"@RG\\tID:{sname}\\tSM:{sname}",
        "-K", "50000000",
        "-v", "1",
        "-t", str(bwa_threads),
        str(reference),
        str(fastqs[0]),
    ]
    if fastqs[1] is not None:
        cmd.append(str(fastqs[1]))
    return cmd


def _structural_filter_cmd(counts_path: Path) -> list[str]:
    """Build the primary mapped-read filter used before BAM writing."""
    return [
        BIN_SAMTOOLS, "view",
        "-b",
        "-u",
        "-F", "0x4",
        "-F", "0x100",
        "-F", "0x200",
        "-F", "0x800",
        "--save-counts", str(counts_path),
        "-o", "-",
    ]


def _paired_same_scaffold_filter_cmd(counts_path: Path) -> list[str]:
    """Keep only mapped paired reads whose mates also map to the same scaffold."""
    return [
        BIN_SAMTOOLS, "view",
        "-b",
        "-u",
        "-e", '((flag&8)==0) && (rnext=="=" || rnext==rname)',
        "--save-counts", str(counts_path),
        "-o", "-",
    ]


def _coord_sort_cmd(
    tmp_prefix: Path,
    threads: int,
    outfile: Path,
) -> list[str]:
    """Build one coordinate-sort command that also writes the CSI index."""
    return [
        BIN_SAMTOOLS, "sort",
        "-m", "256M",
        "-T", str(tmp_prefix),
        "-@", str(threads),
        "--write-index",
        "-O", "bam",
        "-o", str(outfile),
        "-",
    ]


def _namesort_cmd(
    tmp_prefix: Path,
    threads: int,
    outfile: Path,
) -> list[str]:
    """Build one name-grouped sort command for fixmate."""
    return [
        BIN_SAMTOOLS, "sort",
        "-n",
        "-m", "256M",
        "-T", str(tmp_prefix),
        "-@", str(threads),
        "-O", "bam",
        "-o", str(outfile),
        "-",
    ]


def _fixmate_cmd(
    threads: int,
    infile: Path,
    outfile: Path,
) -> list[str]:
    """Build one samtools fixmate command."""
    return [
        BIN_SAMTOOLS, "fixmate",
        "-m",
        "-@", str(threads),
        str(infile),
        str(outfile),
    ]


def _markdup_cmd(
    tmp_prefix: Path,
    threads: int,
    infile: Path,
    outfile: Path,
    stats_path: Path,
    mark_dups_by_umis: bool,
) -> list[str]:
    """Build one samtools markdup duplicate-removal command."""
    cmd = [
        BIN_SAMTOOLS, "markdup",
        "-r",
        "-T", str(tmp_prefix),
        "-s",
        "-f", str(stats_path),
        "-@", str(threads),
        "--write-index",
        str(infile),
        str(outfile),
    ]
    if mark_dups_by_umis:
        cmd.extend([
            "--barcode-rgx", "UMI_([ACGTN]+)",
        ])
    return cmd


def _map_sample(
    sname: str,
    fastqs: Tuple[Path, Path | None],
    reference: Path,
    outdir: Path,
    threads: int,
    is_paired: bool,
    mark_dups_by_coords: bool,
    mark_dups_by_umis: bool,
) -> MappingJobResult:
    """Map one sample and return summary metadata for final-BAM stats."""
    out_bam = _output_bam_path(sname, outdir)
    out_bam_tmp = outdir / f"{sname}.filtered.tmp.bam"
    tmpdir = outdir / "tmpdir"
    tmp_prefix = tmpdir / f"{sname}.tmp.pre"
    structural_counts_path = tmpdir / f"{sname}.tmp.structural_counts.json"
    same_scaffold_counts_path = tmpdir / f"{sname}.tmp.same_scaffold_counts.json"
    dup_stats_path = tmpdir / f"{sname}.tmp.markdup_stats.txt"
    bam_namesort = tmpdir / f"{sname}.tmp.namesort.bam"
    bam_fixmate = tmpdir / f"{sname}.tmp.fixmate.bam"
    bam_coordsort = tmpdir / f"{sname}.tmp.coordsort.bam"
    temp_paths = [
        out_bam_tmp,
        _output_index_path(out_bam_tmp),
        structural_counts_path,
        same_scaffold_counts_path,
        dup_stats_path,
        bam_namesort,
        bam_fixmate,
        bam_coordsort,
        _output_index_path(bam_coordsort),
    ]

    bwa_threads = max(1, int(threads * 0.75))
    sort_threads = max(1, threads - bwa_threads)
    duplicate_stats: dict[str, int] = {}

    try:
        cmd1 = _bwa_mem_cmd(sname, fastqs, reference, bwa_threads)
        cmd2 = _structural_filter_cmd(structural_counts_path)
        filter_pipeline = [cmd1, cmd2]
        if is_paired:
            filter_pipeline.append(_paired_same_scaffold_filter_cmd(same_scaffold_counts_path))

        if is_paired and (mark_dups_by_coords or mark_dups_by_umis):
            cmd3 = _namesort_cmd(tmp_prefix, sort_threads, bam_namesort)
            run_pipeline(filter_pipeline + [cmd3])

            cmd4 = _fixmate_cmd(threads, bam_namesort, bam_fixmate)
            run_pipeline([cmd4])

            cmd5 = [
                BIN_SAMTOOLS, "sort",
                "-m", "256M",
                "-T", str(tmp_prefix),
                "-@", str(threads),
                "-O", "bam",
                "-o", str(bam_coordsort),
                str(bam_fixmate),
            ]
            run_pipeline([cmd5])

            cmd6 = _markdup_cmd(
                tmp_prefix=tmp_prefix,
                threads=threads,
                infile=bam_coordsort,
                outfile=out_bam_tmp,
                stats_path=dup_stats_path,
                mark_dups_by_umis=mark_dups_by_umis,
            )
            run_pipeline([cmd6])
            duplicate_stats = parse_markdup_report(dup_stats_path)
        else:
            cmd3 = _coord_sort_cmd(tmp_prefix, sort_threads, out_bam_tmp)
            run_pipeline(filter_pipeline + [cmd3])

        structural_counts = _load_save_counts(structural_counts_path)
        same_scaffold_counts = (
            _load_save_counts(same_scaffold_counts_path)
            if is_paired
            else {"records_filter_rejected": 0, "records_filter_accepted": structural_counts["records_filter_accepted"]}
        )
        _finalize_indexed_bam(out_bam_tmp, out_bam)
        logger.debug(f"finished mapping: {sname}")
        return MappingJobResult(
            sname=sname,
            bam_path=out_bam,
            is_paired=is_paired,
            nreads_processed=structural_counts["records_processed"],
            nreads_filtered_before_bam_by_unmapped_or_nonprimary=structural_counts["records_filter_rejected"],
            nreads_filtered_before_bam_by_mate_unmapped_or_cross_scaffold=same_scaffold_counts["records_filter_rejected"],
            nreads_written_before_duplicate_removal=same_scaffold_counts["records_filter_accepted"],
            duplicate_stats=duplicate_stats,
        )
    finally:
        for path in temp_paths:
            try:
                path.unlink()
            except FileNotFoundError:
                pass
        _cleanup_prefix(tmp_prefix)


def _select_output_samples(
    fastq_dict: Dict[str, Tuple[Path, Path | None]],
    outdir: Path,
    force: bool,
) -> Dict[str, Tuple[Path, Path | None]]:
    """Apply overwrite policy and return the samples that should be processed."""
    selected: Dict[str, Tuple[Path, Path | None]] = {}
    skipped = []

    for sname, fastq_tuple in fastq_dict.items():
        bam_path = _output_bam_path(sname, outdir)
        index_path = _output_index_path(bam_path)
        if not bam_path.exists() and not index_path.exists():
            selected[sname] = fastq_tuple
            continue

        if force:
            _remove_output_artifacts(bam_path)
            selected[sname] = fastq_tuple
            logger.debug(f"removing existing bam outputs for sample: {sname}")
            continue

        skipped.append(sname)

    if skipped:
        logger.warning(
            "skipping {}/{} samples that already have results (.bam/.bam.csi) in outdir. Use --force to overwrite.",
            len(skipped),
            len(fastq_dict),
        )
    return selected


def _next_stats_path(outdir: Path) -> Path:
    """Return the next free ipyrad_map_stats_N.txt path in the output directory."""
    idx = 0
    while True:
        outstats = outdir / f"ipyrad_map_stats_{idx}.txt"
        if not outstats.exists():
            return outstats
        idx += 1


def run_mapper(
    fastqs,
    outdir: Path,
    reference: Path,
    imap: Path | None,
    cores: int,
    threads: int,
    force: bool,
    mark_dups_by_coords: bool,
    mark_dups_by_umis: bool,
    delim_str: str | None,
    delim_idx: int,
    log_level: str,
):
    """Run the ipyrad2 map workflow."""
    reference = expand_path(reference)
    outdir = expand_path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    tmpdir = outdir / "tmpdir"
    tmpdir.mkdir(parents=True, exist_ok=True)
    imap = expand_path(imap) if imap is not None else None

    fastq_dict, is_paired = prepare_map_samples(
        fastqs=fastqs,
        delim_str=delim_str,
        delim_idx=delim_idx,
        imap=imap,
        tmpdir=tmpdir,
    )

    _check_mapper_dependencies()
    _validate_runtime_settings(
        cores=cores,
        threads=threads,
        is_paired=is_paired,
        mark_dups_by_coords=mark_dups_by_coords,
        mark_dups_by_umis=mark_dups_by_umis,
    )
    _check_duplicate_mode_warnings(
        mark_dups_by_coords=mark_dups_by_coords,
        mark_dups_by_umis=mark_dups_by_umis,
    )

    fastq_dict = _select_output_samples(fastq_dict, outdir, force)
    if not fastq_dict:
        logger.info("all samples are completed.")
        raise SystemExit(0)

    _index_ref_with_bwa(reference)

    workers = max(1, cores // threads)
    logger.info(
        "mapping {} samples to coordinate-sorted BAMs in {}",
        len(fastq_dict),
        outdir,
    )
    logger.info(
        "using up to {} cores (up to {} multi-threaded jobs using {} threads)",
        cores,
        workers,
        threads,
    )

    map_jobs = {
        sname: (
            _map_sample,
            {
                "sname": sname,
                "fastqs": fastq_tuple,
                "reference": reference,
                "outdir": outdir,
                "threads": threads,
                "is_paired": is_paired,
                "mark_dups_by_coords": mark_dups_by_coords,
                "mark_dups_by_umis": mark_dups_by_umis,
            },
        )
        for sname, fastq_tuple in fastq_dict.items()
    }
    map_results = run_with_pool(
        map_jobs,
        log_level,
        max_workers=workers,
        msg="Mapping",
    )

    stats_func = collect_paired_bam_stats if is_paired else collect_single_end_bam_stats
    stats_jobs = {
        sname: (
            stats_func,
            {
                "job_result": job_result,
            },
        )
        for sname, job_result in map_results.items()
    }
    stats = run_with_pool(
        stats_jobs,
        log_level,
        max_workers=workers,
        msg="Gathering mapping stats",
    )

    outstats = _next_stats_path(outdir)
    outstats.write_text(
        render_map_stats_report(stats, is_paired),
        encoding="utf-8",
    )
    logger.info(f"mapping stats written to {outstats}")


if __name__ == "__main__":
    pass
