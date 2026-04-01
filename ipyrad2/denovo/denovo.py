#!/usr/bin/env python

"""Build a denovo reference library by clustering reads and locus consensuses."""

from __future__ import annotations

import os
import shutil
import sys
import csv
from pathlib import Path
from typing import Any

from loguru import logger

from .align import AlignmentRunSummary, write_ordered_consensus_stream_to_file
from .cluster import build_sample_summary, concat_summaries
from .graph import make_global_tables
from ..utils.exceptions import IPyradError
from ..utils.names import get_name_to_fastq_dict
from ..utils.parallel import run_pipeline, run_with_pool


WORKDIR_NAME = "_denovo_work"


def _validate_runtime_args(
    within_similarity: float,
    across_similarity: float,
    min_derep_size: int,
    min_length: int,
    min_merge_overlap: int,
    max_merge_diffs: int,
    cores: int,
    threads: int,
    delim_idx: int,
) -> None:
    """Validate runner arguments for direct library use."""
    if not 0 < within_similarity <= 1:
        raise IPyradError("within_similarity must be > 0 and <= 1")
    if not 0 < across_similarity <= 1:
        raise IPyradError("across_similarity must be > 0 and <= 1")
    if min_derep_size < 1:
        raise IPyradError("min_derep_size must be >= 1")
    if min_length < 1:
        raise IPyradError("min_length must be >= 1")
    if min_merge_overlap < 1:
        raise IPyradError("min_merge_overlap must be >= 1")
    if max_merge_diffs < 0:
        raise IPyradError("max_merge_diffs must be >= 0")
    if cores < 1:
        raise IPyradError("cores must be >= 1")
    if threads < 1:
        raise IPyradError("threads must be >= 1")
    if threads > cores:
        raise IPyradError("threads cannot exceed cores")
    if delim_idx < 1:
        raise IPyradError("delim_idx must be >= 1")


def _is_executable(path: Path) -> bool:
    """Return True if a path exists and is executable."""
    return path.exists() and path.is_file() and os.access(path, os.X_OK)


def _resolve_binary(binary: Path | None, name: str) -> str:
    """Resolve a tool binary from an explicit path, the active env, or PATH."""
    if binary is not None:
        candidate = binary.expanduser().absolute()
        if not _is_executable(candidate):
            raise IPyradError(f"{name} binary is not executable: {candidate}")
        return str(candidate)

    env_candidate = Path(sys.prefix) / "bin" / name
    if _is_executable(env_candidate):
        return str(env_candidate)

    resolved = shutil.which(name)
    if resolved:
        return resolved

    raise IPyradError(
        f"Cannot find the '{name}' executable. Set --{name}-binary explicitly "
        "or install it into the active environment or PATH."
    )


def _validate_fastq_layout(
    fastq_dict: dict[str, tuple[Path, Path | None]],
) -> bool:
    """Validate parsed FASTQ layout and return whether the run is paired-end."""
    if not fastq_dict:
        raise IPyradError("No FASTQ inputs were parsed for denovo.")

    paired_states = {paths[1] is not None for paths in fastq_dict.values()}
    if len(paired_states) > 1:
        raise IPyradError(
            "denovo requires all inputs to be consistently single-end or paired-end; "
            "mixed input layouts were detected."
        )
    return paired_states.pop()


def _iter_denovo_outputs(outdir: Path) -> list[Path]:
    """Return the curated denovo output paths for one output directory."""
    return [
        outdir / "denovo_reference.fa",
        outdir / "loci.mapping.tsv",
        outdir / "loci.stats.tsv",
        outdir / "denovo.stats.txt",
        outdir / "denovo.audit",
    ]


def _prepare_output_paths(
    outdir: Path,
    force: bool,
) -> tuple[Path, dict[str, Path]]:
    """Prepare curated denovo outputs and the internal working directory."""
    outdir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "reference": outdir / "denovo_reference.fa",
        "mapping": outdir / "loci.mapping.tsv",
        "loci_stats": outdir / "loci.stats.tsv",
        "run_stats": outdir / "denovo.stats.txt",
        "audit_dir": outdir / "denovo.audit",
        "workdir": outdir / WORKDIR_NAME,
    }
    existing = [path for path in outputs.values() if path.exists()]
    if existing and not force:
        joined = ", ".join(path.name for path in existing)
        raise IPyradError(
            f"denovo outputs already exist in {outdir}. Use --force to overwrite: {joined}"
        )

    if force:
        for path in _iter_denovo_outputs(outdir):
            if path.exists():
                if path.is_dir():
                    shutil.rmtree(path)
                else:
                    path.unlink()
        if outputs["workdir"].exists():
            shutil.rmtree(outputs["workdir"])

    outputs["workdir"].mkdir(parents=True, exist_ok=True)
    return outputs["workdir"], outputs


def _write_denovo_stats(
    outpath: Path,
    *,
    fastq_dict: dict[str, tuple[Path, Path | None]],
    paired: bool,
    within_similarity: float,
    across_similarity: float,
    min_derep_size: int,
    min_length: int,
    min_merge_overlap: int,
    max_merge_diffs: int,
    allow_reverse_complement: bool,
    cores: int,
    threads: int,
    workers: int,
    graph_splitter: str,
    alignment_summary: AlignmentRunSummary,
    vsearch_binary: str,
    mafft_binary: str,
    keep_intermediates: bool,
    workdir: Path,
    mapping_df: Any,
    stats_df: Any,
    outputs: dict[str, Path],
) -> None:
    """Write a human-readable summary of the denovo run."""
    nloci = int(stats_df.shape[0]) if stats_df is not None else 0
    ncores = int(mapping_df.shape[0]) if mapping_df is not None else 0
    if stats_df is not None and not getattr(stats_df, "empty", False):
        duplicated_components_seen = (
            int(stats_df.loc[stats_df["duplicated_component"], "component_id"].nunique())
            if "duplicated_component" in stats_df.columns
            else 0
        )
        duplicated_components_aligned = (
            int(stats_df.loc[stats_df["aligned_for_reconciliation"], "component_id"].nunique())
            if "aligned_for_reconciliation" in stats_df.columns
            else 0
        )
        components_reconciled = (
            int(stats_df.loc[stats_df["used_reconciliation"], "component_id"].nunique())
            if "used_reconciliation" in stats_df.columns
            else 0
        )
        joined_only_reconciled_loci = int(
            (stats_df["reconcile_mode"] == "joined_only").sum()
        ) if "reconcile_mode" in stats_df.columns else 0
        mixed_reconciled_loci = int(
            (stats_df["reconcile_mode"] == "mixed").sum()
        ) if "reconcile_mode" in stats_df.columns else 0
        mixed_reconciled_groups = int(stats_df["n_reconciled_groups"].sum()) if "n_reconciled_groups" in stats_df.columns else 0
    else:
        duplicated_components_seen = 0
        duplicated_components_aligned = 0
        components_reconciled = 0
        joined_only_reconciled_loci = 0
        mixed_reconciled_loci = 0
        mixed_reconciled_groups = 0
    lines = [
        "Inputs",
        f"  fastq_count: {sum(2 if paths[1] else 1 for paths in fastq_dict.values())}",
        f"  sample_count: {len(fastq_dict)}",
        f"  paired_mode: {'paired-end' if paired else 'single-end'}",
        f"  sample_names: {', '.join(sorted(fastq_dict))}",
        "",
        "Clustering",
        f"  within_similarity: {within_similarity}",
        f"  across_similarity: {across_similarity}",
        f"  min_derep_size: {min_derep_size}",
        f"  min_length: {min_length}",
        f"  min_merge_overlap: {min_merge_overlap}",
        f"  max_merge_diffs: {max_merge_diffs}",
        f"  allow_reverse_complement: {allow_reverse_complement}",
        "",
        "Runtime",
        f"  cores: {cores}",
        f"  vsearch_threads_per_job: {threads}",
        f"  vsearch_worker_processes: {workers}",
        f"  graph_splitter: {graph_splitter}",
        f"  mafft_threads_per_job: {alignment_summary.mafft_threads_per_job}",
        f"  mafft_worker_processes: {alignment_summary.mafft_worker_processes}",
        f"  alignment_mode: {alignment_summary.alignment_mode}",
        f"  mafft_timeout_seconds: {alignment_summary.mafft_timeout_seconds}",
        "  duplicated_component_reconciliation: aligned",
        f"  cluster_spacer_mode: stripped",
        f"  output_spacer_length: {alignment_summary.output_spacer_length}",
        f"  vsearch_binary: {vsearch_binary}",
        f"  mafft_binary: {mafft_binary}",
        f"  keep_intermediates: {keep_intermediates}",
        f"  workdir: {workdir}",
        "",
        "Results",
        f"  consensus_records: {ncores}",
        f"  loci_written: {nloci}",
        f"  single_sequence_loci: {alignment_summary.single_sequence_loci}",
        f"  identical_sequence_loci: {alignment_summary.identical_sequence_loci}",
        f"  mafft_required_loci: {alignment_summary.mafft_required_loci}",
        f"  joined_spacer_loci: {alignment_summary.joined_spacer_loci}",
        f"  mixed_reconciled_spacer_loci: {alignment_summary.mixed_reconciled_spacer_loci}",
        f"  stripped_output_loci: {alignment_summary.stripped_output_loci}",
        f"  duplicated_components_seen: {duplicated_components_seen}",
        f"  duplicated_components_aligned: {duplicated_components_aligned}",
        f"  components_reconciled: {components_reconciled}",
        f"  joined_only_reconciled_loci: {joined_only_reconciled_loci}",
        f"  mixed_reconciled_loci: {mixed_reconciled_loci}",
        f"  mixed_reconciled_groups: {mixed_reconciled_groups}",
        "",
        "Outputs",
        f"  reference: {outputs['reference']}",
        f"  mapping: {outputs['mapping']}",
        f"  loci_stats: {outputs['loci_stats']}",
        f"  run_stats: {outputs['run_stats']}",
        f"  audit_dir: {outputs['audit_dir']}",
    ]
    if not keep_intermediates:
        lines.append("  intermediates: cleaned on success")
    outpath.write_text("\n".join(lines) + "\n", encoding="utf-8")
    logger.info(f"wrote denovo run stats to {outpath}")


def vsearch_pairs(
    sname: str,
    r1: Path,
    r2: Path | None,
    outdir: Path,
    vsearch_binary: str,
    min_derep_size: int,
    min_merge_overlap: int,
    min_length: int,
    max_merge_diffs: int,
    allow_reverse_complement: bool,
    within_similarity: float,
    by_length: bool,
    threads: int,
    paired: bool = False,
) -> None:
    """Run the within-sample vsearch workflow for one sample."""
    unmerged_r1 = outdir / f"{sname}.unmerged_R1.fq"
    unmerged_r2 = outdir / f"{sname}.unmerged_R2.fq"
    merged = outdir / f"{sname}.merged.fa"
    joined = outdir / f"{sname}.joined.fa"
    derep = outdir / f"{sname}.derep.sizesorted.fa"
    consensus = outdir / f"{sname}.consensus.fa"
    clusters = outdir / f"{sname}.clusters.tsv"

    if paired:
        if r2 is None:
            raise IPyradError(f"Missing R2 FASTQ for paired sample: {sname}")
        cmd1 = [
            vsearch_binary, "--fastq_mergepairs", str(r1),
            "--reverse", str(r2),
            "--fastq_minovlen", str(min_merge_overlap),
            "--fastq_maxdiffs", str(max_merge_diffs),
            "--fastq_minlen", str(min_length),
            "--fastq_allowmergestagger",
            "--fasta_width", "0",
            "--fastqout_notmerged_fwd", str(unmerged_r1),
            "--fastqout_notmerged_rev", str(unmerged_r2),
            "--relabel", f"{sname};M",
            "--fastaout", str(merged),
        ]
        logger.debug(" ".join(cmd1))
        run_pipeline([cmd1])

        cmd1 = [
            vsearch_binary, "--fastq_join", str(unmerged_r1),
            "--reverse", str(unmerged_r2),
            "--join_padgap", "N" * 24,
            "--join_padgapq", "I" * 24,
            "--fasta_width", "0",
            "--relabel", f"{sname};J",
            "--fastaout", str(joined),
        ]
        logger.debug(" ".join(cmd1))
        run_pipeline([cmd1])

        cmd1 = ["cat", str(joined), str(merged)]
    else:
        cmd1 = [
            vsearch_binary, "--fastx_subsample", str(r1),
            "--sample_pct", "100",
            "--relabel", f"{sname};S",
            "--fastaout", str(joined),
        ]
        logger.debug(" ".join(cmd1))
        run_pipeline([cmd1])
        cmd1 = ["cat", str(joined)]

    cmd2 = [
        vsearch_binary, "--fastx_uniques", "-",
        "--minuniquesize", str(min_derep_size),
        "--strand", "both" if allow_reverse_complement else "plus",
        "--fasta_width", "0",
        "--sizeout",
        "--relabel_keep",
        "--fastaout", "-",
    ]
    cmd3 = [
        vsearch_binary, "--sortbylength" if by_length else "--sortbysize", "-",
        "--sizein",
        "--sizeout",
        "--fasta_width", "0",
        "--output", str(derep),
    ]
    logger.debug(f"{' '.join(cmd1)} | {' '.join(cmd2)} | {' '.join(cmd3)}")
    run_pipeline([cmd1, cmd2, cmd3])

    cmd1 = [
        vsearch_binary, "--cluster_fast" if by_length else "--cluster_size", str(derep),
        "--id", str(within_similarity),
        "--strand", "both" if allow_reverse_complement else "plus",
        "--maxaccepts", "1",
        "--maxrejects", "0",
        "--query_cov", "0.75",
        "--fasta_width", "0",
        "--qmask", "none",
        "--consout", str(consensus),
        "--uc", str(clusters),
        "--threads", str(threads),
    ]
    logger.debug(" ".join(cmd1))
    run_pipeline([cmd1])


def _write_cluster_sequence_fasta(
    summary_tsv: Path,
    out_fasta: Path,
) -> Path:
    """Write a clustering FASTA from spacer-stripped summary sequences."""
    with open(summary_tsv, "rt", encoding="utf-8", newline="") as infile, open(
        out_fasta,
        "wt",
        encoding="utf-8",
    ) as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        if reader.fieldnames is None or "seed" not in reader.fieldnames:
            raise RuntimeError("concat.summary.tsv is missing required column: seed")
        seq_field = "cluster_sequence" if "cluster_sequence" in reader.fieldnames else "consensus"
        for row in reader:
            seed = str(row["seed"])
            seq = str(row[seq_field]).upper()
            outfile.write(f">{seed}\n{seq}\n")
    return out_fasta


def vsearch_cluster_across(
    outdir: Path,
    summary_tsv: Path,
    across_similarity: float,
    threads: int,
    vsearch_binary: str,
) -> None:
    """Cluster all sample-level consensus sequences across samples."""
    cluster_table = outdir / "global_hits.uc.tsv"
    consensus_concat = outdir / "consensus.concat.fa"
    _write_cluster_sequence_fasta(summary_tsv, consensus_concat)

    cmd1 = [
        vsearch_binary, "--usearch_global", str(consensus_concat),
        "--db", str(consensus_concat),
        "--id", str(across_similarity),
        "--userout", str(cluster_table),
        "--userfields", "query+target+id+qstrand+qcov+ql+tl",
        "--maxaccepts", "0",
        "--maxrejects", "0",
        "--query_cov", "0.75",
        "--self",
        "--qmask", "none",
        "--notmatched", "/dev/null",
        "--fasta_width", "0",
        "--threads", str(threads),
    ]
    run_pipeline([cmd1])


def run_denovo(
    fastqs: list[Path],
    outdir: Path,
    within_similarity: float,
    across_similarity: float,
    min_derep_size: int,
    min_length: int,
    min_merge_overlap: int,
    max_merge_diffs: int,
    delim_str: str | None,
    delim_idx: int,
    allow_reverse_complement: bool,
    cores: int,
    threads: int,
    graph_splitter: str,
    no_alignment: bool,
    force: bool,
    keep_intermediates: bool,
    vsearch_binary: Path | None,
    mafft_binary: Path | None,
    log_level: str,
) -> None:
    """Run the denovo reference construction workflow."""
    _validate_runtime_args(
        within_similarity=within_similarity,
        across_similarity=across_similarity,
        min_derep_size=min_derep_size,
        min_length=min_length,
        min_merge_overlap=min_merge_overlap,
        max_merge_diffs=max_merge_diffs,
        cores=cores,
        threads=threads,
        delim_idx=delim_idx,
    )

    outdir = outdir.expanduser().absolute()
    workdir, outputs = _prepare_output_paths(outdir, force=force)
    vsearch_path = _resolve_binary(vsearch_binary, "vsearch")
    mafft_path = _resolve_binary(mafft_binary, "mafft")
    fastq_dict = get_name_to_fastq_dict(fastqs, delim_str, delim_idx)
    is_paired = _validate_fastq_layout(fastq_dict)
    workers = max(1, cores // threads)

    logger.info(f"paired mode: {'paired-end' if is_paired else 'single-end'}")
    logger.info(f"using vsearch binary: {vsearch_path}")
    logger.info(f"using mafft binary: {mafft_path}")

    msg = "Joining/merging pairs, dereplicating, and clustering" if is_paired else "Dereplicating and clustering"
    jobs: dict[str, tuple[Any, dict[str, Any]]] = {}
    for sname, fastq_tuple in fastq_dict.items():
        kwargs = dict(
            sname=sname,
            r1=fastq_tuple[0],
            r2=fastq_tuple[1],
            outdir=workdir,
            vsearch_binary=vsearch_path,
            min_derep_size=min_derep_size,
            min_length=min_length,
            min_merge_overlap=min_merge_overlap,
            max_merge_diffs=max_merge_diffs,
            allow_reverse_complement=allow_reverse_complement,
            within_similarity=within_similarity,
            by_length=True,
            threads=threads,
            paired=is_paired,
        )
        jobs[sname] = (vsearch_pairs, kwargs)

    run_with_pool(jobs, log_level, workers, msg=msg)

    for sname in fastq_dict:
        build_sample_summary(sname, workdir)
    concat_summaries(workdir)

    logger.info("Clustering consensus sequences across samples")
    vsearch_cluster_across(
        outdir=workdir,
        summary_tsv=workdir / "concat.summary.tsv",
        across_similarity=across_similarity,
        threads=threads,
        vsearch_binary=vsearch_path,
    )

    logger.info("Splitting global clusters and writing locus tables")
    mapping_df, stats_df = make_global_tables(
        workdir,
        graph_splitter=graph_splitter,
        cores=cores,
        log_level=log_level,
        across_similarity=across_similarity,
        mafft_binary=mafft_path,
    )

    if no_alignment:
        logger.info("Selecting longest locus representatives and writing denovo reference")
    else:
        logger.info("Aligning locus consensuses and writing denovo reference")
    alignment_summary = write_ordered_consensus_stream_to_file(
        mapping_tsv=outputs["mapping"],
        summary_tsv=workdir / "concat.summary.tsv",
        out_fa=outputs["reference"],
        mafft_binary=mafft_path,
        log_level=log_level,
        cores=cores,
        alignment_mode="none" if no_alignment else "mafft",
    )

    _write_denovo_stats(
        outputs["run_stats"],
        fastq_dict=fastq_dict,
        paired=is_paired,
        within_similarity=within_similarity,
        across_similarity=across_similarity,
        min_derep_size=min_derep_size,
        min_length=min_length,
        min_merge_overlap=min_merge_overlap,
        max_merge_diffs=max_merge_diffs,
        allow_reverse_complement=allow_reverse_complement,
        cores=cores,
        threads=threads,
        workers=workers,
        graph_splitter=graph_splitter,
        alignment_summary=alignment_summary,
        vsearch_binary=vsearch_path,
        mafft_binary=mafft_path,
        keep_intermediates=keep_intermediates,
        workdir=workdir,
        mapping_df=mapping_df,
        stats_df=stats_df,
        outputs=outputs,
    )

    if not keep_intermediates:
        shutil.rmtree(workdir)


if __name__ == "__main__":
    pass
