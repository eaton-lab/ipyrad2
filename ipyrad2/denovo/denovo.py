#!/usr/bin/env python

"""Dereplicate reads for denovo clustering.

PAIRED END
----------
1. merge or join pairs, concat, dereplicate, sort
2. cluster within
3.

"""

from typing import List
import sys
import shutil
# import tempfile
from pathlib import Path
from loguru import logger
from .cluster import build_sample_summary, concat_summaries
from .graph import make_global_tables
from .align import write_ordered_consensus_stream_to_file
# from .dereplicate import join_pairs_and_merge_derep
from ..utils.exceptions import IPyradError
from ..utils.names import get_name_to_fastq_dict
from ..utils.parallel import run_pipeline, run_with_pool

BIN = Path(sys.prefix) / "bin"
BIN_VSEARCH = str(BIN / "vsearch")


def vsearch_pairs(
    sname: str,
    r1: Path,
    r2: Path,
    outdir: Path,
    min_dereplication_size: int,
    min_merge_overlap: int,
    min_length: int,
    max_merge_diffs: int,
    strand_both: bool,
    similarity_threshold_within: float,
    by_length: bool,
    threads: int,
    paired: bool = False,
):
    unmerged_R1 = outdir / f"{sname}.unmerged_R1.fq"
    unmerged_R2 = outdir / f"{sname}.unmerged_R2.fq"
    merged = outdir / f"{sname}.merged.fa"
    joined = outdir / f"{sname}.joined.fa"
    derep = outdir / f"{sname}.derep.sizesorted.fa"
    consensus = outdir / f"{sname}.consensus.fa"
    clusters = outdir / f"{sname}.clusters.tsv"

    if paired:
        cmd1 = [
            BIN_VSEARCH,
            "--fastq_mergepairs", str(r1),
            "--reverse", str(r2),
            "--fastq_minovlen", str(min_merge_overlap),
            "--fastq_maxdiffs", str(max_merge_diffs),
            "--fastq_minlen", str(min_length),
            "--fastq_allowmergestagger",
            "--fasta_width", "0",
            "--fastqout_notmerged_fwd", str(unmerged_R1),
            "--fastqout_notmerged_rev", str(unmerged_R2),
            "--fasta_width", "0",
            "--relabel", f"{sname};M",
            "--fastaout", str(merged),
        ]
        logger.debug(" ".join(cmd1))
        run_pipeline([cmd1])

        cmd1 = [
            BIN_VSEARCH,
            "--fastq_join", str(unmerged_R1),
            "--reverse", str(unmerged_R2),
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
        # Relabel the R1 data for agreement with PE format in relabeling the reads
        cmd1 = [
            BIN_VSEARCH,
            "--fastx_subsample", str(r1),
            "--sample_pct", "100",
            "--relabel", f"{sname};S",
            "--fastaout", str(joined),
        ]
        logger.debug(" ".join(cmd1))
        run_pipeline([cmd1])

        cmd1 = ["cat", str(joined)]

    cmd2 = [
        BIN_VSEARCH,
        "--fastx_uniques", "-",
        "--minuniquesize", str(min_dereplication_size),
        "--strand", "both" if strand_both else "plus",
        "--fasta_width", "0",
        "--sizeout",
        "--relabel_keep",
        "--fastaout", "-",
    ]
    cmd3 = [
        BIN_VSEARCH,
        "--sortbylength" if by_length else "--sortbysize",
        "-",
        "--sizein",
        "--sizeout",
        "--fasta_width", "0",
        "--output", str(derep),
    ]
    logger.debug(f"{' '.join(cmd1)} | {' '.join(cmd2)} | {' '.join(cmd3)}")
    run_pipeline([cmd1, cmd2, cmd3])

    cmd1 = [
        BIN_VSEARCH,
        "--cluster_fast" if by_length else "--cluster_size", str(derep),
        "--id", str(similarity_threshold_within),
        "--strand", "both" if strand_both else "plus",
        "--maxaccepts", "1",
        "--maxrejects", "0",
        # "--minsl", "0.75",
        "--query_cov", "0.75",
        "--fasta_width", "0",
        "--qmask", "none",
        "--consout", str(consensus),
        "--uc", str(clusters),
        "--threads", str(threads),
    ]
    logger.debug(" ".join(cmd1))
    run_pipeline([cmd1])


def vsearch_cluster_across(
    outdir: Path,
    similarity_threshold_across: float,
    threads: int,
):
    cluster_table = outdir / "global_hits.uc.tsv"
    consensus_concat = outdir / "consensus.concat.fa"

    # create concat.consensus
    cmd1 = ["cat"] + sorted(outdir.glob("*.consensus.fa"))
    run_pipeline([cmd1], consensus_concat)

    # all-vs-all search
    # (keeps edges; cheaper than full clustering if you want graph control)
    # NB: cannot pipe concat into this.
    cmd1 = [
        BIN_VSEARCH,
        "--usearch_global", str(consensus_concat),
        "--db", str(consensus_concat),
        "--id", str(similarity_threshold_across),
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
    fastqs: List[Path],
    outdir: Path,
    similarity_threshold_within: float,
    similarity_threshold_across: float,
    min_dereplication_size: int,
    min_length: int,
    min_merge_overlap: int,
    max_merge_diffs: int,
    delim_str: str | None,
    delim_idx: int,
    strand_both: bool,
    cores: int,
    threads: int,
    force: bool,
    log_level: str,
):
    """..."""
    outdir = outdir.expanduser().absolute()
    outdir.mkdir(exist_ok=True)
    denovo_reference = outdir / "denovo_reference.fa"
    tmpdir = outdir / "tmpdir"
    logger.warning(tmpdir)
    fastq_dict = get_name_to_fastq_dict(fastqs, delim_str, delim_idx)
    is_paired = list(fastq_dict.values())[0][1] is not None
    workers = max(1, cores // threads)

    # -------------------------------------------
    if tmpdir.exists() or denovo_reference.exists():
        if not force:
            raise IPyradError("denovo reference results exist in outdir. Use --force to overwrite.")
        else:
            shutil.rmtree(tmpdir)
    tmpdir.mkdir(exist_ok=True)

    # -------------------------------------------
    msg = "Joining/merging pairs, d" if is_paired else "D"
    logger.info(f"{msg}ereplicating and clustering")
    jobs = {}
    for sname, fastq_tuple in fastq_dict.items():
        kwargs=dict(
            sname=sname,
            r1=fastq_tuple[0],
            r2=fastq_tuple[1],
            outdir=tmpdir,
            min_dereplication_size=min_dereplication_size,
            min_length=min_length,
            min_merge_overlap=min_merge_overlap,
            max_merge_diffs=max_merge_diffs,
            strand_both=strand_both,
            similarity_threshold_within=similarity_threshold_within,
            by_length=True,
            threads=threads,
            paired=is_paired,
        )
        jobs[sname] = (vsearch_pairs, kwargs)
    run_with_pool(jobs, log_level, workers)

    # write sample summary TSVs
    for sname in fastq_dict:
        build_sample_summary(sname, tmpdir)
    concat_summaries(tmpdir)

    logger.info("Clustering consensus sequences across samples")
    vsearch_cluster_across(tmpdir, similarity_threshold_across, threads)

    logger.info("Splitting clusters and writing mapping table")
    mapping_tsv, summary_tsv = make_global_tables(tmpdir)

    logger.info("Aligning and writing denovo consensus reference")
    write_ordered_consensus_stream_to_file(outdir, log_level)

    # -------------------------------------------




if __name__ == "__main__":
    pass
