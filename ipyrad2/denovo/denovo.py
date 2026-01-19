#!/usr/bin/env python

"""Dereplicate reads for denovo clustering.

PAIRED END
----------
1. merge or join pairs, concat, dereplicate, sort
2. cluster within
3.

"""

from typing import List, Dict, Tuple
import numpy as np
import itertools
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
from ..utils.pops import parse_pops_file, parse_imap

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
    msg = "All-by-all clustering"
    run_pipeline([cmd1], msg=msg, quiet=False)


def run_denovo(
    fastqs: List[Path],
    outdir: Path,
    imap: Path | None,
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
    logger.debug(tmpdir)
    fastq_dict = get_name_to_fastq_dict(fastqs, delim_str, delim_idx)
    is_paired = list(fastq_dict.values())[0][1] is not None
    workers = max(1, cores // threads)

    # -------------------------------------------
    # Clean up stale files from previous denovo assembly
    if tmpdir.exists() or denovo_reference.exists():
        if not force:
            raise IPyradError("denovo reference results exist in outdir. Use --force to overwrite.")
        else:
            # Clean up stale bwa-mem2 index files
            suffs = [".pac", ".ann", ".amb", ".0123", ".bwt.2bit.64", ".fai"]  # bwa-mem2
            # don't use Path.with_suffix here b/c '.fa.ann' double suffix is messy.
            paths = [denovo_reference.with_suffix(denovo_reference.suffix + i) for i in suffs]
            for i in paths:
                try:
                    i.unlink()
                except FileNotFoundError:
                    pass
            # Clean up the tmpdir
            shutil.rmtree(tmpdir)
    tmpdir.mkdir(exist_ok=True)

    # -------------------------------------------
    # Use imap for subsetting samples for building denovo reference
    # Potentially return a subset of fastqs determined either by the contents
    # of imap or by randomly selecting 10 samples total
    fastq_dict = _subset_fastqs(imap, fastq_dict)

    # -------------------------------------------
    # vsearch w/in samples (derep/cluster)
    msg = "Joining/merging pairs, d" if is_paired else "D"
    msg = f"{msg}ereplicating and clustering"
    logger.info(msg)
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
    run_with_pool(jobs, log_level, workers, msg=msg)

    logger.success("Building summary tables")
    # write sample summary TSVs
    for sname in fastq_dict:
        build_sample_summary(sname, tmpdir)
    concat_summaries(tmpdir)

    #TODO: Add some logging messages here so people can see progress
    logger.info("Clustering consensus sequences across samples")
    vsearch_cluster_across(tmpdir, similarity_threshold_across, threads)

    logger.info("Splitting clusters and writing mapping table")
    mapping_tsv, summary_tsv = make_global_tables(tmpdir)

    logger.info("Aligning and writing denovo consensus reference")
    write_ordered_consensus_stream_to_file(outdir, log_level)

    # -------------------------------------------


def _subset_fastqs(imap: Path | None,
    fastq_dict: Dict[str, Tuple[Path, Path | None]],
    nsamples: int = 10,
    seed: int | None = None):
    """
    """
    if not seed:
        seed = np.random.randint(0, 1e9)
    rng = np.random.default_rng(np.random.SeedSequence(seed))

    # -------------------------------------------
    # Get imap/minmap for subsetting samples for building denovo reference
    # parse_pops_file is responsible for validating that the minmap pops and
    # imap pops are identical.
    if imap is None:
        imap = {'all': list(fastq_dict.keys())}
        minmap = {'all': nsamples}
    else:
        if not imap.exists():
            raise IPyradError(f"imap file does not exists: {imap}")
        minmap = {}
        try:
            # Favor ipyrad style imap file including sample/pop mapping
            # and trailing minmap line (# pop1:10 Pop2:5 ...)
            imap, minmap = parse_pops_file(imap)
        except IPyradError as e:
            logger.warning(e)
            logger.info("imap file doesn't include minmap info, parsing standard imap file format.")
            imap = parse_imap(imap)
        # Validate names in imap and fastq_dict agree
        # raise error if any imap sample names not in database names
        imapset = set(itertools.chain(*imap.values()))
        badnames = imapset.difference(fastq_dict.keys())
        if badnames:
            raise IPyradError(
                f"Samples {badnames} are not in fastqs list: {fastq_dict.keys()}")

        # Enforce at least one sample per population
        if not minmap:
            if len(imap) > nsamples:
                samples_per_pop = 1
            else:
                # If the # of pops is smaller than nsamples we do a little fudging
                # to get the target number of samples per population, so the sum
                # of samples_per_pop * len(imap) will sometimes be slightly higher
                # or lower than the passed in (hopeful) nsamples value
                samples_per_pop = round(nsamples/len(imap))

            # Have to retain _at_ least the number of samples available, and at
            # most the number determined by dividing nsamples by npops.
            minmap = {pop:min(samples_per_pop, len(samps)) for pop, samps in imap.items()}
    logger.warning(f"sample coverage minmap = {minmap}")
    if sum(minmap.values()) > nsamples:
        logger.error(f"imap file is selecting more than {nsamples} samples. Time to "
            "construct the pseudo-reference increases with increasing numbers of samples.")

    tmp_fastq_dict = {}
    for pop, samps in imap.items():
        # Constrain the number to be sampled from a given population when replace=False
        max_samps = len(imap[pop])
        samps = rng.choice(samps, min(minmap[pop], max_samps), replace=False)
        # Grab the fastq Paths for retained samples
        for samp in samps:
            tmp_fastq_dict[str(samp)] = fastq_dict[samp]

    logger.info("Subsetting populations for construction of pseudo-reference sequence.")
    logger.info(f"Retaining samples: {list(tmp_fastq_dict.keys())}")
    logger.debug(f"Retaining: {tmp_fastq_dict}")

    return tmp_fastq_dict


if __name__ == "__main__":
    pass
