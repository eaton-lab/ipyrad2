#!/usr/bin/env python

"""Map, filter, sort, mark BAMs.

Example
-------
map --fastqs DATA --ref REF -out MAP
"""

from typing import Tuple, Dict
import os
import sys
import json
from collections import defaultdict
from pathlib import Path
from loguru import logger
import numpy as np
import pandas as pd
from ..utils.exceptions import IPyradError
from ..utils.names import get_name_to_fastq_dict
from ..utils.parallel import run_pipeline, run_with_pool


BIN = Path(sys.prefix) / "bin"
BIN_BWA = str(BIN / "bwa-mem2")
BIN_SAMTOOLS = str(BIN / "samtools")  # indexing


def map_filter_sort_mark_pairs(sname: str, fastqs: Tuple[Path, Path], reference: Path, mark_dups_by_umis: bool, min_map_q: int, max_soft_clip: int, max_edit_dist: int, outdir: Path, threads: int) -> Path:
    """Map reads to the reference to get a sorted bam.

    This pipeline is for PE data w duplication marking information,
    meaning either WGS data (coordinates), or ddRAD w/ i5 UMIs stored
    to names (see option in ipyrad2 trim to do this.)

    In the last step w apply the following filters:
        - q = min mapping score
        - (flag&4)==0) && ((flag&8)==0) = both reads mapped
        - (rnext=="=" || rnext==rname) = same scaff
        - (tlen>=-2000 && tlen<=2000) = not more than >2 kb apart.

    NB: some samtools lack abs() so we use (tlen >= -X && tlen <= X)
    NB: -q can filter one pair and leave the other. That's ok for all downstream steps.
    NB: this file will be used in both variants.py and beds.py in the next step.
    """
    # paths
    # logger.debug(f"starting mapping of {sname}: {fastqs}")
    out_bam = outdir / f"{sname}.filtered.bam"
    out_bam_tmp = outdir / f"{sname}.filtered.bam.tmp"
    bam_namesort = outdir / "tmpdir" / f"{sname}.tmp.namesort.bam"
    bam_fixmate = outdir / "tmpdir" / f"{sname}.tmp.fixmate.bam"
    bam_coordsort = outdir / "tmpdir" / f"{sname}.tmp.coordsort.bam"
    bam_markdup = outdir / "tmpdir" / f"{sname}.tmp.markdup.bam"
    tmp_prefix = outdir / "tmpdir" / f"{sname}.tmp.pre"
    tmp_stats1 = outdir / "tmpdir" / f"{sname}.tmp.stats1.json"
    tmp_stats2 = outdir / "tmpdir" / f"{sname}.tmp.stats2.json"
    tmp_stats3 = outdir / "tmpdir" / f"{sname}.tmp.stats3.json"
    tmp_stats4 = outdir / "tmpdir" / f"{sname}.tmp.stats4.json"
    tmp_stats5 = outdir / "tmpdir" / f"{sname}.tmp.stats5.json"
    tmp_stats_dups = outdir / "tmpdir" / f"{sname}.tmp.stats_dups.txt"

    # Split threads between BWA and samtools
    bwa_threads = max(1, int(threads * 0.75))
    sort_threads = max(1, threads - bwa_threads)

    # mapping command
    # additional options to consider for toggle
    # -k 15   # reduce kmer size for small denovo loci
    # -L 1,1  # reduce penalty for clipping
    cmd1 = [
        BIN_BWA, "mem",
        "-Y",                     # mark supplementary with soft-clip. We exclude supplementals anyways.
        "-T", "20",               # minimum alignmnet score to output (default=30). This is different from MAPQ in samtools.
        "-R", f"@RG\\tID:{sname}\\tSM:{sname}",  # store sample name.
        "-K", "50000000",         # stable nbases chunk size. Improves repeatability.
        "-v", "1",                # less verbose.
        "-t", str(bwa_threads),
        str(reference),
        str(fastqs[0]),
        str(fastqs[1]),
    ]

    # drop secondary + supplemental + QCFail is OK here.
    cmd2 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-F", "0x100",    # secondary
        "-F", "0x200",    # qcfail
        "-F", "0x800",    # supplemental
        # "-F", "0xB00", # 0x100, 0x200, 0x800
        "--save-counts", str(tmp_stats1),
        "-o", "-",
    ]

    # name sort is required for fixmate
    cmd3 = [
        BIN_SAMTOOLS, "sort",
        "-n",
        "-m", "256M",
        "-T", str(tmp_prefix),
        "-@", str(sort_threads),
        "-o", str(bam_namesort),
        "-",
    ]
    run_pipeline([cmd1, cmd2, cmd3])

    # fixmate checks and updates tags about pairing
    cmd4 = [
        BIN_SAMTOOLS, "fixmate",
        "-m",          # add mate score tags
        "-@", str(threads),
        str(bam_namesort),
        str(bam_fixmate),
    ]
    run_pipeline([cmd4])
    bam_namesort.unlink()

    # coordinate sort command
    cmd5 = [
        BIN_SAMTOOLS, "sort",
        "-m", "256M",                # tune per-thread memory
        "-T", str(tmp_prefix),
        "-@", str(threads),
        "-o", str(bam_coordsort),
        str(bam_fixmate),
    ]
    run_pipeline([cmd5])
    bam_fixmate.unlink()

    # mark dups in coordinate sorted fixmate bams
    cmd6 = [
        BIN_SAMTOOLS, "markdup",
        "-r",
        "-T", str(tmp_prefix),
        "-s",                       # write stats
        "-f", str(tmp_stats_dups),  # write stats to this file
        "-@", str(threads),
        "--write-index",
        str(bam_coordsort),
        str(bam_markdup),
    ]
    if mark_dups_by_umis:
        cmd6.extend(["--barcode-rgx", "UMI_([ACGTN]+)"])
    run_pipeline([cmd6])
    bam_coordsort.unlink()

    # filter for pairing and quality. Run two views piped to split stats out
    cmd7 = [
        BIN_SAMTOOLS, "view",
        "-q", str(min_map_q),        # apply min map q
        "--save-counts", str(tmp_stats2),
        "-b", "-u",
        "-o", "-",
        str(bam_markdup),
    ]

    # filter to require paired and constrain distance between pairs
    cmd8 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-e", '((flag&4)==0) && ((flag&8)==0) && (rnext=="=" || rnext==rname) && (tlen>=-2000 && tlen<=2000)',
        "--save-counts", str(tmp_stats3),
        "-o", "-",
    ]

    # allow at most this many soft-clipped bases
    cmd9 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-e", f"sclen <= {max_soft_clip}",
        "--save-counts", str(tmp_stats4),
        "-o", "-",
    ]

    # allow at most this many changes relative to the reference
    cmd10 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-e", f"[NM] < {max_edit_dist}",
        "--save-counts", str(tmp_stats5),
        "-o", "-",
    ]
    run_pipeline([cmd7, cmd8, cmd9, cmd10], out_bam_tmp)
    logger.debug(f"finished mapping: {sname}")

    # rename files to finished name (protects for restarting w/o -f)
    os.replace(out_bam_tmp, out_bam)
    cmd = [BIN_SAMTOOLS, "index", "-c", str(out_bam)]
    run_pipeline([cmd])

    # remove tmp files
    bam_markdup.unlink()
    bam_markdup.with_suffix(bam_markdup.suffix + '.csi').unlink()
    return out_bam


def map_filter_sort_pairs(sname: str, fastqs: Tuple[Path, Path], reference: Path, min_map_q: int, max_soft_clip: int, max_edit_dist: int, outdir: Path, threads: int, **kwargs) -> Path:
    """Map reads to the reference to get a sorted bam.

    This pipeline is for PE data w/o duplication marking information.
    """
    # paths
    out_bam = outdir / f"{sname}.filtered.bam"
    out_bam_tmp = outdir / f"{sname}.filtered.bam.tmp"
    tmp_prefix = outdir / "tmpdir" / f"{sname}.tmp.pre"
    tmp_stats1 = outdir / "tmpdir" / f"{sname}.tmp.stats1.json"
    tmp_stats2 = outdir / "tmpdir" / f"{sname}.tmp.stats2.json"
    tmp_stats3 = outdir / "tmpdir" / f"{sname}.tmp.stats3.json"
    tmp_stats4 = outdir / "tmpdir" / f"{sname}.tmp.stats4.json"
    tmp_stats5 = outdir / "tmpdir" / f"{sname}.tmp.stats5.json"

    # Split threads between BWA and samtools
    bwa_threads = max(1, int(threads * 0.75))
    sort_threads = max(1, threads - bwa_threads)

    # mapping command
    cmd1 = [
        BIN_BWA, "mem",
        "-Y",                     # soft-clip supplementary. Wouldn't hurt, but not necessary.
        "-T", "20",               # minimum score to output (default=30). Keep lower scores for now, we filter on MAPQ later.
        "-R", f"@RG\\tID:{sname}\\tSM:{sname}",  # store sample names in bam.
        "-K", "50000000",         # stable chunk size. Improves repeatability.
        "-v", "1",                # less verbose.
        "-t", str(bwa_threads),
        str(reference),
        str(fastqs[0]),
        str(fastqs[1]),
    ]

    # drop secondary + supplemental + QCFail is OK here.
    cmd2 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-F", "0x100",    # secondary
        "-F", "0x200",    # qcfail
        "-F", "0x800",    # supplemental
        # "-F", "0xB00", # 0x100, 0x200, 0x800
        "--save-counts", str(tmp_stats1),
        "-o", "-",
    ]

    # filter for quality
    cmd3 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-q", str(min_map_q),        # apply min map q
        "--save-counts", str(tmp_stats2),
        "-o", "-",
    ]

    # filter to require paired and constrain distance between pairs
    cmd4 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-e", '((flag&4)==0) && ((flag&8)==0) && (rnext=="=" || rnext==rname) && (tlen>=-2000 && tlen<=2000)',
        "--save-counts", str(tmp_stats3),
        "-o", "-",
    ]

    # allow at most this many soft-clipped bases
    cmd5 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-e", f"sclen <= {max_soft_clip}",
        "--save-counts", str(tmp_stats4),
        "-o", "-",
    ]

    # allow at most this many changes relative to the reference
    cmd6 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-e", f"[NM] < {max_edit_dist}",
        "--save-counts", str(tmp_stats5),
        "-o", "-",
    ]

    # coordinate sort
    cmd7 = [
        BIN_SAMTOOLS, "sort",
        "-m", "256M",                # per-thread memory
        "-T", str(tmp_prefix),
        "-@", str(sort_threads),
        "--write-index",
        "-O", "bam",                 # explicity ask for bam format to be safe
        "-o", str(out_bam_tmp),
        "-",
    ]
    cmds = [cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7]
    run_pipeline(cmds)

    # rename files to finished name (protects for restarting w/o -f)
    os.replace(out_bam_tmp, out_bam)
    os.replace(out_bam_tmp.with_suffix(out_bam_tmp.suffix + ".csi"), out_bam.with_suffix(out_bam.suffix + ".csi"))
    logger.debug(f"finished mapping: {sname}")
    return out_bam


def map_filter_sort_single(sname: str, fastqs: Tuple[Path, Path], reference: Path, min_map_q: int, max_soft_clip: int, max_edit_dist: int, outdir: Path, threads: int, **kwargs) -> Path:
    """Map reads to the reference to get a sorted bam.

    This pipeline is for SE data w/o duplication marking information.
    """
    # paths
    out_bam = outdir / f"{sname}.filtered.bam"
    tmp_prefix = outdir / "tmpdir" / f"{sname}.tmp.pre"
    tmp_stats1 = outdir / "tmpdir" / f"{sname}.tmp.stats1.json"
    tmp_stats2 = outdir / "tmpdir" / f"{sname}.tmp.stats2.json"
    tmp_stats3 = outdir / "tmpdir" / f"{sname}.tmp.stats3.json"
    tmp_stats4 = outdir / "tmpdir" / f"{sname}.tmp.stats4.json"
    tmp_stats5 = outdir / "tmpdir" / f"{sname}.tmp.stats5.json"

    # Split threads between BWA and samtools
    bwa_threads = max(1, int(threads * 0.75))
    sort_threads = max(1, threads - bwa_threads)

    # mapping command
    cmd1 = [
        BIN_BWA, "mem",
        # "-Y",                     # soft-clip supplementary. Wouldn't hurt, but not necessary.
        "-T", "20",               # minimum score to output (default=30). Keep lower scores for now, we filter on MAPQ later.
        "-R", f"@RG\\tID:{sname}\\tSM:{sname}",  # store sample names in bam.
        "-K", "50000000",         # stable chunk size. Improves repeatability.
        "-v", "1",                # less verbose.
        "-t", str(bwa_threads),
        str(reference),
        str(fastqs[0]),
    ]

    # drop secondary + supplemental + QCFail is OK here.
    cmd2 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-F", "0x100",    # secondary
        "-F", "0x200",    # qcfail
        "-F", "0x800",    # supplemental
        # "-F", "0xB00", # 0x100, 0x200, 0x800
        "--save-counts", str(tmp_stats1),
        "-o", "-",
    ]
    cmd3 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-q", str(min_map_q),        # apply min map q
        "--save-counts", str(tmp_stats2),
        "-o", "-",
    ]

    cmd4 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        # TODO: Document what this is doing. Lol.
        "-e", '((flag&4)==0) && ((flag&8)==0)',
        "--save-counts", str(tmp_stats3),
        "-o", "-",
    ]

    # allow at most this many soft-clipped bases
    cmd5 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-e", f"sclen <= {max_soft_clip}",
        "--save-counts", str(tmp_stats4),
        "-o", "-",
    ]

    # allow at most this many changes relative to the reference
    cmd6 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-e", f"[NM] < {max_edit_dist}",
        "--save-counts", str(tmp_stats5),
        "-o", "-",
    ]

    # coordinate sort
    cmd7 = [
        BIN_SAMTOOLS, "sort",
        "-m", "256M",                # per-thread memory
        "-T", str(tmp_prefix),
        "-@", str(sort_threads),
        "--write-index",
        "-O", "bam",                 # explicity ask for bam format to be safe
        "-o", str(out_bam),
        "-",
    ]
    cmds = [cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7]
    run_pipeline(cmds)
    logger.debug(f"finished mapping: {sname}")
    return out_bam


def concat_tech_reps_into_tmpdir(imap: Path, tmpdir: Path, fastq_dict: Dict[str, Tuple[Path, Path]]) -> Dict[str, Path]:
    """Return fastq_dict revised by an imap file.

    The imap file is a whitespace delimited table of sample names to
    include in the mapping step, and optionally a second column with
    new names. The names in the imap file should match the names that
    will be parsed by the CLI name-delimiter, i.e., not the file names
    with .trimmed.fastq.gz, etc.

    The imap file can serve three functions here: (1) to subsample
    which samples to process from a glob of potentially many files;
    (2) to rename samples, in which case they are written to
    MAPPED/{newname}.filtered.bam; and (3) to combine technical
    replicates by assigning two or more samples to the same name.

    IMAP file examples
    ------------------
    # one column: only these samples will be processed from glob.
    sample_1A
    sample_2A
    sample_3A

    # two columns: samples with a column 2 value will be renamed to it.
    sample_1A   sample_A
    sample_2B   sample_B
    sample_3C   sample_C

    # two columns: samples with the same column 2 value will be combined and renamed.
    sample_1A   sample_A
    sample_2A   sample_A
    sample_3B   sample_B
    """
    # fill pdict and warn if names don't match any samples.
    snames = set(fastq_dict)
    pname_to_path_tuples = defaultdict(list)
    pname_to_snamelist = defaultdict(list)
    warn_list = []

    # parse the population file.
    with open(imap, 'r') as infile:
        for line in infile:
            if not line.strip():
                continue
            sname, *data = line.strip().split()
            if data:
                pname = data[0]
            else:
                pname = sname

            # if imap sname is in fastq snames then store it and its file paths
            if sname in snames:
                pname_to_path_tuples[pname].append(fastq_dict.pop(sname))
                pname_to_snamelist[pname].append(sname)
            else:
                warn_list.append(sname)

    # warn user about mismatches names and bail if no names were selected
    if warn_list:
        logger.warning(f"One or more names in imap file did not match sample names and will be skipped: {' '.join(warn_list)}")
    if not pname_to_snamelist:
        raise IPyradError("No samples in imap file match parsed sample names. Revise imap file or name parsing args.")

    # report to logger
    logger.info(f"subselecting, renaming, or merging {len(fastq_dict)} samples into {len(pname_to_snamelist)} samples according to: {imap.name}")
    maxlen = max(len(i) for i in pname_to_snamelist)
    for pname, tups in pname_to_path_tuples.items():
        snames = pname_to_snamelist[pname]
        logger.info(f"{pname}{' ' * (maxlen - len(pname))} <- {' + '.join(snames)}")

        # renaming, do not run pipe
        if len(snames) == 1:
            fastq_dict[pname] = tups[0]

        # concating, run pipe
        else:
            out1 = tmpdir / f"{pname}.tmp.R1.fastq.gz"
            cmd = ["cat"] + [str(i[0]) for i in tups]
            run_pipeline([cmd], out1)

            if tups[0][1] is not None:
                out2 = tmpdir / f"{pname}.tmp.R2.fastq.gz"
                cmd = ["cat"] + [str(i[1]) for i in tups]
                run_pipeline([cmd], out2)
                fastq_dict[pname] = (out1, out2)
            else:
                fastq_dict[pname] = (out1, None)
    return fastq_dict


def count_mapped_reads(sname: str, outdir: Path) -> int:
    """Return the number of mapped reads in the filtered/sorted bam.

    Note that for PE data this is the still nreads, so divide by 2 to
    get the n read pairs.
    """
    # tmp.stats1.json contains filter to remove secondary, supplementary, qcfail
    s1 = outdir / "tmpdir" / f"{sname}.tmp.stats1.json"
    with s1.open('r') as indata:
        d1 = json.loads(indata.read())

    # tmp.stats_dups.json (if present) contains filter to remove dups
    sd = outdir / "tmpdir" / f"{sname}.tmp.stats_dups.txt"
    if sd.exists():
        with sd.open('r') as indata:
            a = d1['records_filter_accepted']
            b = [int(i.split()[-1]) for i in indata.readlines() if i.startswith("DUPLICATE TOTAL")][0]
            dd = {
                'records_processed': a,
                'records_filter_accepted': a - b,
                'records_filter_rejected': b,
            }
    else:
        dd = {"records_filter_rejected": 0}

    # tmp.stats2.json contains filter to remove mapping quality < q
    s2 = outdir / "tmpdir" / f"{sname}.tmp.stats2.json"
    with s2.open('r') as indata:
        d2 = json.loads(indata.read())

    # tmp stats3.json contains filter to remove bad pairs.
    s3 = outdir / "tmpdir" / f"{sname}.tmp.stats3.json"
    with s3.open('r') as indata:
        d3 = json.loads(indata.read())

    # tmp stats4.json contains filter to remove max soft clip
    s4 = outdir / "tmpdir" / f"{sname}.tmp.stats4.json"
    with s4.open('r') as indata:
        d4 = json.loads(indata.read())

    # tmp stats5.json contains filter to remove max edit dist
    s5 = outdir / "tmpdir" / f"{sname}.tmp.stats5.json"
    with s5.open('r') as indata:
        d5 = json.loads(indata.read())

    # get mean, std of mapq
    bam_file = outdir / f"{sname}.filtered.bam"
    cmd1 = [BIN_SAMTOOLS, "view", str(bam_file)]
    cmd2 = ["cut", "-f", "5"]
    _, out, _ = run_pipeline([cmd1, cmd2])
    out = np.array(list(map(int, out.decode().strip().split())))
    if out.size:
        mean_mapq = np.mean(out)
        median_mapq = np.median(out)
        stdev_mapq = np.std(out)
    else:
        mean_mapq = float('nan')
        median_mapq = float('nan')
        stdev_mapq = float('nan')

    data = {
        "nreads_processed": d1["records_processed"],
        "nreads_passed_filters": d3["records_filter_accepted"],
        "proportion_passed_filters": float(d5["records_filter_accepted"] / d1["records_processed"]),
        "filtered_by_not_primary": d1["records_filter_rejected"],
        "filtered_by_duplicates": dd["records_filter_rejected"],
        "filtered_by_min_mapq": d2["records_filter_rejected"],
        "filtered_by_bad_pairing": d3["records_filter_rejected"],
        "filtered_by_max_soft_clip": d4["records_filter_rejected"],
        "filtered_by_max_edit_dist": d5["records_filter_rejected"],
        "mapq_mean_after_filters": float(mean_mapq),
        "mapq_median_after_filters": float(median_mapq),
        "mapq_stdev_after_filters": float(stdev_mapq),
    }
    s1.unlink()
    s2.unlink()
    s3.unlink()
    s4.unlink()
    s5.unlink()
    if sd.exists():
        sd.unlink()
    return data


def index_ref_with_bwa(reference: Path) -> None:
    """Index the reference sequence, unless it already exists
    """
    # check that ref file exists
    if not reference.exists():
        raise IPyradError(f"reference path {reference} does not exist.")

    # If reference sequence already exists then bail out of this func
    suffs = [".pac", ".ann", ".amb", ".0123", ".bwt.2bit.64"]  # bwa-mem2
    # don't use Path.with_suffix here b/c '.fa.ann' double suffix is messy.
    paths = [reference.with_suffix(reference.suffix + i) for i in suffs]
    if all(i.exists() for i in paths):
        logger.debug(f"reference is already bwa indexed: {reference}")
        return

    # check that location of reference file is writable before trying to index.
    if not os.access(reference.parent, os.W_OK | os.X_OK):
        raise IPyradError("cannot index reference because you do not have write access to its directory.")

    # bwa index <reference_file>
    logger.info(f"indexing reference: {reference.name}")
    cmd = [str(BIN_BWA), "index", str(reference)]
    logger.debug(f"cmd: {' '.join(cmd)}")
    run_pipeline([cmd])


def run_mapper(
    fastqs: Tuple[Path, Path],
    outdir: Path,
    reference: Path,
    imap: Path,
    min_map_q: int,
    max_soft_clip: int,
    max_edit_dist: int,
    cores: int,
    threads: int,
    force: bool,
    mark_dups_by_coords: bool,
    mark_dups_by_umis: bool,
    delim_str: str | None,
    delim_idx: int,
    log_level: str,
):
    # ------------------------------------------------------------
    # check reference and outdir paths
    reference = reference.expanduser().resolve()
    outdir = outdir.expanduser().resolve()
    tmpdir = outdir / "tmpdir"
    tmpdir.mkdir(exist_ok=True, parents=True)

    # parse dict of {name: (r1, r2)}
    fastq_dict = get_name_to_fastq_dict(fastqs, delim_str, delim_idx)

    # check outdir for existing, warn if any exist, and skip or overwrite if (-f)
    results = [sname for sname in fastq_dict if (outdir / f"{sname}.filtered.bam").exists()]
    if any(results) and not force:
        logger.warning(f"skipping {len(results)}/{len(fastq_dict)} samples that already have results (.bam) in outdir. Use --force to instead overwrite.")
    for sname in results:
        if force:
            rfile = outdir / f"{sname}.filtered.bam"
            rfile.unlink()
            logger.debug(f"removing existing bam file: {rfile}")
        else:
            fastq_dict.pop(sname)
    if not fastq_dict:
        logger.info("all samples are completed.")
        raise SystemExit(0)

    # TODO: auto-tune threads based on nsamples and ncores?
    # 12 cores, 1 sample; threads=12
    # 12 cores, 12 samples; threads=4
    # 24 cores, 2 samples;
    # threads = max(cores, len(fastq_dict))
    workers = max(1, cores // threads)

    # check mark_dups suitability
    if mark_dups_by_coords or mark_dups_by_umis:
        if mark_dups_by_coords and mark_dups_by_umis:
            raise IPyradError("you cannot select both mark_dups_by_coords and mark_dups_by_umis.")
        if list(fastq_dict.values())[0][1] is None:
            raise IPyradError("Data do not appear to be paired. Cannot mark duplicates for SE data.")
        # TODO: check for valid rather than just warn.
        if mark_dups_by_coords:
            logger.warning("marking PCR duplicates by coordinates. Be sure this run includes only WGS samples, not RAD")
        if mark_dups_by_umis:
            logger.warning("marking PCR duplicates by i5 UMIs. Be sure you ran `ipyrad2 trim` with `-U` to store i5 tags for these samples")

    # store whether reads are paired.
    # TODO: this can raise an error when glob catches extras (e.g., not .fq, .gz, etc). Report that their glob might be bad?
    is_paired = False
    try:
        pairs_exist = [(r1.exists() and r2.exists()) for (r1, r2) in fastq_dict.values()]
        if any(pairs_exist):
            if all(pairs_exist):
                is_paired = True
            else:
                raise IPyradError("some but not all files have R1 and R2 pairs. Check inputs.")
    except AttributeError:
        # Single-end: If r2 == None then r2.exists() raises Attribute Error
        is_paired = False

    # index the reference
    index_ref_with_bwa(reference)

    # if tech-reps present concat files into tmpdir and update fastq_dict
    if imap is not None:
        fastq_dict = concat_tech_reps_into_tmpdir(imap, tmpdir, fastq_dict)

    # run map, filter, sort
    logger.info(f"mapping and filtering {len(fastq_dict)} inputs to bams in {outdir}")
    logger.info(f"using up to {cores} cores (up to {workers} multi-threaded jobs using {threads} threads)")
    # logger.info(f"running up to {workers} parallel jobs each using up to {threads} threads")
    jobs = {}
    for sname, fastq_tuple in fastq_dict.items():
        kwargs = dict(
            sname=sname,
            fastqs=fastq_tuple,
            outdir=outdir,
            reference=reference,
            min_map_q=min_map_q,
            max_soft_clip=max_soft_clip,
            max_edit_dist=max_edit_dist,
            mark_dups_by_umis=mark_dups_by_umis,
            threads=threads,
        )
        if is_paired:
            if mark_dups_by_coords or mark_dups_by_umis:
                jobs[sname] = (map_filter_sort_mark_pairs, kwargs)
            else:
                jobs[sname] = (map_filter_sort_pairs, kwargs)
        else:
            if mark_dups_by_coords or mark_dups_by_umis:
                raise NotImplementedError("TODO: SE data mark_duplicates")
            else:
                jobs[sname] = (map_filter_sort_single, kwargs)

    # run mapping jobs in parallel
    run_with_pool(jobs, log_level, workers, msg="Mapping samples to reference")

    # get bam file stats and write to a file
    jobs = {}
    for sname in fastq_dict:
        jobs[sname] = (count_mapped_reads, dict(sname=sname, outdir=outdir))
    stats = run_with_pool(jobs, log_level, workers, msg="Gathering mapping stats")

    # get a new stats outfile path in outdir
    idx = 0
    while 1:
        outstats = outdir / f"ipyrad_map_stats_{idx}.txt"
        if outstats.exists():
            idx += 1
        else:
            break

    # write stats
    df = pd.DataFrame({i: stats[i] for i in sorted(stats)}).T
    df.to_string(
        outstats,
        formatters={
            "nreads_processed": lambda x: f"{int(x)}",
            "nreads_passed_filters": lambda x: f"{int(x)}",
            "proportion_passed_filters": lambda x: f"{x:.3f}",
            "filtered_by_not_primary": lambda x: f"{int(x)}",
            "filtered_by_duplicates": lambda x: f"{int(x)}",
            "filtered_by_min_mapq": lambda x: f"{int(x)}",
            "filtered_by_bad_pairing": lambda x: f"{int(x)}",
            "filtered_by_max_soft_clip": lambda x: f"{int(x)}",
            "filtered_by_max_edit_dist": lambda x: f"{int(x)}",
            "mapq_mean_after_filters": lambda x: f"{x:.3f}",
            "mapq_median_after_filters": lambda x: f"{x:.3f}",
            "mapq_stdev_after_filters": lambda x: f"{x:.3f}",
        },
    )
    logger.info(f"mapping stats written to {outstats}")


if __name__ == "__main__":
    pass

