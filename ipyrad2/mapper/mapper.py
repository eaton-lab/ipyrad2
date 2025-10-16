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


def map_filter_sort_mark_pairs(sname: str, fastqs: Tuple[Path, Path], reference: Path, outdir: Path, mark_dups_by_umis: bool, min_map_q: int, threads: int) -> Path:
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
    out_bam = outdir / f"{sname}.filtered.bam"
    bam_namesort = outdir / "tmpdir" / f"{sname}.tmp.namesort.bam"
    bam_fixmate = outdir / "tmpdir" / f"{sname}.tmp.fixmate.bam"
    bam_coordsort = outdir / "tmpdir" / f"{sname}.tmp.coordsort.bam"
    bam_markdup = outdir / "tmpdir" / f"{sname}.tmp.markdup.bam"
    tmp_prefix = outdir / "tmpdir" / f"{sname}.tmp.pre"
    tmp_stats1 = outdir / "tmpdir" / f"{sname}.tmp.stats1.json"
    tmp_stats2 = outdir / "tmpdir" / f"{sname}.tmp.stats2.json"
    tmp_stats3 = outdir / "tmpdir" / f"{sname}.tmp.stats3.json"
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
        # "-Y",                   # mark supplementary with soft-clip. We exclude supplementals anyways.
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
    cmd8 = [
        BIN_SAMTOOLS, "view",
        "-e", '((flag&4)==0) && ((flag&8)==0) && (rnext=="=" || rnext==rname) && (tlen>=-2000 && tlen<=2000)',
        "--save-counts", str(tmp_stats3),
        "--write-index",
        "-o", str(out_bam),
    ]
    run_pipeline([cmd7, cmd8])
    logger.debug(f"finished mapping: {sname}")
    bam_markdup.unlink()
    bam_markdup.with_suffix(bam_markdup.suffix + '.csi').unlink()
    return out_bam


def map_filter_sort_pairs(sname: str, fastqs: Tuple[Path, Path], reference: Path, min_map_q: int, outdir: Path, threads: int, **kwargs) -> Path:
    """Map reads to the reference to get a sorted bam.

    This pipeline is for PE data w/o duplication marking information.
    """
    # paths
    out_bam = outdir / f"{sname}.filtered.bam"
    # tmp_bam = outdir / f"{sname}.bam.tmp"
    tmp_prefix = outdir / "tmpdir" / f"{sname}.tmp.pre"
    tmp_stats1 = outdir / "tmpdir" / f"{sname}.tmp.stats1.json"
    tmp_stats2 = outdir / "tmpdir" / f"{sname}.tmp.stats2.json"
    tmp_stats3 = outdir / "tmpdir" / f"{sname}.tmp.stats3.json"

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
        "-e", '((flag&4)==0) && ((flag&8)==0) && (rnext=="=" || rnext==rname) && (tlen>=-2000 && tlen<=2000)',
        "--save-counts", str(tmp_stats3),
        "-o", "-",
    ]

    # coordinate sort
    cmd5 = [
        BIN_SAMTOOLS, "sort",
        "-m", "256M",                # per-thread memory
        "-T", str(tmp_prefix),
        "-@", str(sort_threads),
        "--write-index",
        "-O", "bam",                 # explicity ask for bam format to be safe
        "-o", str(out_bam),
        "-",
    ]
    cmds = [cmd1, cmd2, cmd3, cmd4, cmd5]
    run_pipeline(cmds)
    logger.info(f"finished mapping: {sname}")
    return out_bam


def concat_tech_reps_into_tmpdir(imap: Path, tmpdir: Path, fastq_dict: Dict[str, Tuple[Path, Path]]) -> Dict[str, Path]:
    """Return fastq_dict pointing to updated concat paths in tmpdir"""
    # parse the population file
    df = pd.read_csv(imap, header=None, sep=r"\s+")

    # fill pdict and warn if names don't match any samples.
    snames = set(fastq_dict)
    pop2tups = defaultdict(list)
    pop2snames = defaultdict(list)
    for idx in df.index:
        sname, pname, *_ = df.loc[idx]
        if sname in snames:
            pop2tups[pname].append(fastq_dict.pop(sname))
            pop2snames[pname].append(sname)
        else:
            logger.warning(f"sample name '{sname}' from imap file was not found in data. Skipping.")

    # return original dict if nothing to merge/rename
    if not pop2snames:
        return fastq_dict

    # report to logger
    logger.info(f"merging/renaming samples according to imap file: {imap}")
    maxlen = max(len(i) for i in pop2snames)
    for pname, tups in pop2tups.items():
        snames = pop2snames[pname]
        logger.info(f"{pname}{' ' * (maxlen - len(pname))} <- {' + '.join(snames)}")

        # renaming, do not run pipe
        if len(snames) == 1:
            fastq_dict[pname] = tups

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
        "nreads_filtered_by_not_primary": d1["records_filter_rejected"],
        "nreads_filtered_by_duplicates": dd["records_filter_rejected"],
        "nreads_filtered_by_min_mapq": d2["records_filter_rejected"],
        "nreads_filtered_by_bad_pairing": d3["records_filter_rejected"],
        "nreads_passed_filters": d3["records_filter_accepted"],
        "proportion_reads_passed_filters": float(d3["records_filter_accepted"] / d1["records_processed"]),
        "mapq_mean_after_filters": float(mean_mapq),
        "mapq_median_after_filters": float(median_mapq),
        "mapq_stdev_after_filters": float(stdev_mapq),
    }
    s1.unlink()
    s2.unlink()
    s3.unlink()
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
    # run at most this many concurrent jobs
    workers = max(1, cores // threads)

    # check reference and outdir paths
    reference = reference.expanduser().absolute()
    outdir = outdir.expanduser().absolute()
    tmpdir = outdir / "tmpdir"
    tmpdir.mkdir(exist_ok=True, parents=True)

    # parse dict of {name: (r1, r2)}
    fastq_dict = get_name_to_fastq_dict(fastqs, delim_str, delim_idx)

    # check outdir for existing and raise or remove
    result_files = [list(outdir.glob(f"{sname}.*.bam")) for sname in fastq_dict]
    if any(result_files):
        if not force:
            raise IPyradError(".bam exists for >=1 of the selected samples in outdir. Use --force to overwrite.")
        else:
            for bam_list in result_files:
                for bam_file in bam_list:
                    if bam_file.exists():
                        logger.debug(f"removing existing bam file: {bam_file}")

    # check mark_dups suitability
    if mark_dups_by_coords or mark_dups_by_umis:
        if mark_dups_by_coords and mark_dups_by_umis:
            raise IPyradError("you cannot select both mark_dups_by_coords and mark_dups_by_umis.")
        if list(fastq_dict.values())[0][1] is None:
            raise IPyradError("Data do not appear to be paired. Cannot mark duplicates for SE data.")
        # TODO: check for valid rather than just warn.
        if mark_dups_by_coords:
            logger.warning("marking PCR duplicates by coordinates. Data is expected to be WGS, not RAD")
        if mark_dups_by_umis:
            logger.warning("marking PCR duplicates. Data is expected to be RAD with i5 UMIs moved into read names")

    # store whether reads are paired.
    # TODO: this can raise an error when glob catches extras (e.g., not .fq, .gz, etc). Report that their glob might be bad?
    is_paired = False
    pairs_exist = [(r1.exists() and r2.exists()) for (r1, r2) in fastq_dict.values()]
    if any(pairs_exist):
        if all(pairs_exist):
            is_paired = True
        else:
            raise IPyradError("some but not all files have R1 and R2 pairs. Check inputs.")

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
            fastqs=fastq_tuple,
            sname=sname,
            outdir=outdir,
            reference=reference,
            min_map_q=min_map_q,
            mark_dups_by_umis=mark_dups_by_umis,
            threads=threads,
        )
        if is_paired:
            if mark_dups_by_coords or mark_dups_by_umis:
                jobs[sname] = (map_filter_sort_mark_pairs, kwargs)
            else:
                jobs[sname] = (map_filter_sort_pairs, kwargs)
        else:
            raise NotImplementedError("TODO: SE data.")
    # run mapping jobs in parallel
    run_with_pool(jobs, log_level, workers)

    # get bam file stats and write to a file
    jobs = {}
    for sname in fastq_dict:
        jobs[sname] = (count_mapped_reads, dict(sname=sname, outdir=outdir))
    stats = run_with_pool(jobs, log_level, workers)

    # get a new stats outfile path in outdir
    idx = 0
    while 1:
        outstats = outdir / f"ipyrad_map_stats_{idx}.txt"
        if outstats.exists():
            idx += 1
        else:
            break

    # write stats
    df = pd.DataFrame(stats).T
    df.to_string(
        outstats,
        formatters={
            "nreads_processed": lambda x: f"{int(x)}",
            "nreads_filtered_by_not_primary": lambda x: f"{int(x)}",
            "nreads_filtered_by_duplicates": lambda x: f"{int(x)}",
            "nreads_filtered_by_min_mapq": lambda x: f"{int(x)}",
            "nreads_filtered_by_bad_pairing": lambda x: f"{int(x)}",
            "nreads_passed_filters": lambda x: f"{int(x)}",
            "proportion_reads_passed_filters": lambda x: f"{x:.3f}",
            "mapq_mean_after_filters": lambda x: f"{x:.3f}",
            "mapq_median_after_filters": lambda x: f"{x:.3f}",
            "mapq_stdev_after_filters": lambda x: f"{x:.3f}",
        },
    )
    logger.info(f"mapping stats written to {outstats}")


if __name__ == "__main__":
    pass

