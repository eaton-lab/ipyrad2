#!/usr/bin/env python

"""Map, filter, sort, mark BAMs.

Example
-------
map --fastqs DATA --ref REF -out MAP
"""

from typing import Tuple
import os
import sys
from pathlib import Path
from loguru import logger
from ..utils.exceptions import IPyradError
from ..utils.names import get_name_to_fastq_dict
from ..utils.parallel import run_pipeline, run_with_pool


BIN = Path(sys.prefix) / "bin"
BIN_BWA = str(BIN / "bwa-mem2")
BIN_SAMTOOLS = str(BIN / "samtools")  # indexing

LOG_CMD1 = ("bwa-mem2 mem REF R1 R2 "
           "| samtools view -bu -F 0x900 -q 20 "
           "| samtools sort -o - ")
LOG_CMD2 = ("bwa-mem2 mem REF R1 R2 "
           "| samtools view -bu -F 0x900 -q 20 "
           "| samtools sort -n -o - "
           "| samtools fixmate -m - - "
           "| samtools sort -o - "
           "| samtools markdup - - "
           "| samtools view -b -f 0x2 -q 20")
LOG_CMD3 = ("bwa-mem2 mem REF R1 R2 "
           "| samtools view -bu -F 0x900 -q 20 "
           "| samtools sort -n -o - "
           "| samtools fixmate -m - - "
           "| samtools sort -o - "
           "| samtools markdup --barcode-rgx 'UMI_([ACGTN]+)' - - "
           "| samtools view -b -f 0x2 -q 20")


def map_filter_sort_mark_pairs(sname: str, fastqs: Tuple[Path, Path], reference: Path, outdir: Path, mark_dups_by_umis: bool, threads: int) -> Path:
    """Map reads to the reference to get a sorted bam.

    This pipeline is for PE data w duplication marking information,
    meaning either WGS data (coordinates), or ddRAD w/ i5 UMIs stored
    to names (see option in ipyrad2 trim to do this.)
    """
    # paths
    out_bam = outdir / f"{sname}.marked.sorted.bam"
    tmp_bam = outdir / f"{sname}.bam.tmp"
    tmp_prefix = outdir / f"{sname}.sam.tmp"

    # Split threads between BWA and samtools
    bwa_threads = max(1, int(threads * 0.75))
    sort_threads = max(1, threads - bwa_threads)

    # mapping command
    # additional options to consider for toggle
    # -k 15   # reduce kmer size for small denovo loci
    # -L 1,1  # reduce penalty for clipping
    cmd1 = [
        BIN_BWA, "mem",
        "-Y",                     # mark supplementary with soft-clip. Doesn't hurt reference mapping, but helps short loci mapping.
        "-T", "20",               # minimum alignmnet score to output (default=30). This is different from MAPQ in samtools.
        "-R", f"@RG\\tID:{sname}\\tSM:{sname}\\tPL:ILLUMINA",  # store sample name; group can be overriden with -G.
        "-K", "50000000",         # stable nbases chunk size. Improves repeatability.
        "-v", "1",                # less verbose.
        "-t", str(bwa_threads),
        str(reference),
        str(fastqs[0]),
        str(fastqs[1]),
    ]

    # drop unmapped + seconday + supplementary; require proper pair only if paired
    cmd2 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",
        "-o", "-",
    ]

    # [optional] sort by name for fixmate and marking dups
    cmd3 = [
        BIN_SAMTOOLS, "sort",
        "-n",
        "-m", "256M",
        "-T", str(tmp_prefix),
        "-o", "-",
        "-@", str(sort_threads),
    ]

    # [optional] fixmate checks and updates tags about pairing
    cmd4 = [
        BIN_SAMTOOLS, "fixmate",
        "-m",
        "-r",                       # remove unmapped and secondary
        "-", "-",
        "-@", str(bwa_threads),
    ]

    # coordinate sort command
    cmd5 = [
        BIN_SAMTOOLS, "sort",
        "-m", "50M",                # tune per-thread memory
        "-T", str(tmp_prefix),
        "-O", "bam",
        "-o", "-",
        "-@", str(sort_threads),
    ]

    # mark dups in coordinate sorted fixmate bams
    cmd6 = [
        BIN_SAMTOOLS, "markdup",
        "-@", str(bwa_threads),
        "-", "-",
    ]
    if mark_dups_by_umis:
        cmd6.extend(["--barcode-rgx", "UMI_([ACGTN]+)"])

    # final view
    cmd7 = [
        BIN_SAMTOOLS, "view",
        "-F", "4",                   # drop 'unmapped''
        "-F", "8",                   # drop 'mate-unmapped'
        "-F", "2048",                # drop 'secondary''
        "-q", "0",                   # do not yet apply map quality filters
        "-b",
        "-o", str(tmp_bam),
    ]

    cmds = [cmd1, cmd2, cmd3, cmd4, cmd5, cmd6, cmd7]
    run_pipeline(cmds, tmp_bam)
    os.replace(tmp_bam, out_bam)

    # CSI index bam file
    cmd = [BIN_SAMTOOLS, "index", "-c", "--threads", str(threads), str(out_bam)]
    run_pipeline([cmd])
    logger.debug(f"finished mapping: {sname}")
    return out_bam


    # [optional] sort by name for fixmate and marking dups
    cmd3 = [
        BIN_SAMTOOLS, "sort",
        "-n",
        "-m", "256M",
        "-T", str(tmp_prefix),
        "-@", str(sort_threads),
        "-o", "-",
    ]

    # [optional] fixmate checks and updates tags about pairing
    cmd4 = [
        BIN_SAMTOOLS, "fixmate",
        "-m",
        "-r",                       # remove unmapped and secondary
        "-", "-",
    ]


def map_filter_sort_pairs(sname: str, fastqs: Tuple[Path, Path], reference: Path, outdir: Path, threads: int, **kwargs) -> Path:
    """Map reads to the reference to get a sorted bam.

    This pipeline is for PE data w/o duplication marking information.
    """
    # paths
    out_bam = outdir / f"{sname}.sorted.filtered.bam"
    tmp_bam = outdir / f"{sname}.bam.tmp"
    tmp_prefix = outdir / f"{sname}.sam.tmp"

    # Split threads between BWA and samtools
    bwa_threads = max(1, int(threads * 0.75))
    sort_threads = max(1, threads - bwa_threads)

    # mapping command
    cmd1 = [
        BIN_BWA, "mem",
        "-Y",                     # soft-clip supplementary. Wouldn't hurt, but not necessary.
        "-T", "20",               # minimum score to output (default=30). Keep lower scores for now, can filter on higher MAPQ later.
        "-R", f"@RG\\tID:{sname}\\tSM:{sname}\\tPL:ILLUMINA",  # store sample names in bam.
        "-K", "50000000",         # stable chunk size. Improves repeatability.
        "-v", "1",                # less verbose.
        "-t", str(bwa_threads),
        str(reference),
        str(fastqs[0]),
        str(fastqs[1]),
    ]

    # stream converted to bam
    cmd2 = [
        BIN_SAMTOOLS, "view",
        "-u", "-b",                  # stream uncompressed bam
    ]

    # coordinate sort
    cmd3 = [
        BIN_SAMTOOLS, "sort",
        "-m", "256M",                # per-thread memory
        "-T", str(tmp_prefix),
        "-@", str(sort_threads),
        "-O", "bam",                 # explicity ask for bam format to be safe
        "-o", str(tmp_bam),
        "-",
    ]
    cmds = [cmd1, cmd2, cmd3]
    run_pipeline(cmds)

    # apply filters and write to disk
    cmd4 = [
        BIN_SAMTOOLS, "view",
        # "-F", "2060",                # combines three options below
        "-F", "4",                   # drop 'unmapped''
        "-F", "8",                   # drop 'mate-unmapped'
        "-F", "2048",                # drop 'secondary''
        "-q", "0",                   # do not yet apply map quality filters
        "-b",
        "-o", str(out_bam),
        str(tmp_bam),
    ]
    # cmds = [cmd1, cmd2, cmd3, cmd4]
    # cmd_strings = [f"{' '.join(cmd)}" for cmd in cmds]
    # logger.debug(" | ".join(cmd_strings))
    run_pipeline([cmd4])

    # CSI index bam file
    cmd1 = [BIN_SAMTOOLS, "index", "-c", "--threads", str(threads), str(out_bam)]
    run_pipeline([cmd1])
    logger.info(f"finished mapping: {sname}")
    return out_bam


def count_mapped_reads(bam_file: Path, threads: int) -> int:
    """Return the number of mapped reads in the filtered/sorted bam.

    Note that for PE data this is the still nreads, so divide by 2 to
    get the n read pairs.
    """
    # Count number of mapped read pairs
    cmd1 = [BIN_SAMTOOLS, "flagstat", "--threads", str(threads), bam_file]
    _, out, _ = run_pipeline([cmd1])
    lines = out.decode().strip().split("\n")
    for line in lines:
        parts = line.split()
        if parts[-1] == "primary":
            primary_mapped = int(parts[0])
        elif parts[-1] == "duplicates":
            primary_duplicates = int(parts[0])
    return {"mapped_primary": primary_mapped, "mapped_duplicates": primary_duplicates}


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
    outdir.mkdir(exist_ok=True)

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
        if mark_dups_by_coords:
            logger.warning("marking PCR duplicates by coordinates. Data is expected to be WGS, not RAD")
        if mark_dups_by_umis:
            logger.warning("marking PCR duplicates. Data is expected to be RAD with i5 UMIs moved into read names")

    # store whether reads are paired
    is_paired = False
    pairs_exist = [(r1.exists() and r2.exists()) for (r1, r2) in fastq_dict.values()]
    if any(pairs_exist):
        if all(pairs_exist):
            is_paired = True
        else:
            raise IPyradError("some but not all files have R1 and R2 pairs. Check inputs.")

    # index the reference
    index_ref_with_bwa(reference)

    # run map, filter, sort
    logger.info(f"mapping and filtering {len(fastq_dict)} inputs to bams in {outdir}")
    logger.info(f"running up to {workers} parallel jobs each using up to {threads} threads")
    jobs = {}
    for sname, fastq_tuple in fastq_dict.items():
        kwargs = dict(
            fastqs=fastq_tuple,
            sname=sname,
            outdir=outdir,
            reference=reference,
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
    # run jobs in parallel
    bam_dict = run_with_pool(jobs, log_level, workers)

    # get bam file stats and write to a file
    jobs = {}
    for sname, bam_file in bam_dict.items():
        jobs[sname] = (count_mapped_reads, dict(bam_file=bam_file, threads=threads))
    stats = run_with_pool(jobs, log_level, workers)

    # write stats
    handle = outdir / "ipyrad_map_stats.txt"
    with open(handle, 'w') as out:
        out.write("sample\tmapped_primary\tmapped_duplicates\n")
        for key in sorted(stats):
            out.write(f"{key}\t{stats[key]['mapped_primary']}\t{stats[key]['mapped_duplicates']}\n")
        logger.info(f"mapping stats written to {handle}")


if __name__ == "__main__":
    pass
    # PATHS = sorted(Path("/tmp/").glob("test.trimmed.*.gz"))
    # REF = Path("/home/deren/Documents/tools/ipyrad2/examples/LiuLiu-genome/Pcr.genome.1.0.fasta")
    # assert REF.exists()
    # fastq_dict = get_fastq_tuples_dict_from_paths_list(PATHS)
    # fastqs = fastq_dict["test.trimmed.R"]
    # map_filter_sort(
    #     "test",
    #     fastqs,
    #     REF,
    #     "/tmp",
    #     4,
    # )
