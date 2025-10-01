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
import subprocess as sp
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


def map_filter_sort_mark(sname: str, fastqs: Tuple[Path, Path], reference: Path, outdir: Path, mark_dups_by_umis: bool, threads: int) -> Path:
    """Map reads to the reference to get a sorted bam.
    """
    if not reference.exists():
        raise IPyradError(f"reference_sequence not found: {reference}")

    # Inputs: paired if split files exist, otherwise single-end derep
    r1, r2 = fastqs
    paired = r1.exists() and r2.exists()
    if not paired:
        r2 = None
    if not r1.exists():
        raise IPyradError(f"fastq file not found: {r1}")

    # paths
    out_bam = outdir / f"{sname}.marked.sorted.bam"
    tmp_bam = outdir / f"{sname}.bam.tmp"
    tmp_prefix = outdir / f"{sname}.sam.tmp"

    # Split threads between BWA and samtools
    nthreads = max(1, int(threads))
    bwa_threads = max(1, nthreads - 1)
    mid_threads = max(1, int(nthreads / 2))

    # mapping command
    cmd1 = [
        BIN_BWA, "mem",
        "-t", str(bwa_threads),
        "-v", "1",                # less verbose.
        "-K", "50000000",         # stable nbases chunk size. Improves repeatability.
        # "-Y",                   # soft-clip supplementary. Wouldn't hurt, but not necessary.
        # "-M",                   # Picard compatibility. Not necessary, we use samtools fixmate for dups.
        "-R", f"@RG\\tID:{sname}\\tSM:{sname}\\tPL:ILLUMINA",  # store sample name; group can be overriden with -G.
        str(reference),
        str(r1),
    ]
    if r2:
        cmd1.append(str(r2))

    # drop unmapped + seconday + supplementary; require proper pair only if paired
    cmd2 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",         # stream uncompressed bam
        "-F", "0x900",      # exclude secondary and supplemental.
        "-q", "20",         # only MAPQ≥20
        "-@", str(mid_threads),
    ]

    # [optional] sort by name for fixmate and marking dups
    cmd3 = [
        BIN_SAMTOOLS, "sort",
        "-n",
        "-@", str(mid_threads),
        "-m", "50M",
        "-T", str(tmp_prefix),
        "-o", "-",
    ]

    # [optional] fixmate checks and updates tags about pairing
    cmd4 = [
        BIN_SAMTOOLS, "fixmate",
        "-m",
        "-", "-",
        "-@", str(mid_threads),
    ]

    # coordinate sort command
    cmd5 = [
        BIN_SAMTOOLS, "sort",
        "-m", "50M",                # tune per-thread memory
        "-T", str(tmp_prefix),
        "-O", "bam",
        "-o", "-",
        "-@", str(mid_threads),
    ]

    # mark dups in coordinate sorted fixmate bams
    cmd6 = [
        BIN_SAMTOOLS, "markdup",
        "-@", str(mid_threads),
        "-", "-",
    ]
    if mark_dups_by_umis:
        cmd6.extend(["--barcode-rgx", "UMI_([ACGTN]+)"])

    # final view
    cmd7 = [
        BIN_SAMTOOLS, "view",
        "-b",
        "-f", "0x2",              # filter improperly paired
        "-q", "20",
        "-@", str(mid_threads),
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


def map_filter_sort(sname: str, fastqs: Tuple[Path, Path], reference: Path, outdir: Path, threads: int, **kwargs) -> Path:
    """Map reads to the reference to get a sorted bam."""
    if not reference.exists():
        raise IPyradError(f"reference_sequence not found: {reference}")

    # Inputs: paired if split files exist, otherwise single-end derep
    r1, r2 = fastqs
    paired = r1.exists() and r2.exists()
    if not paired:
        r2 = None
    if not r1.exists():
        raise IPyradError(f"fastq file not found: {r1}")

    # paths
    outdir = Path(outdir).expanduser().absolute()
    outdir.mkdir(exist_ok=True)
    out_bam = outdir / f"{sname}.sorted.bam"
    tmp_bam = outdir / f"{sname}.bam.tmp"
    tmp_prefix = outdir / f"{sname}.sam.tmp"

    # Split threads between BWA and samtools
    nthreads = max(1, int(threads))
    bwa_threads = max(1, nthreads - 1)

    # mapping command
    cmd1 = [
        BIN_BWA, "mem",
        "-t", str(bwa_threads),
        "-v", "1",                # less verbose.
        "-K", "50000000",         # stable chunk size. Improves repeatability.
        # "-Y",                   # soft-clip supplementary. Wouldn't hurt, but not necessary.
        # "-M",                   # Picard compatibility. Not necessary, we use samtools fixmate for dups.
        "-R", f"@RG\\tID:{sname}\\tSM:{sname}\\tPL:ILLUMINA",  # not currently used, since we provide custom -G to bcftools.
        str(reference),
        str(r1),
    ]
    if r2:
        cmd1.append(str(r2))

    # drop unmapped + seconday + supplementary; require proper pair only if paired
    cmd2 = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",         # stream uncompressed bam
        # "-F", "0x400",    # exclude optical/dups if marked (nb: bcftools already ignores reads that are marked.)
        "-F", "0x900",      # exclude secondary and supplemental.
        "-q", "20",         # only MAPQ≥20         # TODO: expose as param
        "-@", "1",
    ] + (["-f", "0x2"] if paired else [])

    # coordinate sorted command
    cmd3 = [
        BIN_SAMTOOLS, "sort",
        "-m", "100M",                # tune per-thread memory
        "-T", str(tmp_prefix),
        "-O", "bam",
        "-o", str(tmp_bam),
        "-@", "1",
    ]
    run_pipeline([cmd1, cmd2, cmd3])
    os.replace(tmp_bam, out_bam)

    # CSI index bam file
    cmd1 = [BIN_SAMTOOLS, "index", "-c", "--threads", str(threads), str(out_bam)]
    run_pipeline([cmd1])
    logger.debug(f"finished mapping: {sname}")
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
            raise IPyradError("{prefix}.bam exists for >=1 sample in outdir. Use --force to overwrite.")
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
        if mark_dups_by_coords or mark_dups_by_umis:
            jobs[sname] = (map_filter_sort_mark, kwargs)
        else:
            jobs[sname] = (map_filter_sort, kwargs)
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
