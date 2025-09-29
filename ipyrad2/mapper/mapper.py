#!/usr/bin/env python

"""Map, filter, sort, reads to BAM and delimit RAD locus beds.


Example
-------
map-delim --fastqs-rad ... --fastqs-wgs ... --ref REF --max-reads ... --min-samp-cov 4 --min-read-depth ...
"""

from typing import Tuple
import os
import sys
from pathlib import Path
import subprocess as sp
from loguru import logger
from ..utils.exceptions import IPyradError
from ..utils.parse_names import get_name_to_fastq_dict
# from ..utils.cluster import Cluster
# from ..utils.progress import track_remote_jobs
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


def map_filter_sort_dedup(sname: str, fastqs: Tuple[Path, Path], reference: Path, outdir: Path, umi_tag_in_i5: bool, threads: int) -> Path:
    """Map reads to the reference to get a sorted bam.

    TODO: test
    TODO: test r1 only
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
    log_dir = outdir / "logs"
    log_dir.mkdir(exist_ok=True)
    out_bam = outdir / f"{sname}.marked.sorted.bam"
    tmp_bam = outdir / f"{sname}.bam.tmp"
    tmp_prefix = outdir / f"{sname}.sam.tmp"

    # Split threads between BWA and samtools
    nthreads = max(1, int(threads))
    bwa_threads = max(1, nthreads - 1)
    mid_threads = max(1, int(nthreads / 2))

    # mapping command
    bwa_cmd = [
        BIN_BWA, "mem",
        "-t", str(bwa_threads),
        "-v", "1",                # less verbose.
        "-K", "50000000",         # stable nbases chunk size. Improves repeatability.
        # "-Y",                   # soft-clip supplementary. Wouldn't hurt, but not necessary.
        # "-M",                   # Picard compatibility. Not necessary, we use samtools fixmate for dups.
        # "-R", f"@RG\\tID:{sname}\\tSM:{sname}\\tPL:ILLUMINA",  # not currently used, since we provide custom -G to bcftools.
        str(reference),
        str(r1),
    ]
    if r2:
        bwa_cmd.append(str(r2))

    # drop unmapped + seconday + supplementary; require proper pair only if paired
    view_1_cmd = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",         # stream uncompressed bam
        "-F", "0x900",      # exclude secondary and supplemental.
        "-q", "20",         # only MAPQ≥20
        "-@", str(mid_threads),
    ]

    # [optional] sort by name for fixmate and marking dups
    sort_n_cmd = [
        BIN_SAMTOOLS, "sort",
        "-n",
        "-@", str(mid_threads),
        "-m", "50M",
        "-T", str(tmp_prefix),
        "-o", "-",
    ]

    # [optional] fixmate checks and updates tags about pairing
    fixmate_cmd = [
        BIN_SAMTOOLS, "fixmate",
        "-m",
        "-", "-",
        "-@", str(mid_threads),
    ]

    # coordinate sort command
    sort_c_cmd = [
        BIN_SAMTOOLS, "sort",
        "-m", "50M",                # tune per-thread memory
        "-T", str(tmp_prefix),
        "-O", "bam",
        "-o", "-",
        "-@", str(mid_threads),
    ]

    # mark dups in coordinate sorted fixmate bams
    markdup_cmd = [
        BIN_SAMTOOLS, "markdup",
        "-@", str(mid_threads),
        "-", "-",
    ]
    if umi_tag_in_i5:
        markdup_cmd.extend(["--barcode-rgx", "UMI_([ACGTN]+)"])

    # final view
    view_2_cmd = [
        BIN_SAMTOOLS, "view",
        "-b",
        "-f", "0x2",              # filter improperly paired
        "-q", "20",
        "-@", str(mid_threads),
        "-o", str(tmp_bam),
    ]

    cmds = [bwa_cmd, view_1_cmd, sort_n_cmd, fixmate_cmd, sort_c_cmd, view_2_cmd]
    run_pipeline(cmds, tmp_bam)
    os.replace(tmp_bam, out_bam)

    # CSI index bam file
    cmd = [BIN_SAMTOOLS, "index", "-c", "--threads", str(threads), str(out_bam)]
    res = sp.run(cmd, stdout=sp.DEVNULL, stderr=sp.PIPE)
    if res.returncode:
        raise IPyradError(f"samtools index failed ({res.returncode}).\n{res.stderr.decode(errors='ignore')}")
    # print(f"@@INFO: finished mapping/writing BAM data for {sname}")
    # logger.debug(f"finished mapping/writing BAM data for {sname}")
    return out_bam


def map_filter_sort(sname: str, fastqs: Tuple[Path, Path], reference: Path, outdir: Path, umi_tag_in_i5: bool, threads: int) -> Path:
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
    bwa_cmd = [
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
        bwa_cmd.append(str(r2))

    # drop unmapped + seconday + supplementary; require proper pair only if paired
    view_cmd = [
        BIN_SAMTOOLS, "view",
        "-b", "-u",         # stream uncompressed bam
        # "-F", "0x400",    # exclude optical/dups if marked (nb: bcftools already ignores reads that are marked.)
        "-F", "0x900",      # exclude secondary and supplemental.
        "-q", "20",         # only MAPQ≥20         # TODO: expose as param
        "-@", "1",
    ] + (["-f", "0x2"] if paired else [])

    # coordinate sorted command
    sort_cmd = [
        BIN_SAMTOOLS, "sort",
        "-m", "100M",                # tune per-thread memory
        "-T", str(tmp_prefix),
        "-O", "bam",
        "-o", str(tmp_bam),
        "-@", "1",
    ]

    print(f"@@DEBUG: cmd: {' '.join(map(str, bwa_cmd))}")
    print(f"@@DEBUG: cmd: {' '.join(map(str, view_cmd))}")
    print(f"@@DEBUG: cmd: {' '.join(map(str, sort_cmd))}")

    # Run pipeline and check *all* return codes with real stderr captured
    with sp.Popen(bwa_cmd, stdout=sp.PIPE, stderr=sp.PIPE) as p1:
        with sp.Popen(view_cmd, stdin=p1.stdout, stdout=sp.PIPE, stderr=sp.PIPE) as p2:
            p1.stdout.close()
            with sp.Popen(sort_cmd, stdin=p2.stdout, stdout=sp.DEVNULL, stderr=sp.PIPE) as p3:
                p2.stdout.close()
                _, err3 = p3.communicate()
            _, err2 = p2.communicate()
        _, err1 = p1.communicate()

    # Check in reverse order to surface the first failing stage
    if p3.returncode:
        raise IPyradError(f"samtools sort failed ({p3.returncode}).\n{err3.decode(errors='ignore')}")
    if p2.returncode:
        raise IPyradError(f"samtools view failed ({p2.returncode}).\n{err2.decode(errors='ignore')}")
    if p1.returncode:
        raise IPyradError(f"bwa mem failed ({p1.returncode}).\n{err1.decode(errors='ignore')}")

    # Atomic move
    os.replace(tmp_bam, out_bam)

    # CSI index bam file
    cmd = [BIN_SAMTOOLS, "index", "-c", "--threads", str(threads), str(out_bam)]
    res = sp.run(cmd, stdout=sp.DEVNULL, stderr=sp.PIPE)
    if res.returncode:
        raise IPyradError(f"samtools index failed ({res.returncode}).\n{res.stderr.decode(errors='ignore')}")
    # print(f"@@INFO: finished mapping/writing BAM data for {sname}")
    return out_bam


def count_mapped_reads(bam_file: Path) -> int:
    """Return the number of mapped reads in the filtered/sorted bam.

    Note that for PE data this is the still nreads, so divide by 2 to
    get the n read pairs.
    """
    # Count number of mapped read pairs
    cmd1 = [BIN_SAMTOOLS, "flagstat", bam_file]
    cmd2 = ["grep", "total"]

    # Run pipeline and check *all* return codes with real stderr captured
    with sp.Popen(cmd1, stdout=sp.PIPE, stderr=sp.PIPE) as p1:
        with sp.Popen(cmd2, stdin=p1.stdout, stdout=sp.PIPE, stderr=sp.PIPE) as p2:
            p1.stdout.close()
            line, err2 = p2.communicate()
        _, err1 = p1.communicate()
    # Check in reverse order to surface the first failing stage
    if p2.returncode:
        raise IPyradError(f"grep failed ({p2.returncode}).\n{err2.decode(errors='ignore')}")
    if p1.returncode:
        raise IPyradError(f"flagstat failed ({p1.returncode}).\n{err1.decode(errors='ignore')}")
    nreads_mapped = int(line.decode().strip().split()[0])
    return nreads_mapped


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
    with sp.Popen(cmd, stderr=sp.PIPE, stdout=sp.DEVNULL) as proc:
        error = proc.communicate()[1].decode()

    # error handling for one type of error on stderr
    if proc.returncode:
        if "please use bgzip" in error:
            raise IPyradError(
                "Reference sequence must be de-compressed fasta or bgzip "
                "compressed, your file is probably gzip compressed. The "
                "simplest fix is to gunzip your reference sequence by "
                "running this command: \n"
                f"    gunzip {reference}\n")
        raise IPyradError(error)


def run_mapper(
    fastqs: Tuple[Path, Path],
    outdir: Path,
    reference: Path,
    cores: int,
    threads: int,
    force: bool,
    mark_duplicates: bool,
    umi_tag_in_i5: bool,
    name_parse: Tuple[str, str] | None,
):
    # ------------------------------------------------------------
    # check reference and outdir paths
    reference = reference.expanduser().absolute()
    outdir = outdir.expanduser().absolute()
    outdir.mkdir(exist_ok=True)

    # parse dict of {name: (r1, r2)}
    fastq_dict = get_name_to_fastq_dict(fastqs, name_parse)

    # check outdir for existing and raise or remove
    result_files = [outdir / f"{sname}.sorted.bam" for sname in fastq_dict]
    if any(i.exists() for i in result_files):
        if not force:
            raise IPyradError(f"Bam files exist in outdir: e.g., {result_files[0]}. Use --force to overwrite.")

    # check mark_dups suitability
    if mark_duplicates:
        if list(fastq_dict.values())[0][1] is None:
            raise IPyradError("Data do not appear to be paired. Cannot use mark-duplicates with SE data.")
        if not umi_tag_in_i5:
            logger.warning("marking PCR duplicates by coordinates. Data is expected to be WGS, not RAD")
        if umi_tag_in_i5:
            logger.warning("marking PCR duplicates. Data is expected to be RAD with i5 UMI tags")

    # index the reference
    index_ref_with_bwa(reference)

    # run map, filter, sort
    logger.info(f"mapping and filtering {len(fastq_dict)} inputs to bams in {outdir}")
    logger.info(f"running up to {cores} parallel jobs each using up to {threads} threads")
    jobs = {}
    for sname, fastq_tuple in fastq_dict.items():
        kwargs = dict(
            fastqs=fastq_tuple,
            sname=sname,
            outdir=outdir,
            reference=reference,
            umi_tag_in_i5=umi_tag_in_i5,
            threads=max(1, threads),
        )
        jobs[sname] = kwargs
    bam_dict = run_with_pool(map_filter_sort_dedup, jobs, cores)

    # get bam file stats and write to a file
    jobs = {}
    for sname, bam_file in bam_dict.items():
        jobs[sname] = dict(bam_file=bam_file)
    stats = run_with_pool(count_mapped_reads, jobs, cores)

    # write stats
    handle = outdir / "ipyrad_map_stats.txt"
    with open(handle, 'w') as out:
        out.write("sample\tnreads_mapped\n")
        for key in sorted(stats):
            out.write(f"{key}\t{stats[key]}\n")
        logger.info(f"mapping stats written to {handle}")
    # # ------------------------------------------------------------
    # with Cluster(cores) as ipyclient:
    #     lbview = ipyclient.load_balanced_view()
    #     thview = ipyclient.load_balanced_view(ipyclient.ids[::threads])
    #     jobs = {}

    #     if mark_duplicates:
    #         logger.info(f"mapping/filtering/sorting {len(fastq_dict)} inputs to {outdir}")
    #     else:
    #         logger.info("...")

    #     # run map, filter, sort
    #     for sname, fastq_tuple in fastq_dict.items():
    #         kwargs = dict(
    #             fastqs=fastq_tuple,
    #             sname=sname,
    #             outdir=outdir,
    #             reference=reference,
    #             umi_tag_in_i5=umi_tag_in_i5,
    #             threads=max(1, threads),
    #         )
    #         if mark_duplicates:
    #             jobs[sname] = thview.apply(map_filter_sort_dedup, **kwargs)
    #         else:
    #             jobs[sname] = thview.apply(map_filter_sort, **kwargs)
    #     bam_dict = track_remote_jobs(jobs, ipyclient)

    #     # get bam file stats and write to a file
    #     jobs = {}
    #     for sname in bam_dict:
    #         jobs[sname] = lbview.apply(count_mapped_reads, bam_dict[sname])
    #     stats = track_remote_jobs(jobs, ipyclient)
    #     handle = outdir / "ipyrad_map_stats.txt"
    #     with open(handle, 'w') as out:
    #         out.write("sample\tnreads_mapped\n")
    #         for key in sorted(stats):
    #             out.write(f"{key}\t{stats[key]}\n")
    #     logger.info(f"mapping stats written to {handle}")


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
