#!/usr/bin/env python

"""Methods for quality and adapter trimming rads for RAD-seq analysis.

1. Check file paths and search for sample pairs.
2. Check for restriction overhangs using kmers.
3. Trim reads with fastp.
4. Write stats summary file.

Example
--------
python trim_fastqs.py FQs ...
"""

from typing import Tuple, List, Dict, Any
import sys
import json
from pathlib import Path
from loguru import logger
import pandas as pd
from ..utils.names import get_name_to_fastq_dict
from ..utils.kmers import get_overhang_from_kmers
from ..utils.exceptions import IPyradError
from ..utils.parallel import run_with_pool, run_pipeline


FASTP_BINARY = Path(sys.prefix) / "bin" / "fastp"
ADAPTERS = Path(__file__).absolute().parent / "adapters.fa"


def check_user_re_overhangs(kmer_overhangs: Tuple[str, str], user_overhangs: Tuple[str, str]) -> None:
    """Warn user if found overhang differs from user overhangs.
    """
    if user_overhangs[0] != kmer_overhangs[0]:
        logger.warning(
            "kmer analysis identified the read1 restriction overhang "
            f"as '{kmer_overhangs[0]}', however, you entered {user_overhangs[0]}. "
            "Your entered value is being used, but we recommend comparing "
            "your results with a run using the auto-detected overhangs."
        )
    if user_overhangs[1] != kmer_overhangs[1]:
        logger.warning(
            "kmer analysis identified the read2 restriction overhang "
            f"as '{kmer_overhangs[1]}', however, you entered {user_overhangs[1]}. "
            "Your entered value is being used, but we recommend comparing "
            "your results with a run using the auto-detected overhangs."
        )


def trim_sample_with_fastp(
    fastqs: Tuple[Path, Path],
    sname: str,
    outdir: Path,
    restriction_overhangs: Tuple[str, str],
    max_reads: int,
    min_trimmed_length: int,
    min_quality: int,
    max_low_quality_bases: int,
    phred_qscore_offset: int,
    disable_adapter_trimming: bool,
    disable_quality_filtering: bool,
    umi_tag_in_i5: bool,
    threads: int,
) -> Dict[str, Tuple[Tuple[Path, Path], Dict[str, Any]]]:
    """Run FASTP and return: {sname: ((r1, r2), stats)}
    """
    out1 = outdir / f"{sname}.R1.trimmed.fastq.gz"
    out2 = outdir / f"{sname}.R2.trimmed.fastq.gz"
    stats_html = outdir / f"{sname}.stats.html"
    stats_json = outdir / f"{sname}.stats.json"

    # check if data are paired
    is_paired = fastqs[1].exists() and fastqs[1].name != "-null-"

    # build command
    if is_paired:
        cmd = [
            str(FASTP_BINARY),
            "-i", str(fastqs[0]),
            "-I", str(fastqs[1]),
            "-o", str(out1),
            "-O", str(out2),
            "--detect_adapter_for_pe",
            "-c",                             # paired-end overlap base correction
            "-w", str(threads),   # 2 workers + 2 i/o threads.
            "--trim_front1", str(len(restriction_overhangs[0])),
            "--trim_front2", str(len(restriction_overhangs[1])),
            "--overlap_len_require", "20",    # default is 30
            "--overlap_diff_limit", "5",
        ]
        if umi_tag_in_i5:
            cmd.extend(["-U", "--umi_loc=index2", "--umi_prefix=UMI"])
    else:
        cmd = [
            str(FASTP_BINARY),
            "-i", str(fastqs[0]),
            "-o", str(out1),
            "-w", str(threads),   # 3 workers + 1 i/o threads.
            "--trim_front1", str(len(restriction_overhangs[0])),
        ]

    # common arguments
    cmd.extend([
        "-q", str(min_quality + phred_qscore_offset - 33), # threshold to determine low qual bases
        "-u", "10",                                        # percentage of low qual bases allowed
        "-M", str(min_quality + phred_qscore_offset - 33), # mean quality score used by --cut args
        "--cut_front", "--cut_front_window_size", "4",
        "--cut_tail", "--cut_tail_window_size", "4",
        "-l", str(min_trimmed_length),
        "--trim_poly_g", "--trim_poly_x", "--poly_x_min_len", "10",
        "-y", "-Y", "50",                 # turns on and sets complexity filter to 50
        "--n_base_limit", str(max_low_quality_bases),
        "-j", str(stats_json),
        "-h", str(stats_html),
        "--adapter_fasta", str(ADAPTERS),
    ])

    # mostly for testing, normalize/subsample reads
    if max_reads is not None:
        cmd.extend(["--reads_to_process", str(max_reads)])
    if disable_adapter_trimming:
        cmd.extend("-A")
    if disable_quality_filtering:
        cmd.extend("-Q")

    # run the command in subprocess
    logger.debug(f"CMD: {' '.join(cmd)}")
    run_pipeline([cmd])
    logger.info(f"finished trimming {sname}")
    return None


def write_stats_summary(snames: List[str], outdir: Path):
    """Collect fastp stats from all samples in outdir and write summary.

    If user runs multiple ipyrad trim multiple times with the same
    out directory this will append an int to stats name each time to
    write a new stats file.
    """
    # get a new stats outfile path in outdir
    idx = 0
    while 1:
        outfile = outdir / f"ipyrad_trim_stats_{idx}.txt"
        if outfile.exists():
            idx += 1
        else:
            break

    # load all stats dicts from jsons
    jdata = {}
    snames = sorted(snames)
    for sname in snames:
        stats_file = outdir / f"{sname}.stats.json"
        if stats_file.exists():
            with open(stats_file, 'r', encoding="utf-8") as indata:
                jdata[sname] = json.loads(indata.read())

    # init the dataframe
    df = pd.DataFrame(index=snames, columns=[
        "total_reads_before", "total_bases_before", "q20_rate_before", "q30_rate_before",
        "read1_mean_length_before", "read2_mean_length_before",
        "total_reads_after", "total_bases_after", "q20_rate_after", "q30_rate_after",
        "read1_mean_length_after", "read2_mean_length_after",
        "reads_filtered_by_low_quality",
        "reads_filtered_by_too_many_N",
        "reads_filtered_by_low_complexity",
        "reads_filtered_by_too_short",
        "adapter_trimmed_reads",
        "adapter_trimmed_bases",
    ])

    # fill the dataframe
    for sname in snames:
        j = jdata[sname]
        df.loc[sname, "total_reads_before"] = j["summary"]["before_filtering"]["total_reads"]
        df.loc[sname, "total_bases_before"] = j["summary"]["before_filtering"]["total_bases"]
        df.loc[sname, "q20_rate_before"] = j["summary"]["before_filtering"]["q20_rate"]
        df.loc[sname, "q30_rate_before"] = j["summary"]["before_filtering"]["q30_rate"]
        df.loc[sname, "read1_mean_length_before"] = j["summary"]["before_filtering"]["read1_mean_length"]
        df.loc[sname, "read2_mean_length_before"] = j["summary"]["before_filtering"]["read2_mean_length"]
        df.loc[sname, "total_reads_after"] = j["summary"]["after_filtering"]["total_reads"]
        df.loc[sname, "total_bases_after"] = j["summary"]["after_filtering"]["total_bases"]
        df.loc[sname, "q20_rate_after"] = j["summary"]["after_filtering"]["q20_rate"]
        df.loc[sname, "q30_rate_after"] = j["summary"]["after_filtering"]["q30_rate"]
        df.loc[sname, "read1_mean_length_after"] = j["summary"]["after_filtering"]["read1_mean_length"]
        df.loc[sname, "read2_mean_length_after"] = j["summary"]["after_filtering"]["read2_mean_length"]
        df.loc[sname, "reads_filtered_by_low_quality"] = j["filtering_result"]["low_quality_reads"]
        df.loc[sname, "reads_filtered_by_too_many_N"] = j["filtering_result"]["too_many_N_reads"]
        df.loc[sname, "reads_filtered_by_low_complexity"] = j["filtering_result"]["low_complexity_reads"]
        df.loc[sname, "reads_filtered_by_too_short"] = j["filtering_result"]["too_short_reads"]
        df.loc[sname, "adapter_trimmed_reads"] = j["adapter_cutting"]["adapter_trimmed_reads"]
        df.loc[sname, "adapter_trimmed_bases"] = j["adapter_cutting"]["adapter_trimmed_bases"]

    # write human readable whitespace delimited.
    df.to_string(outfile, float_format=lambda x: f"{x:.6f}")


def run_trimmer(
    fastqs: List[Path],
    outdir: Path,
    restriction_overhangs: Tuple[str, str],
    max_reads: int,
    min_trimmed_length: int,
    min_quality: int,
    max_low_quality_bases: int,
    phred_qscore_offset: int,
    max_reads_kmer: int,
    disable_infer_re_overhangs: bool,
    disable_adapter_trimming: bool,
    disable_quality_filtering: bool,
    cores: int,
    threads: int,
    delim_str: str | None,
    delim_idx: int,
    suffix: str | None,
    umi_tag_in_i5: bool,
    force: bool,
    log_level: str,
):
    # ------------------------------------------------------------
    # parse dict of {name: (r1, r2)}
    # fastq_dict = get_fastq_tuples_dict_from_paths_list(fastqs)
    fastq_dict = get_name_to_fastq_dict(fastqs, delim_str, delim_idx, suffix)

    # check outdir for existing and raise or remove
    result_files = [outdir / f"{sname}.R1.trimmed.fastq.gz" for sname in fastq_dict]
    if any(i.exists() for i in result_files):
        if not force:
            raise IPyradError(f"Trimmed fastqs exist in outdir: e.g., {result_files[0]}. Use --force to overwrite.")
    outdir.mkdir(exist_ok=True)

    # run at most this many concurrent jobs
    workers = max(1, cores // threads)

    # ------------------------------------------------------------
    # infer restriction overhangs by kmer analysis
    if not disable_infer_re_overhangs:
        re1 = get_overhang_from_kmers([i[0] for i in fastq_dict.values()], 20, 100_000, cores, log_level)
        re2 = get_overhang_from_kmers([i[1] for i in fastq_dict.values()], 20, 100_000, cores, log_level)
        logger.info(f"restriction site overhangs inferred by kmer analysis = {re1} {re2}")
        # allow user override but warn if it doesn't match inferred.
        if restriction_overhangs:
            check_user_re_overhangs((re1, re2), restriction_overhangs)
            re1, re2 = restriction_overhangs
    elif restriction_overhangs:
        re1, re2 = restriction_overhangs
    else:
        re1, re2 = ("", "")
    logger.info(f"restriction site overhangs set to {re1} {re2}")

    # ------------------------------------------------------------
    jobs = {}
    logger.info(f"trimming/filtering {len(fastq_dict)} inputs with 'fastp' and writing to {outdir}")
    logger.info(f"running up to {workers} parallel jobs each using up to {threads} threads")
    for sname, fastq_tuple in fastq_dict.items():
        kwargs = dict(
            fastqs=fastq_tuple,
            sname=sname,
            outdir=outdir,
            restriction_overhangs=(re1, re2),
            max_reads=max_reads,
            min_trimmed_length=min_trimmed_length,
            min_quality=min_quality,
            max_low_quality_bases=max_low_quality_bases,
            phred_qscore_offset=phred_qscore_offset,
            disable_adapter_trimming=disable_adapter_trimming,
            disable_quality_filtering=disable_quality_filtering,
            umi_tag_in_i5=umi_tag_in_i5,
            threads=threads,  # recommended >=3 since 2 are used for i/o
        )
        jobs[sname] = (trim_sample_with_fastp, kwargs)
    results = run_with_pool(jobs, log_level, workers)
    write_stats_summary(sorted(results), outdir)


if __name__ == "__main__":

    PATHS = sorted(Path("/home/deren/Documents/tools/ipyrad2/examples/Pedic-PE-ddRAD/").glob("*.gz"))
    # fastq_dict = get_fastq_tuples_dict_from_paths_list(PATHS)
    # max_reads = int(500_000 / len(fastq_dict))
    # logger.warning(max_reads)
    # overhangs = find_re_overhangs(fastq_dict, max_reads)

    # fastqs = fastq_dict["torta-DE758-plate_J2"]

    # x = trim_with_fastp(
    #     fastqs, 'test', "/tmp/", overhangs,
    #     max_reads=10_000,
    #     min_trimmed_length=35,
    #     min_quality=20,
    #     max_low_quality_bases=5,
    #     phred_qscore_offset=33,
    #     disable_adapter_trimming=False,
    #     disable_quality_filtering=False,
    #     threads=2,
    # )

    # print(x)


if __name__ == "__main__":
    pass
