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

from typing import Dict, Tuple, List
import sys
import json
from pathlib import Path
import subprocess as sp
from loguru import logger
from ..utils.parse_names import get_name_to_fastq_dict
from ..utils.kmers import infer_overhang
from ..utils.exceptions import IPyradError
from ..utils.parallel import run_with_pool


FASTP_BINARY = Path(sys.prefix) / "bin" / "fastp"
ADAPTERS = Path(__file__).absolute().parent / "adapters.fa"


def infer_re_overhangs(fastq_dict: Dict[str, Tuple[Path, Path]], max_reads: int) -> None:
    """Infer re overhang from kmer analysis and compare w/ params setting.
    """
    logger.info("inferring restriction site overhangs from kmer analysis")

    # parse first restriction overhang
    read_r1s = [i[0] for i in fastq_dict.values()]
    re1 = infer_overhang(read_r1s, max_reads=max_reads)
    logger.debug(f"inferred R1 restriction overhang as: {re1}")

    # parse second restriction overhang
    read_r2s = [i[1] for i in fastq_dict.values()]
    if all(i.exists() for i in read_r2s):
        re2 = infer_overhang(read_r2s, max_reads=max_reads)
        logger.debug(f"inferred R2 restriction overhang as: {re2}")
    elif not any(i.exists() for i in read_r2s):
        re2 = ""
    else:
        raise ValueError("restriction sites could not be inferred for read2")
    return re1, re2


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
):
    """Run FASTP and return: {stats} (r1, r2)
    """
    outdir = Path(outdir).expanduser().absolute()
    outdir.mkdir(exist_ok=True)
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
        "-r",                                                    # move sliding window front to tail
        "-M", str(min_quality + phred_qscore_offset - 33),       # mean quality in -r window
        "-q", str(min_quality + phred_qscore_offset - 33),       # minqual
        "-l", str(min_trimmed_length),
        "-x",                                                    # trims poly-x tails
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
    with sp.Popen(cmd, stderr=sp.STDOUT, stdout=sp.PIPE) as proc:
        out, _ = proc.communicate()
        if proc.returncode:
            out = out.decode()
            print(
                "@@ERROR: "
                f"\nCMD: {' '.join(cmd)}"
                f"\nerr: {out}",
                flush=True
            )
            raise ValueError(out)

        # callback sent to logger.info on completion
        print(f"@@DEBUG: CMD: {' '.join(cmd)}", flush=True)
        print(f"@@INFO: finished trimming {sname}", flush=True)

    # parse JSON stats
    with open(stats_json, 'r', encoding="utf-8") as indata:
        jdata = json.loads(indata.read())
    return jdata, (out1, out2)


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
    workers: int,
    threads: int,
    name_parse: Tuple[str, str] | None,
    umi_tag_in_i5: bool,
    force: bool,
):
    # ------------------------------------------------------------
    # parse dict of {name: (r1, r2)}
    # fastq_dict = get_fastq_tuples_dict_from_paths_list(fastqs)
    fastq_dict = get_name_to_fastq_dict(fastqs, name_parse)

    # check outdir for existing and raise or remove
    result_files = [outdir / f"{sname}.R1.trimmed.fastq.gz" for sname in fastq_dict]
    if any(i.exists() for i in result_files):
        if not force:
            raise IPyradError(f"Trimmed fastqs exist in outdir: e.g., {result_files[0]}. Use --force to overwrite.")

    # ------------------------------------------------------------
    # infer restriction overhangs by kmer analysis
    if not disable_infer_re_overhangs:
        max_reads_per_sample = int(max_reads_kmer / len(fastq_dict))
        re1, re2 = infer_re_overhangs(fastq_dict, max_reads_per_sample)

        # if user also entered REs then check them here.
        if restriction_overhangs:
            check_user_re_overhangs((re1, re2), restriction_overhangs)
            re1, re2 = restriction_overhangs
    elif restriction_overhangs:
        re1, re2 = restriction_overhangs
    else:
        logger.warning("not trimming restriction site overhangs")
        re1, re2 = ("", "")
    logger.info(f"using restriction site overhangs: {re1} {re2}")

    # ------------------------------------------------------------
    jobs = {}
    logger.info(f"trimming/filtering {len(fastq_dict)} inputs to trimmed fastqs in {outdir}")
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
            threads=max(1, threads - 2),  # uses 2 I/O threads + requested threads
        )
        jobs[sname] = kwargs
    _ = run_with_pool(trim_sample_with_fastp, jobs, workers)



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
