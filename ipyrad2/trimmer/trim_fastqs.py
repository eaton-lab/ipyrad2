#!/usr/bin/env python

"""Methods for quality and adapter trimming RAD-seq reads.

1. Check file paths and search for sample pairs.
2. Check for cutsite motifs using kmers.
3. Trim reads with fastp.
4. Write stats summary file.

Example
--------
python trim_fastqs.py FQs ...
"""

from __future__ import annotations

from typing import Tuple, List, Dict, Any, Sequence
import os
import sys
import json
import gzip
from pathlib import Path
from loguru import logger
import pandas as pd
from ..utils.names import get_name_to_fastq_dict
from ..utils.kmers import InferredJunctionSet, get_overhangs_from_kmers, validate_named_motif_list
from ..utils.exceptions import IPyradError
from ..utils.parallel import run_with_pool, run_pipeline


FASTP_BINARY = Path(sys.prefix) / "bin" / "fastp"
ADAPTERS = Path(__file__).absolute().parent / "adapters.fa"


def _require(condition: bool, message: str) -> None:
    """Raise a user-facing error when a trim precondition is not met."""
    if not condition:
        raise IPyradError(message)


def _validate_trim_config(
    max_reads: int | None,
    min_trimmed_length: int,
    min_quality: int,
    max_unqualified_percent: int,
    min_mean_window_quality: int,
    cut_window_size: int,
    max_reads_kmer: int,
    max_ns: int,
    cores: int,
    threads: int,
) -> None:
    """Validate runtime prerequisites and user-entered numeric arguments."""
    _require(
        FASTP_BINARY.is_file(),
        f"fastp binary was not found at '{FASTP_BINARY}'. Activate the ipyrad2 environment or install fastp.",
    )
    _require(
        os.access(FASTP_BINARY, os.X_OK),
        f"fastp binary is not executable: '{FASTP_BINARY}'.",
    )
    _require(
        ADAPTERS.is_file(),
        f"Adapter FASTA was not found at '{ADAPTERS}'.",
    )
    _require(cores >= 1, "cores must be >= 1.")
    _require(threads >= 1, "threads must be >= 1.")
    _require(threads <= cores, "threads cannot exceed cores.")
    _require(min_trimmed_length >= 1, "min_trimmed_length must be >= 1.")
    _require(min_quality >= 0, "min_quality must be >= 0.")
    _require(
        0 <= max_unqualified_percent <= 100,
        "max_unqualified_percent must be between 0 and 100.",
    )
    _require(
        1 <= min_mean_window_quality <= 36,
        "min_mean_window_quality must be between 1 and 36.",
    )
    _require(
        1 <= cut_window_size <= 1000,
        "cut_window_size must be between 1 and 1000.",
    )
    _require(max_ns >= 0, "max_ns must be >= 0.")
    _require(max_reads_kmer >= 1, "max_reads_kmer must be >= 1.")
    if max_reads is not None:
        _require(max_reads >= 1, "max_reads must be >= 1 when set.")


def _validate_user_cutsite_motifs(
    cutsite_1: str | None,
    cutsite_2: str | None,
) -> Tuple[Tuple[str, ...], Tuple[str, ...]] | None:
    """Validate explicit user-provided cutsite motifs."""
    if cutsite_1 is None and cutsite_2 is None:
        return None
    return (
        validate_named_motif_list(cutsite_1, "R1 cutsite motif"),
        validate_named_motif_list(
            cutsite_2,
            "R2 cutsite motif",
            allow_empty=True,
        ),
    )


def _manual_junction_set(
    motifs: Sequence[str],
    *,
    offset: int = 0,
) -> InferredJunctionSet:
    """Build junction-set metadata for explicit user-entered cutsite motifs."""
    motifs = tuple(motifs)
    return InferredJunctionSet(
        motifs=motifs,
        motif_counts=tuple(0 for _ in motifs),
        offset=offset,
        total_support=0,
        runner_up_offset_support=0,
        candidate_offsets=(offset,),
    )


def _format_motif_set(junction: InferredJunctionSet) -> str:
    """Return a short human-readable description of a motif set."""
    motifs = ", ".join(junction.motifs) if junction.motifs else "<none>"
    return f"[{motifs}] at offset {junction.offset}"


def _warn_multi_motif_inference(
    read_label: str,
    junction: InferredJunctionSet,
    max_reads_kmer: int,
) -> None:
    """Warn when multiple motifs were auto-detected on one read end."""
    if len(junction.motifs) <= 1:
        return
    logger.warning(
        "{} cutsite motif inference found multiple motifs {}. "
        "This can reflect low-quality data or too few sampled reads; consider increasing "
        "--max-reads-kmer from {}. It can also be valid multi-enzyme data such as 3RAD. "
        "Trim will still use the longest inferred junction length. Enter cutsite motifs "
        "manually to suppress this warning.",
        read_label,
        junction.motifs,
        max_reads_kmer,
    )


def _open_trim_input(path: Path):
    """Open one trim input in binary mode after validating its compression."""
    if path.suffix == ".bz2":
        raise IPyradError(
            f"Trim supports only plain FASTQ or .gz-compressed FASTQ inputs: {path}"
        )
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        return opener(path, "rb")
    except OSError as err:
        raise IPyradError(f"Failed to read FASTQ input: {path}") from err


def _fastq_has_complete_first_record(path: Path) -> bool:
    """Return True when a FASTQ has at least one complete record."""
    with _open_trim_input(path) as infile:
        header = infile.readline()
        if header == b"":
            return False
        sequence = infile.readline()
        plus = infile.readline()
        quality = infile.readline()
    if not sequence or not plus or not quality:
        raise IPyradError(f"FASTQ input is truncated or incomplete at the first record: {path}")
    if not header.startswith(b"@"):
        raise IPyradError(f"FASTQ input does not start with a '@' header line: {path}")
    if not plus.startswith(b"+"):
        raise IPyradError(f"FASTQ input is missing the '+' separator line in the first record: {path}")
    return True


def _classify_trim_sample(
    sname: str,
    fastqs: Tuple[Path, Path | None],
) -> Tuple[str, str | None]:
    """Classify one sample as usable or empty."""
    r1_has_reads = _fastq_has_complete_first_record(fastqs[0])
    r2 = fastqs[1]
    if r2 is None or r2.name == "-null-":
        if r1_has_reads:
            return "usable", None
        return "empty_single_end", (
            f"skipping sample '{sname}' because its input FASTQ is empty: {fastqs[0]}"
        )

    r2_has_reads = _fastq_has_complete_first_record(r2)
    if r1_has_reads and r2_has_reads:
        return "usable", None
    if not r1_has_reads and not r2_has_reads:
        return "empty_paired_both", (
            f"skipping sample '{sname}' because both paired FASTQs are empty: "
            f"{fastqs[0]}, {r2}"
        )

    empty_mates = []
    if not r1_has_reads:
        empty_mates.append(f"R1={fastqs[0]}")
    if not r2_has_reads:
        empty_mates.append(f"R2={r2}")
    return "empty_paired_one_mate", (
        f"skipping sample '{sname}' because one paired FASTQ is empty: "
        f"{', '.join(empty_mates)}"
    )


def _partition_usable_samples(
    fastq_dict: Dict[str, Tuple[Path, Path | None]],
) -> Tuple[Dict[str, Tuple[Path, Path | None]], Dict[str, str]]:
    """Split parsed samples into usable and skipped subsets."""
    usable: Dict[str, Tuple[Path, Path | None]] = {}
    skipped: Dict[str, str] = {}
    for sname, fastq_tuple in fastq_dict.items():
        status, reason = _classify_trim_sample(sname, fastq_tuple)
        if status == "usable":
            usable[sname] = fastq_tuple
            continue
        skipped[sname] = reason or f"skipping sample '{sname}' for an unknown reason"
    return usable, skipped


def _sample_output_artifacts(
    sname: str,
    fastq_tuple: Tuple[Path, Path | None],
    outdir: Path,
) -> Tuple[Path, ...]:
    """Return all output artifacts generated for one sample."""
    artifacts = [outdir / f"{sname}.R1.trimmed.fastq.gz"]
    if fastq_tuple[1] is not None and fastq_tuple[1].name != "-null-":
        artifacts.append(outdir / f"{sname}.R2.trimmed.fastq.gz")
    artifacts.extend([
        outdir / f"{sname}.stats.json",
        outdir / f"{sname}.stats.html",
    ])
    return tuple(artifacts)


def _find_existing_trim_artifact(
    fastq_dict: Dict[str, Tuple[Path, Path | None]],
    outdir: Path,
) -> Path | None:
    """Return the first existing trim artifact that would be overwritten."""
    for sname in sorted(fastq_dict):
        for artifact in _sample_output_artifacts(sname, fastq_dict[sname], outdir):
            if artifact.exists():
                return artifact
    return None


def _build_fastp_command(
    fastqs: Tuple[Path, Path | None],
    sname: str,
    outdir: Path,
    cutsite_motifs: Tuple[str, str],
    trim_front_lengths: Tuple[int, int] | None,
    max_reads: int | None,
    min_trimmed_length: int,
    min_quality: int,
    max_unqualified_percent: int,
    min_mean_window_quality: int,
    cut_window_size: int,
    max_ns: int,
    phred64: bool,
    disable_adapter_trimming: bool,
    disable_quality_filtering: bool,
    umi_tag_in_i5: bool,
    threads: int,
) -> List[str]:
    """Build a fastp command for a single sample."""
    out1 = outdir / f"{sname}.R1.trimmed.fastq.gz"
    out2 = outdir / f"{sname}.R2.trimmed.fastq.gz"
    stats_html = outdir / f"{sname}.stats.html"
    stats_json = outdir / f"{sname}.stats.json"
    is_paired = fastqs[1] is not None and fastqs[1].name != "-null-"
    trim_front_lengths = trim_front_lengths or (
        len(cutsite_motifs[0]),
        len(cutsite_motifs[1]),
    )

    if is_paired:
        cmd = [
            str(FASTP_BINARY),
            "-i", str(fastqs[0]),
            "-I", str(fastqs[1]),
            "-o", str(out1),
            "-O", str(out2),
            "--detect_adapter_for_pe",
            "-c",
            "-w", str(threads),
            "--trim_front1", str(trim_front_lengths[0]),
            "--trim_front2", str(trim_front_lengths[1]),
            "--overlap_len_require", "20",
            "--overlap_diff_limit", "5",
        ]
        if umi_tag_in_i5:
            cmd.extend(["-U", "--umi_loc=index2", "--umi_prefix=UMI"])
    else:
        cmd = [
            str(FASTP_BINARY),
            "-i", str(fastqs[0]),
            "-o", str(out1),
            "-w", str(threads),
            "--trim_front1", str(trim_front_lengths[0]),
        ]

    cmd.extend([
        "-q", str(min_quality),
        "-u", str(max_unqualified_percent),
        "-M", str(min_mean_window_quality),
        "-W", str(cut_window_size),
        "--cut_front", "--cut_front_window_size", "5",
        "--cut_tail", "--cut_tail_window_size", "5",
        "--length_required", str(min_trimmed_length),
        "--trim_poly_g",
        "--trim_poly_x",
        "--poly_x_min_len", "10",
        "--low_complexity_filter",
        "--n_base_limit", str(max_ns),
        "--json", str(stats_json),
        "--html", str(stats_html),
        "--adapter_fasta", str(ADAPTERS),
    ])

    if max_reads is not None:
        cmd.extend(["--reads_to_process", str(max_reads)])
    if phred64:
        cmd.append("-6")
    if disable_adapter_trimming:
        cmd.append("-A")
    if disable_quality_filtering:
        cmd.append("-Q")
    return cmd


def _resolve_cutsite_motifs(
    fastq_dict: Dict[str, Tuple[Path, Path | None]],
    cutsite_motifs: Tuple[Tuple[str, ...], Tuple[str, ...]] | None,
    disable_infer_cutsite_motifs: bool,
    max_reads_kmer: int,
    cores: int,
    log_level: str,
) -> Tuple[InferredJunctionSet, InferredJunctionSet]:
    """Infer or accept user-defined cutsite motifs."""
    user_r1 = cutsite_motifs[0] if cutsite_motifs else ()
    user_r2 = cutsite_motifs[1] if cutsite_motifs else ()

    if user_r1:
        re1 = _manual_junction_set(user_r1)
    elif disable_infer_cutsite_motifs:
        re1 = _manual_junction_set(())
    else:
        re1 = get_overhangs_from_kmers(
            [reads[0] for reads in fastq_dict.values()],
            20,
            max_reads_kmer,
            cores,
            log_level,
            candidate_offsets=(0, 1),
            label="R1 cutsite motif inference",
        )
        _warn_multi_motif_inference("R1", re1, max_reads_kmer)

    read2s = [reads[1] for reads in fastq_dict.values() if reads[1] is not None]
    if user_r2:
        re2 = _manual_junction_set(user_r2)
    elif not read2s or disable_infer_cutsite_motifs:
        re2 = _manual_junction_set(())
    else:
        re2 = get_overhangs_from_kmers(
            read2s,
            20,
            max_reads_kmer,
            cores,
            log_level,
            candidate_offsets=(0, 1),
            label="R2 cutsite motif inference",
        )
        _warn_multi_motif_inference("R2", re2, max_reads_kmer)

    return re1, re2


def _load_stats_json(stats_file: Path) -> Dict[str, Any]:
    """Load a fastp stats JSON file or raise a clear user-facing error."""
    if not stats_file.exists():
        raise IPyradError(f"Missing fastp stats report: {stats_file}")
    try:
        with open(stats_file, "r", encoding="utf-8") as indata:
            return json.loads(indata.read())
    except OSError as err:
        raise IPyradError(f"Failed to read fastp stats report: {stats_file}") from err
    except json.JSONDecodeError as err:
        raise IPyradError(f"Failed to parse fastp stats report: {stats_file}") from err


def trim_sample_with_fastp(
    fastqs: Tuple[Path, Path | None],
    sname: str,
    outdir: Path,
    cutsite_motifs: Tuple[str, str],
    trim_front_lengths: Tuple[int, int] | None,
    max_reads: int | None,
    min_trimmed_length: int,
    min_quality: int,
    max_unqualified_percent: int,
    min_mean_window_quality: int,
    cut_window_size: int,
    max_ns: int,
    phred64: bool,
    disable_adapter_trimming: bool,
    disable_quality_filtering: bool,
    umi_tag_in_i5: bool,
    threads: int,
) -> None:
    """Run FASTP and write stats per sample to outdir as json
    """
    cmd = _build_fastp_command(
        fastqs=fastqs,
        sname=sname,
        outdir=outdir,
        cutsite_motifs=cutsite_motifs,
        trim_front_lengths=trim_front_lengths,
        max_reads=max_reads,
        min_trimmed_length=min_trimmed_length,
        min_quality=min_quality,
        max_unqualified_percent=max_unqualified_percent,
        min_mean_window_quality=min_mean_window_quality,
        cut_window_size=cut_window_size,
        max_ns=max_ns,
        phred64=phred64,
        disable_adapter_trimming=disable_adapter_trimming,
        disable_quality_filtering=disable_quality_filtering,
        umi_tag_in_i5=umi_tag_in_i5,
        threads=threads,
    )

    # run the command in subprocess
    logger.debug(f"CMD: {' '.join(cmd)}")
    input_paths = ", ".join(
        str(path) for path in fastqs
        if path is not None and path.name != "-null-"
    )
    try:
        run_pipeline([cmd])
    except Exception as err:
        raise IPyradError(
            f"fastp failed for sample '{sname}' on input(s) {input_paths}: {err}"
        ) from err
    logger.debug(f"finished trimming {sname}")


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
        jdata[sname] = _load_stats_json(stats_file)

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
        df.loc[sname, "total_reads_after"] = j["summary"]["after_filtering"]["total_reads"]
        df.loc[sname, "total_bases_after"] = j["summary"]["after_filtering"]["total_bases"]
        df.loc[sname, "q20_rate_after"] = j["summary"]["after_filtering"]["q20_rate"]
        df.loc[sname, "q30_rate_after"] = j["summary"]["after_filtering"]["q30_rate"]
        df.loc[sname, "read1_mean_length_after"] = j["summary"]["after_filtering"]["read1_mean_length"]
        df.loc[sname, "reads_filtered_by_low_quality"] = j["filtering_result"]["low_quality_reads"]
        df.loc[sname, "reads_filtered_by_too_many_N"] = j["filtering_result"]["too_many_N_reads"]
        df.loc[sname, "reads_filtered_by_low_complexity"] = j["filtering_result"]["low_complexity_reads"]
        df.loc[sname, "reads_filtered_by_too_short"] = j["filtering_result"]["too_short_reads"]
        try:
            df.loc[sname, "adapter_trimmed_reads"] = j["adapter_cutting"]["adapter_trimmed_reads"]
            df.loc[sname, "adapter_trimmed_bases"] = j["adapter_cutting"]["adapter_trimmed_bases"]
        except KeyError:
            # SE data, no adapter trimming info
            pass
        try:
            df.loc[sname, "read2_mean_length_before"] = j["summary"]["before_filtering"]["read2_mean_length"]
            df.loc[sname, "read2_mean_length_after"] = j["summary"]["after_filtering"]["read2_mean_length"]
        except KeyError:
            # If SE, no read2 filtering
            pass

    # Drop na columns for SE data
    df = df.dropna(axis=1)

    # write human readable whitespace delimited.
    df.to_string(outfile, float_format=lambda x: f"{x:.6f}")
    logger.info(f"trimming stats written to {outfile}")


def run_trimmer(
    fastqs: List[Path],
    outdir: Path,
    cutsite_motifs: Sequence[str] | None,
    max_reads: int | None,
    min_trimmed_length: int,
    max_unqualified_percent: int,
    min_quality: int,
    min_mean_window_quality: int,
    cut_window_size: int,
    phred64: bool,
    max_reads_kmer: int,
    max_ns: int,
    disable_infer_cutsite_motifs: bool,
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
    cutsite_motifs = tuple(cutsite_motifs) if cutsite_motifs else None
    cutsite_motifs = _validate_user_cutsite_motifs(
        cutsite_motifs[0] if cutsite_motifs else None,
        cutsite_motifs[1] if cutsite_motifs else None,
    )
    _validate_trim_config(
        max_reads=max_reads,
        min_trimmed_length=min_trimmed_length,
        min_quality=min_quality,
        max_unqualified_percent=max_unqualified_percent,
        min_mean_window_quality=min_mean_window_quality,
        cut_window_size=cut_window_size,
        max_reads_kmer=max_reads_kmer,
        max_ns=max_ns,
        cores=cores,
        threads=threads,
    )

    # ------------------------------------------------------------
    # parse dict of {name: (r1, r2)}
    # fastq_dict = get_fastq_tuples_dict_from_paths_list(fastqs)
    fastq_dict = get_name_to_fastq_dict(fastqs, delim_str, delim_idx, suffix)
    fastq_dict, skipped_samples = _partition_usable_samples(fastq_dict)
    logger.info(
        "trim input preflight found {} usable samples and {} skipped empty samples",
        len(fastq_dict),
        len(skipped_samples),
    )
    for message in skipped_samples.values():
        logger.warning(message)
    if not fastq_dict:
        raise IPyradError("No non-empty FASTQ samples remain after input validation.")

    # check outdir for existing and raise or remove
    outdir = outdir.expanduser().absolute()
    existing_artifact = _find_existing_trim_artifact(fastq_dict, outdir)
    if existing_artifact is not None and not force:
        raise IPyradError(
            f"Trim output artifact exists in outdir: {existing_artifact}. Use --force to overwrite."
        )
    outdir.mkdir(parents=True, exist_ok=True)

    # run at most this many concurrent jobs
    workers = max(1, cores // threads)

    # ------------------------------------------------------------
    # infer cutsite motifs by kmer analysis
    re1, re2 = _resolve_cutsite_motifs(
        fastq_dict=fastq_dict,
        cutsite_motifs=cutsite_motifs,
        disable_infer_cutsite_motifs=disable_infer_cutsite_motifs,
        max_reads_kmer=max_reads_kmer,
        cores=cores,
        log_level=log_level,
    )
    logger.info(
        "cutsite motifs set to R1={} R2={}",
        _format_motif_set(re1),
        _format_motif_set(re2),
    )

    # ------------------------------------------------------------
    jobs = {}
    logger.info(f"trimming/filtering {len(fastq_dict)} samples with 'fastp' and writing to {outdir}")
    logger.info(f"running up to {workers} parallel jobs each using up to {threads} threads")
    for sname, fastq_tuple in fastq_dict.items():
        kwargs = dict(
            fastqs=fastq_tuple,
            sname=sname,
            outdir=outdir,
            cutsite_motifs=(re1.primary_motif, re2.primary_motif),
            trim_front_lengths=(re1.trim_length, re2.trim_length),
            max_reads=max_reads,
            max_ns=max_ns,
            min_trimmed_length=min_trimmed_length,
            min_quality=min_quality,
            max_unqualified_percent=max_unqualified_percent,
            min_mean_window_quality=min_mean_window_quality,
            cut_window_size=cut_window_size,
            phred64=phred64,
            disable_adapter_trimming=disable_adapter_trimming,
            disable_quality_filtering=disable_quality_filtering,
            umi_tag_in_i5=umi_tag_in_i5,
            threads=threads,  # recommended >=3 since 2 are used for i/o
        )
        jobs[sname] = (trim_sample_with_fastp, kwargs)
    results = run_with_pool(jobs, log_level, workers, msg="Trimming")
    write_stats_summary(sorted(results), outdir)


# if __name__ == "__main__":

#     PATHS = sorted(Path("/home/deren/Documents/tools/ipyrad2/examples/Pedic-PE-ddRAD/").glob("*.gz"))
#     # fastq_dict = get_fastq_tuples_dict_from_paths_list(PATHS)
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
