import bz2
import gzip
from pathlib import Path

import pytest

import ipyrad2.trimmer.trim_fastqs as trim_fastqs
from ipyrad2.utils.exceptions import IPyradError
from ipyrad2.utils.kmers import InferredJunction, InferredJunctionSet


def _write_file(path: Path, text: str = "") -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _write_executable(path: Path) -> Path:
    path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    path.chmod(0o755)
    return path


def _write_fastq(path: Path, records: list[tuple[str, str, str]]) -> Path:
    if path.suffix == ".gz":
        opener = gzip.open
    elif path.suffix == ".bz2":
        opener = bz2.open
    else:
        opener = open
    with opener(path, "wt", encoding="utf-8") as out:
        for name, seq, qual in records:
            out.write(f"@{name}\n{seq}\n+\n{qual}\n")
    return path


def _inferred(
    sequence: str,
    *,
    offset: int = 0,
    winner_count: int = 100,
    runner_up_count: int = 0,
    candidate_offsets: tuple[int, ...] = (0, 1),
) -> InferredJunction:
    return InferredJunction(
        sequence=sequence,
        offset=offset,
        k=len(sequence),
        winner_count=winner_count,
        runner_up_count=runner_up_count,
        candidate_offsets=candidate_offsets,
    )


def _junction_set(
    motifs: tuple[str, ...],
    *,
    offset: int = 0,
    counts: tuple[int, ...] | None = None,
    runner_up_offset_support: int = 0,
    candidate_offsets: tuple[int, ...] = (0, 1),
) -> InferredJunctionSet:
    counts = counts or tuple(100 for _ in motifs)
    return InferredJunctionSet(
        motifs=motifs,
        motif_counts=counts,
        offset=offset,
        total_support=sum(counts),
        runner_up_offset_support=runner_up_offset_support,
        candidate_offsets=candidate_offsets,
    )


def test_build_fastp_command_for_paired_end_includes_expected_flags(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastp_binary = _write_executable(tmp_path / "fastp")
    adapters = _write_file(tmp_path / "adapters.fa", ">adapter\nACGT\n")
    monkeypatch.setattr(trim_fastqs, "FASTP_BINARY", fastp_binary)
    monkeypatch.setattr(trim_fastqs, "ADAPTERS", adapters)

    cmd = trim_fastqs._build_fastp_command(
        fastqs=(Path("sample_R1.fastq.gz"), Path("sample_R2.fastq.gz")),
        sname="sample",
        outdir=tmp_path,
        cutsite_motifs=("TGCAG", "CGAT"),
        trim_front_lengths=None,
        max_reads=123,
        min_trimmed_length=35,
        min_quality=20,
        max_unqualified_percent=15,
        min_mean_window_quality=30,
        cut_window_size=5,
        max_ns=5,
        phred64=True,
        disable_adapter_trimming=True,
        disable_quality_filtering=True,
        umi_tag_in_i5=True,
        threads=4,
    )

    assert cmd[0] == str(fastp_binary)
    assert ["-i", "sample_R1.fastq.gz"] == cmd[1:3]
    assert "-I" in cmd
    assert "-O" in cmd
    assert "--detect_adapter_for_pe" in cmd
    assert "--trim_front1" in cmd and cmd[cmd.index("--trim_front1") + 1] == "5"
    assert "--trim_front2" in cmd and cmd[cmd.index("--trim_front2") + 1] == "4"
    assert "--reads_to_process" in cmd and cmd[cmd.index("--reads_to_process") + 1] == "123"
    assert "-6" in cmd
    assert "-A" in cmd
    assert "-Q" in cmd
    assert "-q" not in cmd
    assert "-u" not in cmd
    assert "-M" not in cmd
    assert "-W" not in cmd
    assert "--cut_front" not in cmd
    assert "--cut_tail" not in cmd
    assert "-U" in cmd
    assert "--umi_loc=index2" in cmd
    assert "--umi_prefix=UMI" in cmd
    assert "--adapter_fasta" in cmd
    assert cmd[cmd.index("--adapter_fasta") + 1] == str(adapters)
    assert "--html" in cmd
    assert cmd[cmd.index("--html") + 1] == str(tmp_path / "sample.stats.html")


def test_build_fastp_command_for_single_end_omits_paired_flags(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastp_binary = _write_executable(tmp_path / "fastp")
    adapters = _write_file(tmp_path / "adapters.fa", ">adapter\nACGT\n")
    monkeypatch.setattr(trim_fastqs, "FASTP_BINARY", fastp_binary)
    monkeypatch.setattr(trim_fastqs, "ADAPTERS", adapters)

    cmd = trim_fastqs._build_fastp_command(
        fastqs=(Path("sample.fastq.gz"), None),
        sname="sample",
        outdir=tmp_path,
        cutsite_motifs=("TGCAG", ""),
        trim_front_lengths=None,
        max_reads=None,
        min_trimmed_length=35,
        min_quality=20,
        max_unqualified_percent=15,
        min_mean_window_quality=30,
        cut_window_size=5,
        max_ns=5,
        phred64=False,
        disable_adapter_trimming=False,
        disable_quality_filtering=False,
        umi_tag_in_i5=True,
        threads=3,
    )

    assert cmd[0] == str(fastp_binary)
    assert "-I" not in cmd
    assert "-O" not in cmd
    assert "--detect_adapter_for_pe" not in cmd
    assert "-U" not in cmd
    assert "--trim_front1" in cmd
    assert "--trim_front2" not in cmd
    assert "-6" not in cmd
    assert "--html" in cmd
    assert cmd[cmd.index("--html") + 1] == str(tmp_path / "sample.stats.html")


def test_build_fastp_command_uses_explicit_trim_front_lengths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastp_binary = _write_executable(tmp_path / "fastp")
    adapters = _write_file(tmp_path / "adapters.fa", ">adapter\nACGT\n")
    monkeypatch.setattr(trim_fastqs, "FASTP_BINARY", fastp_binary)
    monkeypatch.setattr(trim_fastqs, "ADAPTERS", adapters)

    cmd = trim_fastqs._build_fastp_command(
        fastqs=(Path("sample.fastq.gz"), None),
        sname="sample",
        outdir=tmp_path,
        cutsite_motifs=("ATCGG", ""),
        trim_front_lengths=(6, 0),
        max_reads=None,
        min_trimmed_length=35,
        min_quality=20,
        max_unqualified_percent=15,
        min_mean_window_quality=30,
        cut_window_size=5,
        max_ns=5,
        phred64=False,
        disable_adapter_trimming=False,
        disable_quality_filtering=False,
        umi_tag_in_i5=False,
        threads=1,
    )

    assert "--trim_front1" in cmd and cmd[cmd.index("--trim_front1") + 1] == "6"


def test_validate_trim_config_requires_fastp_binary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    adapters = _write_file(tmp_path / "adapters.fa", ">adapter\nACGT\n")
    monkeypatch.setattr(trim_fastqs, "FASTP_BINARY", tmp_path / "missing-fastp")
    monkeypatch.setattr(trim_fastqs, "ADAPTERS", adapters)

    with pytest.raises(IPyradError, match="fastp binary was not found"):
        trim_fastqs._validate_trim_config(
            max_reads=None,
            min_trimmed_length=35,
            min_quality=20,
            max_unqualified_percent=15,
            min_mean_window_quality=30,
            cut_window_size=5,
            max_reads_kmer=500_000,
            max_ns=5,
            cores=6,
            threads=3,
        )


def test_write_stats_summary_raises_when_stats_are_missing(tmp_path: Path) -> None:
    with pytest.raises(IPyradError, match="Missing fastp stats report"):
        trim_fastqs.write_stats_summary(["sampleA"], tmp_path)


def test_fastq_has_complete_first_record_detects_empty_plain_fastq(tmp_path: Path) -> None:
    fastq = _write_file(tmp_path / "empty.fastq", "")

    assert trim_fastqs._fastq_has_complete_first_record(fastq) is False


def test_fastq_has_complete_first_record_detects_empty_gzipped_fastq(tmp_path: Path) -> None:
    fastq = tmp_path / "empty.fastq.gz"
    with gzip.open(fastq, "wt", encoding="utf-8") as out:
        out.write("")

    assert trim_fastqs._fastq_has_complete_first_record(fastq) is False


def test_fastq_has_complete_first_record_rejects_incomplete_first_record(tmp_path: Path) -> None:
    fastq = _write_file(tmp_path / "truncated.fastq", "@r1\nACGT\n+\n")

    with pytest.raises(IPyradError, match="truncated or incomplete at the first record"):
        trim_fastqs._fastq_has_complete_first_record(fastq)


def test_fastq_has_complete_first_record_rejects_missing_header_marker(tmp_path: Path) -> None:
    fastq = _write_file(tmp_path / "bad_header.fastq", "r1\nACGT\n+\nIIII\n")

    with pytest.raises(IPyradError, match="does not start with a '@' header line"):
        trim_fastqs._fastq_has_complete_first_record(fastq)


def test_fastq_has_complete_first_record_rejects_missing_plus_marker(tmp_path: Path) -> None:
    fastq = _write_file(tmp_path / "bad_plus.fastq", "@r1\nACGT\n=\nIIII\n")

    with pytest.raises(IPyradError, match="missing the '\\+' separator line"):
        trim_fastqs._fastq_has_complete_first_record(fastq)


def test_fastq_has_complete_first_record_rejects_bz2_inputs(tmp_path: Path) -> None:
    fastq = _write_fastq(
        tmp_path / "sample.fastq.bz2",
        [("r1", "ACGT", "IIII")],
    )

    with pytest.raises(IPyradError, match="plain FASTQ or \\.gz-compressed FASTQ inputs"):
        trim_fastqs._fastq_has_complete_first_record(fastq)


def test_run_trimmer_uses_max_reads_kmer(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastp_binary = _write_executable(tmp_path / "fastp")
    adapters = _write_file(tmp_path / "adapters.fa", ">adapter\nACGT\n")
    monkeypatch.setattr(trim_fastqs, "FASTP_BINARY", fastp_binary)
    monkeypatch.setattr(trim_fastqs, "ADAPTERS", adapters)

    r1 = _write_fastq(tmp_path / "sample_R1.fastq.gz", [("r1", "ACGT", "IIII")])
    r2 = _write_fastq(tmp_path / "sample_R2.fastq.gz", [("r1", "TGCA", "IIII")])
    calls: list[int] = []

    monkeypatch.setattr(
        trim_fastqs,
        "get_name_to_fastq_dict",
        lambda fastqs, delim_str, delim_idx, suffix: {"sample": (r1, r2)},
    )

    def fake_get_overhangs_from_kmers(
        fastqs,
        max_len,
        max_reads,
        workers,
        log_level,
        candidate_offsets=None,
        label=None,
    ):
        calls.append(max_reads)
        if all(path.name.endswith("R2.fastq.gz") for path in fastqs):
            return _junction_set(("CGATC",))
        return _junction_set(("TGCAG",))

    monkeypatch.setattr(trim_fastqs, "get_overhangs_from_kmers", fake_get_overhangs_from_kmers)
    monkeypatch.setattr(
        trim_fastqs,
        "run_with_pool",
        lambda jobs, log_level, max_workers=None, max_inflight=None, msg="Processing": {
            key: None for key in jobs
        },
    )
    monkeypatch.setattr(trim_fastqs, "write_stats_summary", lambda snames, outdir: None)

    trim_fastqs.run_trimmer(
        fastqs=[Path("ignored.fastq.gz")],
        outdir=tmp_path / "trimmed",
        cutsite_motifs=None,
        max_reads=None,
        min_trimmed_length=35,
        max_unqualified_percent=15,
        min_quality=20,
        min_mean_window_quality=30,
        cut_window_size=5,
        phred64=False,
        max_reads_kmer=321_321,
        max_ns=5,
        disable_infer_cutsite_motifs=False,
        disable_adapter_trimming=False,
        disable_quality_filtering=False,
        cores=6,
        threads=3,
        delim_str=None,
        delim_idx=1,
        suffix=None,
        umi_tag_in_i5=False,
        force=False,
        log_level="INFO",
    )

    assert calls == [321_321, 321_321]


def test_partition_usable_samples_skips_empty_inputs(
    tmp_path: Path,
) -> None:
    usable_r1 = _write_fastq(tmp_path / "usable_R1.fastq.gz", [("r1", "ACGT", "IIII")])
    usable_r2 = _write_fastq(tmp_path / "usable_R2.fastq.gz", [("r1", "TGCA", "IIII")])
    empty_single = _write_file(tmp_path / "empty.fastq", "")
    half_empty_r1 = _write_fastq(tmp_path / "half_R1.fastq.gz", [("r1", "ACGT", "IIII")])
    half_empty_r2 = tmp_path / "half_R2.fastq.gz"
    with gzip.open(half_empty_r2, "wt", encoding="utf-8") as out:
        out.write("")

    usable, skipped = trim_fastqs._partition_usable_samples(
        {
            "usable": (usable_r1, usable_r2),
            "empty": (empty_single, None),
            "half": (half_empty_r1, half_empty_r2),
        }
    )

    assert usable == {"usable": (usable_r1, usable_r2)}
    assert "empty" in skipped and "input FASTQ is empty" in skipped["empty"]
    assert "half" in skipped and "one paired FASTQ is empty" in skipped["half"]


def test_run_trimmer_logs_overhang_summary_once(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastp_binary = _write_executable(tmp_path / "fastp")
    adapters = _write_file(tmp_path / "adapters.fa", ">adapter\nACGT\n")
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", [("r1", "ACGT", "IIII")])
    messages: list[str] = []

    stub_logger = type(
        "LoggerStub",
        (),
        {
            "info": staticmethod(lambda *args: messages.append(args[0].format(*args[1:]))),
            "warning": staticmethod(lambda *args: None),
            "debug": staticmethod(lambda *args: None),
        },
    )

    monkeypatch.setattr(trim_fastqs, "FASTP_BINARY", fastp_binary)
    monkeypatch.setattr(trim_fastqs, "ADAPTERS", adapters)
    monkeypatch.setattr(trim_fastqs, "logger", stub_logger)
    monkeypatch.setattr(
        trim_fastqs,
        "run_with_pool",
        lambda jobs, log_level, max_workers=None, max_inflight=None, msg="Processing": {
            key: None for key in jobs
        },
    )
    monkeypatch.setattr(trim_fastqs, "write_stats_summary", lambda snames, outdir: None)
    monkeypatch.setattr(
        trim_fastqs,
        "get_overhangs_from_kmers",
        lambda *args, **kwargs: _junction_set(("TGCAG",), offset=0),
    )

    trim_fastqs.run_trimmer(
        fastqs=[fastq],
        outdir=tmp_path / "trimmed",
        cutsite_motifs=None,
        max_reads=None,
        min_trimmed_length=1,
        max_unqualified_percent=15,
        min_quality=20,
        min_mean_window_quality=30,
        cut_window_size=5,
        phred64=False,
        max_reads_kmer=100,
        max_ns=5,
        disable_infer_cutsite_motifs=False,
        disable_adapter_trimming=False,
        disable_quality_filtering=False,
        cores=1,
        threads=1,
        delim_str=None,
        delim_idx=1,
        suffix=None,
        umi_tag_in_i5=False,
        force=False,
        log_level="INFO",
    )

    assert sum("cutsite motifs set to" in message for message in messages) == 1


def test_run_trimmer_logs_preflight_summary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastp_binary = _write_executable(tmp_path / "fastp")
    adapters = _write_file(tmp_path / "adapters.fa", ">adapter\nACGT\n")
    usable = _write_fastq(tmp_path / "usable.fastq.gz", [("r1", "ACGT", "IIII")])
    empty = _write_file(tmp_path / "empty.fastq", "")
    messages: list[str] = []

    stub_logger = type(
        "LoggerStub",
        (),
        {
            "info": staticmethod(lambda *args: messages.append(args[0].format(*args[1:]))),
            "warning": staticmethod(lambda *args: None),
            "debug": staticmethod(lambda *args: None),
        },
    )

    monkeypatch.setattr(trim_fastqs, "FASTP_BINARY", fastp_binary)
    monkeypatch.setattr(trim_fastqs, "ADAPTERS", adapters)
    monkeypatch.setattr(trim_fastqs, "logger", stub_logger)
    monkeypatch.setattr(
        trim_fastqs,
        "run_with_pool",
        lambda jobs, log_level, max_workers=None, max_inflight=None, msg="Processing": {
            key: None for key in jobs
        },
    )
    monkeypatch.setattr(trim_fastqs, "write_stats_summary", lambda snames, outdir: None)

    trim_fastqs.run_trimmer(
        fastqs=[usable, empty],
        outdir=tmp_path / "trimmed",
        cutsite_motifs=("TGCAG", ""),
        max_reads=None,
        min_trimmed_length=1,
        max_unqualified_percent=15,
        min_quality=20,
        min_mean_window_quality=30,
        cut_window_size=5,
        phred64=False,
        max_reads_kmer=100,
        max_ns=5,
        disable_infer_cutsite_motifs=True,
        disable_adapter_trimming=False,
        disable_quality_filtering=False,
        cores=1,
        threads=1,
        delim_str=None,
        delim_idx=1,
        suffix=None,
        umi_tag_in_i5=False,
        force=False,
        log_level="INFO",
    )

    assert any(
        "trim input preflight found 1 usable samples and 1 skipped empty samples" in message
        for message in messages
    )


def test_run_trimmer_raises_when_stats_artifact_exists_without_force(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastp_binary = _write_executable(tmp_path / "fastp")
    adapters = _write_file(tmp_path / "adapters.fa", ">adapter\nACGT\n")
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", [("r1", "ACGT", "IIII")])
    outdir = tmp_path / "trimmed"
    outdir.mkdir()
    (outdir / "sample.stats.json").write_text("{}", encoding="utf-8")

    monkeypatch.setattr(trim_fastqs, "FASTP_BINARY", fastp_binary)
    monkeypatch.setattr(trim_fastqs, "ADAPTERS", adapters)

    with pytest.raises(IPyradError, match="Trim output artifact exists in outdir: .*sample\\.stats\\.json"):
        trim_fastqs.run_trimmer(
            fastqs=[fastq],
            outdir=outdir,
            cutsite_motifs=("TGCAG", ""),
            max_reads=None,
            min_trimmed_length=1,
            max_unqualified_percent=15,
            min_quality=20,
            min_mean_window_quality=30,
            cut_window_size=5,
            phred64=False,
            max_reads_kmer=100,
            max_ns=5,
            disable_infer_cutsite_motifs=True,
            disable_adapter_trimming=False,
            disable_quality_filtering=False,
            cores=1,
            threads=1,
            delim_str=None,
            delim_idx=1,
            suffix=None,
            umi_tag_in_i5=False,
            force=False,
            log_level="ERROR",
        )


def test_trim_sample_with_fastp_wraps_pipeline_errors_with_sample_context(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastp_binary = _write_executable(tmp_path / "fastp")
    adapters = _write_file(tmp_path / "adapters.fa", ">adapter\nACGT\n")
    r1 = _write_fastq(tmp_path / "sample_R1.fastq.gz", [("r1", "ACGT", "IIII")])
    r2 = _write_fastq(tmp_path / "sample_R2.fastq.gz", [("r1", "TGCA", "IIII")])

    monkeypatch.setattr(trim_fastqs, "FASTP_BINARY", fastp_binary)
    monkeypatch.setattr(trim_fastqs, "ADAPTERS", adapters)
    monkeypatch.setattr(
        trim_fastqs,
        "run_pipeline",
        lambda cmds: (_ for _ in ()).throw(RuntimeError("pipeline failed (rc=1): ['fastp']\nboom")),
    )

    with pytest.raises(
        IPyradError,
        match="fastp failed for sample 'sample' on input\\(s\\).*sample_R1\\.fastq\\.gz, .*sample_R2\\.fastq\\.gz: pipeline failed",
    ):
        trim_fastqs.trim_sample_with_fastp(
            fastqs=(r1, r2),
            sname="sample",
            outdir=tmp_path,
            cutsite_motifs=("TGCAG", "CGATC"),
            trim_front_lengths=None,
            max_reads=None,
            min_trimmed_length=1,
            min_quality=20,
            max_unqualified_percent=15,
            min_mean_window_quality=30,
            cut_window_size=5,
            max_ns=5,
            phred64=False,
            disable_adapter_trimming=False,
            disable_quality_filtering=False,
            umi_tag_in_i5=False,
            threads=1,
        )


def test_resolve_cutsite_motifs_uses_longest_inferred_motif_length(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastq_dict = {"sample": (tmp_path / "sample_R1.fastq.gz", None)}

    def fake_get_overhangs_from_kmers(
        fastqs,
        max_len,
        max_reads,
        workers,
        log_level,
        candidate_offsets=None,
        label=None,
    ):
        return _junction_set(("ATCGG", "ATCGAT"), offset=1, counts=(40, 20))

    monkeypatch.setattr(trim_fastqs, "get_overhangs_from_kmers", fake_get_overhangs_from_kmers)

    re1, re2 = trim_fastqs._resolve_cutsite_motifs(
        fastq_dict=fastq_dict,
        cutsite_motifs=None,
        disable_infer_cutsite_motifs=False,
        max_reads_kmer=100,
        cores=1,
        log_level="ERROR",
    )

    assert re1.motifs == ("ATCGG", "ATCGAT")
    assert re1.offset == 1
    assert re1.trim_length == 7
    assert re2.motifs == ()


def test_resolve_cutsite_motifs_skips_manual_end_inference(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fastq_dict = {
        "sample": (
            tmp_path / "sample_R1.fastq.gz",
            tmp_path / "sample_R2.fastq.gz",
        )
    }
    calls = []

    def fake_get_overhangs_from_kmers(
        fastqs,
        max_len,
        max_reads,
        workers,
        log_level,
        candidate_offsets=None,
        label=None,
    ):
        calls.append(tuple(path.name for path in fastqs))
        return _junction_set(("CGATC",), offset=0)

    monkeypatch.setattr(trim_fastqs, "get_overhangs_from_kmers", fake_get_overhangs_from_kmers)

    re1, re2 = trim_fastqs._resolve_cutsite_motifs(
        fastq_dict=fastq_dict,
        cutsite_motifs=(("ATCGG", "ATCGAT"), ()),
        disable_infer_cutsite_motifs=False,
        max_reads_kmer=100,
        cores=1,
        log_level="ERROR",
    )

    assert calls == [("sample_R2.fastq.gz",)]
    assert re1.motifs == ("ATCGG", "ATCGAT")
    assert re1.trim_length == 6
    assert re2.motifs == ("CGATC",)


def test_validate_user_cutsite_motifs_rejects_monomorphic_motif() -> None:
    with pytest.raises(IPyradError, match="cannot be monomorphic or all N"):
        trim_fastqs._validate_user_cutsite_motifs("AAAAAA", "")


def test_validate_user_cutsite_motifs_accepts_comma_separated_lists() -> None:
    motifs = trim_fastqs._validate_user_cutsite_motifs("ATCGG,ATCGAT", "CGATC")

    assert motifs == (("ATCGG", "ATCGAT"), ("CGATC",))
