import gzip
import json
import os
from pathlib import Path

import pytest

import ipyrad2.trimmer.trim_fastqs as trim_fastqs
from ipyrad2.utils.kmers import InferredJunction, InferredJunctionSet


TRUSEQ_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
TRUSEQ_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


@pytest.fixture(scope="module", autouse=True)
def require_fastp_binary() -> None:
    assert trim_fastqs.FASTP_BINARY.is_file(), (
        f"fastp binary was not found at {trim_fastqs.FASTP_BINARY}"
    )
    assert os.access(trim_fastqs.FASTP_BINARY, os.X_OK), (
        f"fastp binary is not executable: {trim_fastqs.FASTP_BINARY}"
    )
    assert trim_fastqs.ADAPTERS.is_file(), (
        f"adapter FASTA was not found at {trim_fastqs.ADAPTERS}"
    )


@pytest.fixture
def sequential_pool(monkeypatch: pytest.MonkeyPatch) -> None:
    def run_jobs_sequentially(
        jobs,
        log_level,
        max_workers=None,
        max_inflight=None,
        msg="Processing",
    ):
        results = {}
        for key, (func, kwargs) in jobs.items():
            results[key] = func(**kwargs)
        return results

    monkeypatch.setattr(trim_fastqs, "run_with_pool", run_jobs_sequentially)


def _write_fastq(path: Path, records: list[tuple[str, str, str]]) -> Path:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8") as out:
        for name, seq, qual in records:
            out.write(f"@{name}\n{seq}\n+\n{qual}\n")
    return path


def _read_fastq(path: Path) -> list[tuple[str, str, str]]:
    with gzip.open(path, "rt", encoding="utf-8") as infile:
        lines = infile.read().splitlines()
    records = []
    for idx in range(0, len(lines), 4):
        if idx + 3 >= len(lines):
            break
        records.append((lines[idx][1:], lines[idx + 1], lines[idx + 3]))
    return records


def _read_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as infile:
        return json.load(infile)


def _assert_single_record_has_prefix_without_adapter(
    path: Path,
    *,
    prefix: str,
    adapter: str,
) -> None:
    records = _read_fastq(path)
    assert len(records) == 1
    _name, seq, qual = records[0]
    assert seq.startswith(prefix)
    assert adapter not in seq
    assert len(qual) == len(seq)


def _run_trim_sample(
    tmp_path: Path,
    fastqs: tuple[Path, Path | None],
    cutsite_motifs: tuple[str, str],
    **kwargs,
) -> Path:
    defaults = dict(
        sname="sample",
        outdir=tmp_path,
        cutsite_motifs=cutsite_motifs,
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
    defaults.update(kwargs)
    trim_fastqs.trim_sample_with_fastp(fastqs=fastqs, **defaults)
    return tmp_path / "sample.R1.trimmed.fastq.gz"


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


def test_trim_sample_single_end_removes_overhang_and_adapter(tmp_path: Path) -> None:
    seq = "TGCAG" + "ACGTACGTAC" + TRUSEQ_R1
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", [("r1", seq, "I" * len(seq))])

    out1 = _run_trim_sample(
        tmp_path,
        fastqs=(fastq, None),
        cutsite_motifs=("TGCAG", ""),
    )

    _assert_single_record_has_prefix_without_adapter(
        out1,
        prefix="ACGTACGTAC",
        adapter=TRUSEQ_R1,
    )


def test_trim_sample_single_end_supports_plain_fastq_input(tmp_path: Path) -> None:
    seq = "TGCAG" + "ACGTACGTAC" + TRUSEQ_R1
    fastq = _write_fastq(tmp_path / "sample.fastq", [("r1", seq, "I" * len(seq))])

    out1 = _run_trim_sample(
        tmp_path,
        fastqs=(fastq, None),
        cutsite_motifs=("TGCAG", ""),
    )

    _assert_single_record_has_prefix_without_adapter(
        out1,
        prefix="ACGTACGTAC",
        adapter=TRUSEQ_R1,
    )


def test_trim_sample_single_end_disable_adapter_trimming_preserves_adapter_tail(tmp_path: Path) -> None:
    seq = "TGCAG" + "ACGTACGTAC" + TRUSEQ_R1
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", [("r1", seq, "I" * len(seq))])

    out1 = _run_trim_sample(
        tmp_path,
        fastqs=(fastq, None),
        cutsite_motifs=("TGCAG", ""),
        disable_adapter_trimming=True,
    )

    assert _read_fastq(out1) == [("r1", "ACGTACGTAC" + TRUSEQ_R1, "I" * (10 + len(TRUSEQ_R1)))]


def test_trim_sample_paired_end_removes_overhangs_and_adapters_from_both_reads(tmp_path: Path) -> None:
    seq1 = "TGCAG" + "ACGTACGTACGT" + TRUSEQ_R1
    seq2 = "CGATC" + "TGCATGCATGCA" + TRUSEQ_R2
    r1 = _write_fastq(tmp_path / "sample_R1.fastq.gz", [("r1", seq1, "I" * len(seq1))])
    r2 = _write_fastq(tmp_path / "sample_R2.fastq.gz", [("r1", seq2, "I" * len(seq2))])

    _run_trim_sample(
        tmp_path,
        fastqs=(r1, r2),
        cutsite_motifs=("TGCAG", "CGATC"),
    )

    _assert_single_record_has_prefix_without_adapter(
        tmp_path / "sample.R1.trimmed.fastq.gz",
        prefix="ACGTACGTACGT",
        adapter=TRUSEQ_R1,
    )
    _assert_single_record_has_prefix_without_adapter(
        tmp_path / "sample.R2.trimmed.fastq.gz",
        prefix="TGCATGCATGCA",
        adapter=TRUSEQ_R2,
    )


def test_trim_sample_internal_low_quality_read_is_dropped_by_default(tmp_path: Path) -> None:
    seq = "TGCAG" + "ACGTACGTACGTACGTACGT"
    qual = "I" * 5 + "I" * 6 + "!!!!" + "I" * 10
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", [("r1", seq, qual)])

    out1 = _run_trim_sample(
        tmp_path,
        fastqs=(fastq, None),
        cutsite_motifs=("TGCAG", ""),
        disable_adapter_trimming=True,
    )

    assert _read_fastq(out1) == []


def test_trim_sample_internal_low_quality_read_is_retained_when_quality_filtering_is_disabled(tmp_path: Path) -> None:
    seq = "TGCAG" + "ACGTACGTACGTACGTACGT"
    qual = "I" * 5 + "I" * 6 + "!!!!" + "I" * 10
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", [("r1", seq, qual)])

    out1 = _run_trim_sample(
        tmp_path,
        fastqs=(fastq, None),
        cutsite_motifs=("TGCAG", ""),
        disable_adapter_trimming=True,
        disable_quality_filtering=True,
    )

    assert _read_fastq(out1) == [("r1", "ACGTACGTACGTACGTACGT", "IIIIII!!!!IIIIIIIIII")]


def test_trim_sample_phred64_converts_output_qualities_to_phred33(tmp_path: Path) -> None:
    seq = "TGCAGACGTACGTACGT"
    qual = "h" * len(seq)
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", [("r1", seq, qual)])

    out1 = _run_trim_sample(
        tmp_path,
        fastqs=(fastq, None),
        cutsite_motifs=("TGCAG", ""),
        disable_adapter_trimming=True,
    )
    retained_without_phred64 = _read_fastq(out1)

    out2_dir = tmp_path / "phred64"
    out2_dir.mkdir()
    fastq2 = _write_fastq(out2_dir / "sample.fastq.gz", [("r1", seq, qual)])
    out2 = _run_trim_sample(
        out2_dir,
        fastqs=(fastq2, None),
        cutsite_motifs=("TGCAG", ""),
        disable_adapter_trimming=True,
        phred64=True,
    )
    retained_with_phred64 = _read_fastq(out2)

    assert retained_without_phred64 == [("r1", "ACGTACGTACGT", "h" * 12)]
    assert retained_with_phred64 == [("r1", "ACGTACGTACGT", "I" * 12)]


def test_trim_sample_min_trimmed_length_controls_read_retention(tmp_path: Path) -> None:
    seq = "TGCAGACGTAC"
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", [("r1", seq, "I" * len(seq))])

    retained = _run_trim_sample(
        tmp_path,
        fastqs=(fastq, None),
        cutsite_motifs=("TGCAG", ""),
        min_trimmed_length=6,
    )

    drop_dir = tmp_path / "drop"
    drop_dir.mkdir()
    fastq2 = _write_fastq(drop_dir / "sample.fastq.gz", [("r1", seq, "I" * len(seq))])
    dropped = _run_trim_sample(
        drop_dir,
        fastqs=(fastq2, None),
        cutsite_motifs=("TGCAG", ""),
        min_trimmed_length=7,
    )

    assert _read_fastq(retained) == [("r1", "ACGTAC", "IIIIII")]
    assert _read_fastq(dropped) == []


def test_trim_sample_many_ns_are_filtered_out(tmp_path: Path) -> None:
    seq = "TGCAGACGTNNNNNNNNNNACGT"
    fastq = _write_fastq(tmp_path / "sample.fastq.gz", [("r1", seq, "I" * len(seq))])

    out1 = _run_trim_sample(
        tmp_path,
        fastqs=(fastq, None),
        cutsite_motifs=("TGCAG", ""),
        disable_adapter_trimming=True,
        max_ns=5,
    )

    assert _read_fastq(out1) == []


def test_run_trimmer_respects_max_reads_and_writes_single_end_summary(
    tmp_path: Path,
    sequential_pool,
) -> None:
    fastq = _write_fastq(
        tmp_path / "sample.fastq.gz",
        [
            ("r1", "TGCAGACGTACGTACGT", "I" * 17),
            ("r2", "TGCAGTTGCAACGTTGCA", "I" * 18),
        ],
    )
    outdir = tmp_path / "out"

    trim_fastqs.run_trimmer(
        fastqs=[fastq],
        outdir=outdir,
        cutsite_motifs=("TGCAG", ""),
        max_reads=1,
        min_trimmed_length=1,
        max_unqualified_percent=15,
        min_quality=20,
        min_mean_window_quality=30,
        cut_window_size=5,
        phred64=False,
        max_reads_kmer=100,
        max_ns=5,
        disable_infer_cutsite_motifs=True,
        disable_adapter_trimming=True,
        disable_quality_filtering=False,
        cores=1,
        threads=1,
        delim_str=None,
        delim_idx=1,
        suffix=None,
        umi_tag_in_i5=False,
        force=False,
        log_level="ERROR",
        logged_command="ipyrad2 trim -d sample.fastq.gz -o out",
    )

    records = _read_fastq(outdir / "sample.R1.trimmed.fastq.gz")
    stats_json = _read_json(outdir / "sample.stats.json")
    summary_text = (outdir / "ipyrad_trim_stats_0.txt").read_text(encoding="utf-8")
    summary_json = _read_json(outdir / "ipyrad_trim_stats_0.json")

    assert len(records) == 1
    assert stats_json["summary"]["before_filtering"]["total_reads"] == 1
    assert summary_text.startswith("CMD: ipyrad2 trim -d sample.fastq.gz -o out\n\n")
    assert "sample" in summary_text
    assert "read2_mean_length_before" not in summary_text
    assert summary_json["command"] == "ipyrad2 trim -d sample.fastq.gz -o out"
    assert summary_json["sample_summary"][0]["sample"] == "sample"


def test_run_trimmer_writes_paired_end_summary_columns(
    tmp_path: Path,
    sequential_pool,
) -> None:
    seq1 = "TGCAG" + "ACGTACGTACGT" + TRUSEQ_R1
    seq2 = "CGATC" + "TGCATGCATGCA" + TRUSEQ_R2
    r1 = _write_fastq(tmp_path / "sample_R1.fastq.gz", [("r1", seq1, "I" * len(seq1))])
    r2 = _write_fastq(tmp_path / "sample_R2.fastq.gz", [("r1", seq2, "I" * len(seq2))])
    outdir = tmp_path / "out"

    trim_fastqs.run_trimmer(
        fastqs=[r1, r2],
        outdir=outdir,
        cutsite_motifs=("TGCAG", "CGATC"),
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

    summary_text = (outdir / "ipyrad_trim_stats_0.txt").read_text(encoding="utf-8")
    summary_json = _read_json(outdir / "ipyrad_trim_stats_0.json")

    assert (outdir / "sample.R1.trimmed.fastq.gz").exists()
    assert (outdir / "sample.R2.trimmed.fastq.gz").exists()
    assert "sample" in summary_text
    assert "read2_mean_length_before" in summary_text
    assert "command" not in summary_json


def test_run_trimmer_blocks_when_only_r2_output_exists_without_force(
    tmp_path: Path,
    sequential_pool,
) -> None:
    seq1 = "TGCAG" + "ACGTACGTACGT" + TRUSEQ_R1
    seq2 = "CGATC" + "TGCATGCATGCA" + TRUSEQ_R2
    r1 = _write_fastq(tmp_path / "sample_R1.fastq.gz", [("r1", seq1, "I" * len(seq1))])
    r2 = _write_fastq(tmp_path / "sample_R2.fastq.gz", [("r1", seq2, "I" * len(seq2))])
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "sample.R2.trimmed.fastq.gz").write_bytes(b"existing")

    with pytest.raises(
        trim_fastqs.IPyradError,
        match="Trim output artifact exists in outdir: .*sample\\.R2\\.trimmed\\.fastq\\.gz",
    ):
        trim_fastqs.run_trimmer(
            fastqs=[r1, r2],
            outdir=outdir,
            cutsite_motifs=("TGCAG", "CGATC"),
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


def test_run_trimmer_ignores_preexisting_stats_html_file(
    tmp_path: Path,
    sequential_pool,
) -> None:
    fastq = _write_fastq(
        tmp_path / "sample.fastq.gz",
        [("r1", "TGCAGACGTACGTACGT", "I" * 17)],
    )
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "sample.stats.html").write_text("existing", encoding="utf-8")

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
        disable_adapter_trimming=True,
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

    assert (outdir / "sample.stats.json").exists()
    assert not (outdir / "sample.stats.html").exists()


def test_run_trimmer_force_allows_rewrite_when_stats_artifact_exists(
    tmp_path: Path,
    sequential_pool,
) -> None:
    fastq = _write_fastq(
        tmp_path / "sample.fastq.gz",
        [("r1", "TGCAGACGTACGTACGT", "I" * 17)],
    )
    outdir = tmp_path / "out"
    outdir.mkdir()
    (outdir / "sample.stats.html").write_text("existing", encoding="utf-8")

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
        disable_adapter_trimming=True,
        disable_quality_filtering=False,
        cores=1,
        threads=1,
        delim_str=None,
        delim_idx=1,
        suffix=None,
        umi_tag_in_i5=False,
        force=True,
        log_level="ERROR",
    )

    assert _read_fastq(outdir / "sample.R1.trimmed.fastq.gz") == [("r1", "ACGTACGTACGT", "I" * 12)]
    assert (outdir / "sample.stats.json").exists()
    assert not (outdir / "sample.stats.html").exists()


def test_run_trimmer_supports_plain_paired_fastq_inputs(
    tmp_path: Path,
    sequential_pool,
) -> None:
    seq1 = "TGCAGACGTACGT"
    seq2 = "CGATCTGCATGCA"
    r1 = _write_fastq(tmp_path / "sample_R1.fastq", [("r1", seq1, "I" * len(seq1))])
    r2 = _write_fastq(tmp_path / "sample_R2.fastq", [("r1", seq2, "I" * len(seq2))])
    outdir = tmp_path / "out"

    trim_fastqs.run_trimmer(
        fastqs=[r1, r2],
        outdir=outdir,
        cutsite_motifs=("TGCAG", "CGATC"),
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
        disable_adapter_trimming=True,
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

    assert _read_fastq(outdir / "sample.R1.trimmed.fastq.gz") == [("r1", "ACGTACGT", "I" * 8)]
    assert _read_fastq(outdir / "sample.R2.trimmed.fastq.gz") == [("r1", "TGCATGCA", "I" * 8)]


def test_run_trimmer_skips_empty_single_end_samples_in_mixed_batch(
    tmp_path: Path,
    sequential_pool,
) -> None:
    empty_fastq = _write_fastq(tmp_path / "empty.fastq", [])
    full_fastq = _write_fastq(
        tmp_path / "full.fastq",
        [("r1", "TGCAGACGTACGTACGT", "I" * 17)],
    )
    outdir = tmp_path / "out"

    trim_fastqs.run_trimmer(
        fastqs=[empty_fastq, full_fastq],
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
        disable_adapter_trimming=True,
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

    summary_text = (outdir / "ipyrad_trim_stats_0.txt").read_text(encoding="utf-8")

    assert (outdir / "full.R1.trimmed.fastq.gz").exists()
    assert not (outdir / "empty.R1.trimmed.fastq.gz").exists()
    assert "full" in summary_text
    assert "empty" not in summary_text


def test_run_trimmer_skips_paired_sample_when_both_mates_are_empty(
    tmp_path: Path,
    sequential_pool,
) -> None:
    empty_r1 = _write_fastq(tmp_path / "empty_R1.fastq.gz", [])
    empty_r2 = _write_fastq(tmp_path / "empty_R2.fastq.gz", [])
    full_r1 = _write_fastq(
        tmp_path / "full_R1.fastq.gz",
        [("r1", "TGCAGACGTACGTACGT", "I" * 17)],
    )
    full_r2 = _write_fastq(
        tmp_path / "full_R2.fastq.gz",
        [("r1", "CGATCTGCATGCATGCA", "I" * 17)],
    )
    outdir = tmp_path / "out"

    trim_fastqs.run_trimmer(
        fastqs=[empty_r1, empty_r2, full_r1, full_r2],
        outdir=outdir,
        cutsite_motifs=("TGCAG", "CGATC"),
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
        disable_adapter_trimming=True,
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

    assert (outdir / "full.R1.trimmed.fastq.gz").exists()
    assert (outdir / "full.R2.trimmed.fastq.gz").exists()
    assert not (outdir / "empty.R1.trimmed.fastq.gz").exists()
    assert not (outdir / "empty.R2.trimmed.fastq.gz").exists()


def test_run_trimmer_skips_paired_sample_when_one_mate_is_empty(
    tmp_path: Path,
    sequential_pool,
) -> None:
    half_r1 = _write_fastq(
        tmp_path / "half_R1.fastq.gz",
        [("r1", "TGCAGACGTACGTACGT", "I" * 17)],
    )
    half_r2 = _write_fastq(tmp_path / "half_R2.fastq.gz", [])
    full_r1 = _write_fastq(
        tmp_path / "full_R1.fastq.gz",
        [("r1", "TGCAGACGTACGTACGT", "I" * 17)],
    )
    full_r2 = _write_fastq(
        tmp_path / "full_R2.fastq.gz",
        [("r1", "CGATCTGCATGCATGCA", "I" * 17)],
    )
    outdir = tmp_path / "out"

    trim_fastqs.run_trimmer(
        fastqs=[half_r1, half_r2, full_r1, full_r2],
        outdir=outdir,
        cutsite_motifs=("TGCAG", "CGATC"),
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
        disable_adapter_trimming=True,
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

    assert (outdir / "full.R1.trimmed.fastq.gz").exists()
    assert (outdir / "full.R2.trimmed.fastq.gz").exists()
    assert not (outdir / "half.R1.trimmed.fastq.gz").exists()
    assert not (outdir / "half.R2.trimmed.fastq.gz").exists()


def test_run_trimmer_raises_when_all_samples_are_empty(
    tmp_path: Path,
    sequential_pool,
) -> None:
    empty_fastq = _write_fastq(tmp_path / "empty.fastq.gz", [])

    with pytest.raises(
        trim_fastqs.IPyradError,
        match="No non-empty FASTQ samples remain after input validation",
    ):
        trim_fastqs.run_trimmer(
            fastqs=[empty_fastq],
            outdir=tmp_path / "out",
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
            disable_adapter_trimming=True,
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


def test_run_trimmer_uses_inferred_overhangs_to_trim_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    sequential_pool,
) -> None:
    seq1 = "XTGCAGACGTACGT"
    seq2 = "CGATCTGCATGCA"
    r1 = _write_fastq(tmp_path / "sample_R1.fastq.gz", [("r1", seq1, "I" * len(seq1))])
    r2 = _write_fastq(tmp_path / "sample_R2.fastq.gz", [("r1", seq2, "I" * len(seq2))])
    outdir = tmp_path / "out"
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
        calls.append((tuple(path.name for path in fastqs), max_reads, tuple(candidate_offsets or ())))
        if all(path.name.endswith("R2.fastq.gz") for path in fastqs):
            return _junction_set(("CGATC",), offset=0)
        return _junction_set(("TGCAG",), offset=1)

    monkeypatch.setattr(trim_fastqs, "get_overhangs_from_kmers", fake_get_overhangs_from_kmers)

    trim_fastqs.run_trimmer(
        fastqs=[r1, r2],
        outdir=outdir,
        cutsite_motifs=None,
        max_reads=None,
        min_trimmed_length=1,
        max_unqualified_percent=15,
        min_quality=20,
        min_mean_window_quality=30,
        cut_window_size=5,
        phred64=False,
        max_reads_kmer=123,
        max_ns=5,
        disable_infer_cutsite_motifs=False,
        disable_adapter_trimming=True,
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

    assert calls == [
        (("sample_R1.fastq.gz",), 123, (0, 1)),
        (("sample_R2.fastq.gz",), 123, (0, 1)),
    ]
    assert _read_fastq(outdir / "sample.R1.trimmed.fastq.gz") == [("r1", "ACGTACGT", "IIIIIIII")]
    assert _read_fastq(outdir / "sample.R2.trimmed.fastq.gz") == [("r1", "TGCATGCA", "IIIIIIII")]


def test_run_trimmer_uses_longest_inferred_multi_motif_trim_length(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    sequential_pool,
) -> None:
    seq1 = "ATCGATCGTACGT"
    r1 = _write_fastq(tmp_path / "sample.fastq.gz", [("r1", seq1, "I" * len(seq1))])
    outdir = tmp_path / "out"

    monkeypatch.setattr(
        trim_fastqs,
        "get_overhangs_from_kmers",
        lambda *args, **kwargs: _junction_set(("ATCGG", "ATCGAT"), counts=(40, 20)),
    )

    trim_fastqs.run_trimmer(
        fastqs=[r1],
        outdir=outdir,
        cutsite_motifs=None,
        max_reads=None,
        min_trimmed_length=1,
        max_unqualified_percent=15,
        min_quality=20,
        min_mean_window_quality=30,
        cut_window_size=5,
        phred64=False,
        max_reads_kmer=123,
        max_ns=5,
        disable_infer_cutsite_motifs=False,
        disable_adapter_trimming=True,
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

    assert _read_fastq(outdir / "sample.R1.trimmed.fastq.gz") == [("r1", "CGTACGT", "I" * 7)]


def test_run_trimmer_inferred_single_end_trim_matches_detected_motif_length_exactly(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    sequential_pool,
) -> None:
    insert = "GGACTTACGA"
    seq1 = "TGCAG" + insert
    r1 = _write_fastq(tmp_path / "sample.fastq.gz", [("r1", seq1, "I" * len(seq1))])
    outdir = tmp_path / "out"

    monkeypatch.setattr(
        trim_fastqs,
        "get_overhangs_from_kmers",
        lambda *args, **kwargs: _junction_set(("TGCAG",), offset=0),
    )

    trim_fastqs.run_trimmer(
        fastqs=[r1],
        outdir=outdir,
        cutsite_motifs=None,
        max_reads=None,
        min_trimmed_length=1,
        max_unqualified_percent=15,
        min_quality=20,
        min_mean_window_quality=30,
        cut_window_size=5,
        phred64=False,
        max_reads_kmer=123,
        max_ns=5,
        disable_infer_cutsite_motifs=False,
        disable_adapter_trimming=True,
        disable_quality_filtering=True,
        cores=1,
        threads=1,
        delim_str=None,
        delim_idx=1,
        suffix=None,
        umi_tag_in_i5=False,
        force=False,
        log_level="ERROR",
    )

    trimmed = _read_fastq(outdir / "sample.R1.trimmed.fastq.gz")
    assert trimmed == [("r1", insert, "I" * len(insert))]
    assert trimmed[0][1].startswith("GGA")


def test_run_trimmer_inferred_paired_end_trim_does_not_remove_extra_three_bases(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    sequential_pool,
) -> None:
    insert1 = "GGACTTACGA"
    insert2 = "GGATTCAGTA"
    seq1 = "XTGCAG" + insert1
    seq2 = "CGATC" + insert2
    r1 = _write_fastq(tmp_path / "sample_R1.fastq.gz", [("r1", seq1, "I" * len(seq1))])
    r2 = _write_fastq(tmp_path / "sample_R2.fastq.gz", [("r1", seq2, "I" * len(seq2))])
    outdir = tmp_path / "out"

    def fake_get_overhangs_from_kmers(
        fastqs,
        max_len,
        max_reads,
        workers,
        log_level,
        candidate_offsets=None,
        label=None,
    ):
        del max_len, max_reads, workers, log_level, candidate_offsets, label
        if all(path.name.endswith("R2.fastq.gz") for path in fastqs):
            return _junction_set(("CGATC",), offset=0)
        return _junction_set(("TGCAG",), offset=1)

    monkeypatch.setattr(trim_fastqs, "get_overhangs_from_kmers", fake_get_overhangs_from_kmers)

    trim_fastqs.run_trimmer(
        fastqs=[r1, r2],
        outdir=outdir,
        cutsite_motifs=None,
        max_reads=None,
        min_trimmed_length=1,
        max_unqualified_percent=15,
        min_quality=20,
        min_mean_window_quality=30,
        cut_window_size=5,
        phred64=False,
        max_reads_kmer=123,
        max_ns=5,
        disable_infer_cutsite_motifs=False,
        disable_adapter_trimming=True,
        disable_quality_filtering=True,
        cores=1,
        threads=1,
        delim_str=None,
        delim_idx=1,
        suffix=None,
        umi_tag_in_i5=False,
        force=False,
        log_level="ERROR",
    )

    trimmed_r1 = _read_fastq(outdir / "sample.R1.trimmed.fastq.gz")
    trimmed_r2 = _read_fastq(outdir / "sample.R2.trimmed.fastq.gz")
    assert trimmed_r1 == [("r1", insert1, "I" * len(insert1))]
    assert trimmed_r2 == [("r1", insert2, "I" * len(insert2))]
    assert trimmed_r1[0][1].startswith("GGA")
    assert trimmed_r2[0][1].startswith("GGA")
