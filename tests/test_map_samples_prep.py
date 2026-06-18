from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from ipyrad2.mapper.map_samples_prep import apply_imap_to_samples
from ipyrad2.mapper.map_samples_prep import materialize_sample_plan
from ipyrad2.mapper.map_samples_prep import prepare_map_samples
from ipyrad2.mapper.map_samples_prep import unmate_paired_samples
from ipyrad2.utils.exceptions import IPyradError


FASTQ_RECORD = b"@r1\nACGT\n+\n!!!!\n"


def _write_fastq(path: Path, data: bytes = FASTQ_RECORD) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with gzip.open(path, "wb") as out:
            out.write(data)
    else:
        path.write_bytes(data)


def test_prepare_map_samples_merges_plain_and_gz_replicates(tmp_path: Path) -> None:
    fastqs = [
        tmp_path / "sampleA_R1.fastq",
        tmp_path / "sampleA_R2.fastq",
        tmp_path / "sampleB_R1.fastq.gz",
        tmp_path / "sampleB_R2.fastq.gz",
    ]
    _write_fastq(fastqs[0], b"@a1\nAAAA\n+\n!!!!\n")
    _write_fastq(fastqs[1], b"@a1\nTTTT\n+\n!!!!\n")
    _write_fastq(fastqs[2], b"@b1\nCCCC\n+\n!!!!\n")
    _write_fastq(fastqs[3], b"@b1\nGGGG\n+\n!!!!\n")

    imap = tmp_path / "imap.tsv"
    imap.write_text("sampleA merged\nsampleB merged\n", encoding="utf-8")

    fastq_dict, is_paired = prepare_map_samples(
        fastqs=fastqs,
        delim_str=None,
        delim_idx=1,
        imap=imap,
        tmpdir=tmp_path / "tmpdir",
        unmate=False,
    )

    assert is_paired is True
    assert list(fastq_dict) == ["merged"]
    plan = fastq_dict["merged"]
    assert plan.output_name == "merged"
    assert plan.source_names == ("sampleA", "sampleB")
    assert plan.is_paired_input is True
    assert plan.source_fastqs == (
        (fastqs[0], fastqs[1]),
        (fastqs[2], fastqs[3]),
    )


def test_prepare_map_samples_rejects_mixed_se_pe_merge(tmp_path: Path) -> None:
    sample_a_r1 = tmp_path / "sampleA_R1.fastq.gz"
    sample_a_r2 = tmp_path / "sampleA_R2.fastq.gz"
    sample_b = tmp_path / "sampleB.fastq.gz"
    for path in (sample_a_r1, sample_a_r2, sample_b):
        _write_fastq(path)

    imap = tmp_path / "imap.tsv"
    imap.write_text("sampleA merged\nsampleB merged\n", encoding="utf-8")

    with pytest.raises(IPyradError, match="mixed SE and PE"):
        apply_imap_to_samples(
            imap=imap,
            tmpdir=tmp_path / "tmpdir",
            fastq_dict={
                "sampleA": (sample_a_r1, sample_a_r2),
                "sampleB": (sample_b, None),
            },
        )


def test_apply_imap_to_samples_supports_glob_merge_rules(tmp_path: Path) -> None:
    sample_a_r1 = tmp_path / "sampleA_R1.fastq.gz"
    sample_a_r2 = tmp_path / "sampleA_R2.fastq.gz"
    sample_b_r1 = tmp_path / "sampleB_R1.fastq.gz"
    sample_b_r2 = tmp_path / "sampleB_R2.fastq.gz"
    for path, data in (
        (sample_a_r1, b"@a1\nAAAA\n+\n!!!!\n"),
        (sample_a_r2, b"@a1\nTTTT\n+\n!!!!\n"),
        (sample_b_r1, b"@b1\nCCCC\n+\n!!!!\n"),
        (sample_b_r2, b"@b1\nGGGG\n+\n!!!!\n"),
    ):
        _write_fastq(path, data)

    imap = tmp_path / "imap.tsv"
    imap.write_text("sample* merged\n", encoding="utf-8")

    fastq_dict = apply_imap_to_samples(
        imap=imap,
        tmpdir=tmp_path / "tmpdir",
        fastq_dict={
            "sampleA": (sample_a_r1, sample_a_r2),
            "sampleB": (sample_b_r1, sample_b_r2),
        },
    )

    assert list(fastq_dict) == ["merged"]
    merged_r1, merged_r2 = fastq_dict["merged"]
    with gzip.open(merged_r1, "rb") as infile:
        assert infile.read() == b"@a1\nAAAA\n+\n!!!!\n@b1\nCCCC\n+\n!!!!\n"
    with gzip.open(merged_r2, "rb") as infile:
        assert infile.read() == b"@a1\nTTTT\n+\n!!!!\n@b1\nGGGG\n+\n!!!!\n"


def test_apply_imap_to_samples_one_column_glob_keeps_identity_names(tmp_path: Path) -> None:
    sample_a = tmp_path / "sampleA.fastq.gz"
    sample_b = tmp_path / "sampleB.fastq.gz"
    other = tmp_path / "other.fastq.gz"
    for path in (sample_a, sample_b, other):
        _write_fastq(path)

    imap = tmp_path / "imap.tsv"
    imap.write_text("sample*\n", encoding="utf-8")

    fastq_dict = apply_imap_to_samples(
        imap=imap,
        tmpdir=tmp_path / "tmpdir",
        fastq_dict={
            "sampleA": (sample_a, None),
            "sampleB": (sample_b, None),
            "other": (other, None),
        },
    )

    assert list(fastq_dict) == ["sampleA", "sampleB"]
    assert fastq_dict["sampleA"] == (sample_a, None)
    assert fastq_dict["sampleB"] == (sample_b, None)


def test_unmate_paired_samples_concatenates_r1_then_r2(tmp_path: Path) -> None:
    sample_r1 = tmp_path / "sample_R1.fastq.gz"
    sample_r2 = tmp_path / "sample_R2.fastq.gz"
    _write_fastq(sample_r1, b"@r1\nAAAA\n+\n!!!!\n")
    _write_fastq(sample_r2, b"@r2\nTTTT\n+\n!!!!\n")

    result = unmate_paired_samples(
        {"sample": (sample_r1, sample_r2)},
        tmp_path / "tmpdir",
    )

    merged_r1, merged_r2 = result["sample"]
    assert merged_r2 is None
    with gzip.open(merged_r1, "rb") as infile:
        assert infile.read() == b"@r1\nAAAA\n+\n!!!!\n@r2\nTTTT\n+\n!!!!\n"


def test_unmate_paired_samples_reuses_gzip_members_without_recompressing(tmp_path: Path) -> None:
    sample_r1 = tmp_path / "sample_R1.fastq.gz"
    sample_r2 = tmp_path / "sample_R2.fastq.gz"
    _write_fastq(sample_r1, b"@r1\nAAAA\n+\n!!!!\n")
    _write_fastq(sample_r2, b"@r2\nTTTT\n+\n!!!!\n")
    raw_r1 = sample_r1.read_bytes()
    raw_r2 = sample_r2.read_bytes()

    result = unmate_paired_samples(
        {"sample": (sample_r1, sample_r2)},
        tmp_path / "tmpdir",
    )

    merged_r1, merged_r2 = result["sample"]
    assert merged_r2 is None
    assert merged_r1.read_bytes() == raw_r1 + raw_r2


def test_prepare_map_samples_unmates_after_imap_merge(tmp_path: Path) -> None:
    fastqs = [
        tmp_path / "sampleA_R1.fastq.gz",
        tmp_path / "sampleA_R2.fastq.gz",
        tmp_path / "sampleB_R1.fastq.gz",
        tmp_path / "sampleB_R2.fastq.gz",
    ]
    _write_fastq(fastqs[0], b"@a1\nAAAA\n+\n!!!!\n")
    _write_fastq(fastqs[1], b"@a1\nTTTT\n+\n!!!!\n")
    _write_fastq(fastqs[2], b"@b1\nCCCC\n+\n!!!!\n")
    _write_fastq(fastqs[3], b"@b1\nGGGG\n+\n!!!!\n")
    imap = tmp_path / "imap.tsv"
    imap.write_text("sampleA merged\nsampleB merged\n", encoding="utf-8")

    fastq_dict, is_paired = prepare_map_samples(
        fastqs=fastqs,
        delim_str=None,
        delim_idx=1,
        imap=imap,
        tmpdir=tmp_path / "tmpdir",
        unmate=True,
    )

    assert is_paired is True
    plan = fastq_dict["merged"]
    materialized, temp_paths = materialize_sample_plan(
        "merged",
        plan,
        tmp_path / "tmpdir",
        unmate=True,
    )
    merged_fastq, merged_r2 = materialized
    assert merged_r2 is None
    assert temp_paths == [merged_fastq]
    with gzip.open(merged_fastq, "rb") as infile:
        assert infile.read() == (
            b"@a1\nAAAA\n+\n!!!!\n"
            b"@b1\nCCCC\n+\n!!!!\n"
            b"@a1\nTTTT\n+\n!!!!\n"
            b"@b1\nGGGG\n+\n!!!!\n"
        )


def test_prepare_map_samples_rejects_unmate_on_single_end_inputs(tmp_path: Path) -> None:
    fastq = tmp_path / "sample.fastq.gz"
    _write_fastq(fastq)

    with pytest.raises(IPyradError, match="--unmate can only be used with paired-end FASTQ inputs"):
        prepare_map_samples(
            fastqs=[fastq],
            delim_str=None,
            delim_idx=1,
            imap=None,
            tmpdir=tmp_path / "tmpdir",
            unmate=True,
        )


def test_prepare_map_samples_normalizes_trimmed_single_end_name(tmp_path: Path) -> None:
    fastq = tmp_path / "sample.trimmed.fastq.gz"
    _write_fastq(fastq)

    fastq_dict, is_paired = prepare_map_samples(
        fastqs=[fastq],
        delim_str=None,
        delim_idx=1,
        imap=None,
        tmpdir=tmp_path / "tmpdir",
        unmate=False,
    )

    assert is_paired is False
    plan = fastq_dict["sample"]
    assert plan.output_name == "sample"
    assert plan.source_names == ("sample",)
    assert plan.source_fastqs == ((fastq, None),)
    assert plan.is_paired_input is False


def test_prepare_map_samples_matches_imap_against_normalized_trimmed_name(tmp_path: Path) -> None:
    fastq = tmp_path / "sample.trimmed.fastq.gz"
    _write_fastq(fastq)
    imap = tmp_path / "imap.tsv"
    imap.write_text("sample\n", encoding="utf-8")

    fastq_dict, is_paired = prepare_map_samples(
        fastqs=[fastq],
        delim_str=None,
        delim_idx=1,
        imap=imap,
        tmpdir=tmp_path / "tmpdir",
        unmate=False,
    )

    assert is_paired is False
    plan = fastq_dict["sample"]
    assert plan.output_name == "sample"
    assert plan.source_names == ("sample",)
    assert plan.source_fastqs == ((fastq, None),)
    assert plan.is_paired_input is False


def test_prepare_map_samples_rejects_trimmed_name_collision(tmp_path: Path) -> None:
    fastq_a = tmp_path / "sample.fastq.gz"
    fastq_b = tmp_path / "sample.trimmed.fastq.gz"
    _write_fastq(fastq_a)
    _write_fastq(fastq_b)

    with pytest.raises(IPyradError, match="collide after internal workflow-suffix normalization"):
        prepare_map_samples(
            fastqs=[fastq_a, fastq_b],
            delim_str=None,
            delim_idx=1,
            imap=None,
            tmpdir=tmp_path / "tmpdir",
            unmate=False,
        )
