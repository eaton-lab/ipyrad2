from pathlib import Path

import pytest

from ipyrad2.utils.exceptions import IPyradError
from ipyrad2.utils.names import get_name_to_fastq_dict
from ipyrad2.utils.names import get_paths_list_from_fastq_str


def test_fastq_globs_expand_env_vars(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fastq = tmp_path / "sample.fastq.gz"
    fastq.touch()
    monkeypatch.setenv("IPYRAD2_FASTQ_DIR", str(tmp_path))

    result = get_paths_list_from_fastq_str("$IPYRAD2_FASTQ_DIR/*.fastq.gz")

    assert result == [fastq]


def test_common_paired_end_names_are_grouped_by_mate_token(tmp_path: Path) -> None:
    fastq2 = tmp_path / "sample_R2_001.fastq.gz"
    fastq1 = tmp_path / "sample_R1_001.fastq.gz"
    fastq2.touch()
    fastq1.touch()

    result = get_name_to_fastq_dict([fastq2, fastq1], None, None)

    assert result == {"sample": (fastq1, fastq2)}


def test_shared_trailing_suffix_after_mate_token_is_still_paired(tmp_path: Path) -> None:
    fastq2 = tmp_path / "iTru7_301_01b_S13_L008_R2_001-sub.fastq.gz"
    fastq1 = tmp_path / "iTru7_301_01b_S13_L008_R1_001-sub.fastq.gz"
    fastq2.touch()
    fastq1.touch()

    result = get_name_to_fastq_dict([fastq2, fastq1], None, None)

    assert result == {"iTru7_301_01b_S13_L008": (fastq1, fastq2)}


def test_split_pair_suffixes_are_grouped_by_secondary_mate_fallback(tmp_path: Path) -> None:
    fastq1a = tmp_path / "iTru7_301_01b_S13_L008_R1_001-sub.fastq.gz"
    fastq2a = tmp_path / "iTru7_301_01b_S13_L008_R2_001-sub.fastq.gz"
    fastq1b = tmp_path / "iTru7_301_01b_S13_L008_R1_002-sub.fastq.gz"
    fastq2b = tmp_path / "iTru7_301_01b_S13_L008_R2_002-sub.fastq.gz"
    for path in (fastq1a, fastq2a, fastq1b, fastq2b):
        path.touch()

    result = get_name_to_fastq_dict([fastq2b, fastq1a, fastq2a, fastq1b], None, None)

    assert result == {
        "iTru7_301_01b_S13_L008_001-sub": (fastq1a, fastq2a),
        "iTru7_301_01b_S13_L008_002-sub": (fastq1b, fastq2b),
    }


def test_secondary_mate_fallback_supports_dotted_suffix_patterns(tmp_path: Path) -> None:
    fastq1a = tmp_path / "sample.R1.partA.fastq.gz"
    fastq2a = tmp_path / "sample.R2.partA.fastq.gz"
    fastq1b = tmp_path / "sample.R1.partB.fastq.gz"
    fastq2b = tmp_path / "sample.R2.partB.fastq.gz"
    for path in (fastq1a, fastq2a, fastq1b, fastq2b):
        path.touch()

    result = get_name_to_fastq_dict([fastq2b, fastq1a, fastq2a, fastq1b], None, None)

    assert result == {
        "sample.partA": (fastq1a, fastq2a),
        "sample.partB": (fastq1b, fastq2b),
    }


def test_mismatched_trailing_suffix_after_mate_token_raises(tmp_path: Path) -> None:
    fastq1 = tmp_path / "sample_L001_R1_001-sub.fastq.gz"
    fastq2 = tmp_path / "sample_L001_R2_001-trim.fastq.gz"
    fastq1.touch()
    fastq2.touch()

    with pytest.raises(IPyradError, match="same trailing suffix"):
        get_name_to_fastq_dict([fastq1, fastq2], None, None)


def test_ambiguous_single_end_fastqs_are_not_collapsed_into_a_pair(tmp_path: Path) -> None:
    fastq1 = tmp_path / "sampleA_rep.fastq.gz"
    fastq2 = tmp_path / "sampleB_rep.fastq.gz"
    fastq1.touch()
    fastq2.touch()

    result = get_name_to_fastq_dict([fastq1, fastq2], None, None)

    assert result == {
        "sampleA_rep": (fastq1, None),
        "sampleB_rep": (fastq2, None),
    }


def test_incomplete_split_pairs_fall_back_to_single_end_names(tmp_path: Path) -> None:
    fastq1a = tmp_path / "sample_R1_001-sub.fastq.gz"
    fastq2a = tmp_path / "sample_R2_001-sub.fastq.gz"
    fastq1b = tmp_path / "sample_R1_002-sub.fastq.gz"
    for path in (fastq1a, fastq2a, fastq1b):
        path.touch()

    with pytest.raises(IPyradError, match="do not form complete consistent R1/R2 pairs"):
        get_name_to_fastq_dict([fastq1a, fastq2a, fastq1b], None, None)


def test_mixed_paired_end_and_unrecognized_layout_raises(tmp_path: Path) -> None:
    fastq1 = tmp_path / "sample_R1.fastq.gz"
    fastq2 = tmp_path / "sample_R2.fastq.gz"
    weird = tmp_path / "sample_extra.fastq.gz"
    for path in (fastq1, fastq2, weird):
        path.touch()

    with pytest.raises(IPyradError, match="unrecognized alongside paired-end-looking files"):
        get_name_to_fastq_dict([fastq1, fastq2, weird], None, None)


def test_delim_parsing_raises_on_ambiguous_non_paired_groups(tmp_path: Path) -> None:
    fastq1 = tmp_path / "sample_laneA.fastq.gz"
    fastq2 = tmp_path / "sample_laneB.fastq.gz"
    fastq1.touch()
    fastq2.touch()

    with pytest.raises(IPyradError, match="ambiguous non-paired sample names"):
        get_name_to_fastq_dict([fastq1, fastq2], "_lane", 1)


def test_negative_delim_index_raises(tmp_path: Path) -> None:
    fastq = tmp_path / "sample.fastq.gz"
    fastq.touch()

    with pytest.raises(IPyradError, match="delim_index must be >= 1"):
        get_name_to_fastq_dict([fastq], "_", -1)


def test_empty_inputs_fail_fast() -> None:
    with pytest.raises(IPyradError, match="No fastq data were provided"):
        get_name_to_fastq_dict([], None, None)


def test_overlapping_input_patterns_raise_instead_of_deduping(tmp_path: Path) -> None:
    fastq = tmp_path / "sample_R1.fastq.gz"
    fastq.touch()

    with pytest.raises(IPyradError, match="matched more than once"):
        get_paths_list_from_fastq_str(
            [tmp_path / "sample_R1.fastq.gz", tmp_path / "sample_R*.fastq.gz"]
        )
