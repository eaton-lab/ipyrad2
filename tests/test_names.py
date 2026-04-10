from pathlib import Path

import pytest

import ipyrad2.utils.names as names_module
from ipyrad2.utils.exceptions import IPyradError
from ipyrad2.utils.names import get_name_to_fastq_dict
from ipyrad2.utils.names import get_paths_list_from_fastq_str


class _StubLogger:
    def __init__(self) -> None:
        self.warnings = []

    def info(self, *args, **kwargs) -> None:
        pass

    def warning(self, message, *args, **kwargs) -> None:
        if args:
            message = message.format(*args)
        self.warnings.append(message)


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


def test_mismatched_trailing_suffix_after_mate_token_is_quiet_single_end(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    logger = _StubLogger()
    monkeypatch.setattr(names_module, "logger", logger)
    fastq1 = tmp_path / "sample_L001_R1_001-sub.fastq.gz"
    fastq2 = tmp_path / "sample_L001_R2_001-trim.fastq.gz"
    fastq1.touch()
    fastq2.touch()

    result = get_name_to_fastq_dict([fastq1, fastq2], None, None)

    assert result == {
        "sample_L001_R1_001-sub": (fastq1, None),
        "sample_L001_R2_001-trim": (fastq2, None),
    }
    assert logger.warnings == []


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


def test_terminal_r2_in_single_end_name_is_quiet_single_end(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    logger = _StubLogger()
    monkeypatch.setattr(names_module, "logger", logger)
    fastq = tmp_path / "Larkspur2.fastq.gz"
    fastq.touch()

    result = names_module.get_name_to_fastq_dict([fastq], None, None)

    assert result == {"Larkspur2": (fastq, None)}
    assert logger.warnings == []


def test_r1_only_single_end_outputs_are_quiet_single_end(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    logger = _StubLogger()
    monkeypatch.setattr(names_module, "logger", logger)
    fastqs = []
    for idx in range(1, 5):
        fastq = tmp_path / f"barbeyi-CO1-{idx:02d}_R1.fastq.gz"
        fastq.touch()
        fastqs.append(fastq)

    result = names_module.get_name_to_fastq_dict(fastqs, None, None)

    assert result == {
        f"barbeyi-CO1-{idx:02d}_R1": (fastqs[idx - 1], None)
        for idx in range(1, 5)
    }
    assert logger.warnings == []


def test_unmatched_r1_r2_tokens_are_quiet_single_end(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    logger = _StubLogger()
    monkeypatch.setattr(names_module, "logger", logger)
    fastq1 = tmp_path / "alpha_R1.fastq.gz"
    fastq2 = tmp_path / "beta_R2.fastq.gz"
    fastq1.touch()
    fastq2.touch()

    result = names_module.get_name_to_fastq_dict([fastq1, fastq2], None, None)

    assert result == {
        "alpha_R1": (fastq1, None),
        "beta_R2": (fastq2, None),
    }
    assert logger.warnings == []


def test_low_complete_pair_ratio_warns_and_falls_back(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    logger = _StubLogger()
    monkeypatch.setattr(names_module, "logger", logger)
    fastqs = [
        tmp_path / "sample_R1.fastq.gz",
        tmp_path / "sample_R2.fastq.gz",
        tmp_path / "alpha.fastq.gz",
        tmp_path / "beta.fastq.gz",
        tmp_path / "gamma.fastq.gz",
    ]
    for fastq in fastqs:
        fastq.touch()

    result = names_module.get_name_to_fastq_dict(fastqs, None, None)

    assert result == {
        "alpha": (fastqs[2], None),
        "beta": (fastqs[3], None),
        "gamma": (fastqs[4], None),
        "sample_R1": (fastqs[0], None),
        "sample_R2": (fastqs[1], None),
    }
    assert any("complete auto-detected pairs cover 2/5" in warning for warning in logger.warnings)


def test_majority_complete_split_pairs_warn_and_fall_back_to_single_end(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    logger = _StubLogger()
    monkeypatch.setattr(names_module, "logger", logger)
    fastq1a = tmp_path / "sample_R1_001-sub.fastq.gz"
    fastq2a = tmp_path / "sample_R2_001-sub.fastq.gz"
    fastq1b = tmp_path / "sample_R1_002-sub.fastq.gz"
    for path in (fastq1a, fastq2a, fastq1b):
        path.touch()

    result = get_name_to_fastq_dict([fastq1a, fastq2a, fastq1b], None, None)

    assert result == {
        "sample_R1_001-sub": (fastq1a, None),
        "sample_R1_002-sub": (fastq1b, None),
        "sample_R2_001-sub": (fastq2a, None),
    }
    assert any("complete auto-detected pairs cover 2/3" in warning for warning in logger.warnings)


def test_mixed_paired_end_and_unrecognized_layout_warns_and_falls_back(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    logger = _StubLogger()
    monkeypatch.setattr(names_module, "logger", logger)
    fastq1 = tmp_path / "sample_R1.fastq.gz"
    fastq2 = tmp_path / "sample_R2.fastq.gz"
    weird = tmp_path / "sample_extra.fastq.gz"
    for path in (fastq1, fastq2, weird):
        path.touch()

    result = get_name_to_fastq_dict([fastq1, fastq2, weird], None, None)

    assert result == {
        "sample_R1": (fastq1, None),
        "sample_R2": (fastq2, None),
        "sample_extra": (weird, None),
    }
    assert any("unrecognized alongside paired-end-looking files" in warning for warning in logger.warnings)
    assert any("complete auto-detected pairs cover 2/3" in warning for warning in logger.warnings)


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
