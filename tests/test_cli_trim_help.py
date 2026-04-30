import argparse
from pathlib import Path

import pytest

from ipyrad2.cli.cli_main import setup_parsers
from ipyrad2.cli.cli_trim import validate_trim_args


def _get_trim_parser() -> argparse.ArgumentParser:
    parser = setup_parsers()
    subparsers = next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    return subparsers.choices["trim"]


def test_trim_help_groups_examples_and_phred64_are_updated() -> None:
    help_text = _get_trim_parser().format_help()

    expected_sections = [
        "Core inputs:",
        "Filtering and trimming:",
        "Cutsite motifs and adapters:",
        "Performance and compatibility:",
        "Sample naming and UMI:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    expected_order = [
        "-d, --fastqs",
        "-o, --out",
        "-f, --force",
        "-q, --min-quality",
        "-u, --max-unqualified-percent",
        "-M, --min-mean-window-quality",
        "-W, --cut-window-size",
        "-n, --max-ns",
        "-e, --min-trimmed-length",
        "-Q, --disable-quality-filtering",
        "-e1, --cutsite-1",
        "-e2, --cutsite-2",
        "-k, --max-reads-kmer",
        "-E, --disable-infer-cutsite-motifs",
        "-A, --disable-adapter-trimming",
        "-x, --max-reads",
        "-c, --cores",
        "-t, --threads",
        "--phred64",
        "-dx, --delim-str",
        "-di, --delim-idx",
        "-s, --suffix",
        "-U, --umi-tag-in-i5",
        "-l, --log-level",
        "-h, --help",
    ]
    start = help_text.index("Core inputs:")
    indices = []
    for item in expected_order:
        idx = help_text.index(item, start)
        indices.append(idx)
        start = idx + 1
    assert indices == sorted(indices)

    assert "ipyrad2 trim: trim reads for quality, adapters, and cutsite motifs" in help_text
    assert "$ ipyrad2 trim -d DATA/*.gz -o OUT --phred64 -U" in help_text
    assert "--phred64" in help_text
    assert "Use commas for multiple motifs" in help_text
    assert "positive from left, negative from right" in help_text
    assert "ipyrad trim" not in help_text
    assert "--phred-qscore-offset" not in help_text
    assert "-nx" not in help_text
    assert "-ni" not in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_trim_parser_defaults_are_updated() -> None:
    args = setup_parsers().parse_args(["trim", "-d", "a.fastq.gz"])

    assert args.subcommand == "trim"
    assert args.fastqs == [Path("a.fastq.gz")]
    assert args.out == Path("TRIMMED")
    assert args.force is False
    assert args.max_unqualified_percent == 15
    assert args.min_quality == 20
    assert args.min_mean_window_quality == 30
    assert args.cut_window_size == 5
    assert args.max_ns == 5
    assert args.min_trimmed_length == 35
    assert args.max_reads is None
    assert args.cutsite_1 is None
    assert args.cutsite_2 is None
    assert args.max_reads_kmer == 500_000
    assert args.disable_infer_cutsite_motifs is False
    assert args.disable_adapter_trimming is False
    assert args.disable_quality_filtering is False
    assert args.cores == 6
    assert args.threads == 3
    assert args.phred64 is False
    assert args.delim_str is None
    assert args.delim_idx == 1
    assert args.suffix is None
    assert args.umi_tag_in_i5 is False
    assert args.log_level == "INFO"
    assert not hasattr(args, "log_file")
    assert not hasattr(args, "phred_qscore_offset")


def test_trim_parser_accepts_negative_delim_idx() -> None:
    args = setup_parsers().parse_args(
        ["trim", "-d", "a.fastq.gz", "-dx", "_R", "-di", "-1"]
    )

    assert args.delim_idx == -1


def test_trim_parser_accepts_comma_separated_manual_overhangs() -> None:
    args = setup_parsers().parse_args(
        ["trim", "-d", "a.fastq.gz", "-e1", "ATCGG,ATCGAT", "-e2", "CGATC"]
    )

    assert args.cutsite_1 == "ATCGG,ATCGAT"
    assert args.cutsite_2 == "CGATC"


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (["-d", "a.fastq.gz", "-c", "0"], "--cores must be >= 1"),
        (["-d", "a.fastq.gz", "-t", "0"], "--threads must be >= 1"),
        (["-d", "a.fastq.gz", "-c", "2", "-t", "3"], "--threads cannot exceed --cores"),
        (["-d", "a.fastq.gz", "-u", "101"], "--max-unqualified-percent must be between 0 and 100"),
        (["-d", "a.fastq.gz", "-M", "0"], "--min-mean-window-quality must be between 1 and 36"),
        (["-d", "a.fastq.gz", "-W", "1001"], "--cut-window-size must be between 1 and 1000"),
        (["-d", "a.fastq.gz", "-e", "0"], "--min-trimmed-length must be >= 1"),
        (["-d", "a.fastq.gz", "-k", "0"], "--max-reads-kmer must be >= 1"),
    ],
)
def test_validate_trim_args_rejects_invalid_values(
    argv: list[str],
    message: str,
    capsys: pytest.CaptureFixture[str],
) -> None:
    parser = _get_trim_parser()
    args = parser.parse_args(argv)

    with pytest.raises(SystemExit):
        validate_trim_args(args, parser)

    assert message in capsys.readouterr().err
