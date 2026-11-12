import argparse
import json
from pathlib import Path

import pytest

from ipyrad2.cli.cli_demux import validate_demux_args
from ipyrad2.cli.cli_main import command_line, setup_parsers


def _get_demux_parser() -> argparse.ArgumentParser:
    parser = setup_parsers()
    subparsers = next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    return subparsers.choices["demux"]


def test_demux_schema_matches_parser(capsys) -> None:
    command_line(["schema", "demux"])
    output = capsys.readouterr()
    schema = json.loads(output.out)
    assert not output.err
    assert schema["schema_version"] == 1
    assert schema["command"] == ["ipyrad2", "demux"]

    parser = _get_demux_parser()
    options = schema["options"]
    assert set(options) == {action.dest for action in parser._actions}
    for action in parser._actions:
        option = options[action.dest]
        assert option["flags"] == action.option_strings
        assert option["default"] == action.default
        assert option["required"] == action.required
        if action.choices is not None:
            assert option["choices"] == list(action.choices)

    assert options["tmpdir"]["type"] == "path"
    assert options["pigz"]["type"] == "boolean"
    assert options["max_reads"]["type"] == "integer"
    assert options["max_reads"]["minimum"] == 1
    assert "first N input read pairs" in options["max_reads"]["description"]
    assert options["max_mismatch"]["maximum"] == 2
    assert options["delim_idx"]["excluded_values"] == [0]
    assert any(
        rule.get("requires") == "pigz" and rule["option"] == "tmpdir"
        for rule in schema["dependencies"]
    )
    assert any(
        rule.get("conflicts_with") == "i7"
        and rule["option"] == "allow_leading_barcode_deletion"
        for rule in schema["dependencies"]
    )


def test_schema_help_and_invalid_tool(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        command_line(["--help"])
    assert exc.value.code == 0
    help_text = capsys.readouterr().out
    assert help_text.index("treeslider") < help_text.index("    schema ")

    with pytest.raises(SystemExit) as exc:
        command_line(["schema", "map"])
    assert exc.value.code == 2
    assert "invalid choice" in capsys.readouterr().err


def test_demux_help_groups_and_examples_are_updated() -> None:
    help_text = _get_demux_parser().format_help()

    expected_sections = [
        "Core inputs:",
        "Demultiplexing mode:",
        "Sample naming and pairing:",
        "Cutsite motifs:",
        "Performance and sampling:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    expected_order = [
        "-d, --fastqs",
        "-b, --barcodes",
        "-o, --out",
        "-f, --force",
        "--i7",
        "-m, --max_mismatch",
        "-M, --merge-technical-replicates",
        "--barcode-boundary-slack",
        "--allow-leading-barcode-deletion",
        "-dx, --delim-str",
        "-di, --delim-idx",
        "-e1, --cutsite-1",
        "-e2, --cutsite-2",
        "-E, --disable-infer-cutsite-motifs",
        "-c, --cores",
        "-k, --chunksize",
        "-x, --max_reads",
        "--max-reads-kmer",
        "--pigz",
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

    assert "ipyrad2 demux: demultiplex pooled reads to sample files by barcode or index" in help_text
    assert "$ ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.csv --log-level DEBUG" in help_text
    assert "Use commas for multiple motifs" in help_text
    assert "positive from left, negative from right" in help_text
    assert "compared to inferred motifs unless -E is set" in help_text
    assert "overrides inference" not in help_text
    assert "ipyrad demux" not in help_text
    assert "--logger" not in help_text
    assert "*.fastqs.gz" not in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_demux_parser_defaults_are_unchanged() -> None:
    args = setup_parsers().parse_args(["demux", "-d", "a.fastq.gz", "-b", "bars.tsv"])

    assert args.subcommand == "demux"
    assert args.fastqs == [Path("a.fastq.gz")]
    assert args.out == Path("DEMUX")
    assert args.force is False
    assert args.barcodes == Path("bars.tsv")
    assert args.cutsite_1 is None
    assert args.cutsite_2 is None
    assert args.max_mismatch == 0
    assert args.max_reads is None
    assert args.max_reads_kmer == 100_000
    assert args.pigz is False
    assert args.chunksize == 10_000_000
    assert args.delim_str is None
    assert args.delim_idx == 1
    assert args.disable_infer_cutsite_motifs is False
    assert args.merge_technical_replicates is False
    assert args.barcode_boundary_slack == 1
    assert args.allow_leading_barcode_deletion is False
    assert args.i7 is False
    assert args.cores == 4
    assert args.log_level == "INFO"
    assert not hasattr(args, "log_file")


def test_demux_parser_accepts_negative_delim_idx() -> None:
    args = setup_parsers().parse_args(
        ["demux", "-d", "a.fastq.gz", "-b", "bars.tsv", "-dx", "_R", "-di", "-1"]
    )

    assert args.delim_str == "_R"
    assert args.delim_idx == -1


def test_demux_parser_accepts_comma_separated_manual_overhangs() -> None:
    args = setup_parsers().parse_args(
        [
            "demux",
            "-d",
            "a.fastq.gz",
            "-b",
            "bars.tsv",
            "-e1",
            "ATCGG,ATCGAT",
        ]
    )

    assert args.cutsite_1 == "ATCGG,ATCGAT"


def test_demux_parser_accepts_pigz_flag() -> None:
    args = setup_parsers().parse_args(
        ["demux", "-d", "a.fastq.gz", "-b", "bars.tsv", "--pigz"]
    )

    assert args.pigz is True


def test_demux_parser_accepts_barcode_boundary_slack_zero() -> None:
    args = setup_parsers().parse_args(
        [
            "demux",
            "-d",
            "a.fastq.gz",
            "-b",
            "bars.tsv",
            "--barcode-boundary-slack",
            "0",
        ]
    )

    assert args.barcode_boundary_slack == 0


def test_demux_parser_accepts_leading_barcode_deletion_flag() -> None:
    args = setup_parsers().parse_args(
        [
            "demux",
            "-d",
            "a.fastq.gz",
            "-b",
            "bars.tsv",
            "--allow-leading-barcode-deletion",
        ]
    )

    assert args.allow_leading_barcode_deletion is True


def test_demux_parser_rejects_invalid_barcode_boundary_slack() -> None:
    with pytest.raises(SystemExit):
        setup_parsers().parse_args(
            [
                "demux",
                "-d",
                "a.fastq.gz",
                "-b",
                "bars.tsv",
                "--barcode-boundary-slack",
                "2",
            ]
        )


@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (["demux", "-d", "a.fastq.gz", "-b", "bars.tsv", "-m", "-1"], "--max-mismatch must be between 0 and 2"),
        (["demux", "-d", "a.fastq.gz", "-b", "bars.tsv", "-m", "3"], "--max-mismatch must be between 0 and 2"),
        (["demux", "-d", "a.fastq.gz", "-b", "bars.tsv", "--max-reads-kmer", "0"], "--max-reads-kmer must be >= 1"),
        (["demux", "-d", "a.fastq.gz", "-b", "bars.tsv", "-di", "0"], "--delim-idx cannot be 0"),
        (
            [
                "demux",
                "-d",
                "a.fastq.gz",
                "-b",
                "bars.tsv",
                "--i7",
                "--allow-leading-barcode-deletion",
            ],
            "--allow-leading-barcode-deletion applies only to inline barcode demux",
        ),
    ],
)
def test_demux_parser_validation(argv: list[str], message: str, capsys: pytest.CaptureFixture[str]) -> None:
    parser = setup_parsers()
    subparsers = next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    args = parser.parse_args(argv)

    with pytest.raises(SystemExit):
        validate_demux_args(args, subparsers.choices["demux"])
    assert message in capsys.readouterr().err


def test_demux_tmpdir_parsing_and_pigz_requirement():
    parser = setup_parsers()
    args = parser.parse_args(["demux", "-d", "a.gz", "-b", "bars.tsv", "--pigz", "--tmpdir", "/tmp"])
    assert str(args.tmpdir) == "/tmp"
    assert args.pigz
    args.pigz = False
    with pytest.raises(SystemExit):
        validate_demux_args(args, parser)
