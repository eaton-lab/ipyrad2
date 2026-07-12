import argparse
from pathlib import Path

import pytest

from ipyrad2.cli.cli_inspect import validate_inspect_args
from ipyrad2.cli.cli_main import setup_parsers


def _get_inspect_parser() -> argparse.ArgumentParser:
    parser = setup_parsers()
    subparsers = next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    return subparsers.choices["inspect"]


def test_inspect_help_describes_assembly_directory() -> None:
    help_text = _get_inspect_parser().format_help()

    assert "ipyrad2 inspect: launch the interactive assembly browser" in help_text
    assert "Directory containing ipyrad2 assembly output files." in help_text
    assert "$ ipyrad2 inspect OUT/ip2-pe_outfiles" in help_text


def test_inspect_parser_accepts_assembly_directory() -> None:
    args = setup_parsers().parse_args(["inspect", "outfiles"])

    assert args.subcommand == "inspect"
    assert args.assembly_dir == Path("outfiles")


def test_validate_inspect_args_rejects_missing_directory(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    parser = _get_inspect_parser()
    missing = tmp_path / "missing"
    args = parser.parse_args([str(missing)])

    with pytest.raises(SystemExit):
        validate_inspect_args(args, parser)

    assert "assembly output directory does not exist" in capsys.readouterr().err
