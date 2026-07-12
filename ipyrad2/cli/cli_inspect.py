#!/usr/bin/env python

"""Inspect command-line parser."""

import argparse
from argparse import ArgumentParser, Namespace
from pathlib import Path

from .common import RAW_HELP_FORMATTER


EPILOG = """\
Examples
--------
$ ipyrad2 inspect ./myassembly_outfiles
"""


def validate_inspect_args(args: Namespace, parser: ArgumentParser) -> None:
    """Validate inspect CLI arguments after parsing."""
    if not args.assembly_dir.exists():
        parser.error(f"assembly output directory does not exist: {args.assembly_dir}")
    if not args.assembly_dir.is_dir():
        parser.error(f"assembly output path is not a directory: {args.assembly_dir}")


def _setup_inspect_subparser(
    subparsers: argparse._SubParsersAction,
    header: str = None,
) -> None:
    """Add `ipyrad2 inspect` subcommand parser."""
    tool = subparsers.add_parser(
        "inspect",
        description=header,
        help="Launch the interactive assembly browser.",
        epilog=EPILOG,
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )
    core = tool.add_argument_group("Core inputs")

    core.add_argument(
        "assembly_dir",
        metavar="Path",
        type=Path,
        help="Directory containing ipyrad2 assembly output files.",
    )
    tool.add_argument(
        "-h", "--help",
        action="help",
        help="show this help message and exit",
    )
