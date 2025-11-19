#!/usr/env/bin python

"""Window extracter command line
"""

import argparse
import ipyrad2 as ip
import sys
from loguru import logger
from pathlib import Path
from .make_wide import make_wide
from .cli_wex import _setup_wex_subparser
from .cli_lex import _setup_lex_subparser
from ..analysis.window_extracter import run_window_extracter
from ..analysis.locus_extracter import run_locus_extracter

VERSION = str(ip.__version__)

HEADER = f"""
-------------------------------------------------------------
ipyrad [v.{VERSION}]
Interactive assembly and analysis of RAD-seq data
-------------------------------------------------------------\
"""

DESCRIPTION = "ipyrad analysis command line. Select a positional subcommand:"


EPILOG = r"""
Examples
--------
$ ipyrad analysis wex -d seqs.hdf5 --print-scaffold-table
$ ipyrad analysis wex -d seqs.hdf5 -o OUT/ -n TEST -m 10
"""


def _setup_analysis_subparser(subparser: argparse._SubParsersAction, header: str = None) -> None:
    """Add `ipyrad analysis` subcommand parser.

    """
    analysis_parser = subparser.add_parser(
        "analysis",
        #description=f"{HEADER}\n{DESCRIPTION}",
        description=header,
        help="Utilities for downstream analysis",
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter, max_help_position=60, width=140))

    analysis_subparser = analysis_parser.add_subparsers(
        dest="tool",
        required=True)
    _setup_wex_subparser(analysis_subparser, f"{HEADER}\nipyrad analysis wex: window extracter to filter and write concatenated alignments")
    _setup_lex_subparser(analysis_subparser, f"{HEADER}\nipyrad analysis lex: locus extracter to select and write random alignments")


def run_analysis_tool(args):

    # WEX: --------------------------------------------------------
    if args.tool == "wex":
        logger.info("-------------------------------------------------------")
        logger.info("----- ipyrad wex: extract alignments from windows -----")
        logger.info("-------------------------------------------------------")
        logger.info(f"CMD: ipyrad {' '.join(sys.argv[1:])}")
        run_window_extracter(
            data=args.data,
            name=args.name,
            outdir=args.out,
            out_format=args.out_format,
            windows=args.windows,
            min_sample_coverage=args.min_sample_coverage,
            max_sample_missing=args.max_sample_missing,
            imap=args.imap,
            minmap=args.minmap,
            exclude=args.exclude,
            print_scaffold_table=args.print_scaffold_table,
            stdout=args.stdout,
            force=args.force,
        )
        sys.exit(0)

    if args.tool == "lex":
        logger.info("-------------------------------------------------------")
        logger.info("--- ipyrad lex: extract random alignments from hdf5 ---")
        logger.info("-------------------------------------------------------")
        logger.info(f"CMD: ipyrad {' '.join(sys.argv[1:])}")
        run_locus_extracter(
            data=args.data,
            name=args.name,
            outdir=args.out,
            nloci=args.nloci,
            length=args.length,
            windows=args.windows,
            min_sample_coverage=args.min_sample_coverage,
            max_sample_missing=args.max_sample_missing,
            imap=args.imap,
            minmap=args.minmap,
            exclude=args.exclude,
            print_scaffold_table=args.print_scaffold_table,
            stdout=args.stdout,
            force=args.force,
        )

        sys.exit(0)

