#!/usr/bin/env python

"""VCF to HDF5 converter command line."""

import argparse
from pathlib import Path

from .common import RAW_HELP_FORMATTER


EPILOG = r"""
Examples
--------
$ ipyrad2 analysis vcf-to-hdf5 -d variants.vcf.gz -o OUT/
$ ipyrad2 analysis vcf-to-hdf5 -d variants.vcf.gz -o OUT/ -n snps -b 20000
"""


def _setup_vcf_to_hdf5_subparser(
    subparsers: argparse._SubParsersAction,
    header: str = None,
) -> None:
    """Add `ipyrad2 analysis vcf-to-hdf5` subcommand parser."""
    tool = subparsers.add_parser(
        "vcf-to-hdf5",
        description=header,
        help="Convert a VCF file into an SNP-capable HDF5 database.",
        epilog=EPILOG,
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )

    core = tool.add_argument_group("Core inputs")
    core.add_argument(
        "-d", "--data", metavar="Path", type=Path, required=True,
        help="Path to a VCF or bgzipped VCF file.",
    )
    core.add_argument(
        "-n", "--name", metavar="str", type=str, default="snps",
        help="Prefix name for the converted HDF5 database. [default=snps]",
    )
    core.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="output-vcf-to-hdf5",
        help="Directory to write the converted HDF5 database. [default=output-vcf-to-hdf5]",
    )

    conversion = tool.add_argument_group("Conversion")
    conversion.add_argument(
        "-b", "--ld-block-size", metavar="int", type=int, default=10000,
        help="Linkage-block size used for generic VCF inputs. [default=10000]",
    )
    conversion.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite an existing converted HDF5 database with the same name.",
    )

    logging = tool.add_argument_group("Logging")
    logging.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Log level (TRACE, DEBUG, INFO, WARNING, ERROR) [default=INFO]",
    )
    logging.add_argument(
        "-h", "--help", action="help",
        help="Show this help message and exit.",
    )
