#!/usr/bin/env python

"""
"""

import argparse
import glob
import os
import sys
import traceback
from ipyrad2.cli.make_wide import make_wide
from ..demuxer import run_demuxer
from ..trimmer import run_trimmer
from ..denovo import run_denovo
from ..mapper import run_mapper
from ..assembler import run_assembler
from ..utils.logger import set_log_level
from ..utils.exceptions import IPyradError
from ..utils.params import read_params, new_params
from argparse import Namespace
from loguru import logger
from pathlib import Path
import ipyrad2 as ip

VERSION = str(ip.__version__)

HEADER = f"""
-------------------------------------------------------------
ipyrad [v.{VERSION}]
Interactive assembly and analysis of RAD-seq data
-------------------------------------------------------------\
"""

DESCRIPTION = "ipyrad2 classic command line tool."

EPILOG = """\

Examples
--------
ipyrad2-classic -n data                       ## create new file called params-data.txt 
ipyrad2-classic -p params-data.txt -s 123     ## run only steps 1-3 of assembly.
ipyrad2-classic -p params-data.txt -s 3 -f    ## run step 3, overwrite existing data.
"""


def setup_parsers() -> argparse.ArgumentParser:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = argparse.ArgumentParser(
        prog="ipyrad",
        description=f"{HEADER}\n{DESCRIPTION}",
        epilog=EPILOG,
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter),
        add_help=False,
    )
    parser.add_argument("-n", action='store', dest='new', help="create new file 'params-{new}.txt' in current directory")
    parser.add_argument("-p", action='store', dest='params', help="path to params file for Assembly")
    parser.add_argument("-s", action='store', dest='steps', help="Set of assembly steps to run, e.g., -s 123")
    parser.add_argument("-c", action='store', dest='cores', type=int, default=8, help="number of CPU cores to use (Default=8)")
    parser.add_argument("-t", action='store', dest='threads', type=int, default=2, help="tune threading of multi-threaded binaries (Default=2)")
    parser.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Log level (DEBUG, INFO, WARN, ERROR) [default=INFO]",
    )
    parser.add_argument(
        "-L", "--log-file", metavar="Path", type=Path,
        help="Log file. Logging to stdout is also appended to this file. [default=None]."
    )
    parser.add_argument("-f", "--force", action="store_true", help="force overwrite of existing data")
    parser.add_argument("-d", "--debug", action="store_true", help="Print debug information")
    parser.add_argument("-v", "--version", action='version', version=f"ipyrad {VERSION}")
    parser.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)

    return parser


def main():
    try:
        command_line()
    except KeyboardInterrupt:
        logger.error("interrupted by user. Shutting down.")
        sys.exit(1)
    # expected error, only report message no traceback
    except IPyradError as exc:
        logger.error(f"Error: {exc}")
        logger.error("see error message above. Shutting down.")
        sys.exit(1)
    # raise with traceback
    except Exception as exc:
        logger.exception("Unexpected error: see traceback below.")
        raise exc


def command_line():
    parser = setup_parsers()
    args = parser.parse_args()

    if args.new:
        _flagnew(args.new)
        sys.exit(0)

    elif args.params is not None:
        params = read_params(args.params)
        if not os.path.exists(params.main.project_dir):
            os.mkdir(params.main.project_dir)
    else:
        sys.exit("Classic mode requires either -n or -p")

    # LOGGING: -----------------------------------------------------
    if hasattr(args, "log_level"):
        set_log_level(args.log_level, args.log_file)

    # DEMUX: -------------------------------------------------------
    if "1" in args.steps:
        s1_args = params.demux
        s1_args.subcommand = "demux"
        # Black magic to merge s1 specific args with the few useful ones
        # we read from the cli, e.g. cores, force, and logging info
        s1_args = Namespace(**{**vars(s1_args), **vars(args)})
        # Update demux params from the params file
        s1_args.fastqs = params.main.raw_fastq_path
        s1_args.barcodes = params.main.barcodes_path

        s1_args.out = Path(params.main.project_dir) / (params.main.name + "_fastqs")
        ip.cli.cli_main.run_subcommand(s1_args, _exit=False)

    # TRIM: -------------------------------------------------------
    if "2" in args.steps:
        s2_args = params.trim
        s2_args.subcommand = "trim"
        s2_args = Namespace(**{**vars(s2_args), **vars(args)})
        s2_args.fastqs = Path(params.main.project_dir) / (params.main.name + "_fastqs/*.gz")
        s2_args.out = Path(params.main.project_dir) / (params.main.name + "_edits")
        ip.cli.cli_main.run_subcommand(s2_args, _exit=False)

    # DENOVO: --------------------------------------------------------
    if "3" in args.steps:
        if not os.path.exists(params.main.reference_sequence):
            s3_args = params.denovo
            s3_args.subcommand = "denovo"
            s3_args = Namespace(**{**vars(s3_args), **vars(args)})
            s3_args.fastqs = Path(params.main.project_dir) / (params.main.name + "_edits/*.gz")
            s3_args.out = Path(params.main.project_dir) / (params.main.name + "_reference")
            ip.cli.cli_main.run_subcommand(s3_args, _exit=False)
        else:
            logger.info("Reference sequence exists, skipping denovo reference assembly.")

    # MAP: --------------------------------------------------------
    if "4" in args.steps:
        s4_args = params.map
        s4_args.subcommand = "map"
        s4_args = Namespace(**{**vars(s4_args), **vars(args)})
        s4_args.fastqs = Path(params.main.project_dir) / (params.main.name + "_edits/*.gz")
        # If user passed in reference then use this else use the default ref from step 3
        if os.path.exists(params.main.reference_sequence):
            s4_args.reference = Path(params.main.reference_sequence)
        else:
            s4_args.reference = Path(params.main.project_dir) / (params.main.name + "_reference/denovo_reference.fa")
        s4_args.out = Path(params.main.project_dir) / (params.main.name + "_mapped")
        ip.cli.cli_main.run_subcommand(s4_args, _exit=False)

    # ASSEMBLE: ---------------------------------------------------
    if "5" in args.steps:
        s5_args = params.assemble
        s5_args.subcommand = "assemble"
        s5_args = Namespace(**{**vars(s5_args), **vars(args)})
        s5_args.name = params.main.name
        bams = glob.glob(str(Path(params.main.project_dir) / (params.main.name + "_mapped/*.bam")))
        s5_args.rad_bams = [Path(x) for x in bams]
        # TODO: Handle wgs_bams in classic mode
        s5_args.wgs_bams = None
        # Toggle whether to use the passed in or denovo constructed reference sequence
        if os.path.exists(params.main.reference_sequence):
            s5_args.reference = Path(params.main.reference_sequence)
        else:
            s5_args.reference = Path(params.main.project_dir) / (params.main.name + "_reference/denovo_reference.fa")

        s5_args.out = Path(params.main.project_dir) / (params.main.name + "_outfiles")
        ip.cli.cli_main.run_subcommand(s5_args, _exit=False)

    sys.exit(0)


def _flagnew(name):

        new_params(name)
        # print log to screen
        print(f"\n  New file 'params-{name}.txt' created\n")



if __name__ == "__main__":

    main()
