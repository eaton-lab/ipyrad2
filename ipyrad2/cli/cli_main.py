#!/usr/bin/env python

"""
"""

import sys
import argparse
from .make_wide import make_wide
from .cli_demux import _setup_demux_subparser
from .cli_trim import _setup_trim_subparser
from .cli_denovo import _setup_denovo_subparser
from .cli_map import _setup_map_subparser
from .cli_assemble import _setup_assemble_subparser
from .cli_wex import _setup_wex_subparser
from ..demuxer import run_demuxer
from ..trimmer import run_trimmer
from ..denovo import run_denovo
from ..mapper import run_mapper
from ..assembler import run_assembler
from ..analysis.window_extracter import run_window_extracter
from ..utils.logger import set_log_level
from ..utils.exceptions import IPyradError
from loguru import logger
import ipyrad2 as ip

VERSION = str(ip.__version__)

HEADER = f"""
-------------------------------------------------------------
ipyrad2 [v.{VERSION}]
Interactive assembly and analysis of RAD-seq data
-------------------------------------------------------------\
"""

DESCRIPTION = "ipyrad2 command line tool. Select a positional subcommand:"

EPILOG = """\
Note
----
Each subcommand has its own help screen, e.g.,:
$ ipyrad2 demux -h

Examples
--------
# demux: demultiplexing data to samples by index or barcode
$ ipyrad2 demux -d RAW/*.fastq.gz -b BARCODES.csv -m 1 -c 10 -o ./demux

# trim: trim reads for quality, adapters, and restriction overhangs
$ ipyrad2 trim -d DATA/*.fastq.gz -o TRIMMED/ -q 20 -n 5 -c 10

# map: map reads to a reference genome and filter and sort BAMs
$ ipyrad2 map -d DATA/*.fastq.gz -o BAMs -c 10

# assemble: delimit rad loci, call variants, filter, and write outputs
$ ipyrad2 assemble -d BAMS/*.bam -o OUT -p TEST -m 4 -q 20 -c 10
"""


def setup_parsers() -> argparse.ArgumentParser:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = argparse.ArgumentParser(
        prog="ipyrad2",
        description=f"{HEADER}\n{DESCRIPTION}",
        epilog=EPILOG,
        formatter_class=make_wide(argparse.RawDescriptionHelpFormatter),
        add_help=False,
    )
    parser.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)
    parser.add_argument("-v", "--version", action='version', version=f"ipyrad {VERSION}")
    subparser = parser.add_subparsers(help="sub-commands", dest="subcommand")

    # add subcommands: these messages are subcommand headers
    _setup_demux_subparser(subparser, f"{HEADER}\nipyrad2 demux: demultiplex pooled reads to sample files by index/barcode")
    _setup_trim_subparser(subparser, f"{HEADER}\nipyrad2 trim: trim for quality, adapters, and restriction overhangs")
    _setup_denovo_subparser(subparser, f"{HEADER}\nipyrad2 denovo: construct a reference locus library")
    _setup_map_subparser(subparser, f"{HEADER}\nipyrad2 map: reference map, filter, and sort reads to bam files")
    _setup_assemble_subparser(subparser, f"{HEADER}\nipyrad2 assemble: delimit loci, call variants, and write outputs")
    _setup_wex_subparser(subparser, f"{HEADER}\nipyrad2 wex: window extracter to filter and write concatenated alignments")
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

    # LOGGING: -----------------------------------------------------
    if hasattr(args, "log_level"):
        set_log_level(args.log_level, args.log_file)

    if args.subcommand not in ["demux", "trim", "denovo", "map", "assemble", "wex"]:
        # NO SUBCOMMAND: print help
        parser.print_help()
        sys.exit(0)
    else:
        run_subcommand(args)



def run_subcommand(args, _exit=True):
    # DEMUX: -------------------------------------------------------
    if args.subcommand == "demux":
        logger.info("---------------------------------------------------------")
        logger.info("----- ipyrad demux: demultiplexing reads to samples -----")
        logger.info("---------------------------------------------------------")
        logger.info(f"CMD: ipyrad {' '.join(sys.argv[1:])}")
        run_demuxer(
            fastqs=args.fastqs,
            outdir=args.out,
            barcodes=args.barcodes,
            re1=args.restriction_overhang_1,
            re2=args.restriction_overhang_2,
            max_mismatch=args.max_mismatch,
            chunksize=args.chunksize,
            i7=args.i7,
            disable_infer_re_overhangs=args.disable_infer_re_overhangs,
            merge_technical_replicates=args.merge_technical_replicates,
            cores=args.cores,
            max_reads=args.max_reads,
            log_level=args.log_level,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # TRIM: -------------------------------------------------------
    if args.subcommand == "trim":
        logger.info("----------------------------------------------------------")
        logger.info("----- ipyrad trim: quality, adapter, and RE trimming -----")
        logger.info("----------------------------------------------------------")
        logger.info(f"CMD: ipyrad {' '.join(sys.argv[1:])}")
        run_trimmer(
            fastqs=args.fastqs,
            outdir=args.out,
            max_unqualified_percent=args.max_unqualified_percent,
            min_quality=args.min_quality,
            min_mean_window_quality=args.min_mean_window_quality,
            cut_window_size=args.cut_window_size,
            max_ns=args.max_ns,
            min_trimmed_length=args.min_trimmed_length,
            max_reads=args.max_reads,
            restriction_overhangs=args.restriction_overhangs,
            max_reads_kmer=args.max_reads_kmer,
            phred_qscore_offset=args.phred_qscore_offset,
            disable_infer_re_overhangs=args.disable_infer_re_overhangs,
            disable_adapter_trimming=args.disable_adapter_trimming,
            disable_quality_filtering=args.disable_quality_filtering,
            cores=args.cores,
            threads=args.threads,
            delim_str=args.delim_str,
            delim_idx=args.delim_idx,
            suffix=args.suffix,
            umi_tag_in_i5=args.umi_tag_in_i5,
            force=args.force,
            log_level=args.log_level,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # DENOVO: --------------------------------------------------------
    if args.subcommand == "denovo":
        logger.info("------------------------------------------------------------")
        logger.info("----- ipyrad denovo: construct locus reference library -----")
        logger.info("------------------------------------------------------------")
        logger.info(f"CMD: ipyrad {' '.join(sys.argv[1:])}")
        run_denovo(
            fastqs=args.fastqs,
            outdir=args.out,
            similarity_threshold_within=args.similarity_threshold_within,
            similarity_threshold_across=args.similarity_threshold_across,
            min_dereplication_size=args.min_dereplication_size,
            min_length=args.min_length,
            min_merge_overlap=args.min_merge_overlap,
            max_merge_diffs=args.max_merge_diffs,
            cores=args.cores,
            threads=args.threads,
            force=args.force,
            strand_both=args.strand_both,
            delim_str=args.delim_str,
            delim_idx=args.delim_idx,
            log_level=args.log_level,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # MAP: --------------------------------------------------------
    if args.subcommand == "map":
        logger.info("----------------------------------------------------------")
        logger.info("----- ipyrad map: map, filter and sort reads to bams -----")
        logger.info("----------------------------------------------------------")
        logger.info(f"CMD: ipyrad {' '.join(sys.argv[1:])}")
        run_mapper(
            fastqs=args.fastqs,
            reference=args.reference,
            outdir=args.out,
            imap=args.imap,
            min_map_q=args.min_map_q,
            max_edit_dist=args.max_edit_dist,
            max_soft_clip=args.max_soft_clip,
            mark_dups_by_coords=args.mark_dups_by_coords,
            mark_dups_by_umis=args.mark_dups_by_umis,
            cores=args.cores,
            threads=args.threads,
            force=args.force,
            delim_str=args.delim_str,
            delim_idx=args.delim_idx,
            log_level=args.log_level,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # ASSEMBLE: ---------------------------------------------------
    if args.subcommand == "assemble":
        logger.info("-----------------------------------------------------------")
        logger.info("----- ipyrad assemble: delimit loci and call variants -----")
        logger.info("-----------------------------------------------------------")
        logger.info(f"CMD: ipyrad {' '.join(sys.argv[1:])}")
        run_assembler(
            rad_bams=args.rad_bams,
            wgs_bams=args.wgs_bams,
            reference=args.reference,
            outdir=args.out,
            name=args.name,
            loci_bed=args.loci_bed,
            min_map_q=args.min_map_q,
            min_site_q=args.min_site_q,
            min_geno_q=args.min_geno_q,
            min_base_q=args.min_base_q,
            min_sample_depth=args.min_sample_depth,
            min_locus_sample_coverage=args.min_locus_sample_coverage,
            min_locus_trim_sample_coverage=args.min_locus_trim_sample_coverage,
            min_locus_length=args.min_locus_length,
            min_locus_merge_distance=args.min_locus_merge_distance,
            max_locus_hetero_frequency=args.max_locus_hetero_frequency,
            max_locus_variant_frequency=args.max_locus_variant_frequency,
            populations=args.populations,
            masks=args.masks,
            exclude_reference=args.exclude_reference,
            cores=args.cores,
            threads=args.threads,
            force=args.force,
            log_level=args.log_level,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # WEX: --------------------------------------------------------
    if args.subcommand == "wex":
        logger.info("-------------------------------------------------------")
        logger.info("----- ipyrad wex: extract alignments from windows -----")
        logger.info("-------------------------------------------------------")
        logger.info(f"CMD: ipyrad {' '.join(sys.argv[1:])}")
        run_window_extracter(
            data=args.data,
            name=args.name,
            outdir=args.out,
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
        if _exit: sys.exit(0)  # noqa: E701


if __name__ == "__main__":

    main()
