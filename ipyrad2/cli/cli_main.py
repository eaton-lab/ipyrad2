#!/usr/bin/env python

"""Top-level CLI parser and subcommand dispatcher."""

import sys
import argparse
from .common import RAW_HELP_FORMATTER
from .command_log import format_logged_command
from .cli_demux import _setup_demux_subparser, validate_demux_args
from .cli_trim import _setup_trim_subparser, validate_trim_args
from .cli_denovo import _setup_denovo_subparser, validate_denovo_args
from .cli_map import _setup_map_subparser
from .cli_assemble import _setup_assemble_subparser
from .cli_analysis import (
    ANALYSIS_TOOL_NAMES,
    RESERVED_TOOL_NAMES,
    _setup_analysis_tool_subparsers,
    run_analysis_tool,
)

import importlib
# from ..demuxer import run_demuxer
# from ..trimmer import run_trimmer
# from ..denovo import run_denovo
# from ..mapper import run_mapper
# from ..assembler import run_assembler
from ..utils.logger import set_log_level
from ..utils.exceptions import IPyradError
from loguru import logger
from ipyrad2 import __version__ as VERSION

HEADER = f"""
-------------------------------------------------------------
ipyrad2 [v.{VERSION}]
Interactive assembly and analysis of RAD-seq data
-------------------------------------------------------------\
"""

DESCRIPTION = "ipyrad2 command line tool. Select a positional subcommand:"

TOP_LEVEL_HELP = f"""{HEADER}
{DESCRIPTION}

assembly subcommands
    demux                                    Demultiplex pooled data to samples by index or barcode.
    trim                                     Trim reads for quality, adapters, and cutsite motifs using 'fastp'.
    denovo                                   Optionally construct a reference library by de novo clustering reads.
    map                                      Map reads to a reference with 'bwa-mem2' and 'samtools' to write coordinate-sorted BAMs.
    assemble                                 Delimit loci, call variants, and write assembled outputs, stats, and database (HDF5)

data export/conversion subcommands
    wex                                      Extract loci from HDF5 file, filter, and write as concatenated matrix to various formats
    lex                                      Extract loci from HDF5 file, filter, and write as multi-locus data to various formats
    snpex                                    Extract SNPs from HDF5 file, filter, optionally impute, and write to various formats
    vcf2hdf5                                 Convert an external VCF to HDF5 for use in analysis tools below.

analysis subcommands
    pca                                      Infer population structure from pca, tsne, or umap on filtered SNPs
    dapc                                     Infer population genetic clustering by discriminant analysis of principal components
    snmf                                     Infer population genetic clustering by non-negative matrix factorization
    admixture                                Infer population genetic clustering with external ADMIXTURE
    popgen                                   Infer population genetic statistics for one or more populations
    bpp                                      Infer species tree; species delim; or MSC+ model from multi-locus data
    baba                                     Infer admixture metrics from ABBA/BABA and related SNP patterns
    treeslider                               Infer gene trees for each qualified locus or refmapped genomic window of loci

options:
  -v, --version                              show program's version number and exit

Note
----
Each subcommand has its own help screen, e.g.,:
$ ipyrad2 demux -h

Assembly pipeline
-----------------
$ ipyrad2 demux    -d RAW/*.fastq.gz     -o DATA/    -b BARCODES.csv -m 1 -c 10
$ ipyrad2 trim     -d DATA/*.fastq.gz    -o TRIMMED/ -q 20 -n 5 -c 10
$ ipyrad2 map      -d TRIMMED/*.fastq.gz -o MAPPED   -r REF.fa -c 10
$ ipyrad2 assemble -d MAPPED/*.bam       -o OUT      -r REF.fa -m 4 -qm 20 -c 10

Export/Conversion examples
--------------------------
$ ipyrad2 wex -d OUT/HDF5 -m 10
$ ipyrad2 lex -d OUT/HDF5 -N 1000 -L 100 -m 10

Analysis examples
-----------------
$ ipyrad2 popgen -d OUT/HDF5 -i IMAP -g MINMAP
$ ipyrad2 pca -d OUT/HDF5 -i IMAP -g MINMAP -I sample --plot
"""


class TopLevelParser(argparse.ArgumentParser):
    """Root parser with grouped top-level help text."""

    def format_help(self) -> str:
        return TOP_LEVEL_HELP


def setup_parsers() -> argparse.ArgumentParser:
    """Setup and return an ArgumentParser w/ subcommands."""
    parser = TopLevelParser(
        prog="ipyrad2",
        description=f"{HEADER}\n{DESCRIPTION}",
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )
    parser.add_argument('-h', '--help', action='help', help=argparse.SUPPRESS)
    parser.add_argument("-v", "--version", action='version', version=f"ipyrad2 {VERSION}")
    subparser = parser.add_subparsers(
        help="sub-commands",
        dest="subcommand",
        parser_class=argparse.ArgumentParser,
    )

    # add subcommands: these messages are subcommand headers
    _setup_demux_subparser(subparser, f"{HEADER}\nipyrad2 demux: demultiplex pooled reads to sample files by barcode or index")
    _setup_trim_subparser(subparser, f"{HEADER}\nipyrad2 trim: trim reads for quality, adapters, and cutsite motifs")
    _setup_denovo_subparser(subparser, f"{HEADER}\nipyrad2 denovo: construct a reference locus library")
    _setup_map_subparser(subparser, f"{HEADER}\nipyrad2 map: map reads and write coordinate-sorted BAM files")
    _setup_assemble_subparser(subparser, f"{HEADER}\nipyrad2 assemble: delimit loci, call variants, and write outputs")
    _setup_analysis_tool_subparsers(subparser, HEADER)
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

    if args.subcommand == "demux":
        subparsers = next(
            action for action in parser._actions
            if isinstance(action, argparse._SubParsersAction)
        )
        validate_demux_args(args, subparsers.choices["demux"])

    if args.subcommand == "trim":
        subparsers = next(
            action for action in parser._actions
            if isinstance(action, argparse._SubParsersAction)
        )
        validate_trim_args(args, subparsers.choices["trim"])

    if args.subcommand == "denovo":
        subparsers = next(
            action for action in parser._actions
            if isinstance(action, argparse._SubParsersAction)
        )
        validate_denovo_args(args, subparsers.choices["denovo"])

    # LOGGING: -----------------------------------------------------
    if hasattr(args, "log_level"):
        set_log_level(args.log_level)

    if args.subcommand is None:
        # NO SUBCOMMAND: print help
        parser.print_help()
        sys.exit(0)
    else:
        run_subcommand(args)



def run_subcommand(args, _exit=True):
    # DEMUX: -------------------------------------------------------
    if args.subcommand == "demux":

        module = importlib.import_module("..demuxer", package=__package__)
        run_demuxer = getattr(module, "run_demuxer")

        logger.info("---------------------------------------------------------")
        logger.info("----- ipyrad2 demux: demultiplexing reads to samples -----")
        logger.info("---------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_demuxer(
            fastqs=args.fastqs,
            outdir=args.out,
            barcodes=args.barcodes,
            cutsite_1=args.cutsite_1,
            cutsite_2=args.cutsite_2,
            max_mismatch=args.max_mismatch,
            barcode_boundary_slack=args.barcode_boundary_slack,
            chunksize=args.chunksize,
            i7=args.i7,
            disable_infer_cutsite_motifs=args.disable_infer_cutsite_motifs,
            merge_technical_replicates=args.merge_technical_replicates,
            cores=args.cores,
            max_reads=args.max_reads,
            max_reads_kmer=args.max_reads_kmer,
            log_level=args.log_level,
            pigz=args.pigz,
            force=args.force,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # TRIM: -------------------------------------------------------
    if args.subcommand == "trim":

        module = importlib.import_module("..trimmer", package=__package__)
        run_trimmer = getattr(module, "run_trimmer")

        logger.info("----------------------------------------------------------")
        logger.info("----- ipyrad2 trim: quality, adapter, and cutsite motif trimming -----")
        logger.info("----------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
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
            cutsite_motifs=(args.cutsite_1, args.cutsite_2),
            max_reads_kmer=args.max_reads_kmer,
            phred64=args.phred64,
            disable_infer_cutsite_motifs=args.disable_infer_cutsite_motifs,
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

        module = importlib.import_module("..denovo", package=__package__)
        run_denovo = getattr(module, "run_denovo")

        logger.info("------------------------------------------------------------")
        logger.info("----- ipyrad2 denovo: construct locus reference library -----")
        logger.info("------------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_denovo(
            fastqs=args.fastqs,
            outdir=args.out,
            within_similarity=args.within_similarity,
            across_similarity=args.across_similarity,
            min_derep_size=args.min_derep_size,
            min_length=args.min_length,
            min_merge_overlap=args.min_merge_overlap,
            max_merge_diffs=args.max_merge_diffs,
            cores=args.cores,
            threads=args.threads,
            no_alignment=args.no_alignment,
            force=args.force,
            imap=args.imap,
            use_all_samples=args.use_all_samples,
            allow_reverse_complement=args.allow_reverse_complement,
            delim_str=args.delim_str,
            delim_idx=args.delim_idx,
            keep_intermediates=args.keep_intermediates,
            log_level=args.log_level,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # MAP: --------------------------------------------------------
    if args.subcommand == "map":

        module = importlib.import_module("..mapper", package=__package__)
        run_mapper = getattr(module, "run_mapper")

        logger.info("--------------------------------------------------------------")
        logger.info("----- ipyrad2 map: map reads and write coordinate-sorted BAMs -----")
        logger.info("--------------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_mapper(
            fastqs=args.fastqs,
            reference=args.reference,
            outdir=args.out,
            imap=args.imap,
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

        module = importlib.import_module("..assembler", package=__package__)
        run_assembler = getattr(module, "run_assembler")

        logger.info("-----------------------------------------------------------")
        logger.info("----- ipyrad2 assemble: delimit loci and call variants -----")
        logger.info("-----------------------------------------------------------")
        logger.info(f"CMD: {format_logged_command(sys.argv[1:])}")
        run_assembler(
            rad_bams=args.rad_bams,
            wgs_bams=args.wgs_bams,
            reference=args.reference,
            outdir=args.out,
            name=args.name,
            loci_bed=args.loci_bed,
            min_map_q=args.min_map_q,
            max_tlen=args.max_tlen,
            max_softclip=args.max_softclip,
            max_nm=args.max_nm,
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
            max_sample_hetero_frequency=args.max_sample_hetero_frequency,
            softclip_len_threshold=args.softclip_len_threshold,
            softclip_frac_max=args.softclip_frac_max,
            depth_z_max=args.depth_z_max,
            third_frac_cut=args.third_frac_cut,
            min_3allele_sites=args.min_3allele_sites,
            maf_threshold=args.maf_threshold,
            max_sites_above_maf=args.max_sites_above_maf,
            paralog_fail_frac_max=args.paralog_fail_frac_max,
            populations=args.populations,
            rename_bams=args.rename_bams,
            masks=args.masks,
            cores=args.cores,
            threads=args.threads,
            force=args.force,
            log_level=args.log_level,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # EXPORT / ANALYSIS: -----------------------------------------
    if args.subcommand in ANALYSIS_TOOL_NAMES:
        run_analysis_tool(args, _exit=_exit)
        return

    if args.subcommand in RESERVED_TOOL_NAMES:
        run_analysis_tool(args, _exit=_exit)
        return


if __name__ == "__main__":

    main()
