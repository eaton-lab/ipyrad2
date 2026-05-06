#!/usr/bin/env python

"""Top-level CLI parser and subcommand dispatcher."""

import argparse
import importlib
import sys
from typing import Optional, Sequence, Tuple

from loguru import logger

from .common import RAW_HELP_FORMATTER
from .command_log import format_logged_command
from ..utils.logger import set_log_level
from ..utils.exceptions import IPyradError
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


_CORE_SUBCOMMAND_SPECS = {
    "demux": {
        "module": ".cli_demux",
        "setup": "_setup_demux_subparser",
        "validator": "validate_demux_args",
        "header": "ipyrad2 demux: demultiplex pooled reads to sample files by barcode or index",
    },
    "trim": {
        "module": ".cli_trim",
        "setup": "_setup_trim_subparser",
        "validator": "validate_trim_args",
        "header": "ipyrad2 trim: trim reads for quality, adapters, and cutsite motifs",
    },
    "denovo": {
        "module": ".cli_denovo",
        "setup": "_setup_denovo_subparser",
        "validator": "validate_denovo_args",
        "header": "ipyrad2 denovo: construct a reference locus library",
    },
    "map": {
        "module": ".cli_map",
        "setup": "_setup_map_subparser",
        "validator": None,
        "header": "ipyrad2 map: map reads and write coordinate-sorted BAM files",
    },
    "assemble": {
        "module": ".cli_assemble",
        "setup": "_setup_assemble_subparser",
        "validator": None,
        "header": "ipyrad2 assemble: delimit loci, call variants, and write outputs",
    },
}

_CORE_TOOL_NAMES = tuple(_CORE_SUBCOMMAND_SPECS)
ANALYSIS_TOOL_NAMES = (
    "wex",
    "lex",
    "snpex",
    "vcf2hdf5",
    "pca",
    "snmf",
    "dapc",
    "admixture",
    "popgen",
    "baba",
    "bpp",
)
RESERVED_TOOL_NAMES = (
    "treeslider",
)
ALL_TOP_LEVEL_COMMANDS = _CORE_TOOL_NAMES + ANALYSIS_TOOL_NAMES + RESERVED_TOOL_NAMES


class TopLevelParser(argparse.ArgumentParser):
    """Root parser with grouped top-level help text."""

    def format_help(self) -> str:
        return TOP_LEVEL_HELP


def setup_parsers() -> argparse.ArgumentParser:
    """Setup and return an ArgumentParser w/ subcommands."""
    return _setup_parsers()


def _load_cli_module(module_name: str):
    """Import one CLI helper module lazily."""
    return importlib.import_module(module_name, package=__package__)


def _setup_parsers(
    selected_subcommands: Optional[Sequence[str]] = None,
) -> argparse.ArgumentParser:
    """Setup and return an ArgumentParser with lazily imported subcommands."""
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

    selected = None if selected_subcommands is None else set(selected_subcommands)

    for name in _CORE_TOOL_NAMES:
        if selected is not None and name not in selected:
            continue
        spec = _CORE_SUBCOMMAND_SPECS[name]
        module = _load_cli_module(spec["module"])
        setup = getattr(module, spec["setup"])
        setup(subparser, f"{HEADER}\n{spec['header']}")

    if selected is None or selected.intersection(ANALYSIS_TOOL_NAMES + RESERVED_TOOL_NAMES):
        analysis_module = _load_cli_module(".cli_analysis")
        analysis_module._setup_analysis_tool_subparsers(
            subparser,
            HEADER,
            selected_tools=None if selected is None else selected,
        )
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


def _get_subparsers_action(parser: argparse.ArgumentParser) -> argparse._SubParsersAction:
    """Return the top-level subparsers action."""
    return next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )


def _validate_subcommand_args(args, parser: argparse.ArgumentParser) -> None:
    """Run any subcommand-specific argparse validation hooks."""
    if args.subcommand not in _CORE_SUBCOMMAND_SPECS:
        return
    validator_name = _CORE_SUBCOMMAND_SPECS[args.subcommand]["validator"]
    if validator_name is None:
        return
    module = _load_cli_module(_CORE_SUBCOMMAND_SPECS[args.subcommand]["module"])
    validator = getattr(module, validator_name)
    subparsers = _get_subparsers_action(parser)
    validator(args, subparsers.choices[args.subcommand])


def _print_top_level_help() -> None:
    """Print the curated top-level help text."""
    print(TOP_LEVEL_HELP)


def _is_top_level_help_argv(argv: Sequence[str]) -> bool:
    """Return True when argv requests top-level help or no command."""
    return not argv or argv[0] in {"-h", "--help"}


def _help_only_selected_subcommands(argv: Sequence[str]) -> Optional[Tuple[str, ...]]:
    """Return a single-command selection when argv only requests one help screen."""
    if len(argv) < 2:
        return None
    if argv[0] not in ALL_TOP_LEVEL_COMMANDS:
        return None
    if not any(token in {"-h", "--help"} for token in argv[1:]):
        return None
    return (argv[0],)


def command_line(argv: Optional[Sequence[str]] = None):
    argv = list(sys.argv[1:] if argv is None else argv)

    if _is_top_level_help_argv(argv):
        _print_top_level_help()
        sys.exit(0)

    parser = _setup_parsers(selected_subcommands=_help_only_selected_subcommands(argv))
    args = parser.parse_args(argv)
    _validate_subcommand_args(args, parser)

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
        logged_command = format_logged_command(sys.argv[1:])

        module = importlib.import_module("..demuxer", package=__package__)
        run_demuxer = getattr(module, "run_demuxer")

        logger.info("---------------------------------------------------------")
        logger.info("----- ipyrad2 demux: demultiplexing reads to samples -----")
        logger.info("---------------------------------------------------------")
        logger.info(f"CMD: {logged_command}")
        run_demuxer(
            fastqs=args.fastqs,
            outdir=args.out,
            barcodes=args.barcodes,
            delim_str=args.delim_str,
            delim_idx=args.delim_idx,
            cutsite_1=args.cutsite_1,
            cutsite_2=args.cutsite_2,
            max_mismatch=args.max_mismatch,
            barcode_boundary_slack=args.barcode_boundary_slack,
            allow_leading_barcode_deletion=args.allow_leading_barcode_deletion,
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
            logged_command=logged_command,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # TRIM: -------------------------------------------------------
    if args.subcommand == "trim":
        logged_command = format_logged_command(sys.argv[1:])

        module = importlib.import_module("..trimmer", package=__package__)
        run_trimmer = getattr(module, "run_trimmer")

        logger.info("----------------------------------------------------------")
        logger.info("----- ipyrad2 trim: quality, adapter, and cutsite motif trimming -----")
        logger.info("----------------------------------------------------------")
        logger.info(f"CMD: {logged_command}")
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
            logged_command=logged_command,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # DENOVO: --------------------------------------------------------
    if args.subcommand == "denovo":
        logged_command = format_logged_command(sys.argv[1:])

        module = importlib.import_module("..denovo", package=__package__)
        run_denovo = getattr(module, "run_denovo")

        logger.info("------------------------------------------------------------")
        logger.info("----- ipyrad2 denovo: construct locus reference library -----")
        logger.info("------------------------------------------------------------")
        logger.info(f"CMD: {logged_command}")
        run_denovo(
            fastqs=args.fastqs,
            outdir=args.out,
            within_similarity=args.within_similarity,
            across_similarity=args.across_similarity,
            query_cov=args.query_cov,
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
            logged_command=logged_command,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # MAP: --------------------------------------------------------
    if args.subcommand == "map":
        logged_command = format_logged_command(sys.argv[1:])

        module = importlib.import_module("..mapper", package=__package__)
        run_mapper = getattr(module, "run_mapper")

        logger.info("--------------------------------------------------------------")
        logger.info("----- ipyrad2 map: map reads and write coordinate-sorted BAMs -----")
        logger.info("--------------------------------------------------------------")
        logger.info(f"CMD: {logged_command}")
        run_mapper(
            fastqs=args.fastqs,
            reference=args.reference,
            outdir=args.out,
            imap=args.imap,
            unmate=args.unmate,
            mark_dups_by_coords=args.mark_dups_by_coords,
            mark_dups_by_umis=args.mark_dups_by_umis,
            cores=args.cores,
            threads=args.threads,
            force=args.force,
            reindex_reference=args.reindex_reference,
            delim_str=args.delim_str,
            delim_idx=args.delim_idx,
            log_level=args.log_level,
            logged_command=logged_command,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # ASSEMBLE: ---------------------------------------------------
    if args.subcommand == "assemble":
        logged_command = format_logged_command(sys.argv[1:])

        module = importlib.import_module("..assembler", package=__package__)
        run_assembler = getattr(module, "run_assembler")

        logger.info("-----------------------------------------------------------")
        logger.info("----- ipyrad2 assemble: delimit loci and call variants -----")
        logger.info("-----------------------------------------------------------")
        logger.info(f"CMD: {logged_command}")
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
            min_aligned_len=args.min_aligned_len,
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
            subsample=args.subsample,
            populations=args.populations,
            rename_bams=args.rename_bams,
            masks=args.masks,
            cores=args.cores,
            threads=args.threads,
            force=args.force,
            log_level=args.log_level,
            logged_command=logged_command,
        )
        if _exit: sys.exit(0)  # noqa: E701

    # EXPORT / ANALYSIS: -----------------------------------------
    if args.subcommand in ANALYSIS_TOOL_NAMES:
        analysis_module = _load_cli_module(".cli_analysis")
        analysis_module.run_analysis_tool(args, _exit=_exit)
        return

    if args.subcommand in RESERVED_TOOL_NAMES:
        analysis_module = _load_cli_module(".cli_analysis")
        analysis_module.run_analysis_tool(args, _exit=_exit)
        return


if __name__ == "__main__":

    main()
