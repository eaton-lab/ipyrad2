#!/usr/bin/env python

"""Assemble command-line parser."""

import argparse
from pathlib import Path
from .common import RAW_HELP_FORMATTER


EPILOG = """\
Examples
--------
$ ipyrad2 assemble -d BAMS/RAD/*.bam -r REF.fa -o OUT -m 4 -qm 20
$ ipyrad2 assemble -d BAMS/RAD/*.bam -w BAMS/WGS/*.bam -r REF.fa -o OUT -m 4 -qm 20
$ ipyrad2 assemble -d BAMS/RAD/*.bam -r REF.fa -b loci.bed -o OUT --max-tlen 2000
$ ipyrad2 assemble -d BAMS/RAD/*.bam -w BAMS/WGS/*.bam -r REF.fa --subsample keep.tsv -o OUT
$ ipyrad2 assemble -d BAMS/RAD/*.bam -r REF.fa -p pops.tsv -o OUT
$ ipyrad2 assemble -d BAMS/RAD/*.bam -r REF.fa --rename rename.tsv -o OUT
"""


def _setup_assemble_subparser(subparsers: argparse._SubParsersAction, header: str = None) -> None:
    """Add `ipyrad2 assemble` subcommand parser."""
    tool = subparsers.add_parser(
        "assemble",
        description=header,
        help="Delimit loci, call variants, and write assembled outputs.",
        epilog=EPILOG,
        formatter_class=RAW_HELP_FORMATTER,
        add_help=False,
    )
    core = tool.add_argument_group("Core inputs")
    mapped = tool.add_argument_group("Mapped-read filters")
    bed = tool.add_argument_group("Locus BED delimiting")
    locus = tool.add_argument_group("Locus and variant filters")
    paralogs = tool.add_argument_group("Paralog filters")
    naming = tool.add_argument_group("Sample naming, grouping, and masks")
    performance = tool.add_argument_group("Performance and overwrite")
    logging = tool.add_argument_group("Logging")

    core.add_argument(
        "-d", "--rad-bams", metavar="Path", type=Path, nargs="*",
        help="RAD BAM inputs that delimit loci unless --loci-bed is provided; also assembled.",
    )
    core.add_argument(
        "-w", "--wgs-bams", metavar="Path", type=Path, nargs="*",
        help="Optional WGS BAM inputs assembled only within loci defined by RAD samples or by --loci-bed.",
    )
    core.add_argument(
        "-r", "--reference", metavar="Path", type=Path, required=True,
        help="Reference FASTA that was used for mapping and is reused here for locus extraction and calling.",
    )
    core.add_argument(
        "-b", "--loci-bed", metavar="Path", type=Path,
        help="BED of loci to assemble instead of delimiting shared loci from RAD samples.",
    )
    core.add_argument(
        "-n", "--name", metavar="str", type=str, default="assembly",
        help="Prefix for assembled output files. [default=%(default)s]",
    )
    core.add_argument(
        "-o", "--out", metavar="Path", type=Path, default="./OUT",
        help="Output directory for assembled loci, VCFs, and stats. [default=%(default)s]",
    )

    mapped.add_argument(
        "-qm", "--min-map-q", metavar="int", type=int, default=30,
        help="Discard mapped reads with MAPQ below this threshold. [default=%(default)s]",
    )
    mapped.add_argument(
        "-ms", "--max-softclip", metavar="int", type=int, default=None,
        help="Discard mapped reads with more than this many soft-clipped bases. [default=%(default)s]",
    )
    mapped.add_argument(
        "-ml", "--min-aligned-len", metavar="int", type=int, default=None,
        help="Discard mapped reads with qlen-sclen below this threshold. [default=%(default)s]",
    )
    mapped.add_argument(
        "-me", "--max-nm", metavar="int", type=int, default=None,
        help="Discard mapped reads with NM above this threshold. [default=%(default)s]",
    )
    mapped.add_argument(
        "-mt", "--max-tlen", metavar="int", type=int, default=1000,
        help="Discard pairs with absolute TLEN above this threshold (PE data only). [default=%(default)s]",
    )
    # mapped.add_argument(
    #     "--require-same-scaffold", action="store_true",
    #     help="For paired data, require both mates map to the same scaffold.",
    # )

    bed.add_argument(
        "-m", "--min-locus-sample-coverage", metavar="int", type=int, default=4,
        help="Min RAD samples with data required to retain a locus. [default=%(default)s]",
    )
    bed.add_argument(
        "-z", "--min-locus-length", metavar="int", type=int, default=25,
        help="Min locus length after edge trimming and filtering. [default=%(default)s]",
    )
    bed.add_argument(
        "-g", "--min-locus-merge-distance", metavar="int", type=int, default=300,
        help="Merge nearby locus intervals within this distance when delimiting shared loci. [default=%(default)s]",
    )


    locus.add_argument(
        "-qb", "--min-base-q", metavar="int", type=int, default=13,
        help="Min base quality used by mpileup and related downstream calling steps. [default=%(default)s]",
    )
    locus.add_argument(
        "-qs", "--min-site-q", metavar="int", type=int, default=13,
        help="Min site QUAL required for variant sites retained after joint calling. [default=%(default)s]",
    )
    locus.add_argument(
        "-qg", "--min-geno-q", metavar="int", type=int, default=13,
        help="Min per-sample genotype quality retained in the filtered VCF. [default=%(default)s]",
    )
    locus.add_argument(
        "-s", "--min-sample-depth", metavar="int", type=int, default=5,
        help="Min within-sample read depth to keep a genotype call instead of masking it. [default=%(default)s]",
    )
    locus.add_argument(
        "-u", "--max-locus-hetero-frequency", metavar="float", type=float, default=0.3,
        help="Max fraction of samples heterozygous at the same site before marking as a paralog [default=%(default)s]",
    )
    locus.add_argument(
        "-y", "--max-locus-variant-frequency", metavar="float", type=float, default=1.0,
        help="Max fraction of sites in a locus that can be variant before the locus is filtered. [default=%(default)s]",
    )
    locus.add_argument(
        "-a", "--min-locus-trim-sample-coverage", metavar="int", type=int, default=4,
        help="Min number of samples with non-N calls required to keep positions at locus edges. [default=%(default)s]",
    )

    paralogs.add_argument(
        "--depth-z-max", metavar="float", type=float, default=7.0,
        help="Max per-sample read-depth z-score to tag a locus as a high-depth outlier. [default=%(default)s]",
    )
    paralogs.add_argument(
        "--softclip-len-threshold", metavar="int", type=int, default=20,
        help="Max soft-clipped bases to label as read as 'highly clipped'. [default=%(default)s]",
    )
    paralogs.add_argument(
        "--softclip-frac-max", metavar="float", type=float, default=0.5,
        help="Max per-sample frac. reads that are highly clipped to tag a locus as paralog-like. [default=%(default)s]",
    )
    paralogs.add_argument(
        "--third-frac-cut", metavar="float", type=float, default=0.10,
        help="Min third-allele fraction at a SNP site to count as strong multi-allelic evidence. [default=%(default)s]",
    )
    paralogs.add_argument(
        "--min-3allele-sites", metavar="int", type=int, default=2,
        help="Min strong 3-allele SNP sites before a locus is tagged as paralog-like. [default=%(default)s]",
    )
    paralogs.add_argument(
        "--maf-threshold", metavar="float", type=float, default=0.20,
        help="Min minor-allele frequency counted as excess allelic variation in a locus. [default=%(default)s]",
    )
    paralogs.add_argument(
        "--max-sites-above-maf", metavar="int", type=int, default=8,
        help="Max SNP sites above --maf-threshold before a locus is tagged as paralog-like. [default=%(default)s]",
    )
    paralogs.add_argument(
        "--paralog-fail-frac-max", metavar="float", type=float, default=0.10,
        help="Max frac. samples with data allowed to fail before a locus is dropped globally. [default=%(default)s]",
    )
    paralogs.add_argument(
        "--max-sample-hetero-frequency", metavar="float", type=float, default=0.10,
        help="Max frac. of heterozygous or QUAL-masked SNP sites in a sample at a locus. [default=%(default)s]",
    )

    naming.add_argument(
        "--subsample", metavar="Path", type=Path,
        help="File whose first column selects BAM filenames or sample names; extra columns ignored",
    )
    naming.add_argument(
        "--populations", metavar="Path", type=Path,
        help="File mapping BAM basenames to group names for population-level variant calls",
    )
    naming.add_argument(
        "--rename", metavar="Path", type=Path,
        help="File mapping BAM basenames to new names for outputs; overrides BAM headers",
    )
    naming.add_argument(
        "--masks", metavar="str", nargs="*", type=str,
        help="Optional site patterns to mask in final assembled sequences. [default=None]",
    )

    performance.add_argument(
        "-c", "--cores", metavar="int", type=int, default=6,
        help="Maximum total cores to use. [default=%(default)s]",
    )
    performance.add_argument(
        "-t", "--threads", metavar="int", type=int, default=3,
        help="Threads per multithreaded job; larger values reduce parallel job count. [default=%(default)s]",
    )
    performance.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite assemble outputs for this run.",
    )
    performance.add_argument(
        "--keep-tmpdir", action="store_true",
        help="Keep the assemble tmpdir after success; useful for testing.",
    )


    logging.add_argument(
        "-l", "--log-level", metavar="str", type=str, default="INFO",
        help="Logging verbosity. [default=%(default)s]",
    )
    logging.add_argument(
        "-h", "--help", action="help",
        help="Show this help message and exit.",
    )
