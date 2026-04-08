import argparse
from pathlib import Path

from ipyrad2.cli.cli_main import setup_parsers


def _get_assemble_parser() -> argparse.ArgumentParser:
    parser = setup_parsers()
    subparsers = next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    return subparsers.choices["assemble"]


def _assert_help_entry_is_single_line(help_text: str, option_label: str) -> None:
    lines = help_text.splitlines()
    idx = next(i for i, line in enumerate(lines) if line.strip().startswith(option_label))
    if idx + 1 >= len(lines):
        return
    next_line = lines[idx + 1].strip()
    is_continuation = bool(next_line) and not next_line.startswith("-") and not next_line.endswith(":")
    assert not is_continuation, f"{option_label} help wrapped to another line: {lines[idx + 1]!r}"


def test_assemble_help_groups_examples_and_current_descriptions() -> None:
    help_text = _get_assemble_parser().format_help()

    expected_sections = [
        "Core inputs:",
        "Mapped-read filters:",
        "Locus BED delimiting:",
        "Locus and variant filters:",
        "Paralog filters:",
        "Sample naming, grouping, and masks:",
        "Performance and overwrite:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    expected_order = [
        "-d, --rad-bams",
        "-w, --wgs-bams",
        "-r, --reference",
        "-b, --loci-bed",
        "-n, --name",
        "-o, --out",
        "-qm, --min-map-q",
        "-ms, --max-softclip",
        "-me, --max-nm",
        "-mt, --max-tlen",
        "-m, --min-locus-sample-coverage",
        "-z, --min-locus-length",
        "-g, --min-locus-merge-distance",
        "-qb, --min-base-q",
        "-qs, --min-site-q",
        "-qg, --min-geno-q",
        "-s, --min-sample-depth",
        "-u, --max-locus-hetero-frequency",
        "-y, --max-locus-variant-frequency",
        "-a, --min-locus-trim-sample-coverage",
        "--depth-z-max",
        "--softclip-len-threshold",
        "--softclip-frac-max",
        "--third-frac-cut",
        "--min-3allele-sites",
        "--maf-threshold",
        "--max-sites-above-maf",
        "--paralog-fail-frac-max",
        "--max-sample-hetero-frequency",
        "-p, --populations",
        "--rename-bams",
        "-x, --masks",
        "-c, --cores",
        "-t, --threads",
        "-f, --force",
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

    assert "ipyrad2 assemble: delimit loci, call variants, and write outputs" in help_text
    assert "$ ipyrad2 assemble -d BAMS/RAD/*.bam -r REF.fa -o OUT -m 4 -qm 20" in help_text
    assert "$ ipyrad2 assemble -d BAMS/RAD/*.bam -r REF.fa -p pops.tsv -o OUT" in help_text
    assert "$ ipyrad2 assemble -d BAMS/RAD/*.bam -r REF.fa --rename-bams rename.tsv -o OUT" in help_text
    assert "RAD BAM inputs that delimit loci unless --loci-bed is provided" in help_text
    assert "Discard mapped reads with MAPQ below this threshold." in help_text
    assert "Locus BED delimiting:" in help_text
    assert "Min third-allele fraction at a SNP site" in help_text
    assert "heterozygous/IUPAC plus masked-N-at-variable-site" in help_text
    assert "BED of loci to assemble instead of delimiting shared loci from RAD samples." in help_text
    assert "Population file for grouped calling; sample/group or classic pop_assign format." in help_text
    assert "mapping BAM basenames to final sample names" in help_text
    assert "ipyrad assemble -d" not in help_text
    assert "--ref REF" not in help_text
    assert "--out OUT" not in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")

    single_line_options = [
        "-d, --rad-bams",
        "-w, --wgs-bams",
        "-r, --reference",
        "-b, --loci-bed",
        "-o, --out",
        "-qm, --min-map-q",
        "-ms, --max-softclip",
        "-me, --max-nm",
        "-mt, --max-tlen",
        "-m, --min-locus-sample-coverage",
        "-z, --min-locus-length",
        "-g, --min-locus-merge-distance",
        "-qb, --min-base-q",
        "-qs, --min-site-q",
        "-qg, --min-geno-q",
        "-s, --min-sample-depth",
        "-u, --max-locus-hetero-frequency",
        "-y, --max-locus-variant-frequency",
        "-a, --min-locus-trim-sample-coverage",
        "--depth-z-max",
        "--softclip-len-threshold",
        "--softclip-frac-max",
        "--third-frac-cut",
        "--min-3allele-sites",
        "--maf-threshold",
        "--max-sites-above-maf",
        "--paralog-fail-frac-max",
        "-x, --masks",
        "-t, --threads",
    ]
    for option_label in single_line_options:
        _assert_help_entry_is_single_line(help_text, option_label)


def test_assemble_parser_defaults_match_current_cli() -> None:
    args = setup_parsers().parse_args(["assemble", "-d", "a.bam", "-r", "ref.fa"])

    assert args.subcommand == "assemble"
    assert args.rad_bams == [Path("a.bam")]
    assert args.wgs_bams is None
    assert args.reference == Path("ref.fa")
    assert args.loci_bed is None
    assert args.name == "assembly"
    assert args.out == Path("OUT")
    assert args.min_map_q == 10
    assert args.max_softclip is None
    assert args.max_nm is None
    assert args.max_tlen is None
    assert not hasattr(args, "require_same_scaffold")
    assert args.min_locus_sample_coverage == 4
    assert args.min_locus_length == 25
    assert args.min_locus_merge_distance == 300
    assert args.min_base_q == 13
    assert args.min_site_q == 13
    assert args.min_geno_q == 13
    assert args.min_sample_depth == 1
    assert args.max_locus_hetero_frequency == 0.3
    assert args.max_locus_variant_frequency == 1.0
    assert args.min_locus_trim_sample_coverage == 4
    assert args.depth_z_max == 7.0
    assert args.softclip_len_threshold == 20
    assert args.softclip_frac_max == 0.5
    assert args.third_frac_cut == 0.10
    assert args.min_3allele_sites == 2
    assert args.maf_threshold == 0.20
    assert args.max_sites_above_maf == 8
    assert args.paralog_fail_frac_max == 0.10
    assert args.max_sample_hetero_frequency == 0.10
    assert args.populations is None
    assert args.rename_bams is None
    assert args.masks is None
    assert args.cores == 6
    assert args.threads == 3
    assert args.force is False
    assert args.log_level == "INFO"


def test_assemble_parser_accepts_loci_bed_with_only_wgs_bams() -> None:
    args = setup_parsers().parse_args(
        ["assemble", "-w", "wgs.bam", "-b", "loci.bed", "-r", "ref.fa"]
    )

    assert args.subcommand == "assemble"
    assert args.rad_bams is None
    assert args.wgs_bams == [Path("wgs.bam")]
    assert args.loci_bed == Path("loci.bed")
    assert args.reference == Path("ref.fa")
