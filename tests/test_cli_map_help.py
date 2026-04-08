import argparse
from pathlib import Path

from ipyrad2.cli.cli_main import setup_parsers


def _get_map_parser() -> argparse.ArgumentParser:
    parser = setup_parsers()
    subparsers = next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    return subparsers.choices["map"]


def test_map_help_groups_and_examples_are_updated() -> None:
    help_text = _get_map_parser().format_help()

    expected_sections = [
        "Core inputs:",
        "Duplicate removal:",
        "Sample naming and grouping:",
        "Performance and overwrite:",
        "Logging:",
    ]
    positions = [help_text.index(section) for section in expected_sections]
    assert positions == sorted(positions)

    expected_order = [
        "-d, --fastqs",
        "-r, --reference",
        "-o, --out",
        "-m, --mark-dups-by-coords",
        "-u, --mark-dups-by-umis",
        "-i, --imap",
        "-dx, --delim-str",
        "-di, --delim-idx",
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

    assert "ipyrad2 map: map reads and write coordinate-sorted BAM files" in help_text
    assert "$ ipyrad2 map -d DATA/*.fastq.gz -r REF.fa -o BAMS -i IMAP.tsv" in help_text
    assert "$ ipyrad2 map -d DATA/*.fastq.gz -r REF.fa -o BAMS -u" in help_text
    assert "Reference FASTA to index and map against." in help_text
    assert "Output directory for coordinate-sorted BAMs and map stats." in help_text
    assert "ipyrad map" not in help_text
    assert "trimmed read files" not in help_text
    assert "--min-map-q" not in help_text
    assert "--max-edit-dist" not in help_text
    assert "--max-soft-clip" not in help_text
    assert help_text.index("-h, --help") > help_text.index("Logging:")


def test_map_parser_defaults_are_unchanged() -> None:
    args = setup_parsers().parse_args(["map", "-d", "a.fastq.gz", "-r", "ref.fa"])

    assert args.subcommand == "map"
    assert args.fastqs == [Path("a.fastq.gz")]
    assert args.reference == Path("ref.fa")
    assert args.out == Path("MAPPED")
    assert args.imap is None
    assert args.mark_dups_by_coords is False
    assert args.mark_dups_by_umis is False
    assert args.cores == 6
    assert args.threads == 3
    assert args.force is False
    assert args.delim_str is None
    assert args.delim_idx == 1
    assert args.log_level == "INFO"
