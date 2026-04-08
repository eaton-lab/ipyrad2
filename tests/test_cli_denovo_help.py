import argparse
from pathlib import Path

import pytest

from ipyrad2.cli.cli_main import setup_parsers


def _get_denovo_parser() -> argparse.ArgumentParser:
    parser = setup_parsers()
    subparsers = next(
        action for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    return subparsers.choices["denovo"]


def test_denovo_help_uses_grouped_layout_and_updated_examples() -> None:
    help_text = _get_denovo_parser().format_help()

    assert "ipyrad2 denovo: construct a reference locus library" in help_text
    assert "Core inputs:" in help_text
    assert "Sample selection:" in help_text
    assert "Clustering and consensus:" in help_text
    assert "Sample naming and library type:" in help_text
    assert "Runtime:" in help_text
    assert "Logging:" in help_text
    assert "$ ipyrad2 denovo -d DATA/*.fastq.gz -o output-denovo" in help_text
    assert "$ ipyrad2 denovo -d DATA/*.fastq.gz -o OUT --imap denovo.imap.tsv" in help_text
    assert "$ ipyrad2 denovo -d DATA/*.fastq.gz -o OUT --no-alignment" in help_text
    assert "--keep-intermediates" in help_text
    assert "--imap" in help_text
    assert "--use-all-samples" in help_text
    assert "--no-alignment" in help_text
    assert "ipyrad denovo" not in help_text


def test_denovo_parser_defaults_are_updated() -> None:
    args = setup_parsers().parse_args(["denovo", "-d", "a.fastq.gz"])

    assert args.subcommand == "denovo"
    assert args.fastqs == [Path("a.fastq.gz")]
    assert args.out == Path("output-denovo")
    assert args.within_similarity == 0.95
    assert args.across_similarity == 0.85
    assert args.min_derep_size == 5
    assert args.min_length == 35
    assert args.min_merge_overlap == 20
    assert args.max_merge_diffs == 4
    assert args.allow_reverse_complement is False
    assert args.cores == 6
    assert args.threads == 3
    assert args.no_alignment is False
    assert args.keep_intermediates is False
    assert args.force is False
    assert args.imap is None
    assert args.use_all_samples is False
    assert args.delim_str is None
    assert args.delim_idx == 1
    assert args.log_level == "INFO"


def test_denovo_parser_rejects_removed_graph_splitter_flag() -> None:
    with pytest.raises(SystemExit):
        setup_parsers().parse_args(
            ["denovo", "-d", "a.fastq.gz", "--graph-splitter", "constrained"]
        )
