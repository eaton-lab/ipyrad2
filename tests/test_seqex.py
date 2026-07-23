import json
from pathlib import Path

import h5py
import numpy as np
import pytest
from loguru import logger

import ipyrad2.analysis as analysis_api
import ipyrad2.analysis.extracters as extracters_api
from ipyrad2.analysis.extracters.seqex import SeqexEngine
from ipyrad2.analysis.extracters.seqex import run_seqex
from ipyrad2.cli.cli_main import setup_parsers
from ipyrad2.utils.exceptions import IPyradError


def _write_h5(
    path: Path,
    sequences: list[str],
    rows: list[tuple[int, int, int, int, int]],
    names: list[str] | None = None,
    scaffold_names: list[str] | None = None,
    scaffold_lengths: list[int] | None = None,
) -> Path:
    dtype = h5py.string_dtype(encoding="utf-8")
    phy = np.vstack(
        [np.frombuffer(sequence.encode(), dtype=np.uint8) for sequence in sequences]
    )
    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["scaffold_names"] = np.array(
            scaffold_names or ["chr1"],
            dtype=dtype,
        )
        io5.attrs["scaffold_lengths"] = np.array(
            scaffold_lengths or [len(sequences[0])]
        )
        io5.attrs["names"] = np.array(names or ["s1", "s2"], dtype=dtype)
        io5.create_dataset("phy", data=phy)
        phymap = io5.create_dataset("phymap", data=np.array(rows, dtype=np.uint64))
        phymap.attrs["columns"] = np.array(
            ["scaff", "phy0", "phy1", "pos0", "pos1"],
            dtype=dtype,
        )
    return path


def _run(tmp_path: Path, h5: Path, **kwargs):
    params = {
        "data": h5,
        "name": "alignment",
        "outdir": tmp_path / "OUT",
        "out_format": "phy",
        "windows": None,
        "max_loci": None,
        "random_seed": None,
        "min_length": None,
        "min_sample_coverage": 1,
        "max_sample_missing": 1.0,
        "exclude": None,
        "include_reference": False,
        "imap": None,
        "minmap": None,
        "concatenate": False,
        "split": False,
        "append_population": False,
        "print_scaffold_table": False,
        "stdout": False,
        "force": True,
        "logged_command": "ipyrad2 seqex",
    }
    params.update(kwargs)
    return run_seqex(**params)


def test_seqex_is_the_public_complete_locus_engine() -> None:
    assert analysis_api.SeqexEngine is SeqexEngine
    assert analysis_api.run_seqex is run_seqex
    assert extracters_api.SeqexEngine is SeqexEngine
    assert extracters_api.run_seqex is run_seqex
    with pytest.raises(AttributeError):
        getattr(analysis_api, "LocusExtracter")
    with pytest.raises(AttributeError):
        getattr(extracters_api, "LocusExtracter")
    with pytest.raises(AttributeError):
        getattr(analysis_api, "WindowExtracter")
    with pytest.raises(AttributeError):
        getattr(extracters_api, "WindowExtracter")


def test_seqex_cli_has_simplified_locus_interface() -> None:
    parser = setup_parsers()
    args = parser.parse_args(
        [
            "seqex",
            "-d",
            "data.hdf5",
            "-N",
            "10",
            "-s",
            "2",
            "-X",
            "-O",
            "fa",
            "-c",
            "3",
            "-a",
            "-i",
            "imap.tsv",
        ]
    )
    assert args.max_loci == 10
    assert args.random_seed == 2
    assert args.split is True
    assert args.out_format == "fa"
    assert args.cores == 3
    assert args.append_population is True
    help_text = parser._subparsers._group_actions[0].choices["seqex"].format_help()
    assert "--unit" not in help_text
    assert "--filter-scope" not in help_text
    assert "--clip" not in help_text
    assert "bpp" not in help_text


def test_seqex_filters_locus_then_sites_and_writes_multilocus_fasta(
    tmp_path: Path,
) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["ACGTACGT", "ACNNNCGT"],
        [(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
    )
    loci = _run(
        tmp_path,
        h5,
        out_format="fa",
        min_sample_coverage=2,
    )
    assert [locus.spec.label for locus in loci] == ["chr1:1-4", "chr1:5-8"]
    text = (tmp_path / "OUT" / "alignment.fa").read_text()
    assert ">s1|chr1:1-4\nAC\n" in text
    assert ">s2|chr1:1-4\nAC\n" in text
    assert ">s1|chr1:5-8\nCGT\n" in text
    assert ">s2|chr1:5-8\nCGT\n" in text
    stats = (tmp_path / "OUT" / "alignment.stats.txt").read_text()
    assert "total_sites_written: 5" in stats
    assert "total_bases_written: 10" in stats
    assert "full_matrix_bases: 10" in stats
    assert "non_missing_occupancy: 1.000000" in stats
    assert "max_samples: 2" in stats
    assert "mean_samples: 2.000000" in stats
    assert "selected_windows: none\n" in stats
    assert "\t" not in stats
    assert "sample  population  written_final" in stats
    assert "locus_index  locus" in stats
    stats_json = json.loads((tmp_path / "OUT" / "alignment.stats.json").read_text())
    assert stats_json["seqex_summary"]["output_layout"] == "multi-locus"
    assert stats_json["seqex_summary"]["max_loci"] is None
    assert stats_json["seqex_summary"]["selected_windows"] is None
    assert stats_json["output_summary"]["max_samples"] == 2
    assert stats_json["output_summary"]["mean_samples"] == 2.0
    assert stats_json["sample_occupancy"][0]["written_final"] is True
    assert stats_json["written_loci"][0] == {
        "locus_index": 1,
        "locus": "chr1:1-4",
        "source_locus": "chr1:1-4",
        "selected_window": "chr1:1-8",
        "clipped": False,
        "raw_samples": 2,
        "raw_sites": 4,
        "filtered_samples": 2,
        "filtered_sites": 2,
        "concat_start": None,
        "concat_end": None,
    }


def test_seqex_rejects_locus_when_no_site_meets_coverage(tmp_path: Path) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["AANNGGGG", "NNCCGGGG"],
        [(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
    )
    loci = _run(tmp_path, h5, min_sample_coverage=2)
    assert [locus.spec.label for locus in loci] == ["chr1:5-8"]
    stats = (tmp_path / "OUT" / "alignment.stats.txt").read_text()
    assert "rejected_site_coverage: 1" in stats


def test_seqex_coordinate_window_clips_overlapping_locus_automatically(
    tmp_path: Path,
) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["AAAACCCC", "AAAACCCC"],
        [(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
    )

    loci = _run(tmp_path, h5, windows=["chr1:6-6"])

    assert [locus.spec.label for locus in loci] == ["chr1:6-6"]
    assert loci[0].spec.clipped is True
    text = (tmp_path / "OUT" / "alignment.phy").read_text(encoding="utf-8")
    assert text.count("2 1\n") == 1
    stats = json.loads((tmp_path / "OUT" / "alignment.stats.json").read_text())
    assert stats["seqex_summary"]["clipping_mode"] == "automatic"
    assert stats["seqex_summary"]["coordinate_clipping_applied"] is True
    assert stats["written_loci"][0]["source_locus"] == "chr1:5-8"
    assert stats["written_loci"][0]["selected_window"] == "chr1:6-6"


def test_seqex_scaffold_selector_keeps_complete_loci(tmp_path: Path) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["AAAACCCC", "AAAACCCC"],
        [(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
    )

    loci = _run(tmp_path, h5, windows=["chr1"])

    assert [locus.spec.label for locus in loci] == ["chr1:1-4", "chr1:5-8"]
    assert not any(locus.spec.clipped for locus in loci)


def test_seqex_bed_and_region_apply_identical_clipping(tmp_path: Path) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["AACCGGTT", "AACCGGTT"],
        [(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
    )
    bed = tmp_path / "window.bed"
    bed.write_text("chr1\t2\t6\n", encoding="utf-8")

    region = _run(
        tmp_path,
        h5,
        outdir=tmp_path / "REGION",
        windows=["chr1:3-6"],
    )
    from_bed = _run(
        tmp_path,
        h5,
        outdir=tmp_path / "BED",
        windows=[str(bed)],
    )

    assert [locus.spec.label for locus in region] == ["chr1:3-4", "chr1:5-6"]
    assert [locus.spec.label for locus in from_bed] == [
        locus.spec.label for locus in region
    ]
    assert (tmp_path / "REGION" / "alignment.phy").read_text() == (
        tmp_path / "BED" / "alignment.phy"
    ).read_text()


def test_seqex_internal_clip_false_keeps_complete_overlapping_locus(
    tmp_path: Path,
) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["AAAACCCC", "AAAACCCC"],
        [(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
    )

    loci = _run(tmp_path, h5, windows=["chr1:6-6"], clip=False)

    assert [locus.spec.label for locus in loci] == ["chr1:5-8"]
    assert loci[0].spec.clipped is False


def test_seqex_disjoint_coordinate_windows_split_one_source_locus(
    tmp_path: Path,
) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["AACCGGTT", "AACCGGTT"],
        [(0, 0, 8, 1, 8)],
    )

    loci = _run(tmp_path, h5, windows=["chr1:2-3", "chr1:6-7"])

    assert [locus.spec.label for locus in loci] == ["chr1:2-3", "chr1:6-7"]
    text = (tmp_path / "OUT" / "alignment.phy").read_text()
    assert "AC" in text
    assert "GT" in text


def test_seqex_min_length_is_applied_to_clipped_fragment(tmp_path: Path) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["AACCGGTT", "AACCGGTT"],
        [(0, 0, 8, 1, 8)],
    )

    with pytest.raises(IPyradError, match="No loci passed"):
        _run(tmp_path, h5, windows=["chr1:2-3"], min_length=3)


def test_seqex_stdout_and_existing_output_protection(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["ACGT", "ACGT"],
        [(0, 0, 4, 1, 4)],
    )
    stdout_dir = tmp_path / "STDOUT"

    _run(tmp_path, h5, outdir=stdout_dir, stdout=True)

    captured = capsys.readouterr()
    assert captured.out.startswith("2 4\n")
    assert not (stdout_dir / "alignment.phy").exists()
    assert (stdout_dir / "alignment.stats.txt").exists()
    assert (stdout_dir / "alignment.stats.json").exists()

    file_dir = tmp_path / "FILES"
    _run(tmp_path, h5, outdir=file_dir)
    with pytest.raises(IPyradError, match="Output file already exists"):
        _run(tmp_path, h5, outdir=file_dir, force=False)


def test_seqex_sampling_is_reproducible_and_written_in_genomic_order(
    tmp_path: Path,
) -> None:
    rows = [(0, i, i + 1, i + 1, i + 1) for i in range(8)]
    h5 = _write_h5(tmp_path / "data.hdf5", ["ACGTACGT", "ACGTACGT"], rows)
    first = _run(
        tmp_path,
        h5,
        outdir=tmp_path / "ONE",
        max_loci=3,
        random_seed=42,
    )
    second = _run(
        tmp_path,
        h5,
        outdir=tmp_path / "TWO",
        max_loci=3,
        random_seed=42,
    )
    first_indices = [locus.spec.index for locus in first]
    assert first_indices == [locus.spec.index for locus in second]
    assert first_indices == sorted(first_indices)
    assert len(first_indices) == 3


def test_seqex_parallel_filtering_matches_serial_output(tmp_path: Path) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["ACNNACNN", "ACGTACGT"],
        [(0, 0, 2, 1, 2), (0, 4, 6, 5, 6)],
    )
    serial = _run(
        tmp_path,
        h5,
        outdir=tmp_path / "SERIAL",
        windows=["chr1:1-1", "chr1:5-5"],
        max_loci=1,
        random_seed=9,
        cores=1,
    )
    parallel = _run(
        tmp_path,
        h5,
        outdir=tmp_path / "PARALLEL",
        windows=["chr1:1-1", "chr1:5-5"],
        max_loci=1,
        random_seed=9,
        cores=2,
    )
    assert [locus.spec.label for locus in serial] == [
        locus.spec.label for locus in parallel
    ]
    assert serial[0].spec.clipped is True
    assert (tmp_path / "SERIAL" / "alignment.phy").read_text() == (
        tmp_path / "PARALLEL" / "alignment.phy"
    ).read_text()


def test_seqex_concatenate_applies_global_missingness_and_tracks_partitions(
    tmp_path: Path,
) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["AAAANNNNNNNN", "AAAACCCCGGGG"],
        [(0, 0, 4, 1, 4), (0, 4, 8, 5, 8), (0, 8, 12, 9, 12)],
    )
    messages: list[str] = []
    sink = logger.add(
        messages.append,
        format="{level.name}:{message}",
        level="DEBUG",
    )
    try:
        _run(
            tmp_path,
            h5,
            concatenate=True,
            max_sample_missing=0.5,
        )
    finally:
        logger.remove(sink)
    text = (tmp_path / "OUT" / "alignment.phy").read_text()
    assert text.splitlines()[0] == "1 12"
    assert "s1" not in text
    assert "AAAACCCCGGGG" in text
    stats = (tmp_path / "OUT" / "alignment.stats.txt").read_text()
    assert "chr1:1-4" in stats
    assert "chr1:9-12" in stats
    assert "total_sites_written: 12" in stats
    assert "total_bases_written: 12" in stats
    assert "full_matrix_bases: 12" in stats
    assert "non_missing_bases: 12" in stats
    assert "non_missing_occupancy: 1.000000" in stats
    assert "max_samples: 1" in stats
    assert "mean_samples: 1.000000" in stats
    stats_json = json.loads((tmp_path / "OUT" / "alignment.stats.json").read_text())
    assert stats_json["output_summary"]["max_samples"] == 1
    assert stats_json["output_summary"]["mean_samples"] == 1.0
    assert stats_json["sample_occupancy"][0]["written_final"] is False
    assert stats_json["sample_occupancy"][1]["written_final"] is True
    assert stats_json["written_loci"][0]["concat_start"] == 1
    assert stats_json["written_loci"][0]["concat_end"] == 4
    assert any(
        "DEBUG:-r dropped s1 from 2 written locus/loci" in msg for msg in messages
    )
    assert any(
        "INFO:-r dropped s1 from the concatenated alignment" in msg for msg in messages
    )
    assert any(
        "wrote 3 filtered loci concatenated into one PHYLIP alignment" in msg
        for msg in messages
    )
    assert any("wrote stats report to:" in msg for msg in messages)
    assert not any(".stats.json" in msg for msg in messages)


def test_seqex_split_nexus_and_append_population(tmp_path: Path) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["ACGT", "ACGT"],
        [(0, 0, 4, 1, 4)],
    )
    imap = tmp_path / "imap.tsv"
    imap.write_text("s1 popA\ns2 popB\n")
    _run(
        tmp_path,
        h5,
        out_format="nex",
        split=True,
        imap=imap,
        append_population=True,
    )
    path = tmp_path / "OUT" / "alignment.chr1_1-4.nex"
    text = path.read_text()
    assert text.startswith("#nexus")
    assert "popA^s1" in text
    assert "popB^s2" in text


@pytest.mark.parametrize(
    ("split", "description"),
    [
        (False, "as independent records in one PHYLIP file"),
        (True, "as separate PHYLIP files"),
    ],
)
def test_seqex_logger_describes_nonconcatenated_layout(
    tmp_path: Path,
    split: bool,
    description: str,
) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["ACGT", "ACGT"],
        [(0, 0, 4, 1, 4)],
    )
    messages: list[str] = []
    sink = logger.add(messages.append, format="{message}", level="INFO")
    try:
        _run(tmp_path, h5, split=split)
    finally:
        logger.remove(sink)
    assert any(
        f"wrote 1 filtered loci {description}" in message for message in messages
    )
    assert any("wrote stats report to:" in message for message in messages)


def test_seqex_validates_incompatible_and_dependent_options(tmp_path: Path) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["ACGT", "ACGT"],
        [(0, 0, 4, 1, 4)],
    )
    with pytest.raises(IPyradError, match="random-seed requires"):
        _run(tmp_path, h5, random_seed=1)
    with pytest.raises(IPyradError, match="append-population requires"):
        _run(tmp_path, h5, append_population=True)
    with pytest.raises(IPyradError, match="split cannot"):
        _run(tmp_path, h5, split=True, stdout=True)
    with pytest.raises(IPyradError, match="cores must"):
        _run(tmp_path, h5, cores=0)


@pytest.mark.parametrize(
    ("names", "imap_text", "message"),
    [
        (["s^1", "s2"], "s^1 popA\ns2 popB\n", "samples=s\\^1"),
        (["s1", "s2"], "s1 pop^A\ns2 popB\n", "populations=pop\\^A"),
    ],
)
def test_seqex_append_population_rejects_caret_delimiter_collisions(
    tmp_path: Path,
    names: list[str],
    imap_text: str,
    message: str,
) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["ACGT", "ACGT"],
        [(0, 0, 4, 1, 4)],
        names=names,
    )
    imap = tmp_path / "imap.tsv"
    imap.write_text(imap_text, encoding="utf-8")
    with pytest.raises(IPyradError, match=message):
        _run(tmp_path, h5, imap=imap, append_population=True)


def test_seqex_mixed_scaffold_and_coordinate_selectors_clip_independently(
    tmp_path: Path,
) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["AAAACCCC", "AAAACCCC"],
        [(0, 0, 4, 1, 4), (1, 4, 8, 1, 4)],
        scaffold_names=["chr1", "chr2"],
        scaffold_lengths=[4, 4],
    )

    loci = _run(tmp_path, h5, windows=["chr1", "chr2:2-3"])

    assert [locus.spec.label for locus in loci] == ["chr1:1-4", "chr2:2-3"]
    assert [locus.spec.clipped for locus in loci] == [False, True]
    stats_text = (tmp_path / "OUT" / "alignment.stats.txt").read_text()
    stats_json = json.loads((tmp_path / "OUT" / "alignment.stats.json").read_text())
    assert "selected_windows: chr1, chr2:2-3\n" in stats_text
    assert "selected_windows: chr1:1-4, chr2:2-3\n" not in stats_text
    assert stats_json["seqex_summary"]["selected_windows"] == [
        "chr1",
        "chr2:2-3",
    ]


def test_seqex_expands_imap_globs_and_applies_minmap(tmp_path: Path) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["ACGT", "NNNN", "ACGT"],
        [(0, 0, 4, 1, 4)],
        names=["barbeyi-01", "barbeyi-02", "geyeri-01"],
    )
    imap = tmp_path / "imap.tsv"
    imap.write_text("barbeyi*\tbarbeyi\ngeyeri*\tgeyeri\n", encoding="utf-8")
    minmap = tmp_path / "minmap.tsv"
    minmap.write_text("barbeyi\t1\ngeyeri\t1\n", encoding="utf-8")

    loci = _run(tmp_path, h5, imap=imap, minmap=minmap)

    assert loci[0].names == ["barbeyi-01", "barbeyi-02", "geyeri-01"]


def test_seqex_reference_selection_and_explicit_exclude_precedence(
    tmp_path: Path,
) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["ACGT", "ACGT", "ACGT"],
        [(0, 0, 4, 1, 4)],
        names=["s1", "assembly_reference_sequence", "s2"],
    )

    default = _run(tmp_path, h5, outdir=tmp_path / "DEFAULT")
    included = _run(
        tmp_path,
        h5,
        outdir=tmp_path / "INCLUDED",
        include_reference=True,
    )
    excluded = _run(
        tmp_path,
        h5,
        outdir=tmp_path / "EXCLUDED",
        include_reference=True,
        exclude=["assembly_reference_sequence"],
    )

    assert default[0].names == ["s1", "s2"]
    assert included[0].names == ["s1", "assembly_reference_sequence", "s2"]
    assert excluded[0].names == ["s1", "s2"]


def test_seqex_print_scaffold_table_and_rejects_overlapping_windows(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["ACGT", "ACGT"],
        [(0, 0, 4, 1, 4)],
    )

    result = _run(tmp_path, h5, print_scaffold_table=True)

    assert result is None
    assert "scaffold_name\tscaffold_length" in capsys.readouterr().out
    with pytest.raises(IPyradError, match="windows cannot overlap"):
        _run(tmp_path, h5, windows=["chr1:1-3", "chr1:3-4"])


def test_seqex_rejects_invalid_internal_clip_value(tmp_path: Path) -> None:
    h5 = _write_h5(
        tmp_path / "data.hdf5",
        ["ACGT", "ACGT"],
        [(0, 0, 4, 1, 4)],
    )

    with pytest.raises(IPyradError, match="clip must"):
        _run(tmp_path, h5, clip="yes")
