import json
from pathlib import Path

import h5py
import numpy as np
import pytest
from loguru import logger

from ipyrad2.analysis.extracters.seqex import run_seqex
from ipyrad2.cli.cli_main import setup_parsers
from ipyrad2.utils.exceptions import IPyradError


def _write_h5(
    path: Path,
    sequences: list[str],
    rows: list[tuple[int, int, int, int, int]],
    names: list[str] | None = None,
) -> Path:
    dtype = h5py.string_dtype(encoding="utf-8")
    phy = np.vstack(
        [np.frombuffer(sequence.encode(), dtype=np.uint8) for sequence in sequences]
    )
    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["scaffold_names"] = np.array(["chr1"], dtype=dtype)
        io5.attrs["scaffold_lengths"] = np.array([len(sequences[0])])
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
    assert "\t" not in stats
    assert "sample  population  written_final" in stats
    assert "locus_index  locus" in stats
    stats_json = json.loads((tmp_path / "OUT" / "alignment.stats.json").read_text())
    assert stats_json["seqex_summary"]["output_layout"] == "multi-locus"
    assert stats_json["seqex_summary"]["max_loci"] is None
    assert stats_json["output_summary"]["max_samples"] == 2
    assert stats_json["output_summary"]["mean_samples"] == 2.0
    assert stats_json["sample_occupancy"][0]["written_final"] is True
    assert stats_json["written_loci"][0] == {
        "locus_index": 1,
        "locus": "chr1:1-4",
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
        max_loci=1,
        random_seed=9,
        cores=1,
    )
    parallel = _run(
        tmp_path,
        h5,
        outdir=tmp_path / "PARALLEL",
        max_loci=1,
        random_seed=9,
        cores=2,
    )
    assert [locus.spec.label for locus in serial] == [
        locus.spec.label for locus in parallel
    ]
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
    assert "s1^popA" in text
    assert "s2^popB" in text


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
