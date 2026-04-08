from pathlib import Path
import re

import h5py
import numpy as np
import pandas as pd
import pytest
from loguru import logger

from ipyrad2.analysis.extracters.locus_extracter import LocusExtracter
from ipyrad2.analysis.extracters.window_extracter import WindowExtracter
from ipyrad2.utils.exceptions import IPyradError


def _write_test_h5(
    path: Path,
    sequences: list[str],
    *,
    rows: list[tuple[int, int, int, int, int]],
    sample_names: list[str] | None = None,
    scaffold_length: int | None = None,
    scaffold_names: list[str] | None = None,
    scaffold_lengths: list[int] | None = None,
) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    phy = np.vstack([
        np.frombuffer(sequence.encode("utf-8"), dtype=np.uint8)
        for sequence in sequences
    ])
    if scaffold_names is None:
        scaffold_names = ["chr1"]
    if scaffold_lengths is None:
        if scaffold_length is None:
            raise ValueError("scaffold_length is required when scaffold_lengths is not provided")
        scaffold_lengths = [scaffold_length]
    if sample_names is None:
        sample_names = ["s1", "s2", "s3"]

    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["scaffold_names"] = np.array(scaffold_names, dtype=string_dtype)
        io5.attrs["scaffold_lengths"] = np.array(scaffold_lengths, dtype=np.uint64)
        io5.attrs["names"] = np.array(sample_names, dtype=string_dtype)
        io5.create_dataset("phy", data=phy)
        phymap = io5.create_dataset("phymap", data=np.array(rows, dtype=np.uint64))
        phymap.attrs["columns"] = np.array(
            ["scaff", "phy0", "phy1", "pos0", "pos1"],
            dtype=string_dtype,
        )

    return path


def _make_wex(
    tmp_path: Path,
    h5: Path,
    *,
    windows: list[str] | None = None,
    out_format: str = "phy",
    stdout: bool = False,
    force: bool = False,
    exclude=None,
    include_reference: bool = False,
    imap=None,
    minmap=None,
    max_sample_missing: float = 1.0,
) -> WindowExtracter:
    return WindowExtracter(
        data=h5,
        name="alignment",
        outdir=tmp_path / "OUT",
        out_format=out_format,
        windows=windows,
        min_sample_coverage=1,
        max_sample_missing=max_sample_missing,
        exclude=exclude,
        include_reference=include_reference,
        imap=imap,
        minmap=minmap,
        stdout=stdout,
        force=force,
    )


def _make_lex(
    tmp_path: Path,
    h5: Path,
    *,
    out_format: str = "bpp",
    nloci: int = 1,
    min_length: int = 4,
    stdout: bool = False,
    force: bool = False,
    exclude=None,
    include_reference: bool = False,
    imap=None,
    minmap=None,
) -> LocusExtracter:
    return LocusExtracter(
        data=h5,
        name="alignment",
        outdir=tmp_path / "LEX",
        out_format=out_format,
        nloci=nloci,
        min_length=min_length,
        windows=["chr1"],
        min_sample_coverage=1,
        max_sample_missing=1.0,
        exclude=exclude,
        include_reference=include_reference,
        imap=imap,
        minmap=minmap,
        stdout=stdout,
        force=force,
    )


def _parse_lex_stats(stats_path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    header: dict[str, str] = {}
    rows: list[dict[str, str]] = []
    columns: list[str] | None = None
    in_summary = False
    in_table = False
    for line in stats_path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        if line == "Summary":
            in_summary = True
            in_table = False
            continue
        if line == "Accepted loci":
            in_summary = False
            in_table = True
            columns = None
            continue
        if set(line) == {"-"}:
            continue
        if in_summary:
            key, value = line.split(": ", 1)
            header[key] = value
            continue
        if in_table:
            parts = re.split(r" {2,}", line.strip())
            if columns is None:
                columns = parts
                continue
            rows.append(dict(zip(columns, parts, strict=True)))
    return header, rows


def test_wex_phy_output_writes_stats_and_keeps_filtered_rows_aligned(
    tmp_path: Path,
) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "NNNNTCGT", "TCGTACGA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        scaffold_length=8,
    )
    tool = _make_wex(
        tmp_path,
        h5,
        windows=["chr1:1-4"],
        force=True,
        max_sample_missing=0.5,
    )

    tool._write_to_phy()

    phy = (tmp_path / "OUT" / "alignment.phy").read_text(encoding="utf-8")
    stats = (tmp_path / "OUT" / "alignment.stats.tsv").read_text(encoding="utf-8")

    lines = phy.splitlines()
    assert lines[0] == "2 4"
    assert lines[1].strip().endswith("ACGT")
    assert lines[2].strip().endswith("TCGT")
    assert "nvariants_in_windows_after_filtering" in stats
    assert "nvariants_in_windows_afater_filtering" not in stats


def test_wex_phy_stdout_writes_stats_without_crashing(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "NNNNTCGT", "TCGTACGA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        scaffold_length=8,
    )
    tool = _make_wex(
        tmp_path,
        h5,
        windows=["chr1:1-4"],
        stdout=True,
        force=True,
    )

    tool._write_to_phy()

    captured = capsys.readouterr()
    stats = (tmp_path / "OUT" / "alignment.stats.tsv").read_text(encoding="utf-8")

    assert captured.out.startswith("3 4\n")
    assert "outfile\tSTDOUT" in stats


def test_wex_fasta_output_and_force_handling(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "NNNNTCGT", "TCGTACGA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        scaffold_length=8,
    )
    tool = _make_wex(
        tmp_path,
        h5,
        windows=["chr1:1-4"],
        out_format="fa",
        force=False,
    )

    tool._write_to_fa()

    fasta = (tmp_path / "OUT" / "alignment.fa").read_text(encoding="utf-8")
    assert fasta == ">s1\nACGT\n>s2\nNNNN\n>s3\nTCGT\n"

    with pytest.raises(IPyradError, match="already exists"):
        tool._write_to_fa()

    forced = _make_wex(
        tmp_path,
        h5,
        windows=["chr1:1-4"],
        out_format="fa",
        force=True,
    )
    forced._write_to_fa()


def test_wex_bed_windows_match_cli_region_selection(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "NNNNTCGT", "TCGTACGA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        scaffold_length=8,
    )
    bed = tmp_path / "windows.bed"
    bed.write_text("chr1\t0\t4\n", encoding="utf-8")

    cli_tool = _make_wex(tmp_path, h5, windows=["chr1:1-4"], force=True)
    bed_tool = _make_wex(tmp_path, h5, windows=[str(bed)], force=True)

    cli_names, cli_seqs = cli_tool._run()
    bed_names, bed_seqs = bed_tool._run()

    assert cli_names == bed_names
    assert np.array_equal(cli_seqs, bed_seqs)
    assert bed_tool.selected_windows == ["chr1:1-4"]


def test_wex_without_windows_selects_all_scaffolds_and_logs_hint(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAACCCC", "TTTTGGGG", "CCCCAAAA"],
        rows=[(0, 0, 4, 1, 4), (1, 4, 8, 1, 4)],
        scaffold_names=["chr1", "chr2"],
        scaffold_lengths=[4, 4],
    )
    tool = _make_wex(tmp_path, h5, force=True)

    messages: list[str] = []
    sink_id = logger.add(messages.append, format="{message}")
    try:
        names, seqs = tool._run()
    finally:
        logger.remove(sink_id)

    assert names == ["s1", "s2", "s3"]
    assert seqs.shape == (3, 8)
    assert tool.selected_windows == ["chr1:1-4", "chr2:1-4"]
    assert seqs[0].tobytes().decode("utf-8") == "AAAACCCC"
    assert any(
        "No windows specified; selecting the full length of all scaffolds." in str(msg)
        and "-w" in str(msg)
        and "-P" in str(msg)
        for msg in messages
    )


def test_wex_accepts_imap_and_minmap_files(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "NNNNTCGT", "TCGTACGA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        scaffold_length=8,
    )
    imap = tmp_path / "imap.tsv"
    imap.write_text("s1\tpop1\ns3\tpop2\n", encoding="utf-8")
    minmap = tmp_path / "minmap.tsv"
    minmap.write_text("pop1\t1\npop2\t1\n", encoding="utf-8")

    tool = _make_wex(
        tmp_path,
        h5,
        windows=["chr1:1-4"],
        force=True,
        imap=imap,
        minmap=minmap,
    )

    names, seqs = tool._run()

    assert tool.imap == {"pop1": ["s1"], "pop2": ["s3"]}
    assert tool.minmap == {"pop1": 1, "pop2": 1}
    assert names == ["s1", "s3"]
    assert seqs.shape == (2, 4)


def test_wex_excludes_reference_by_default(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "NNNNTCGT", "TCGTACGA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        sample_names=["s1", "assembly_reference_sequence", "s3"],
        scaffold_length=8,
    )
    tool = _make_wex(tmp_path, h5, windows=["chr1:1-4"], force=True)

    names, seqs = tool._run()

    assert names == ["s1", "s3"]
    assert seqs.shape == (2, 4)


def test_lex_excludes_reference_by_default(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "NNNNTCGT", "TCGTACGA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        sample_names=["s1", "assembly_reference_sequence", "s3"],
        scaffold_length=8,
    )
    lex = _make_lex(tmp_path, h5, out_format="phy", force=True)

    assert lex.wex.snames == ["s1", "s3"]


def test_wex_include_reference_flag_restores_reference_without_imap(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "NNNNTCGT", "TCGTACGA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        sample_names=["s1", "assembly_reference_sequence", "s3"],
        scaffold_length=8,
    )
    tool = _make_wex(
        tmp_path,
        h5,
        windows=["chr1:1-4"],
        include_reference=True,
        force=True,
    )

    names, seqs = tool._run()

    assert names == ["s1", "assembly_reference_sequence", "s3"]
    assert seqs.shape == (3, 4)


def test_lex_imap_membership_includes_reference_without_flag(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "NNNNTCGT", "TCGTACGA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        sample_names=["s1", "assembly_reference_sequence", "s3"],
        scaffold_length=8,
    )
    imap = {"pop1": ["s1", "assembly_reference_sequence"], "pop2": ["s3"]}
    minmap = {"pop1": 1, "pop2": 1}

    lex = _make_lex(tmp_path, h5, out_format="phy", imap=imap, minmap=minmap, force=True)

    assert lex.wex.snames == ["s1", "assembly_reference_sequence", "s3"]
    assert lex.wex.imap == imap


def test_wex_explicit_exclude_wins_over_include_reference_and_imap(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "NNNNTCGT", "TCGTACGA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        sample_names=["s1", "assembly_reference_sequence", "s3"],
        scaffold_length=8,
    )
    imap = {"pop1": ["s1", "assembly_reference_sequence"], "pop2": ["s3"]}
    minmap = {"pop1": 1, "pop2": 1}
    tool = _make_wex(
        tmp_path,
        h5,
        windows=["chr1:1-4"],
        exclude=["assembly_reference_sequence"],
        include_reference=True,
        imap=imap,
        minmap=minmap,
        force=True,
    )

    names, _ = tool._run()

    assert names == ["s1", "s3"]
    assert tool.imap == {"pop1": ["s1"], "pop2": ["s3"]}


def test_lex_include_reference_with_imap_requires_reference_group_assignment(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "NNNNTCGT", "TCGTACGA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        sample_names=["s1", "assembly_reference_sequence", "s3"],
        scaffold_length=8,
    )
    imap = {"pop1": ["s1"], "pop2": ["s3"]}
    minmap = {"pop1": 1, "pop2": 1}

    with pytest.raises(
        IPyradError,
        match="assembly_reference_sequence was requested with -R, but it must also be assigned to an IMAP group.",
    ):
        _make_lex(
            tmp_path,
            h5,
            out_format="phy",
            include_reference=True,
            imap=imap,
            minmap=minmap,
            force=True,
        )


def test_lex_bpp_output_still_works_and_respects_force(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "ACGTTCGT", "TCGTACGA"],
        rows=[(0, 0, 8, 1, 8)],
        scaffold_length=8,
    )
    lex = _make_lex(tmp_path, h5, force=False)
    lex.loci = pd.DataFrame([{"chrom": "0", "startpos": "0", "endpos": "4"}])

    lex._write_loci()

    outfile = tmp_path / "LEX" / "alignment.phy"
    stats_path = tmp_path / "LEX" / "alignment.stats.txt"
    assert outfile.exists()
    assert stats_path.exists()
    assert outfile.read_text(encoding="utf-8").startswith("3 4")
    assert sorted(path.name for path in (tmp_path / "LEX").glob("*.stats.txt")) == [
        "alignment.stats.txt"
    ]
    stats_text = stats_path.read_text(encoding="utf-8")
    assert stats_text.startswith("Summary\n-------\n")
    assert "\nAccepted loci\n-------------\n" in stats_text
    header, rows = _parse_lex_stats(stats_path)
    assert header["tool"] == "lex"
    assert header["nloci_written"] == "1"
    assert len(rows) == 1
    assert rows[0]["outfile"].endswith("alignment.phy")

    with pytest.raises(IPyradError, match="already exists"):
        lex._write_loci()

    forced = _make_lex(tmp_path, h5, force=True)
    forced.loci = pd.DataFrame([{"chrom": "0", "startpos": "0", "endpos": "4"}])
    forced._write_loci()


def test_lex_write_loci_uses_phymap_rows_not_window_index(monkeypatch, tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAACCCC", "TTTTGGGG", "CCCCAAAA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        scaffold_length=8,
    )
    lex = LocusExtracter(
        data=h5,
        name="alignment",
        outdir=tmp_path / "LEX",
        out_format="phy",
        nloci=1,
        min_length=4,
        windows=["chr1"],
        min_sample_coverage=1,
        max_sample_missing=1.0,
        exclude=None,
        imap=None,
        minmap=None,
        stdout=True,
        force=True,
    )
    lex.loci = pd.DataFrame([{"chrom": "1", "startpos": "0", "endpos": "4"}])

    windows: list[str] = []

    def _fake_write_to_phy() -> None:
        windows.append(lex.wex.windows[0])

    def _fake_stats() -> dict[str, object]:
        return {
            "outfile": "STDOUT",
            "nsamples_before_filtering": 3,
            "nsites_in_windows_before_filtering": 4,
            "nvariants_in_windows_before_filtering": 1,
            "nsamples_after_filtering": 3,
            "nsites_in_windows_after_filtering": 4,
            "nvariants_in_windows_after_filtering": 1,
        }

    def _fake_write_to_phy_with_stats(*, write_stats=True, return_alignment=False, return_stats=False, **kwargs):
        _fake_write_to_phy()
        alignment = "3 4\ns1    AAAA\ns2    AAAA\ns3    AAAA\n"
        if return_alignment and return_stats:
            return alignment, _fake_stats()
        if return_alignment:
            return alignment
        if return_stats:
            return _fake_stats()
        return None

    monkeypatch.setattr(lex.wex, "_write_to_phy", _fake_write_to_phy_with_stats)

    lex._write_loci()

    assert windows == ["chr1:5-8"]


def test_lex_get_loci_respects_selected_windows_and_min_length(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAACCCC", "TTTTGGGG", "CCCCAAAA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        scaffold_length=8,
    )
    lex = LocusExtracter(
        data=h5,
        name="alignment",
        outdir=tmp_path / "LEX",
        out_format="phy",
        nloci=1,
        min_length=4,
        windows=["chr1:5-8"],
        min_sample_coverage=1,
        max_sample_missing=1.0,
        exclude=None,
        imap=None,
        minmap=None,
        stdout=False,
        force=True,
    )

    lex._get_loci()

    assert lex.eligible_loci_before_filtering == 1
    assert lex.loci.to_dict("records") == [
        {"chrom": 1, "startpos": 0, "endpos": 4, "raw_length": 4}
    ]


def test_lex_get_loci_raises_when_no_raw_loci_meet_min_length(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAACCCC", "TTTTGGGG", "CCCCAAAA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        scaffold_length=8,
    )
    lex = _make_lex(tmp_path, h5, out_format="phy", nloci=1, min_length=5, force=True)

    with pytest.raises(IPyradError, match="No loci met the minimum length requirement before filtering"):
        lex._get_loci()


@pytest.mark.parametrize(
    ("out_format", "suffix"),
    [("phy", "phy"), ("nex", "nex")],
)
def test_lex_writes_single_summary_stats_file_for_sequence_outputs(
    tmp_path: Path,
    out_format: str,
    suffix: str,
) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAACCCC", "TTTTGGGG", "CCCCAAAA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        scaffold_length=8,
    )
    lex = _make_lex(tmp_path, h5, out_format=out_format, nloci=2, force=True)
    lex.loci = pd.DataFrame(
        [
            {"chrom": "0", "startpos": "0", "endpos": "4"},
            {"chrom": "1", "startpos": "0", "endpos": "4"},
        ]
    )

    lex._write_loci()

    stats_path = tmp_path / "LEX" / "alignment.stats.txt"
    assert stats_path.exists()
    assert sorted(path.name for path in (tmp_path / "LEX").glob("*.stats.txt")) == [
        "alignment.stats.txt"
    ]
    assert (tmp_path / "LEX" / f"chr1:1-4.{suffix}").exists()
    assert (tmp_path / "LEX" / f"chr1:5-8.{suffix}").exists()

    stats_text = stats_path.read_text(encoding="utf-8")
    assert stats_text.startswith("Summary\n-------\n")
    assert "\nAccepted loci\n-------------\n" in stats_text
    header, rows = _parse_lex_stats(stats_path)
    assert header["tool"] == "lex"
    assert header["out_format"] == out_format
    assert header["nloci_requested"] == "2"
    assert header["nloci_written"] == "2"
    assert header["min_length_requested"] == "4"
    assert header["eligible_loci_before_filtering"] == "2"
    assert header["loci_rejected_after_filtering"] == "0"
    assert len(rows) == 2
    assert rows[0]["locus_name"] == "chr1:1-4"
    assert rows[0]["start"] == "1"
    assert rows[0]["end"] == "4"
    assert rows[0]["outfile"].endswith(f"chr1:1-4.{suffix}")
    assert rows[1]["outfile"].endswith(f"chr1:5-8.{suffix}")


def test_lex_stdout_writes_single_summary_stats_file(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAACCCC", "TTTTGGGG", "CCCCAAAA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        scaffold_length=8,
    )
    lex = _make_lex(tmp_path, h5, out_format="phy", nloci=2, stdout=True, force=True)
    lex.loci = pd.DataFrame(
        [
            {"chrom": "0", "startpos": "0", "endpos": "4"},
            {"chrom": "1", "startpos": "0", "endpos": "4"},
        ]
    )

    lex._write_loci()

    captured = capsys.readouterr()
    stats_path = tmp_path / "LEX" / "alignment.stats.txt"
    assert captured.out.startswith("3 4\n")
    assert stats_path.exists()
    assert sorted(path.name for path in (tmp_path / "LEX").glob("*.stats.txt")) == [
        "alignment.stats.txt"
    ]
    header, rows = _parse_lex_stats(stats_path)
    assert header["nloci_written"] == "2"
    assert [row["outfile"] for row in rows] == ["STDOUT", "STDOUT"]


def test_lex_bpp_stdout_does_not_duplicate_locus_output(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAACCCC", "TTTTGGGG", "CCCCAAAA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        scaffold_length=8,
    )
    lex = _make_lex(tmp_path, h5, out_format="bpp", nloci=2, stdout=True, force=True)
    lex.loci = pd.DataFrame(
        [
            {"chrom": "0", "startpos": "0", "endpos": "4"},
            {"chrom": "1", "startpos": "0", "endpos": "4"},
        ]
    )

    lex._write_loci()

    captured = capsys.readouterr()
    stats_path = tmp_path / "LEX" / "alignment.stats.txt"
    assert captured.out.count("3 4") == 2
    assert stats_path.exists()
    header, rows = _parse_lex_stats(stats_path)
    assert header["out_format"] == "bpp"
    assert [row["outfile"] for row in rows] == ["STDOUT", "STDOUT"]


def test_lex_warns_and_writes_fewer_loci_when_post_filter_length_is_too_short(
    tmp_path: Path,
) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAACCCC", "TTTTGGGG", "CCCCAAAA"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8), (0, 8, 12, 9, 12)],
        scaffold_length=12,
    )
    lex = _make_lex(tmp_path, h5, out_format="phy", nloci=3, min_length=4, force=True)
    lex.loci = pd.DataFrame(
        [
            {"chrom": 0, "startpos": 0, "endpos": 4, "raw_length": 4},
            {"chrom": 1, "startpos": 0, "endpos": 4, "raw_length": 4},
            {"chrom": 2, "startpos": 0, "endpos": 4, "raw_length": 4},
        ]
    )
    lex.eligible_loci_before_filtering = 3

    stats_queue = iter(
        [
            {
                "outfile": tmp_path / "LEX" / "chr1:1-4.phy",
                "nsamples_before_filtering": 3,
                "nsites_in_windows_before_filtering": 4,
                "nvariants_in_windows_before_filtering": 1,
                "nsamples_after_filtering": 3,
                "nsites_in_windows_after_filtering": 3,
                "nvariants_in_windows_after_filtering": 1,
            },
            {
                "outfile": tmp_path / "LEX" / "chr1:5-8.phy",
                "nsamples_before_filtering": 3,
                "nsites_in_windows_before_filtering": 4,
                "nvariants_in_windows_before_filtering": 1,
                "nsamples_after_filtering": 3,
                "nsites_in_windows_after_filtering": 4,
                "nvariants_in_windows_after_filtering": 1,
            },
            {
                "outfile": tmp_path / "LEX" / "chr1:9-12.phy",
                "nsamples_before_filtering": 3,
                "nsites_in_windows_before_filtering": 4,
                "nvariants_in_windows_before_filtering": 1,
                "nsamples_after_filtering": 3,
                "nsites_in_windows_after_filtering": 2,
                "nvariants_in_windows_after_filtering": 1,
            },
        ]
    )

    def _fake_write_to_phy(*, return_alignment=False, return_stats=False, **kwargs):
        stats = next(stats_queue)
        alignment = "3 4\ns1    AAAA\ns2    AAAA\ns3    AAAA\n"
        if return_alignment and return_stats:
            return alignment, stats
        if return_alignment:
            return alignment
        if return_stats:
            return stats
        return None

    messages: list[str] = []
    sink_id = logger.add(messages.append, format="{message}", level="WARNING")
    try:
        monkeypatch = pytest.MonkeyPatch()
        monkeypatch.setattr(lex.wex, "_write_to_phy", _fake_write_to_phy)
        lex._write_loci()
    finally:
        logger.remove(sink_id)
        monkeypatch.undo()

    stats_path = tmp_path / "LEX" / "alignment.stats.txt"
    header, rows = _parse_lex_stats(stats_path)
    assert header["nloci_requested"] == "3"
    assert header["nloci_written"] == "1"
    assert header["eligible_loci_before_filtering"] == "3"
    assert header["loci_rejected_after_filtering"] == "2"
    assert len(rows) == 1
    assert rows[0]["locus_name"] == "chr1:5-8"
    assert any("only 1 met the minimum length requirement" in str(msg) for msg in messages)
