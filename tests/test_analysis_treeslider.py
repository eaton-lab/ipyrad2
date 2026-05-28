from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

import ipyrad2.analysis.methods.treeslider as treeslider_mod


def _write_sequence_h5(
    path: Path,
    sequences: list[str],
    *,
    rows: list[tuple[int, int, int, int, int]],
    scaffold_names: list[str] | None = None,
    scaffold_lengths: list[int] | None = None,
    sample_names: list[str] | None = None,
) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    if scaffold_names is None:
        scaffold_names = ["chr1"]
    if scaffold_lengths is None:
        scaffold_lengths = [len(sequences[0])]
    if sample_names is None:
        sample_names = ["s1", "s2", "s3"]
    phy = np.vstack(
        [np.frombuffer(sequence.encode("utf-8"), dtype=np.uint8) for sequence in sequences]
    )

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


def _mock_raxml_run_success(cmd, cwd, capture_output, text, check):
    workdir = Path(cwd)
    prefix = cmd[cmd.index("--prefix") + 1]
    msa = Path(cmd[cmd.index("--msa") + 1])
    names = [
        line[1:].strip()
        for line in msa.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]
    tree = f"({','.join(names)});"
    (workdir / f"{prefix}.raxml.bestTree").write_text(tree, encoding="utf-8")
    if "--all" in cmd:
        (workdir / f"{prefix}.raxml.support").write_text(tree, encoding="utf-8")
    return treeslider_mod.subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")


def test_run_treeslider_locus_mode_writes_manifest_trees_and_cleans_stage(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    h5 = _write_sequence_h5(
        tmp_path / "assembly.hdf5",
        [
            "AAAACCCCGGGGTTTT",
            "AATACCCCNNNNTTTN",
            "AACACCCCNNNNNNNN",
        ],
        rows=[
            (0, 0, 4, 1, 4),
            (0, 4, 8, 5, 8),
            (0, 8, 12, 9, 12),
            (0, 12, 16, 13, 16),
        ],
    )
    monkeypatch.setattr(treeslider_mod, "_resolve_binary", lambda _binary: "/usr/bin/raxml-ng")
    monkeypatch.setattr(treeslider_mod.subprocess, "run", _mock_raxml_run_success)

    treeslider_mod.run_treeslider_method(
        data=h5,
        name="slider",
        outdir=tmp_path / "OUT",
        window_size=None,
        slide_size=None,
        print_scaffold_table=False,
        scaffolds=None,
        min_sample_coverage=2,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        min_sample_alignment_length=3,
        min_alignment_length=1,
        threads="auto",
        workers="auto",
        bs_trees=0,
        model="GTR+G",
        raxml_ng_binary=None,
        seed=11,
        force=True,
        redo=False,
        log_level="INFO",
    )

    manifest = pd.read_csv(tmp_path / "OUT" / "slider.stats.tsv", sep="\t")
    stats = (tmp_path / "OUT" / "slider.stats.txt").read_text(encoding="utf-8")
    trees = (tmp_path / "OUT" / "slider.trees.nex").read_text(encoding="utf-8")

    assert manifest["status"].tolist() == [
        "tree_completed",
        "polytomy_written",
        "skipped_no_data",
        "skipped_few_samples",
    ]
    assert manifest.loc[0, "tree_source"] == "raxml-ng"
    assert manifest.loc[1, "tree_source"] == "polytomy"
    assert manifest.loc[2, "status_detail"].startswith("No sites remained")
    assert manifest.loc[3, "samples_dropped_by_min_sample_alignment_length"] == "s3"
    assert "tree_completed: 1" in stats
    assert "polytomy_written: 1" in stats
    assert "skipped_no_data: 1" in stats
    assert "skipped_few_samples: 1" in stats
    assert "Tree window_000001" in trees
    assert "Tree window_000002" in trees
    assert "window_000003" not in trees
    assert "window_000004" not in trees
    assert not (tmp_path / "OUT" / ".slider.stage").exists()


def test_run_treeslider_redo_retries_failed_windows_only(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    h5 = _write_sequence_h5(
        tmp_path / "assembly.hdf5",
        [
            "AAAAGGGG",
            "AATAGGGA",
            "AACAGGGT",
        ],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
    )
    monkeypatch.setattr(treeslider_mod, "_resolve_binary", lambda _binary: "/usr/bin/raxml-ng")

    calls: list[str] = []
    fail_once = {"window_000001"}

    def _mock_run(cmd, cwd, capture_output, text, check):
        prefix = cmd[cmd.index("--prefix") + 1]
        calls.append(prefix)
        workdir = Path(cwd)
        msa = Path(cmd[cmd.index("--msa") + 1])
        names = [
            line[1:].strip()
            for line in msa.read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        ]
        if prefix in fail_once:
            fail_once.remove(prefix)
            return treeslider_mod.subprocess.CompletedProcess(
                cmd,
                1,
                stdout="",
                stderr="resource failure",
            )
        tree = f"({','.join(names)});"
        (workdir / f"{prefix}.raxml.bestTree").write_text(tree, encoding="utf-8")
        return treeslider_mod.subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(treeslider_mod.subprocess, "run", _mock_run)

    common_kwargs = dict(
        data=h5,
        name="redo",
        outdir=tmp_path / "OUT",
        window_size=None,
        slide_size=None,
        print_scaffold_table=False,
        scaffolds=None,
        min_sample_coverage=2,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        min_sample_alignment_length=1,
        min_alignment_length=1,
        threads=1,
        workers=1,
        bs_trees=0,
        model="GTR+G",
        raxml_ng_binary=None,
        seed=5,
        force=True,
        redo=False,
        log_level="INFO",
    )

    treeslider_mod.run_treeslider_method(**common_kwargs)
    manifest = pd.read_csv(tmp_path / "OUT" / "redo.stats.tsv", sep="\t")
    assert manifest["status"].tolist() == ["tree_failed", "tree_completed"]
    assert calls == ["window_000001", "window_000002"]

    calls.clear()
    common_kwargs["force"] = False
    common_kwargs["redo"] = True
    treeslider_mod.run_treeslider_method(**common_kwargs)

    manifest = pd.read_csv(tmp_path / "OUT" / "redo.stats.tsv", sep="\t")
    trees = (tmp_path / "OUT" / "redo.trees.nex").read_text(encoding="utf-8")
    assert manifest["status"].tolist() == ["tree_completed", "tree_completed"]
    assert calls == ["window_000001"]
    assert "Tree window_000001" in trees
    assert "Tree window_000002" in trees


def test_run_treeslider_prints_scaffold_table_and_returns(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    h5 = _write_sequence_h5(
        tmp_path / "assembly.hdf5",
        ["AAAATTTT", "AAAATTTT", "AAAATTTT"],
        rows=[(0, 0, 4, 1, 4), (1, 4, 8, 1, 4)],
        scaffold_names=["chr1", "chr2"],
        scaffold_lengths=[4, 4],
    )

    treeslider_mod.run_treeslider_method(
        data=h5,
        name="print",
        outdir=tmp_path / "OUT",
        window_size=None,
        slide_size=None,
        print_scaffold_table=True,
        scaffolds=None,
        min_sample_coverage=2,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        min_sample_alignment_length=1,
        min_alignment_length=1,
        threads="auto",
        workers="auto",
        bs_trees=0,
        model="GTR+G",
        raxml_ng_binary=None,
        seed=None,
        force=False,
        redo=False,
        log_level="INFO",
    )

    captured = capsys.readouterr()
    assert "scaffold_name\tscaffold_length" in captured.out
    assert "chr1\t4" in captured.out
    assert not (tmp_path / "OUT" / "print.stats.tsv").exists()
