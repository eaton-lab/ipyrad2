from __future__ import annotations

import multiprocessing as mp
import os
from pathlib import Path
import signal
import time

import h5py
import numpy as np
import pandas as pd
import pytest
from loguru import logger

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


def _mock_raxml_pipeline_success(cmds, **_kwargs):
    cmd = cmds[0]
    prefix = Path(cmd[cmd.index("--prefix") + 1])
    msa = Path(cmd[cmd.index("--msa") + 1])
    names = [
        line[1:].strip()
        for line in msa.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]
    tree = f"({','.join(names)});"
    prefix.parent.mkdir(parents=True, exist_ok=True)
    (prefix.parent / f"{prefix.name}.raxml.bestTree").write_text(tree, encoding="utf-8")
    if "--all" in cmd:
        (prefix.parent / f"{prefix.name}.raxml.support").write_text(tree, encoding="utf-8")
    return 0, b"", b""


def _run_pool_iter_sequential(job_iter, _log_level, **_kwargs):
    for key, (func, kwargs) in job_iter:
        yield key, func(**kwargs)


def _assert_interrupted_exit(code: int) -> None:
    assert code == -signal.SIGINT or code >= 128


def _pid_is_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _run_treeslider_until_interrupted(
    data: Path,
    outdir: Path,
    raxml_binary: Path,
    pid_file: Path,
) -> None:
    os.environ["IPYRAD2_TEST_TREESLIDER_PID_FILE"] = str(pid_file)
    treeslider_mod.run_treeslider_method(
        data=data,
        name="interrupt",
        outdir=outdir,
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
        workers=2,
        bs_trees=0,
        model="GTR+G",
        raxml_ng_binary=raxml_binary,
        seed=17,
        force=True,
        redo=False,
        log_level="WARNING",
    )


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
    monkeypatch.setattr(treeslider_mod, "run_pipeline", _mock_raxml_pipeline_success)
    monkeypatch.setattr(treeslider_mod, "run_with_pool_iter", _run_pool_iter_sequential)

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

    def _mock_run_pipeline(cmds, **_kwargs):
        cmd = cmds[0]
        prefix = Path(cmd[cmd.index("--prefix") + 1])
        calls.append(prefix.name)
        msa = Path(cmd[cmd.index("--msa") + 1])
        names = [
            line[1:].strip()
            for line in msa.read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        ]
        if prefix.name in fail_once:
            fail_once.remove(prefix.name)
            raise RuntimeError("resource failure")
        tree = f"({','.join(names)});"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        (prefix.parent / f"{prefix.name}.raxml.bestTree").write_text(tree, encoding="utf-8")
        return 0, b"", b""

    monkeypatch.setattr(treeslider_mod, "run_pipeline", _mock_run_pipeline)
    monkeypatch.setattr(treeslider_mod, "run_with_pool_iter", _run_pool_iter_sequential)

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


def test_run_treeslider_reports_filter_and_tree_progress(
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
    monkeypatch.setattr(treeslider_mod, "run_pipeline", _mock_raxml_pipeline_success)
    monkeypatch.setattr(treeslider_mod, "run_with_pool_iter", _run_pool_iter_sequential)

    events: list[tuple[str, int, int, str] | tuple[str, int, int] | tuple[str, int]] = []

    class _ProgressStub:
        def __init__(self, njobs, start=None, message="") -> None:
            self.njobs = njobs
            self.finished = 0
            self.message = message
            events.append(("init", njobs, self.finished, message))

        def update(self) -> None:
            events.append(("update", self.njobs, self.finished, self.message))

        def close(self) -> None:
            events.append(("close", self.njobs, self.finished))

    monkeypatch.setattr(treeslider_mod, "ProgressBar", _ProgressStub)
    messages: list[str] = []
    sink_id = logger.add(messages.append, format="{message}")
    try:
        treeslider_mod.run_treeslider_method(
            data=h5,
            name="progress",
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
    finally:
        logger.remove(sink_id)

    assert ("init", 4, 0, "Filtering windows - total jobs: 4") in events
    assert ("update", 4, 4, "Filtering windows - total jobs: 4") in events
    assert ("close", 4, 4) in events
    assert ("init", 2, 0, "Inferring trees - total jobs: 2") in events
    assert ("update", 2, 2, "Inferring trees - total jobs: 2") in events
    assert ("close", 2, 2) in events
    assert any("filtering windows and writing alignment files" in message for message in messages)
    assert any("inferring trees for accepted windows" in message for message in messages)


def test_run_treeslider_skips_implicit_binary_lookup_when_no_tree_jobs_remain(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    h5 = _write_sequence_h5(
        tmp_path / "assembly.hdf5",
        ["AAAA", "AAAA", "AAAA"],
        rows=[(0, 0, 4, 1, 4)],
    )

    def _unexpected_resolve(_binary):
        raise AssertionError("implicit binary lookup should not be called")

    monkeypatch.setattr(treeslider_mod, "_resolve_binary", _unexpected_resolve)

    treeslider_mod.run_treeslider_method(
        data=h5,
        name="no-raxml-needed",
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
        threads="auto",
        workers="auto",
        bs_trees=0,
        model="GTR+G",
        raxml_ng_binary=None,
        seed=None,
        force=True,
        redo=False,
        log_level="INFO",
    )

    manifest = pd.read_csv(tmp_path / "OUT" / "no-raxml-needed.stats.tsv", sep="\t")
    stats = (tmp_path / "OUT" / "no-raxml-needed.stats.txt").read_text(encoding="utf-8")
    assert manifest["status"].tolist() == ["polytomy_written"]
    assert "raxml_ng_binary: not_used" in stats


def test_run_treeslider_terminal_resume_skips_implicit_binary_lookup(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    h5 = _write_sequence_h5(
        tmp_path / "assembly.hdf5",
        ["AAAA", "AAAA", "AAAA"],
        rows=[(0, 0, 4, 1, 4)],
    )
    outdir = tmp_path / "OUT"
    name = "terminal"
    outdir.mkdir(parents=True, exist_ok=True)
    phymap = treeslider_mod._load_phymap(h5)
    specs = treeslider_mod._plan_locus_windows(phymap, ["chr1"], ["chr1"])
    manifest = treeslider_mod._initialize_manifest(specs)
    manifest.at[1, "status"] = "tree_completed"
    manifest.at[1, "status_detail"] = "Tree inference completed."
    manifest.at[1, "nsamples_before_filtering"] = 3
    manifest.at[1, "nsites_before_filtering"] = 4
    manifest.at[1, "nvariants_before_filtering"] = 1
    manifest.at[1, "nsites_after_site_filter"] = 4
    manifest.at[1, "nsamples_after_sample_length_filter"] = 3
    manifest.at[1, "nsites_after_sample_length_filter"] = 4
    manifest.at[1, "nsamples_after_filtering"] = 3
    manifest.at[1, "nsites_after_filtering"] = 4
    manifest.at[1, "nvariants_after_filtering"] = 1
    manifest.at[1, "retained_sample_names"] = "s1,s2,s3"
    manifest.at[1, "tree_newick"] = "(s1,s2,s3);"
    manifest.at[1, "tree_source"] = "raxml-ng"
    treeslider_mod._write_manifest(manifest, outdir / f"{name}.stats.tsv")

    def _unexpected_resolve(_binary):
        raise AssertionError("implicit binary lookup should not be called")

    monkeypatch.setattr(treeslider_mod, "_resolve_binary", _unexpected_resolve)

    treeslider_mod.run_treeslider_method(
        data=h5,
        name=name,
        outdir=outdir,
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

    stats = (outdir / f"{name}.stats.txt").read_text(encoding="utf-8")
    trees = (outdir / f"{name}.trees.nex").read_text(encoding="utf-8")
    assert "raxml_ng_binary: not_used" in stats
    assert "Tree window_000001" in trees


def test_run_treeslider_still_validates_explicit_binary_when_unused(
    tmp_path: Path,
) -> None:
    h5 = _write_sequence_h5(
        tmp_path / "assembly.hdf5",
        ["AAAA", "AAAA", "AAAA"],
        rows=[(0, 0, 4, 1, 4)],
    )

    with pytest.raises(treeslider_mod.IPyradError, match="Could not find the requested raxml-ng binary"):
        treeslider_mod.run_treeslider_method(
            data=h5,
            name="explicit-invalid",
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
            threads="auto",
            workers="auto",
            bs_trees=0,
            model="GTR+G",
            raxml_ng_binary=tmp_path / "missing-raxml-ng",
            seed=None,
            force=True,
            redo=False,
            log_level="INFO",
        )


def test_run_treeslider_cleans_stale_pending_tree_workdirs_before_resume(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    h5 = _write_sequence_h5(
        tmp_path / "assembly.hdf5",
        ["AAAA", "AATA", "AACA"],
        rows=[(0, 0, 4, 1, 4)],
    )
    outdir = tmp_path / "OUT"
    name = "resume"
    stage_dir = outdir / f".{name}.stage"
    align_dir = stage_dir / "alignments"
    align_dir.mkdir(parents=True, exist_ok=True)
    alignment_path = align_dir / "window_000001.fa"
    alignment_path.write_text(">s1\nAAAA\n>s2\nAATA\n>s3\nAACA\n", encoding="utf-8")

    phymap = treeslider_mod._load_phymap(h5)
    specs = treeslider_mod._plan_locus_windows(phymap, ["chr1"], ["chr1"])
    manifest = treeslider_mod._initialize_manifest(specs)
    manifest.at[1, "status"] = "accepted_pending_tree"
    manifest.at[1, "status_detail"] = "Alignment staged for tree inference."
    manifest.at[1, "nsamples_before_filtering"] = 3
    manifest.at[1, "nsites_before_filtering"] = 4
    manifest.at[1, "nvariants_before_filtering"] = 1
    manifest.at[1, "nsites_after_site_filter"] = 4
    manifest.at[1, "nsamples_after_sample_length_filter"] = 3
    manifest.at[1, "nsites_after_sample_length_filter"] = 4
    manifest.at[1, "nsamples_after_filtering"] = 3
    manifest.at[1, "nsites_after_filtering"] = 4
    manifest.at[1, "nvariants_after_filtering"] = 1
    manifest.at[1, "retained_sample_names"] = "s1,s2,s3"
    manifest.at[1, "alignment_path"] = str(alignment_path)
    treeslider_mod._write_manifest(manifest, outdir / f"{name}.stats.tsv")

    stale_dir = stage_dir / "raxml" / "window_000001"
    stale_dir.mkdir(parents=True, exist_ok=True)
    stale_marker = stale_dir / "stale.txt"
    stale_marker.write_text("stale", encoding="utf-8")

    def _mock_run_pipeline_assert_clean(cmds, **_kwargs):
        prefix = Path(cmds[0][cmds[0].index("--prefix") + 1])
        assert not (prefix.parent / "stale.txt").exists()
        (prefix.parent / f"{prefix.name}.raxml.bestTree").write_text("(s1,s2,s3);", encoding="utf-8")
        return 0, b"", b""

    monkeypatch.setattr(treeslider_mod, "_resolve_binary", lambda _binary: "/usr/bin/raxml-ng")
    monkeypatch.setattr(treeslider_mod, "run_pipeline", _mock_run_pipeline_assert_clean)
    monkeypatch.setattr(treeslider_mod, "run_with_pool_iter", _run_pool_iter_sequential)

    treeslider_mod.run_treeslider_method(
        data=h5,
        name=name,
        outdir=outdir,
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
        seed=9,
        force=False,
        redo=False,
        log_level="INFO",
    )

    manifest = pd.read_csv(outdir / f"{name}.stats.tsv", sep="\t")
    assert manifest["status"].tolist() == ["tree_completed"]
    assert not stage_dir.exists()


@pytest.mark.skipif(os.name != "posix", reason="Ctrl-C tests require POSIX")
def test_run_treeslider_interrupt_cleans_children_and_leaves_resumable_manifest(
    tmp_path: Path,
) -> None:
    h5 = _write_sequence_h5(
        tmp_path / "assembly.hdf5",
        ["ACGTACGT", "AGGAAGGA", "ATGCATGC"],
        rows=[
            (0, 0, 2, 1, 2),
            (0, 2, 4, 3, 4),
            (0, 4, 6, 5, 6),
            (0, 6, 8, 7, 8),
        ],
    )
    outdir = tmp_path / "OUT"
    pid_file = tmp_path / "raxml_pids.txt"
    raxml = tmp_path / "raxml-ng"
    raxml.write_text(
        "#!/bin/sh\n"
        "echo $$ >> \"$IPYRAD2_TEST_TREESLIDER_PID_FILE\"\n"
        "sleep 10\n"
        "prefix=\"\"\n"
        "while [ $# -gt 0 ]; do\n"
        "  if [ \"$1\" = \"--prefix\" ]; then\n"
        "    prefix=\"$2\"\n"
        "    shift 2\n"
        "    continue\n"
        "  fi\n"
        "  shift\n"
        "done\n"
        "printf '(s1,s2,s3);\\n' > \"${prefix}.raxml.bestTree\"\n",
        encoding="utf-8",
    )
    raxml.chmod(0o755)

    ctx = mp.get_context("spawn")
    proc = ctx.Process(
        target=_run_treeslider_until_interrupted,
        args=(h5, outdir, raxml, pid_file),
    )
    proc.start()

    deadline = time.monotonic() + 5.0
    while time.monotonic() < deadline and not pid_file.exists():
        time.sleep(0.05)
    assert pid_file.exists()

    try:
        os.kill(proc.pid, signal.SIGINT)
        time.sleep(0.05)
        os.kill(proc.pid, signal.SIGINT)
    except ProcessLookupError:
        pass

    start = time.time()
    proc.join(timeout=2.15)
    elapsed = time.time() - start
    if proc.is_alive():
        try:
            os.kill(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.join(0.5)
        raise AssertionError(
            f"child did not exit within 2.00s (+0.15s epsilon); elapsed={elapsed:.3f}s"
        )

    _assert_interrupted_exit(proc.exitcode)

    for line in pid_file.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        pid = int(line.strip())
        deadline = time.monotonic() + 1.0
        while time.monotonic() < deadline and _pid_is_alive(pid):
            time.sleep(0.05)
        assert not _pid_is_alive(pid)

    manifest = pd.read_csv(outdir / "interrupt.stats.tsv", sep="\t")
    stage_dir = outdir / ".interrupt.stage"
    align_dir = stage_dir / "alignments"
    raxml_root = stage_dir / "raxml"

    assert "accepted_pending_tree" in manifest["status"].tolist()
    assert "tree_completed" not in manifest["status"].tolist()
    assert align_dir.exists()
    assert any(align_dir.glob("window_*.fa"))
    if raxml_root.exists():
        assert not any(raxml_root.glob("window_*"))
    assert (outdir / "interrupt.stats.txt").exists()
    assert (outdir / "interrupt.trees.nex").exists()


def test_parallel_filter_jobs_match_serial_results(tmp_path: Path) -> None:
    h5 = _write_sequence_h5(
        tmp_path / "assembly.hdf5",
        [
            "AAAACCCCGGGGTTTT",
            "AAAACCCCGGGGTTTT",
            "AAAACCCCGGGGTTTT",
        ],
        rows=[
            (0, 0, 4, 1, 4),
            (0, 4, 8, 5, 8),
            (0, 8, 12, 9, 12),
            (0, 12, 16, 13, 16),
        ],
    )
    common_kwargs = dict(
        data=h5,
        name="slider",
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
        threads="auto",
        workers="auto",
        bs_trees=0,
        model="GTR+G",
        raxml_ng_binary=None,
        seed=17,
        force=True,
        redo=False,
        log_level="WARNING",
    )

    with pytest.raises(treeslider_mod.IPyradError, match="--jobs must be at least 1"):
        treeslider_mod.run_treeslider_method(
            outdir=tmp_path / "INVALID",
            jobs=0,
            **common_kwargs,
        )

    serial_out = tmp_path / "SERIAL"
    parallel_out = tmp_path / "PARALLEL"
    treeslider_mod.run_treeslider_method(
        outdir=serial_out,
        jobs=1,
        **common_kwargs,
    )
    treeslider_mod.run_treeslider_method(
        outdir=parallel_out,
        jobs=2,
        **common_kwargs,
    )

    serial_manifest = pd.read_csv(
        serial_out / "slider.stats.tsv",
        sep="\t",
        keep_default_na=False,
    )
    parallel_manifest = pd.read_csv(
        parallel_out / "slider.stats.tsv",
        sep="\t",
        keep_default_na=False,
    )
    pd.testing.assert_frame_equal(serial_manifest, parallel_manifest)
    assert serial_manifest["status"].tolist() == ["polytomy_written"] * 4
    assert (serial_out / "slider.trees.nex").read_text(encoding="utf-8") == (
        parallel_out / "slider.trees.nex"
    ).read_text(encoding="utf-8")

    serial_stats = (serial_out / "slider.stats.txt").read_text(encoding="utf-8")
    parallel_stats = (parallel_out / "slider.stats.txt").read_text(encoding="utf-8")
    assert "filter_jobs_requested: 1" in serial_stats
    assert "filter_jobs_resolved: 1" in serial_stats
    assert "filter_jobs_requested: 2" in parallel_stats
    assert "filter_jobs_resolved: 2" in parallel_stats
    assert not (serial_out / ".slider.stage").exists()
    assert not (parallel_out / ".slider.stage").exists()


def test_filter_manifest_writes_are_batched(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    h5 = _write_sequence_h5(
        tmp_path / "assembly.hdf5",
        ["AAAACCCCGGGGTTTT"],
        rows=[
            (0, 0, 4, 1, 4),
            (0, 4, 8, 5, 8),
            (0, 8, 12, 9, 12),
            (0, 12, 16, 13, 16),
        ],
        sample_names=["s1"],
    )
    monkeypatch.setattr(treeslider_mod, "run_with_pool_iter", _run_pool_iter_sequential)
    monkeypatch.setattr(treeslider_mod, "FILTER_MANIFEST_CHECKPOINT_BATCH_SIZE", 2)
    monkeypatch.setattr(
        treeslider_mod,
        "FILTER_MANIFEST_CHECKPOINT_SECONDS",
        float("inf"),
    )
    real_write_manifest = treeslider_mod._write_manifest
    manifest_writes: list[Path] = []

    def _record_manifest_write(manifest: pd.DataFrame, path: Path) -> None:
        real_write_manifest(manifest, path)
        manifest_writes.append(Path(path))

    monkeypatch.setattr(treeslider_mod, "_write_manifest", _record_manifest_write)

    outdir = tmp_path / "OUT"
    treeslider_mod.run_treeslider_method(
        data=h5,
        name="slider",
        outdir=outdir,
        window_size=None,
        slide_size=None,
        print_scaffold_table=False,
        scaffolds=None,
        min_sample_coverage=1,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        min_sample_alignment_length=1,
        min_alignment_length=1,
        jobs=2,
        threads="auto",
        workers="auto",
        bs_trees=0,
        model="GTR+G",
        raxml_ng_binary=None,
        seed=17,
        force=True,
        redo=False,
        log_level="WARNING",
    )

    manifest = pd.read_csv(outdir / "slider.stats.tsv", sep="\t")
    assert manifest["status"].tolist() == ["skipped_few_samples"] * 4
    assert manifest_writes == [outdir / "slider.stats.tsv"] * 3


def test_force_and_redo_are_rejected_before_output_cleanup(tmp_path: Path) -> None:
    outdir = tmp_path / "OUT"
    outdir.mkdir()
    manifest_path = outdir / "slider.stats.tsv"
    manifest_path.write_text("preserve me\n", encoding="utf-8")

    with pytest.raises(
        treeslider_mod.IPyradError,
        match="--force and --redo cannot be used together",
    ):
        treeslider_mod.run_treeslider_method(
            data=tmp_path / "missing.hdf5",
            name="slider",
            outdir=outdir,
            window_size=None,
            slide_size=None,
            print_scaffold_table=False,
            scaffolds=None,
            min_sample_coverage=4,
            imap=None,
            minmap=None,
            exclude=None,
            include_reference=False,
            min_sample_alignment_length=1,
            min_alignment_length=1,
            jobs=1,
            threads="auto",
            workers="auto",
            bs_trees=0,
            model="GTR+G",
            raxml_ng_binary=None,
            seed=None,
            force=True,
            redo=True,
            log_level="WARNING",
        )

    assert manifest_path.read_text(encoding="utf-8") == "preserve me\n"
