from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import ipyrad2.analysis.methods.bpp as bpp_mod
from ipyrad2.analysis.methods.bpp import Bpp
from ipyrad2.analysis.methods.bpp import Transformer
from ipyrad2.analysis.methods.bpp import _call_bpp
from ipyrad2.utils.exceptions import IPyradError


def _make_fake_bpp_binary(path: Path, body: str = "exit 0\n") -> Path:
    path.write_text(f"#!/bin/sh\n{body}", encoding="utf-8")
    path.chmod(0o755)
    return path


def test_bpp_honors_explicit_binary_path(tmp_path: Path) -> None:
    binary = _make_fake_bpp_binary(tmp_path / "bpp")
    tool = Bpp(name="demo", binary=str(binary))
    assert tool.kwargs["binary"] == str(binary.resolve())


def test_bpp_finds_binary_on_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    binary = _make_fake_bpp_binary(tmp_path / "bpp")
    monkeypatch.setattr(bpp_mod.shutil, "which", lambda exe: str(binary) if exe == "bpp" else None)
    tool = Bpp(name="demo")
    assert tool.kwargs["binary"] == str(binary.resolve())


def test_write_ctlfile_uses_current_bpp_prior_syntax(tmp_path: Path) -> None:
    binary = _make_fake_bpp_binary(tmp_path / "bpp")
    tool = Bpp(
        name="demo",
        workdir=tmp_path,
        guidetree="(a,b);",
        imap={"a": ["a"], "b": ["b"]},
        binary=str(binary),
    )
    tool.tree = bpp_mod.toytree.tree(tool.guidetree)
    tool.lex = SimpleNamespace(outfile=str(tmp_path / "demo.seqfile.txt"))
    tool._name = "demo_r0"
    tool._seed = 123

    tool._write_mapfile()
    ctlfile = Path(tool._write_ctlfile())
    text = ctlfile.read_text(encoding="utf-8")

    assert "thetaprior = invgamma 3 0.002 E" in text
    assert "tauprior = invgamma 3 0.002" in text
    assert "phiprior = beta 1 1" in text
    assert "print = 1 0 0 0 0" in text


def test_transformer_uses_empirical_parameter_samples() -> None:
    df = pd.DataFrame({
        "tau_demo": [0.5, 1.0],
        "theta_demo": [0.04, 0.08],
    })
    tx = Transformer(df, 1.0, 2.0, 1.0, 2.0, seed=7)
    tx.gentime_rvs = np.array([10.0, 20.0])
    tx.mutrate_rvs = np.array([0.01, 0.02])

    np.testing.assert_allclose(tx.transform("tau_demo"), np.array([500.0, 1000.0]))
    np.testing.assert_allclose(tx.transform("theta_demo"), np.array([1.0, 1.0]))


def test_bpp_transform_uses_transformed_tau_for_multitree(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    class FakeTransformer:
        def __init__(self, df, gentime_min, gentime_max, mutrate_min, mutrate_max, seed=123):
            del gentime_min, gentime_max, mutrate_min, mutrate_max
            self.df = df
            self.seed = seed

        def transform(self, colname):
            if colname.startswith("tau_"):
                return np.array([100.0, 200.0])
            if colname.startswith("theta_"):
                return np.array([1000.0, 2000.0])
            raise IPyradError(f"unexpected parameter: {colname}")

    monkeypatch.setattr(bpp_mod, "Transformer", FakeTransformer)

    binary = _make_fake_bpp_binary(tmp_path / "bpp")
    tool = Bpp(
        name="demo",
        guidetree="(a,b);",
        imap={"a": ["a"], "b": ["b"]},
        binary=str(binary),
    )
    mcmc = pd.DataFrame({
        "tau___a___b": [1.0, 2.0],
        "theta___a": [0.01, 0.02],
    })

    divs, popsize, newtree, mtree = tool.transform(mcmc, 1.0, 2.0, 1.0, 2.0, nsamp=2)

    assert divs.loc["median"].iloc[0] == 150.0
    assert popsize.loc["median"].iloc[0] == 1500.0
    assert all(node.dist == 150.0 for node in newtree.treenode.traverse() if node.is_leaf())
    leaf_dists = [
        [node.dist for node in tree.treenode.traverse() if node.is_leaf()]
        for tree in mtree.treelist
    ]
    assert leaf_dists[0] == [100.0, 100.0]
    assert leaf_dists[1] == [200.0, 200.0]


def test_run_skips_existing_ctl_without_submitting_job(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    class DummyLex:
        def __init__(self, *args, **kwargs):
            del args, kwargs
            self.outfile = str(tmp_path / "loci.bpp.txt")

        def _run(self, postfix: str) -> None:
            self.outfile = str(tmp_path / f"loci_{postfix}.bpp.txt")

    observed = {}

    def _fake_run_with_pool(jobs, log_level, cores, msg):
        observed["jobs"] = jobs
        observed["log_level"] = log_level
        observed["cores"] = cores
        observed["msg"] = msg
        return {}

    binary = _make_fake_bpp_binary(tmp_path / "bpp")
    data = tmp_path / "assembly.hdf5"
    data.write_text("", encoding="utf-8")
    (tmp_path / "demo_r0.ctl.txt").write_text("existing\n", encoding="utf-8")

    monkeypatch.setattr(bpp_mod, "LocusExtracter", DummyLex)
    monkeypatch.setattr(bpp_mod, "run_with_pool", _fake_run_with_pool)

    tool = Bpp(
        name="demo",
        data=str(data),
        workdir=tmp_path,
        guidetree="(a,b);",
        imap={"a": ["a"], "b": ["b"]},
        binary=str(binary),
    )
    tool._run(force=False, nreps=1, dry_run=False)

    assert observed.get("jobs", {}) == {}


def test_summarize_00_returns_concatenated_mcmc_table(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    binary = _make_fake_bpp_binary(tmp_path / "bpp")
    tool = Bpp(
        name="demo",
        workdir=tmp_path,
        guidetree="(a,b);",
        imap={"a": ["a"], "b": ["b"]},
        binary=str(binary),
    )

    mcmc1 = tmp_path / "demo_r0.mcmc.txt"
    mcmc2 = tmp_path / "demo_r1.mcmc.txt"
    cond1 = tmp_path / "demo_r0.conditional_a1b1.txt"
    cond2 = tmp_path / "demo_r1.conditional_a1b1.txt"
    ctl = tmp_path / "demo_r0.ctl.txt"
    ctl.write_text("jobname = demo\nprint = 1 0 0 0 0\n", encoding="utf-8")

    pd.DataFrame({"theta___a": [1.0, 2.0]}).to_csv(mcmc1, sep="\t")
    pd.DataFrame({"theta___a": [3.0, 4.0]}).to_csv(mcmc2, sep="\t")
    pd.DataFrame({"conditional": [10.0]}).to_csv(cond1, sep="\t")
    pd.DataFrame({"conditional": [20.0]}).to_csv(cond2, sep="\t")

    tool.files.mcmcfiles = [str(mcmc1), str(mcmc2)]
    monkeypatch.setattr(bpp_mod, "_call_bpp", lambda binary, ctlfile, alg: None)
    monkeypatch.setattr(tool, "_parse_A00_out", lambda ofile: pd.DataFrame({"ok": [1.0]}))

    _, concat = tool._summarize_00(False)

    assert list(concat.columns) == ["theta___a"]
    assert concat.shape[0] == 4


def test_call_bpp_runs_in_ctl_directory_and_cleans_side_effects(tmp_path: Path) -> None:
    binary = _make_fake_bpp_binary(
        tmp_path / "bpp",
        body="pwd > called_from.txt\n: > FigTree.tre\n: > SeedUsed\nexit 0\n",
    )
    ctlfile = tmp_path / "job.ctl.txt"
    ctlfile.write_text("jobname = job\n", encoding="utf-8")

    _call_bpp(str(binary), str(ctlfile), "00")

    assert (tmp_path / "called_from.txt").read_text(encoding="utf-8").strip() == str(tmp_path)
    assert (tmp_path / "job.figtree.nex").exists()
    assert not (tmp_path / "SeedUsed").exists()
