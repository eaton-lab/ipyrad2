from __future__ import annotations

import os
from pathlib import Path

import h5py
import numpy as np
import pytest

import ipyrad2.analysis.methods.bpp as bpp_mod
from ipyrad2.analysis.methods.bpp import Bpp
from ipyrad2.analysis.methods.bpp import _call_bpp
from ipyrad2.utils.exceptions import IPyradError


def _write_test_h5(
    path: Path,
    sequences: list[str],
    *,
    rows: list[tuple[int, int, int, int, int]],
    sample_names: list[str],
    scaffold_length: int,
) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    phy = np.vstack(
        [np.frombuffer(sequence.encode("utf-8"), dtype=np.uint8) for sequence in sequences]
    )
    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["scaffold_names"] = np.array(["chr1"], dtype=string_dtype)
        io5.attrs["scaffold_lengths"] = np.array([scaffold_length], dtype=np.uint64)
        io5.attrs["names"] = np.array(sample_names, dtype=string_dtype)
        io5.create_dataset("phy", data=phy)
        phymap = io5.create_dataset("phymap", data=np.array(rows, dtype=np.uint64))
        phymap.attrs["columns"] = np.array(
            ["scaff", "phy0", "phy1", "pos0", "pos1"],
            dtype=string_dtype,
        )
    return path


def _make_fake_bpp_binary(path: Path, body: str = "exit 0\n") -> Path:
    path.write_text(f"#!/bin/sh\n{body}", encoding="utf-8")
    path.chmod(0o755)
    return path


class _DummySeqex:
    def __init__(self, **kwargs) -> None:
        self.alignment_path = Path(kwargs["outdir"]) / f"{kwargs['name']}.phy"
        self.stats_path = Path(kwargs["outdir"]) / f"{kwargs['name']}.stats.txt"
        self.stats_json_path = (
            Path(kwargs["outdir"]) / f"{kwargs['name']}.stats.json"
        )

    def run(self) -> list[object]:
        self.alignment_path.write_text(
            "2 4\na^a     AAAA\nb^b     CCCC\n",
            encoding="utf-8",
        )
        self.stats_path.write_text("# Seqex Summary\n", encoding="utf-8")
        self.stats_json_path.write_text("{}\n", encoding="utf-8")
        return [object(), object()]


def test_bpp_write_inputs_does_not_resolve_binary(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(bpp_mod, "SeqexEngine", _DummySeqex)
    monkeypatch.setattr(
        bpp_mod,
        "_resolve_bpp_binary",
        lambda binary: (_ for _ in ()).throw(AssertionError(binary)),
    )

    data = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAA", "CCCC"],
        rows=[(0, 0, 4, 1, 4)],
        sample_names=["a", "b"],
        scaffold_length=4,
    )
    tool = Bpp(
        data=data,
        name="demo",
        outdir=tmp_path,
        tree="(a,b);",
        imap={"a": ["a"], "b": ["b"]},
        minmap={"a": 1, "b": 1},
        max_loci=2,
        min_length=4,
    )

    paths = tool.write_inputs()

    assert paths.seqfile.exists()
    assert paths.mapfile.exists()
    assert paths.ctlfile.exists()


def test_bpp_write_ctlfile_uses_new_defaults(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(bpp_mod, "SeqexEngine", _DummySeqex)
    data = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAA", "CCCC"],
        rows=[(0, 0, 4, 1, 4)],
        sample_names=["a", "b"],
        scaffold_length=4,
    )
    tool = Bpp(
        data=data,
        name="demo",
        outdir=tmp_path,
        tree="(a,b);",
        imap={"a": ["a"], "b": ["b"]},
        minmap={"a": 1, "b": 1},
        max_loci=2,
        min_length=4,
        seed=123,
    )

    tool.write_inputs()
    text = tool.paths.ctlfile.read_text(encoding="utf-8")

    assert "thetaprior = invgamma 3 0.03 E" in text
    assert "tauprior = invgamma 3 0.03" in text
    assert "alphaprior = 1 1 4" in text
    assert "locusrate = 1 2 3 2 iid" in text
    assert "clock = 2 10.0 100.0 5.0 iid LN" in text
    assert "phase = 1 1" in text
    assert "print = 1 0 0 1 0" in text
    assert "thetamodel = linked-none" in text
    assert "geneflow = 0" in text
    assert "seed = 123" in text


def test_bpp_msc_m_supports_parenthesized_references(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setattr(bpp_mod, "SeqexEngine", _DummySeqex)
    data = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAA", "CCCC", "GGGG"],
        rows=[(0, 0, 4, 1, 4)],
        sample_names=["a", "b", "c"],
        scaffold_length=4,
    )
    tool = Bpp(
        data=data,
        name="demo",
        outdir=tmp_path,
        tree="((a,b),c);",
        imap={"a": ["a"], "b": ["b"], "c": ["c"]},
        minmap={"a": 1, "b": 1, "c": 1},
        max_loci=2,
        min_length=4,
        msc_m=["a,b", "(a,b),c"],
    )

    tool.write_inputs()
    text = tool.paths.ctlfile.read_text(encoding="utf-8")

    assert "geneflow = 1" in text
    assert "wprior = 2 200" in text
    assert "migration = 2" in text
    assert "  ___a ___b" in text
    assert "  (___a,___b) ___c" in text


def test_bpp_speciestree_rejects_speciesmodelprior_two_or_three(tmp_path: Path) -> None:
    data = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAA", "CCCC"],
        rows=[(0, 0, 4, 1, 4)],
        sample_names=["a", "b"],
        scaffold_length=4,
    )
    with pytest.raises(IPyradError, match="speciesmodelprior 0 or 1"):
        Bpp(
            data=data,
            name="demo",
            outdir=tmp_path,
            tree="(a,b);",
            imap={"a": ["a"], "b": ["b"]},
            minmap={"a": 1, "b": 1},
            max_loci=2,
            min_length=4,
            speciestree=True,
            speciesmodelprior=2,
        )


def test_bpp_threads_must_have_one_or_three_integers(tmp_path: Path) -> None:
    data = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAA", "CCCC"],
        rows=[(0, 0, 4, 1, 4)],
        sample_names=["a", "b"],
        scaffold_length=4,
    )
    with pytest.raises(IPyradError, match="1 or 3 positive integers"):
        Bpp(
            data=data,
            name="demo",
            outdir=tmp_path,
            tree="(a,b);",
            imap={"a": ["a"], "b": ["b"]},
            minmap={"a": 1, "b": 1},
            max_loci=2,
            min_length=4,
            threads=[1, 2],
        )


def test_bpp_locus_sampling_is_deterministic_with_seed(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAACCCCGGGG", "TTTTGGGGCCCC"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8), (0, 8, 12, 9, 12)],
        sample_names=["a", "b"],
        scaffold_length=12,
    )
    common = dict(
        data=h5,
        tree="(sp1,sp2);",
        imap={"sp1": ["a"], "sp2": ["b"]},
        minmap={"sp1": 1, "sp2": 1},
        max_loci=2,
        min_length=4,
        seed=11,
        force=True,
    )
    tool1 = Bpp(name="demo", outdir=tmp_path / "run1", **common)
    tool2 = Bpp(name="demo", outdir=tmp_path / "run2", **common)

    tool1.write_inputs()
    tool2.write_inputs()

    assert tool1.paths.seqfile.read_text(encoding="utf-8") == tool2.paths.seqfile.read_text(
        encoding="utf-8"
    )


def test_bpp_stages_seqex_names_actual_locus_count_and_effective_threads(
    tmp_path: Path,
) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAACCCC", "TTTTGGGG"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        sample_names=["a", "b"],
        scaffold_length=8,
    )
    tool = Bpp(
        data=h5,
        name="demo",
        outdir=tmp_path / "run",
        tree="(sp1,sp2);",
        imap={"sp1": ["a"], "sp2": ["b"]},
        minmap={"sp1": 1, "sp2": 1},
        max_loci=3,
        min_length=4,
        threads=(3, 4, 1),
        seed=11,
    )

    paths = tool.write_inputs()

    sequence_text = paths.seqfile.read_text(encoding="utf-8")
    assert sequence_text.count("2 4\n") == 2
    assert "sp1^a" in sequence_text
    assert "sp2^b" in sequence_text
    assert "^___" not in sequence_text
    map_rows = {
        tuple(line.split())
        for line in paths.mapfile.read_text(encoding="utf-8").splitlines()
    }
    assert map_rows == {("a", "___sp1"), ("b", "___sp2")}
    ctl_text = paths.ctlfile.read_text(encoding="utf-8")
    assert "nloci = 2" in ctl_text
    assert "threads = 2 4 1" in ctl_text
    assert paths.statsfile.exists()
    assert paths.statsjson.exists()
    assert tool.seqex.cores == 3


@pytest.mark.skipif(
    not os.environ.get("IPYRAD2_BPP_BINARY"),
    reason="set IPYRAD2_BPP_BINARY to run the external BPP smoke test",
)
def test_bpp_real_binary_accepts_seqex_staged_inputs(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAACCCC", "TTTTGGGG"],
        rows=[(0, 0, 4, 1, 4), (0, 4, 8, 5, 8)],
        sample_names=["a", "b"],
        scaffold_length=8,
    )
    tool = Bpp(
        data=h5,
        name="smoke",
        outdir=tmp_path / "run",
        tree="(sp1,sp2);",
        imap={"sp1": ["a"], "sp2": ["b"]},
        minmap={"sp1": 1, "sp2": 1},
        max_loci=2,
        min_length=4,
        burnin=1,
        samplefreq=1,
        nsample=1,
        threads=1,
        seed=11,
        binary=os.environ["IPYRAD2_BPP_BINARY"],
    )

    paths = tool.run()

    assert paths.outfile.exists()
    assert paths.mcmcfile.exists()


def test_bpp_expands_glob_imap_entries_against_hdf5_sample_names(tmp_path: Path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAA", "AAAT", "CCCC"],
        rows=[(0, 0, 4, 1, 4)],
        sample_names=["barbeyi-01", "barbeyi-02", "geyeri-01"],
        scaffold_length=4,
    )
    imap = tmp_path / "imap.tsv"
    imap.write_text(
        "barbeyi*\tbarbeyi\n"
        "geyeri*\tgeyeri\n",
        encoding="utf-8",
    )
    minmap = tmp_path / "minmap.tsv"
    minmap.write_text("barbeyi\t1\ngeyeri\t1\n", encoding="utf-8")

    tool = Bpp(
        data=h5,
        name="demo",
        outdir=tmp_path / "out",
        tree="(barbeyi,geyeri);",
        imap=imap,
        minmap=minmap,
        max_loci=1,
        min_length=4,
    )

    assert tool.imap == {
        "barbeyi": ["barbeyi-01", "barbeyi-02"],
        "geyeri": ["geyeri-01"],
    }
    assert tool.minmap == {"barbeyi": 1, "geyeri": 1}


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
