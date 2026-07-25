from __future__ import annotations

import os
from pathlib import Path

import h5py
import numpy as np
import pytest

import ipyrad2.analysis.methods.bpp as bpp_mod
from ipyrad2.analysis.methods.bpp import Bpp
from ipyrad2.analysis.methods.bpp import _BPP_BINARY_SPECS
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
            ["scaff", "phy0", "phy1", "pos0", "pos1"], dtype=string_dtype
        )
    return path


def _make_fake_bpp_binary(path: Path, version="4.8.7", body="exit 0\n") -> Path:
    path.write_text(
        "#!/bin/sh\n"
        f'if [ "$1" = "--help" ]; then echo "BPP {version}"; exit 0; fi\n'
        f"{body}",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


class _DummySeqex:
    def __init__(self, **kwargs) -> None:
        self.alignment_path = Path(kwargs["outdir"]) / f"{kwargs['name']}.phy"
        self.stats_path = Path(kwargs["outdir"]) / f"{kwargs['name']}.stats.txt"
        self.stats_json_path = Path(kwargs["outdir"]) / f"{kwargs['name']}.stats.json"
        self.cores = kwargs["cores"]
        self.random_seed = kwargs["random_seed"]

    def run(self) -> list[object]:
        self.alignment_path.write_text(
            "2 4\na^a     AAAA\nb^b     CCCC\n", encoding="utf-8"
        )
        self.stats_path.write_text("# Seqex Summary\n", encoding="utf-8")
        self.stats_json_path.write_text("{}\n", encoding="utf-8")
        return [object(), object()]


def _data(tmp_path: Path, names=("a", "b")) -> Path:
    sequences = ["AAAA", "CCCC", "GGGG"][: len(names)]
    return _write_test_h5(
        tmp_path / "assembly.hdf5",
        sequences,
        rows=[(0, 0, 4, 1, 4)],
        sample_names=list(names),
        scaffold_length=4,
    )


def _kwargs(tmp_path: Path, **overrides):
    names = overrides.pop("names", ("a", "b"))
    values = dict(
        data=_data(tmp_path, names),
        name="demo",
        outdir=tmp_path / "out",
        tree="(" + ",".join(names) + ");",
        imap={name: [name] for name in names},
        minmap={name: 1 for name in names},
        max_loci=2,
        min_length=4,
    )
    values.update(overrides)
    return values


def test_write_only_does_not_resolve_binary(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(bpp_mod, "SeqexEngine", _DummySeqex)
    monkeypatch.setattr(
        bpp_mod,
        "_resolve_bpp_binary",
        lambda binary: (_ for _ in ()).throw(AssertionError(binary)),
    )
    tool = Bpp(**_kwargs(tmp_path))
    paths = tool.write_inputs()
    assert paths.seqfile.exists()
    assert paths.mapfile.exists()
    assert paths.ctlfile.exists()


def test_control_file_uses_conservative_defaults(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(bpp_mod, "SeqexEngine", _DummySeqex)
    tool = Bpp(**_kwargs(tmp_path, locus_seed=11, mcmc_seed=123))
    tool.write_inputs()
    text = tool.paths.ctlfile.read_text(encoding="utf-8")

    assert "speciesdelimitation = 0" in text
    assert "speciestree = 0" in text
    assert "speciesmodelprior" not in text
    assert "thetaprior = invgamma 3 0.002" in text
    assert "tauprior = invgamma 3 0.03" in text
    assert "alphaprior" not in text
    assert "locusrate = 0" in text
    assert "clock = 1" in text
    assert "model = JC69" in text
    assert "phase = 1 1" in text
    assert "print = 1 0 0 0 0" in text
    assert "burnin = 10000" in text
    assert "sampfreq = 10" in text
    assert "nsample = 100000" in text
    assert "seed = 123" in text
    assert tool.seqex.random_seed == 11


@pytest.mark.parametrize(
    ("analysis", "delimitation", "tree", "prior"),
    [
        ("a00", "0", "0", None),
        ("A01", "0", "1", "1"),
        ("a10", "1 0 2", "0", "1"),
        ("A11", "1 0 2", "1", "1"),
    ],
)
def test_analysis_modes_render_correct_controls(
    monkeypatch, tmp_path, analysis, delimitation, tree, prior
) -> None:
    monkeypatch.setattr(bpp_mod, "SeqexEngine", _DummySeqex)
    tool = Bpp(**_kwargs(tmp_path, analysis=analysis))
    tool.write_inputs()
    text = tool.paths.ctlfile.read_text(encoding="utf-8")
    assert tool.analysis == analysis.casefold()
    assert f"speciesdelimitation = {delimitation}" in text
    assert f"speciestree = {tree}" in text
    if prior is None:
        assert "speciesmodelprior" not in text
    else:
        assert f"speciesmodelprior = {prior}" in text


def test_msc_i_and_msc_m_render_gene_flow(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(bpp_mod, "SeqexEngine", _DummySeqex)
    msci = Bpp(
        **_kwargs(
            tmp_path,
            outdir=tmp_path / "msci",
            analysis="MSC-I",
            tree="(a,b)[&phi=0.5];",
            theta_model="linked-msci",
        )
    )
    msci.write_inputs()
    msci_text = msci.paths.ctlfile.read_text(encoding="utf-8")
    assert "geneflow = 1" in msci_text
    assert "phiprior = 1 1" in msci_text
    assert "thetamodel = linked-msci" in msci_text

    mscm = Bpp(
        **_kwargs(
            tmp_path,
            outdir=tmp_path / "mscm",
            analysis="MSC-M",
            migration=["a,b", "b,a"],
            theta_model="linked-mscm",
        )
    )
    mscm.write_inputs()
    mscm_text = mscm.paths.ctlfile.read_text(encoding="utf-8")
    assert "geneflow = 1" in mscm_text
    assert "wprior = 2 200" in mscm_text
    assert "migration = 2" in mscm_text
    assert "  ___a ___b" in mscm_text
    assert "  ___b ___a" in mscm_text


@pytest.mark.parametrize(
    ("analysis", "prior"),
    [("a01", 2), ("a10", 3), ("a00", 1), ("msc-i", 0)],
)
def test_species_model_prior_compatibility(tmp_path, analysis, prior) -> None:
    kwargs = _kwargs(tmp_path, analysis=analysis, species_model_prior=prior)
    if analysis == "msc-i":
        kwargs["tree"] = "(a,b)[&phi=0.5];"
    with pytest.raises(IPyradError, match="species_model_prior"):
        Bpp(**kwargs)


def test_a11_accepts_all_species_model_priors(tmp_path) -> None:
    for prior in range(4):
        tool = Bpp(
            **_kwargs(
                tmp_path,
                outdir=tmp_path / f"p{prior}",
                analysis="a11",
                species_model_prior=prior,
            )
        )
        assert tool.species_model_prior == prior


@pytest.mark.parametrize(
    "theta_prior,expected",
    [
        (("invgamma", 3, 0.002, "integrate"), "invgamma 3 0.002 int"),
        (("gamma", 2, 2000), "gamma 2 2000"),
        (("beta", 2, 3, 0.001, 0.1), "beta 2 3 0.001 0.1"),
    ],
)
def test_theta_prior_forms(monkeypatch, tmp_path, theta_prior, expected) -> None:
    monkeypatch.setattr(bpp_mod, "SeqexEngine", _DummySeqex)
    tool = Bpp(**_kwargs(tmp_path, theta_prior=theta_prior))
    tool.write_inputs()
    assert f"thetaprior = {expected}" in tool.paths.ctlfile.read_text(encoding="utf-8")


def test_optional_models_checkpoint_and_print_controls(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(bpp_mod, "SeqexEngine", _DummySeqex)
    tool = Bpp(
        **_kwargs(
            tmp_path,
            substitution_model="hky",
            alpha_prior=(1, 1, 4),
            locus_rate=("1", 2, 3, 2, "IID"),
            clock=("2", 10, 100, 5, "iid", "ln"),
            clean_data=True,
            checkpoint=(100000, 10000),
            print_locus_rates=True,
            print_gene_trees=True,
            print_substitution_parameters=True,
        )
    )
    tool.write_inputs()
    text = tool.paths.ctlfile.read_text(encoding="utf-8")
    assert "model = HKY" in text
    assert "alphaprior = 1 1 4" in text
    assert "locusrate = 1 2 3 2 iid" in text
    assert "clock = 2 10 100 5 iid LN" in text
    assert "cleandata = 1" in text
    assert "checkpoint = 100000 10000" in text
    assert "print = 1 1 0 1 1" in text


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"burnin": 200}, "greater than 200"),
        ({"locus_rate": ("2", "rates.txt")}, "not supported"),
        ({"locus_rate": ("3", 2, 3)}, "not supported"),
        (
            {
                "locus_rate": ("1", 2, 3, 2, "iid"),
                "clock": ("2", 10, 100, 5, "dir", "LN"),
            },
            "same iid or dir",
        ),
        ({"theta_model": "linked-msci"}, "only valid for MSC-I"),
    ],
)
def test_invalid_model_settings_fail_early(tmp_path, changes, message) -> None:
    with pytest.raises(IPyradError, match=message):
        Bpp(**_kwargs(tmp_path, **changes))


def test_constraint_file_is_validated_prefixed_and_staged(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(bpp_mod, "SeqexEngine", _DummySeqex)
    constraint = tmp_path / "constraints.txt"
    constraint.write_text(
        "define pair as (a,b); # keep comment\nconstraint pair;\noutgroup a;\n",
        encoding="utf-8",
    )
    tool = Bpp(
        **_kwargs(tmp_path, analysis="a11", constraint_file=constraint)
    )
    tool.write_inputs()
    staged = tool.paths.constraintfile.read_text(encoding="utf-8")
    assert "define pair as (___a,___b); # keep comment" in staged
    assert "constraint pair;" in staged
    assert "outgroup ___a;" in staged
    ctl = tool.paths.ctlfile.read_text(encoding="utf-8")
    assert f"constraintfile = {tool.paths.constraintfile}" in ctl


def test_constraint_file_rejects_unknown_population(tmp_path) -> None:
    constraint = tmp_path / "constraints.txt"
    constraint.write_text("constraint (a,c);\n", encoding="utf-8")
    with pytest.raises(IPyradError, match="unknown names: c"):
        Bpp(**_kwargs(tmp_path, analysis="a01", constraint_file=constraint))


def test_locus_sampling_is_independent_of_mcmc_seed(tmp_path) -> None:
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
        locus_seed=11,
    )
    first = Bpp(name="demo", outdir=tmp_path / "run1", mcmc_seed=1, **common)
    second = Bpp(name="demo", outdir=tmp_path / "run2", mcmc_seed=2, **common)
    first.write_inputs()
    second.write_inputs()
    assert first.paths.seqfile.read_text() == second.paths.seqfile.read_text()
    assert "seed = 1" in first.paths.ctlfile.read_text()
    assert "seed = 2" in second.paths.ctlfile.read_text()


def test_actual_locus_count_and_threads_are_staged(tmp_path) -> None:
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
        locus_seed=11,
    )
    paths = tool.write_inputs()
    assert paths.seqfile.read_text().count("2 4\n") == 2
    assert "sp1^a" in paths.seqfile.read_text()
    assert {tuple(line.split()) for line in paths.mapfile.read_text().splitlines()} == {
        ("a", "___sp1"), ("b", "___sp2")
    }
    ctl = paths.ctlfile.read_text()
    assert "nloci = 2" in ctl
    assert "threads = 2 4 1" in ctl
    assert tool.seqex.cores == 3


def test_glob_imap_entries_expand_against_hdf5_samples(tmp_path) -> None:
    h5 = _write_test_h5(
        tmp_path / "assembly.hdf5",
        ["AAAA", "AAAT", "CCCC"],
        rows=[(0, 0, 4, 1, 4)],
        sample_names=["barbeyi-01", "barbeyi-02", "geyeri-01"],
        scaffold_length=4,
    )
    imap = tmp_path / "imap.tsv"
    imap.write_text("barbeyi*\tbarbeyi\ngeyeri*\tgeyeri\n")
    minmap = tmp_path / "minmap.tsv"
    minmap.write_text("barbeyi\t1\ngeyeri\t1\n")
    tool = Bpp(
        data=h5,
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


def test_run_requires_supported_parseable_bpp_version(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(bpp_mod, "SeqexEngine", _DummySeqex)
    old = _make_fake_bpp_binary(tmp_path / "old-bpp", version="4.8.6")
    with pytest.raises(IPyradError, match="4.8.7 is required"):
        Bpp(**_kwargs(tmp_path, binary=str(old))).run()

    unknown = tmp_path / "unknown-bpp"
    unknown.write_text("#!/bin/sh\nexit 0\n")
    unknown.chmod(0o755)
    with pytest.raises(IPyradError, match="Could not determine"):
        Bpp(**_kwargs(tmp_path, outdir=tmp_path / "unknown", binary=str(unknown))).run()


def test_run_logs_mcmc_output_path(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(bpp_mod, "SeqexEngine", _DummySeqex)
    binary = _make_fake_bpp_binary(tmp_path / "bpp")
    messages = []
    monkeypatch.setattr(
        bpp_mod.logger,
        "info",
        lambda message, *args: messages.append(message.format(*args)),
    )
    monkeypatch.setattr(bpp_mod, "_call_bpp", lambda *args: None)

    tool = Bpp(**_kwargs(tmp_path, binary=str(binary)))
    tool.run()

    assert any(
        message == f"running BPP and writing MCMC samples to {tool.paths.mcmcfile}"
        for message in messages
    )


def test_all_download_specs_target_bpp_487() -> None:
    assert {spec[0] for spec in _BPP_BINARY_SPECS.values()} == {"4.8.7"}


@pytest.mark.skipif(
    not os.environ.get("IPYRAD2_BPP_BINARY"),
    reason="set IPYRAD2_BPP_BINARY to run the external BPP smoke test",
)
def test_real_binary_accepts_staged_inputs(tmp_path) -> None:
    tool = Bpp(
        **_kwargs(
            tmp_path,
            burnin=201,
            samplefreq=1,
            nsample=1,
            threads=1,
            locus_seed=11,
            mcmc_seed=11,
            binary=os.environ["IPYRAD2_BPP_BINARY"],
        )
    )
    paths = tool.run()
    assert paths.outfile.exists()
    assert paths.mcmcfile.exists()


def test_call_bpp_runs_in_ctl_directory_and_cleans_side_effects(tmp_path) -> None:
    binary = _make_fake_bpp_binary(
        tmp_path / "bpp",
        body="pwd > called_from.txt\n: > FigTree.tre\n: > SeedUsed\nexit 0\n",
    )
    ctlfile = tmp_path / "job.ctl.txt"
    ctlfile.write_text("jobname = job\n", encoding="utf-8")
    _call_bpp(str(binary), str(ctlfile), "a11")
    assert (tmp_path / "called_from.txt").read_text().strip() == str(tmp_path)
    assert (tmp_path / "job.figtree.nex").exists()
    assert not (tmp_path / "SeedUsed").exists()
