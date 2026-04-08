from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

from ipyrad2.analysis.methods.popgen import run_popgen_method
from ipyrad2.analysis.methods.popgen import runner as popgen_runner
from ipyrad2.analysis.methods.popgen import seq_backend as popgen_seq_backend
from ipyrad2.analysis.methods.popgen.estimators import summarize_genotype_block
from ipyrad2.analysis.methods.popgen.estimators import summarize_genotype_site
from ipyrad2.analysis.methods.popgen.estimators import summarize_sequence_block
from ipyrad2.analysis.methods.popgen.estimators import summarize_sequence_site
from ipyrad2.utils.exceptions import IPyradError


def _geno_triplet(a: int, b: int, char: str) -> np.ndarray:
    return np.array([a, b, ord(char)], dtype=np.uint8)


def _write_combined_popgen_h5(path: Path, *, include_phy: bool = True, include_genos: bool = True) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    names = ["s1", "s2", "s3", "s4"]

    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(names, dtype=string_dtype)
        io5.attrs["nsnps"] = 4

        if include_phy:
            io5.attrs["scaffold_names"] = np.array(["chr1"], dtype=string_dtype)
            io5.attrs["scaffold_lengths"] = np.array([4], dtype=np.uint64)
            phy = np.vstack(
                [
                    np.frombuffer(b"AAAA", dtype=np.uint8),
                    np.frombuffer(b"ARAA", dtype=np.uint8),
                    np.frombuffer(b"GGAA", dtype=np.uint8),
                    np.frombuffer(b"GGAA", dtype=np.uint8),
                ]
            )
            io5.create_dataset("phy", data=phy)
            phymap = io5.create_dataset(
                "phymap",
                data=np.array([(0, 0, 4, 1, 4)], dtype=np.uint64),
            )
            phymap.attrs["columns"] = np.array(
                ["scaff", "phy0", "phy1", "pos0", "pos1"],
                dtype=string_dtype,
            )

        if include_genos:
            genos = np.array(
                [
                    [
                        _geno_triplet(0, 0, "A"),
                        _geno_triplet(0, 0, "A"),
                        _geno_triplet(0, 0, "A"),
                        _geno_triplet(0, 0, "A"),
                    ],
                    [
                        _geno_triplet(0, 0, "A"),
                        _geno_triplet(0, 1, "R"),
                        _geno_triplet(0, 0, "A"),
                        _geno_triplet(0, 0, "A"),
                    ],
                    [
                        _geno_triplet(1, 1, "G"),
                        _geno_triplet(1, 1, "G"),
                        _geno_triplet(0, 0, "A"),
                        _geno_triplet(0, 0, "A"),
                    ],
                    [
                        _geno_triplet(1, 1, "G"),
                        _geno_triplet(1, 1, "G"),
                        _geno_triplet(0, 0, "A"),
                        _geno_triplet(0, 0, "A"),
                    ],
                ],
                dtype=np.uint8,
            )
            io5.create_dataset("genos", data=genos)
            io5["genos"].attrs["names"] = np.array(names, dtype=string_dtype)
            io5.create_dataset(
                "reference",
                data=np.array([ord("A"), ord("A"), ord("A"), ord("A")], dtype=np.uint8),
            )
            snpsmap = io5.create_dataset(
                "snpsmap",
                data=np.array(
                    [
                        [0, 0, 0, 0, 0],
                        [0, 1, 1, 0, 1],
                        [1, 0, 0, 0, 2],
                        [1, 1, 1, 0, 3],
                    ],
                    dtype=np.uint32,
                ),
            )
            snpsmap.attrs["columns"] = np.array(
                ["loc", "loc_idx", "loc_pos", "scaff", "pos"],
                dtype=string_dtype,
            )

    return path


def _write_imap(path: Path) -> Path:
    path.write_text(
        "s1\tpop1\ns2\tpop1\ns3\tpop2\ns4\tpop2\n",
        encoding="utf-8",
    )
    return path


def _write_large_sequence_popgen_h5(path: Path) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    names = ["s1", "s2", "s3", "s4"]
    nloci = 3
    loclen = 2500
    nsites = nloci * loclen
    rng = np.random.default_rng(7)
    phy = rng.choice(
        np.array([ord("A"), ord("C"), ord("G"), ord("T"), ord("N")], dtype=np.uint8),
        size=(len(names), nsites),
        p=[0.24, 0.24, 0.24, 0.24, 0.04],
    )
    phymap = np.zeros((nloci, 5), dtype=np.uint64)
    for idx in range(nloci):
        start = idx * loclen
        end = start + loclen
        phymap[idx] = [0, start, end, start + 1, end]

    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(names, dtype=string_dtype)
        io5.attrs["scaffold_names"] = np.array(["chr1"], dtype=string_dtype)
        io5.attrs["scaffold_lengths"] = np.array([nsites], dtype=np.uint64)
        io5.create_dataset("phy", data=phy)
        ds = io5.create_dataset("phymap", data=phymap)
        ds.attrs["columns"] = np.array(
            ["scaff", "phy0", "phy1", "pos0", "pos1"],
            dtype=string_dtype,
        )
    return path


def _write_windowed_sequence_popgen_h5(path: Path) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    names = ["s1", "s2", "s3", "s4"]
    phy = np.vstack(
        [
            np.frombuffer(b"AAAACCCCAAAA", dtype=np.uint8),
            np.frombuffer(b"ARAACCCCAAAA", dtype=np.uint8),
            np.frombuffer(b"GGAATTTTAAAA", dtype=np.uint8),
            np.frombuffer(b"GGAATTTTGGGG", dtype=np.uint8),
        ]
    )
    phymap = np.array(
        [
            (0, 0, 4, 1, 4),
            (0, 4, 8, 5, 8),
            (1, 8, 12, 1, 4),
        ],
        dtype=np.uint64,
    )

    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(names, dtype=string_dtype)
        io5.attrs["scaffold_names"] = np.array(["chr1", "chr2"], dtype=string_dtype)
        io5.attrs["scaffold_lengths"] = np.array([8, 4], dtype=np.uint64)
        io5.create_dataset("phy", data=phy)
        ds = io5.create_dataset("phymap", data=phymap)
        ds.attrs["columns"] = np.array(
            ["scaff", "phy0", "phy1", "pos0", "pos1"],
            dtype=string_dtype,
        )
    return path


def _write_missing_sequence_popgen_h5(path: Path) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    names = ["s1", "s2", "s3", "s4"]
    phy = np.vstack(
        [
            np.frombuffer(b"AN-A", dtype=np.uint8),
            np.frombuffer(b"ARNA", dtype=np.uint8),
            np.frombuffer(b"GGAA", dtype=np.uint8),
            np.frombuffer(b"GG-A", dtype=np.uint8),
        ]
    )

    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(names, dtype=string_dtype)
        io5.attrs["scaffold_names"] = np.array(["chr1"], dtype=string_dtype)
        io5.attrs["scaffold_lengths"] = np.array([4], dtype=np.uint64)
        io5.create_dataset("phy", data=phy)
        ds = io5.create_dataset(
            "phymap",
            data=np.array([(0, 0, 4, 1, 4)], dtype=np.uint64),
        )
        ds.attrs["columns"] = np.array(
            ["scaff", "phy0", "phy1", "pos0", "pos1"],
            dtype=string_dtype,
        )
    return path


def test_run_popgen_sequence_backend_writes_full_panel_from_combined_hdf5(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    h5 = _write_combined_popgen_h5(tmp_path / "assembly.hdf5")
    imap = _write_imap(tmp_path / "imap.tsv")
    logged_messages: list[str] = []

    def _capture_info(message, *args) -> None:
        logged_messages.append(str(message).format(*args))

    monkeypatch.setattr(popgen_runner.logger, "info", _capture_info)

    run_popgen_method(
        data=h5,
        name="pop",
        outdir=tmp_path / "OUT",
        stats="all",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=None,
        exclude=None,
        include_reference=False,
        subsample_unlinked=False,
        random_seed=None,
        cores=1,
        force=True,
        log_level="INFO",
    )

    outdir = tmp_path / "OUT"
    manifest_text = (outdir / "pop.manifest.txt").read_text(encoding="utf-8")
    sample_stats_text = (outdir / "pop.sample_stats.tsv").read_text(encoding="utf-8")
    population_text = (outdir / "pop.population_stats.tsv").read_text(encoding="utf-8")
    sample_stats = pd.read_csv(outdir / "pop.sample_stats.tsv", sep="\t")
    global_stats = pd.read_csv(outdir / "pop.global_stats.tsv", sep="\t")
    population = pd.read_csv(outdir / "pop.population_stats.tsv", sep="\t")
    pairwise = pd.read_csv(outdir / "pop.pairwise_stats.tsv", sep="\t")
    sfs = pd.read_csv(outdir / "pop.sfs.tsv", sep="\t")

    assert not (outdir / "pop.stats.txt").exists()
    assert not (outdir / "pop.sample_data_summary.tsv").exists()
    assert "backend_used: sequence" in manifest_text
    assert any("pop.manifest.txt" in message for message in logged_messages)
    assert any("pop.sample_stats.tsv" in message for message in logged_messages)
    assert any("pop.global_stats.tsv" in message for message in logged_messages)
    assert "Inputs\n------" in manifest_text
    assert "Requested Stats\n---------------" in manifest_text
    assert "requested_stats: ['pi', 'dxy', 'fst', 'tajima_d', 'theta_w', 'heterozygosity', 'fis', 'fit', 'sfs']" in manifest_text
    assert "fis_formula: fis = 1 - Ho/He" in manifest_text
    assert "fit_formula: fit = 1 - Ho/Ht_total" in manifest_text
    assert "sample_stats_rows: 4" in manifest_text
    assert "global_stats_rows: 1" in manifest_text
    assert "sample_data_summary_rows_in_manifest: 4" in manifest_text
    assert "Sample Data Summary\n-------------------" in manifest_text
    assert (
        "sample\tmissing_fraction\tpost_imputation_missing_fraction\timputation_algorithm\timputed_genotype_fraction"
        in manifest_text
    )
    assert "s1\t0.00000000\t0.00000000\tnot-imputed\t0.00000000" in manifest_text
    assert list(sample_stats["sample"]) == ["s1", "s2", "s3", "s4"]
    assert list(sample_stats["population"]) == ["pop1", "pop1", "pop2", "pop2"]
    assert np.allclose(sample_stats["sites_total"], 4)
    assert np.allclose(sample_stats["sites_called"], [4, 4, 4, 4])
    assert np.allclose(sample_stats["called_fraction"], [1.0, 1.0, 1.0, 1.0])
    assert np.allclose(sample_stats["sites_missing"], [0, 0, 0, 0])
    assert np.allclose(sample_stats["missing_fraction"], [0.0, 0.0, 0.0, 0.0])
    assert np.allclose(sample_stats["homozygous_sites"], [4, 3, 4, 4])
    assert np.allclose(sample_stats["heterozygous_sites"], [0, 1, 0, 0])
    assert np.allclose(sample_stats["observed_heterozygosity"], [0.0, 0.25, 0.0, 0.0])
    assert list(global_stats.columns) == [
        "sites_used_heterozygosity",
        "observed_heterozygosity",
        "expected_heterozygosity_total",
        "fit",
    ]
    assert int(global_stats.iloc[0]["sites_used_heterozygosity"]) == 4
    assert np.isclose(float(global_stats.iloc[0]["observed_heterozygosity"]), 0.0625)
    assert np.isclose(
        float(global_stats.iloc[0]["fit"]),
        1.0
        - (
            float(global_stats.iloc[0]["observed_heterozygosity"])
            / float(global_stats.iloc[0]["expected_heterozygosity_total"])
        ),
    )
    assert "called_fraction\t" in sample_stats_text
    assert "1.00000000\t0\t0.00000000\t4\t0\t0.00000000" in sample_stats_text
    assert "1.00000000\t0\t0.00000000\t3\t1\t0.25000000" in sample_stats_text
    assert set(population["population"]) == {"pop1", "pop2"}
    assert {"pi", "theta_w", "tajima_d", "observed_heterozygosity", "expected_heterozygosity", "fis"}.issubset(population.columns)
    assert "0.00000000" in population_text
    pop1 = population.loc[population["population"] == "pop1"].iloc[0]
    pop2 = population.loc[population["population"] == "pop2"].iloc[0]
    assert float(pop1["pi"]) > 0.0
    assert float(pop2["pi"]) == 0.0
    assert float(pop1["observed_heterozygosity"]) > 0.0
    assert np.isnan(float(pop2["fis"]))
    assert set(pairwise.columns) >= {"population1", "population2", "dxy", "fst"}
    assert float(pairwise.iloc[0]["fst"]) > 0.5
    assert set(sfs["population"]) == {"pop1"}


def test_run_popgen_snp_backend_supports_supported_subset_and_unlinked_sampling(tmp_path: Path) -> None:
    h5 = _write_combined_popgen_h5(tmp_path / "snps.hdf5", include_phy=False, include_genos=True)
    imap = _write_imap(tmp_path / "imap.tsv")

    run_popgen_method(
        data=h5,
        name="pop",
        outdir=tmp_path / "OUT",
        stats="all",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=None,
        exclude=None,
        include_reference=False,
        subsample_unlinked=True,
        random_seed=7,
        cores=1,
        force=True,
        log_level="INFO",
    )

    outdir = tmp_path / "OUT"
    manifest_text = (outdir / "pop.manifest.txt").read_text(encoding="utf-8")
    sample_stats_text = (outdir / "pop.sample_stats.tsv").read_text(encoding="utf-8")
    sample_stats = pd.read_csv(outdir / "pop.sample_stats.tsv", sep="\t")
    global_stats = pd.read_csv(outdir / "pop.global_stats.tsv", sep="\t")
    population = pd.read_csv(outdir / "pop.population_stats.tsv", sep="\t")
    pairwise = pd.read_csv(outdir / "pop.pairwise_stats.tsv", sep="\t")
    sfs = pd.read_csv(outdir / "pop.sfs.tsv", sep="\t")

    assert not (outdir / "pop.stats.txt").exists()
    assert not (outdir / "pop.sample_data_summary.tsv").exists()
    assert "backend_used: snp" in manifest_text
    assert "subsample_unlinked: True" in manifest_text
    assert "requested_stats: ['fst', 'heterozygosity', 'fis', 'fit', 'sfs']" in manifest_text
    assert "fit_formula: fit = 1 - Ho/Ht_total" in manifest_text
    assert "sample_stats_rows: 4" in manifest_text
    assert "global_stats_rows: 1" in manifest_text
    assert "sample_data_summary_rows_in_manifest: 4" in manifest_text
    assert "s1\t0.00000000\t0.00000000\tnot-imputed\t0.00000000" in manifest_text
    assert list(sample_stats["population"]) == ["pop1", "pop1", "pop2", "pop2"]
    assert np.allclose(sample_stats["sites_total"], 1)
    assert np.allclose(sample_stats["sites_called"], [1, 1, 1, 1])
    assert np.allclose(sample_stats["called_fraction"], [1.0, 1.0, 1.0, 1.0])
    assert np.allclose(sample_stats["sites_missing"], [0, 0, 0, 0])
    assert np.allclose(sample_stats["missing_fraction"], [0.0, 0.0, 0.0, 0.0])
    assert np.allclose(sample_stats["homozygous_sites"], [1, 0, 1, 1])
    assert np.allclose(sample_stats["heterozygous_sites"], [0, 1, 0, 0])
    assert np.allclose(sample_stats["observed_heterozygosity"], [0.0, 1.0, 0.0, 0.0])
    assert int(global_stats.iloc[0]["sites_used_heterozygosity"]) == 1
    assert np.isclose(float(global_stats.iloc[0]["observed_heterozygosity"]), 0.25)
    assert np.isclose(
        float(global_stats.iloc[0]["fit"]),
        1.0
        - (
            float(global_stats.iloc[0]["observed_heterozygosity"])
            / float(global_stats.iloc[0]["expected_heterozygosity_total"])
        ),
    )
    assert "1.00000000\t0\t0.00000000\t0\t1\t1.00000000" in sample_stats_text
    assert "pi" not in population.columns
    assert {"observed_heterozygosity", "expected_heterozygosity", "fis"}.issubset(population.columns)
    assert list(pairwise.columns) == ["population1", "population2", "sites_used", "fst"]
    assert float(pairwise.iloc[0]["fst"]) > 0.5
    assert np.isnan(float(population.loc[population["population"] == "pop2", "fis"].iloc[0]))
    assert not sfs.empty


def test_run_popgen_sequence_sample_stats_counts_N_and_gap_calls_as_missing(
    tmp_path: Path,
) -> None:
    h5 = _write_missing_sequence_popgen_h5(tmp_path / "assembly.hdf5")
    imap = _write_imap(tmp_path / "imap.tsv")

    run_popgen_method(
        data=h5,
        name="pop",
        outdir=tmp_path / "OUT",
        stats="heterozygosity",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=None,
        exclude=None,
        include_reference=False,
        subsample_unlinked=False,
        random_seed=None,
        cores=1,
        force=True,
        log_level="INFO",
    )

    outdir = tmp_path / "OUT"
    sample_stats = pd.read_csv(outdir / "pop.sample_stats.tsv", sep="\t")
    sample_stats_text = (outdir / "pop.sample_stats.tsv").read_text(encoding="utf-8")
    manifest_text = (outdir / "pop.manifest.txt").read_text(encoding="utf-8")

    assert not (outdir / "pop.global_stats.tsv").exists()
    assert "s1\t0.50000000\t0.50000000\tnot-imputed\t0.00000000" in manifest_text
    assert "s2\t0.25000000\t0.25000000\tnot-imputed\t0.00000000" in manifest_text
    assert np.allclose(sample_stats["sites_total"], [4, 4, 4, 4])
    assert np.allclose(sample_stats["sites_called"], [2, 3, 4, 3])
    assert np.allclose(sample_stats["sites_missing"], [2, 1, 0, 1])
    assert np.allclose(sample_stats["called_fraction"], [0.5, 0.75, 1.0, 0.75])
    assert np.allclose(sample_stats["missing_fraction"], [0.5, 0.25, 0.0, 0.25])
    assert np.allclose(sample_stats["homozygous_sites"], [2, 2, 4, 3])
    assert np.allclose(sample_stats["heterozygous_sites"], [0, 1, 0, 0])
    assert np.allclose(sample_stats["observed_heterozygosity"], [0.0, 1.0 / 3.0, 0.0, 0.0])
    assert "0.50000000" in sample_stats_text
    assert "0.33333333" in sample_stats_text


def test_run_popgen_rejects_sequence_only_stats_on_snp_only_hdf5(tmp_path: Path) -> None:
    h5 = _write_combined_popgen_h5(tmp_path / "snps.hdf5", include_phy=False, include_genos=True)

    with pytest.raises(IPyradError, match="require sequence HDF5"):
        run_popgen_method(
            data=h5,
            name="pop",
            outdir=tmp_path / "OUT",
            stats="pi",
            min_sample_coverage=1,
            max_sample_missing=1.0,
            min_minor_allele_frequency=0.0,
            imap=None,
            minmap=None,
            exclude=None,
            include_reference=False,
            subsample_unlinked=False,
            random_seed=None,
            cores=1,
            force=True,
            log_level="INFO",
        )


def test_run_popgen_without_imap_uses_single_all_population(tmp_path: Path) -> None:
    h5 = _write_combined_popgen_h5(tmp_path / "snps.hdf5", include_phy=False, include_genos=True)

    run_popgen_method(
        data=h5,
        name="pop",
        outdir=tmp_path / "OUT",
        stats="heterozygosity,sfs",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        subsample_unlinked=False,
        random_seed=None,
        cores=1,
        force=True,
        log_level="INFO",
    )

    population = pd.read_csv(tmp_path / "OUT" / "pop.population_stats.tsv", sep="\t")
    manifest_text = (tmp_path / "OUT" / "pop.manifest.txt").read_text(encoding="utf-8")
    sample_stats = pd.read_csv(tmp_path / "OUT" / "pop.sample_stats.tsv", sep="\t")

    assert list(population["population"]) == ["all"]
    assert "samples_final: ['s1', 's2', 's3', 's4']" in manifest_text
    assert sample_stats["population"].eq("all").all()
    assert not (tmp_path / "OUT" / "pop.global_stats.tsv").exists()
    assert not (tmp_path / "OUT" / "pop.pairwise_stats.tsv").exists()


def test_run_popgen_sequence_backend_writes_genomic_windows_and_keeps_genomewide_outputs(
    tmp_path: Path,
) -> None:
    h5 = _write_windowed_sequence_popgen_h5(tmp_path / "assembly.hdf5")
    imap = _write_imap(tmp_path / "imap.tsv")

    run_popgen_method(
        data=h5,
        name="pop",
        outdir=tmp_path / "OUT",
        stats="pi,fst,fis,sfs",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=None,
        exclude=None,
        include_reference=False,
        subsample_unlinked=False,
        random_seed=None,
        window_size=4,
        step_size=4,
        cores=1,
        force=True,
        log_level="INFO",
    )

    manifest_text = (tmp_path / "OUT" / "pop.manifest.txt").read_text(encoding="utf-8")
    genomewide = pd.read_csv(tmp_path / "OUT" / "pop.population_stats.tsv", sep="\t")
    sample_stats = pd.read_csv(tmp_path / "OUT" / "pop.sample_stats.tsv", sep="\t")
    windowed = pd.read_csv(tmp_path / "OUT" / "pop.window_population_stats.tsv", sep="\t")
    window_pairwise = pd.read_csv(tmp_path / "OUT" / "pop.window_pairwise_stats.tsv", sep="\t")

    assert "window_mode: genomic" in manifest_text
    assert "windows_planned: 3" in manifest_text
    assert "windows_written: 3" in manifest_text
    assert "window_sfs_note: Windowed SFS is not written in this phase." in manifest_text
    assert not (tmp_path / "OUT" / "pop.global_stats.tsv").exists()
    assert (tmp_path / "OUT" / "pop.sfs.tsv").exists()
    assert list(sample_stats.columns) == [
        "sample",
        "population",
        "sites_total",
        "sites_called",
        "called_fraction",
        "sites_missing",
        "missing_fraction",
        "homozygous_sites",
        "heterozygous_sites",
        "observed_heterozygosity",
    ]
    assert not windowed.empty
    assert not window_pairwise.empty
    assert set(genomewide["population"]) == {"pop1", "pop2"}
    assert list(windowed.columns[:9]) == [
        "window_id",
        "window_mode",
        "scaffold",
        "start",
        "end",
        "first_locus",
        "last_locus",
        "nloci",
        "sites_total",
    ]
    assert set(windowed["window_mode"]) == {"genomic"}
    assert set(windowed["scaffold"]) == {"chr1", "chr2"}
    assert set(zip(windowed["scaffold"], windowed["start"], windowed["end"])) == {
        ("chr1", 1, 4),
        ("chr1", 5, 8),
        ("chr2", 1, 4),
    }
    assert "fis" in windowed.columns
    assert set(window_pairwise["population1"]) == {"pop1"}
    assert set(window_pairwise["population2"]) == {"pop2"}


def test_run_popgen_sequence_backend_writes_locus_windows(tmp_path: Path) -> None:
    h5 = _write_windowed_sequence_popgen_h5(tmp_path / "assembly.hdf5")
    imap = _write_imap(tmp_path / "imap.tsv")

    run_popgen_method(
        data=h5,
        name="pop",
        outdir=tmp_path / "OUT",
        stats="pi,fst",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=None,
        exclude=None,
        include_reference=False,
        subsample_unlinked=False,
        random_seed=None,
        loci_per_window=2,
        locus_step=1,
        cores=1,
        force=True,
        log_level="INFO",
    )

    windowed = pd.read_csv(tmp_path / "OUT" / "pop.window_population_stats.tsv", sep="\t")
    pairwise = pd.read_csv(tmp_path / "OUT" / "pop.window_pairwise_stats.tsv", sep="\t")

    assert set(windowed["window_mode"]) == {"locus"}
    assert set(windowed["window_id"]) == {1, 2, 3}
    second_window = windowed.loc[windowed["window_id"] == 2].iloc[0]
    assert second_window["scaffold"] == "multiple"
    assert int(second_window["first_locus"]) == 2
    assert int(second_window["last_locus"]) == 3
    assert int(second_window["nloci"]) == 2
    assert set(pairwise["window_id"]) == {1, 2, 3}


def test_run_popgen_rejects_windowing_on_snp_only_hdf5(tmp_path: Path) -> None:
    h5 = _write_combined_popgen_h5(tmp_path / "snps.hdf5", include_phy=False, include_genos=True)

    with pytest.raises(IPyradError, match="Windowed popgen statistics currently require sequence HDF5"):
        run_popgen_method(
            data=h5,
            name="pop",
            outdir=tmp_path / "OUT",
            stats="fst",
            min_sample_coverage=1,
            max_sample_missing=1.0,
            min_minor_allele_frequency=0.0,
            imap=None,
            minmap=None,
            exclude=None,
            include_reference=False,
            subsample_unlinked=False,
            random_seed=None,
            window_size=2,
            cores=1,
            force=True,
            log_level="INFO",
        )


def test_run_popgen_windowing_with_sfs_only_writes_no_window_tables(tmp_path: Path) -> None:
    h5 = _write_windowed_sequence_popgen_h5(tmp_path / "assembly.hdf5")

    run_popgen_method(
        data=h5,
        name="pop",
        outdir=tmp_path / "OUT",
        stats="sfs",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        subsample_unlinked=False,
        random_seed=None,
        window_size=4,
        cores=1,
        force=True,
        log_level="INFO",
    )

    manifest_text = (tmp_path / "OUT" / "pop.manifest.txt").read_text(encoding="utf-8")
    assert "window_mode: genomic" in manifest_text
    assert "window_sfs_note: Windowed SFS is not written in this phase." in manifest_text
    assert not (tmp_path / "OUT" / "pop.global_stats.tsv").exists()
    assert (tmp_path / "OUT" / "pop.sfs.tsv").exists()
    assert (tmp_path / "OUT" / "pop.sample_stats.tsv").exists()
    assert not (tmp_path / "OUT" / "pop.window_population_stats.tsv").exists()
    assert not (tmp_path / "OUT" / "pop.window_pairwise_stats.tsv").exists()


def test_vectorized_sequence_and_genotype_summaries_match_sitewise_helpers() -> None:
    seq_block = np.array(
        [
            list(map(ord, "ARGN")),
            list(map(ord, "AGGN")),
            list(map(ord, "GGAN")),
        ],
        dtype=np.uint8,
    )
    seq_summary = summarize_sequence_block(seq_block, include_minor_allele_count=True)
    for idx in range(seq_block.shape[1]):
        site = summarize_sequence_site(seq_block[:, idx])
        assert int(seq_summary.called_samples[idx]) == site.called_samples
        assert int(seq_summary.chromosome_count[idx]) == site.chromosome_count
        assert np.array_equal(seq_summary.allele_counts[idx], site.allele_counts)
        assert seq_summary.segregating[idx] == site.segregating
        assert seq_summary.biallelic[idx] == site.biallelic
        if np.isnan(site.pi):
            assert np.isnan(seq_summary.pi[idx])
        else:
            assert np.isclose(seq_summary.pi[idx], site.pi)
        if np.isnan(site.observed_heterozygosity):
            assert np.isnan(seq_summary.observed_heterozygosity[idx])
        else:
            assert np.isclose(
                seq_summary.observed_heterozygosity[idx],
                site.observed_heterozygosity,
            )
        if np.isnan(site.expected_heterozygosity):
            assert np.isnan(seq_summary.expected_heterozygosity[idx])
        else:
            assert np.isclose(
                seq_summary.expected_heterozygosity[idx],
                site.expected_heterozygosity,
            )
        expected_minor = 0 if site.minor_allele_count is None else site.minor_allele_count
        assert int(seq_summary.minor_allele_count[idx]) == expected_minor

    geno_block = np.array(
        [
            [0, 1, 2, 255],
            [0, 1, 2, 255],
            [0, 2, 2, 255],
        ],
        dtype=np.uint8,
    )
    geno_summary = summarize_genotype_block(geno_block, include_minor_allele_count=True)
    for idx in range(geno_block.shape[1]):
        site = summarize_genotype_site(geno_block[:, idx])
        assert int(geno_summary.called_samples[idx]) == site.called_samples
        assert int(geno_summary.chromosome_count[idx]) == site.chromosome_count
        assert np.array_equal(geno_summary.allele_counts[idx], site.allele_counts)
        assert geno_summary.segregating[idx] == site.segregating
        assert geno_summary.biallelic[idx] == site.biallelic
        if np.isnan(site.pi):
            assert np.isnan(geno_summary.pi[idx])
        else:
            assert np.isclose(geno_summary.pi[idx], site.pi)
        expected_minor = 0 if site.minor_allele_count is None else site.minor_allele_count
        assert int(geno_summary.minor_allele_count[idx]) == expected_minor


def test_run_popgen_sequence_backend_multicore_matches_serial(tmp_path: Path) -> None:
    h5 = _write_large_sequence_popgen_h5(tmp_path / "assembly.hdf5")
    imap = _write_imap(tmp_path / "imap.tsv")

    run_popgen_method(
        data=h5,
        name="pop",
        outdir=tmp_path / "OUT1",
        stats="all",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=None,
        exclude=None,
        include_reference=False,
        subsample_unlinked=False,
        random_seed=None,
        cores=1,
        force=True,
        log_level="INFO",
    )
    run_popgen_method(
        data=h5,
        name="pop",
        outdir=tmp_path / "OUT2",
        stats="all",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=None,
        exclude=None,
        include_reference=False,
        subsample_unlinked=False,
        random_seed=None,
        cores=2,
        force=True,
        log_level="INFO",
    )

    for suffix in (
        "sample_stats.tsv",
        "global_stats.tsv",
        "population_stats.tsv",
        "pairwise_stats.tsv",
        "sfs.tsv",
    ):
        left = pd.read_csv(tmp_path / "OUT1" / f"pop.{suffix}", sep="\t")
        right = pd.read_csv(tmp_path / "OUT2" / f"pop.{suffix}", sep="\t")
        pd.testing.assert_frame_equal(left, right)


def test_run_popgen_sequence_no_window_reads_chunks_once_without_missing_prefilter(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    h5 = _write_combined_popgen_h5(tmp_path / "assembly.hdf5", include_phy=True, include_genos=False)
    imap = _write_imap(tmp_path / "imap.tsv")
    original = popgen_seq_backend._load_sequence_chunk_from_phy
    calls = {"count": 0}

    def _counting_loader(*args, **kwargs):
        calls["count"] += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(popgen_seq_backend, "_load_sequence_chunk_from_phy", _counting_loader)

    run_popgen_method(
        data=h5,
        name="pop",
        outdir=tmp_path / "OUT",
        stats="pi,fst",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=None,
        exclude=None,
        include_reference=False,
        subsample_unlinked=False,
        random_seed=None,
        cores=1,
        force=True,
        log_level="INFO",
    )

    assert calls["count"] == 1


def test_run_popgen_sequence_windowed_reads_chunks_once_without_missing_prefilter(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    h5 = _write_windowed_sequence_popgen_h5(tmp_path / "assembly.hdf5")
    imap = _write_imap(tmp_path / "imap.tsv")
    original = popgen_seq_backend._load_sequence_chunk_from_phy
    calls = {"count": 0}

    def _counting_loader(*args, **kwargs):
        calls["count"] += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(popgen_seq_backend, "_load_sequence_chunk_from_phy", _counting_loader)

    run_popgen_method(
        data=h5,
        name="pop",
        outdir=tmp_path / "OUT",
        stats="pi,fst",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=None,
        exclude=None,
        include_reference=False,
        subsample_unlinked=False,
        random_seed=None,
        window_size=4,
        step_size=4,
        cores=1,
        force=True,
        log_level="INFO",
    )

    assert calls["count"] == 1


def test_run_popgen_sequence_windowed_with_missing_prefilter_reads_chunks_twice(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    h5 = _write_windowed_sequence_popgen_h5(tmp_path / "assembly.hdf5")
    imap = _write_imap(tmp_path / "imap.tsv")
    original = popgen_seq_backend._load_sequence_chunk_from_phy
    calls = {"count": 0}

    def _counting_loader(*args, **kwargs):
        calls["count"] += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(popgen_seq_backend, "_load_sequence_chunk_from_phy", _counting_loader)

    run_popgen_method(
        data=h5,
        name="pop",
        outdir=tmp_path / "OUT",
        stats="pi,fst",
        min_sample_coverage=1,
        max_sample_missing=0.9,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=None,
        exclude=None,
        include_reference=False,
        subsample_unlinked=False,
        random_seed=None,
        window_size=4,
        step_size=4,
        cores=1,
        force=True,
        log_level="INFO",
    )

    assert calls["count"] == 2
