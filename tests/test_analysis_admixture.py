from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest

import ipyrad2.analysis.methods.admixture as admixture_mod
from ipyrad2.analysis.extracters.snps_extracter import SNPsExtracter
from ipyrad2.utils.exceptions import IPyradError


def _write_phase2_snps_h5(path: Path) -> Path:
    string_dtype = h5py.string_dtype(encoding="utf-8")
    names = ["a1", "a2", "a3", "b1", "b2", "b3"]
    doses = np.array(
        [
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 255, 0],
            [2, 2, 2, 2, 2, 1],
            [2, 2, 2, 1, 2, 2],
            [2, 2, 2, 2, 255, 2],
        ],
        dtype=np.uint8,
    )
    genos = np.zeros((doses.shape[0], doses.shape[1], 3), dtype=np.uint8)
    for row in range(doses.shape[0]):
        for col in range(doses.shape[1]):
            dose = int(doses[row, col])
            if dose == 255:
                genos[row, col, :2] = 255
                genos[row, col, 2] = ord("N")
            elif dose == 0:
                genos[row, col, :2] = [0, 0]
                genos[row, col, 2] = ord("A")
            elif dose == 1:
                genos[row, col, :2] = [0, 1]
                genos[row, col, 2] = ord("R")
            else:
                genos[row, col, :2] = [1, 1]
                genos[row, col, 2] = ord("G")
    snpsmap = np.array(
        [
            [0, 0, 0, 0, 100],
            [0, 1, 1, 0, 120],
            [1, 0, 0, 0, 200],
            [1, 1, 1, 0, 220],
            [2, 0, 0, 1, 50],
            [2, 1, 1, 1, 70],
        ],
        dtype=np.uint32,
    )

    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(names, dtype=string_dtype)
        io5.attrs["nsnps"] = int(snpsmap.shape[0])
        io5.create_dataset("genos", data=genos)
        io5.create_dataset(
            "reference",
            data=np.array([ord("A")] * snpsmap.shape[0], dtype=np.uint8),
        )
        ds = io5.create_dataset("snpsmap", data=snpsmap)
        ds.attrs["columns"] = np.array(
            ["loc", "loc_idx", "loc_pos", "scaff", "pos"],
            dtype=string_dtype,
        )
    return path


def _mock_admixture_run(cmd, cwd, capture_output, text, check):
    stage_dir = Path(cwd)
    bed_name = Path(cmd[-2]).name
    stem = Path(bed_name).stem
    k = int(cmd[-1])

    nsamples = sum(1 for _ in open(stage_dir / f"{stem}.fam", encoding="utf-8"))
    nsnps = sum(1 for _ in open(stage_dir / f"{stem}.bim", encoding="utf-8"))

    if k == 2:
        q = np.array(
            [
                [0.95, 0.05],
                [0.90, 0.10],
                [0.88, 0.12],
                [0.10, 0.90],
                [0.08, 0.92],
                [0.05, 0.95],
            ],
            dtype=float,
        )
        p = np.vstack(
            [
                np.linspace(0.1, 0.3, nsnps),
                np.linspace(0.9, 0.7, nsnps),
            ]
        )
    else:
        q = np.array(
            [
                [0.80, 0.10, 0.10],
                [0.78, 0.12, 0.10],
                [0.30, 0.55, 0.15],
                [0.10, 0.15, 0.75],
                [0.10, 0.20, 0.70],
                [0.12, 0.18, 0.70],
            ],
            dtype=float,
        )
        p = np.vstack(
            [
                np.linspace(0.1, 0.2, nsnps),
                np.linspace(0.4, 0.5, nsnps),
                np.linspace(0.8, 0.9, nsnps),
            ]
        )
    assert q.shape[0] == nsamples
    np.savetxt(stage_dir / f"{stem}.{k}.Q", q, fmt="%.6f")
    np.savetxt(stage_dir / f"{stem}.{k}.P", p, fmt="%.6f")

    stdout = ""
    if "--cv=5" in cmd:
        cv = {2: 0.33, 3: 0.21}[k]
        stdout = f"CV error (K={k}): {cv}\n"

    return admixture_mod.subprocess.CompletedProcess(
        cmd,
        0,
        stdout=stdout,
        stderr="",
    )


def test_run_admixture_method_writes_curated_outputs_and_cleans_intermediates(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")
    observed = {}
    original_write_plink = SNPsExtracter.write_plink

    def _wrapped_write_plink(self, prefix, view, *, impute_method=None):
        observed["impute_method"] = impute_method
        return original_write_plink(self, prefix, view, impute_method=impute_method)

    monkeypatch.setattr(admixture_mod.shutil, "which", lambda name: "/usr/bin/admixture")
    monkeypatch.setattr(admixture_mod.subprocess, "run", _mock_admixture_run)
    monkeypatch.setattr(SNPsExtracter, "write_plink", _wrapped_write_plink)

    admixture_mod.run_admixture_method(
        data=h5,
        name="admix",
        outdir=tmp_path / "OUT",
        k=2,
        k_range=None,
        binary=None,
        min_sample_coverage=2,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        impute_method="sample",
        subsample=True,
        random_seed=7,
        keep_intermediates=False,
        cores=2,
        force=True,
        log_level="INFO",
    )

    membership = pd.read_csv(tmp_path / "OUT" / "admix.membership.tsv", sep="\t")
    allele_freqs = pd.read_csv(tmp_path / "OUT" / "admix.allele_frequencies.tsv", sep="\t")
    assignments = pd.read_csv(tmp_path / "OUT" / "admix.assignments.tsv", sep="\t")
    k_scan = pd.read_csv(tmp_path / "OUT" / "admix.k_scan.tsv", sep="\t")
    sample_summary = pd.read_csv(
        tmp_path / "OUT" / "admix.sample_data_summary.tsv",
        sep="\t",
    )
    stats = (tmp_path / "OUT" / "admix.stats.txt").read_text(encoding="utf-8")

    assert observed["impute_method"] == "sample"
    assert {"sample", "cluster1", "cluster2"} == set(membership.columns)
    assert {"marker_id", "cluster1", "cluster2"} == set(allele_freqs.columns)
    assert {"sample", "assigned_cluster", "assignment_score"} == set(assignments.columns)
    assert {"k", "cv_error", "selected"} == set(k_scan.columns)
    assert {
        "sample",
        "missing_fraction",
        "post_imputation_missing_fraction",
        "imputation_algorithm",
        "imputed_genotype_fraction",
    }.issubset(sample_summary.columns)
    assert sample_summary["imputation_algorithm"].eq("sample").all()
    assert np.allclose(sample_summary["post_imputation_missing_fraction"], 0.0)
    assert np.allclose(
        sample_summary["imputed_genotype_fraction"],
        sample_summary["missing_fraction"],
    )
    assert not (tmp_path / "OUT" / "admix.sample_missing.tsv").exists()
    assert "tool: admixture" in stats
    assert "impute_method: sample" in stats
    assert "keep_intermediates: False" in stats
    assert "k_selected: 2" in stats
    assert not (tmp_path / "OUT" / "admix.intermediates").exists()


def test_run_admixture_method_scans_k_and_keeps_requested_intermediates(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")
    observed = {}
    original_write_plink = SNPsExtracter.write_plink

    def _wrapped_write_plink(self, prefix, view, *, impute_method=None):
        observed["impute_method"] = impute_method
        return original_write_plink(self, prefix, view, impute_method=impute_method)

    monkeypatch.setattr(admixture_mod.shutil, "which", lambda name: "/usr/bin/admixture")
    monkeypatch.setattr(admixture_mod.subprocess, "run", _mock_admixture_run)
    monkeypatch.setattr(SNPsExtracter, "write_plink", _wrapped_write_plink)

    admixture_mod.run_admixture_method(
        data=h5,
        name="admixscan",
        outdir=tmp_path / "OUT",
        k=None,
        k_range="2:3",
        binary=None,
        min_sample_coverage=2,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        impute_method="none",
        subsample=False,
        random_seed=11,
        keep_intermediates=True,
        cores=1,
        force=True,
        log_level="INFO",
    )

    membership = pd.read_csv(tmp_path / "OUT" / "admixscan.membership.tsv", sep="\t")
    k_scan = pd.read_csv(tmp_path / "OUT" / "admixscan.k_scan.tsv", sep="\t")
    stats = (tmp_path / "OUT" / "admixscan.stats.txt").read_text(encoding="utf-8")
    stage_dir = tmp_path / "OUT" / "admixscan.intermediates"

    assert observed["impute_method"] is None
    assert {"sample", "cluster1", "cluster2", "cluster3"} == set(membership.columns)
    assert k_scan["selected"].sum() == 1
    assert int(k_scan.loc[k_scan["selected"], "k"].iloc[0]) == 3
    assert "selected_cv_error: 0.21" in stats
    assert "impute_method: None" in stats
    assert "keep_intermediates: True" in stats
    assert stage_dir.exists()
    assert (stage_dir / "admixscan.bed").exists()
    assert (stage_dir / "admixscan.bim").exists()
    assert (stage_dir / "admixscan.fam").exists()
    assert (stage_dir / "admixscan.3.P").exists()
    assert (stage_dir / "admixscan.3.Q").exists()
    assert (stage_dir / "admixscan.3.log").exists()


def test_run_admixture_method_rejects_vcf_and_missing_binary(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    with pytest.raises(IPyradError, match="requires an SNP-capable HDF5 input"):
        admixture_mod.run_admixture_method(
            data=tmp_path / "variants.vcf.gz",
            name="bad",
            outdir=tmp_path / "OUT",
            k=2,
            k_range=None,
            binary=None,
            min_sample_coverage=2,
            max_sample_missing=1.0,
            min_minor_allele_frequency=0.0,
            imap=None,
            minmap=None,
            exclude=None,
            include_reference=False,
            impute_method="sample",
            subsample=True,
            random_seed=1,
            keep_intermediates=False,
            cores=1,
            force=True,
            log_level="INFO",
        )

    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")
    monkeypatch.setattr(admixture_mod.shutil, "which", lambda name: None)
    with pytest.raises(IPyradError, match="Could not find the `admixture` binary"):
        admixture_mod.run_admixture_method(
            data=h5,
            name="bad",
            outdir=tmp_path / "OUT2",
            k=2,
            k_range=None,
            binary=None,
            min_sample_coverage=2,
            max_sample_missing=1.0,
            min_minor_allele_frequency=0.0,
            imap=None,
            minmap=None,
            exclude=None,
            include_reference=False,
            impute_method="sample",
            subsample=True,
            random_seed=1,
            keep_intermediates=False,
            cores=1,
            force=True,
            log_level="INFO",
        )
