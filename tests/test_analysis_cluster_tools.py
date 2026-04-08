from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from loguru import logger

from ipyrad2.analysis.methods.dapc import run_dapc_method
from ipyrad2.analysis.methods.snmf import run_snmf_method


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


def test_run_snmf_method_writes_membership_assignment_and_stats(tmp_path: Path) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")

    run_snmf_method(
        data=h5,
        name="snmf",
        outdir=tmp_path / "OUT",
        k=2,
        k_range=None,
        min_sample_coverage=2,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        impute_method="sample",
        subsample=True,
        random_seed=3,
        cores=1,
        force=True,
        log_level="INFO",
    )

    membership = pd.read_csv(tmp_path / "OUT" / "snmf.membership.tsv", sep="\t")
    allele_freqs = pd.read_csv(tmp_path / "OUT" / "snmf.allele_frequencies.tsv", sep="\t")
    assignments = pd.read_csv(tmp_path / "OUT" / "snmf.assignments.tsv", sep="\t")
    k_scan = pd.read_csv(tmp_path / "OUT" / "snmf.k_scan.tsv", sep="\t")
    sample_summary = pd.read_csv(
        tmp_path / "OUT" / "snmf.sample_data_summary.tsv",
        sep="\t",
    )
    stats = (tmp_path / "OUT" / "snmf.stats.txt").read_text(encoding="utf-8")

    assert {"sample", "cluster1", "cluster2"} == set(membership.columns)
    assert {"marker_id", "cluster1", "cluster2"} == set(allele_freqs.columns)
    assert {"sample", "assigned_cluster", "assignment_score"} == set(assignments.columns)
    assert {
        "k",
        "mean_cross_entropy",
        "sd_cross_entropy",
        "best_reconstruction_err",
        "best_n_iter",
        "selected",
    } == set(k_scan.columns)
    assert {
        "sample",
        "missing_fraction",
        "post_imputation_missing_fraction",
        "imputation_algorithm",
        "imputed_genotype_fraction",
    }.issubset(sample_summary.columns)
    assert sample_summary["imputation_algorithm"].eq("sample").all()
    assert "tool: snmf" in stats
    assert "impute_method: sample" in stats
    assert "selected_cross_entropy:" in stats
    assert "alpha_W: 0.0001" in stats
    assert "alpha_H: same" in stats
    assert "l1_ratio: 1.0" in stats
    assert "n_init: 10" in stats
    assert "cv_replicates: 5" in stats
    assert "cv_holdout: 0.1" in stats
    assert "nmf_tol: 0.001" in stats
    assert "nmf_max_iter: 3000" in stats
    assert "selected_fit_hit_max_iter:" in stats
    assert "capped_fits_for_selected_k:" in stats
    assert "total_fits_for_selected_k:" in stats
    assert "exported_snps: 3" in stats


def test_run_snmf_method_supports_k_scan_and_none_imputation(tmp_path: Path) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")

    run_snmf_method(
        data=h5,
        name="scan",
        outdir=tmp_path / "OUT",
        k=None,
        k_range="2:3",
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
        cores=1,
        force=True,
        log_level="INFO",
    )

    k_scan = pd.read_csv(tmp_path / "OUT" / "scan.k_scan.tsv", sep="\t")
    stats = (tmp_path / "OUT" / "scan.stats.txt").read_text(encoding="utf-8")
    assert len(k_scan) == 2
    assert k_scan["selected"].sum() == 1
    assert "impute_method: None" in stats
    assert "subsample: False" in stats
    assert "selected_cross_entropy:" in stats


def test_run_snmf_method_logs_output_directory_once(tmp_path: Path) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")
    messages: list[str] = []
    sink_id = logger.add(messages.append, format="{message}", level="INFO")
    try:
        run_snmf_method(
            data=h5,
            name="logged",
            outdir=tmp_path / "OUT",
            k=2,
            k_range=None,
            min_sample_coverage=2,
            max_sample_missing=1.0,
            min_minor_allele_frequency=0.0,
            imap=None,
            minmap=None,
            exclude=None,
            include_reference=False,
            impute_method="sample",
            subsample=True,
            random_seed=3,
            cores=1,
            force=True,
            log_level="INFO",
        )
    finally:
        logger.remove(sink_id)

    extraction_idx = next(
        idx for idx, msg in enumerate(messages) if "SNP extraction summary" in msg
    )
    imputation_idx = next(
        idx for idx, msg in enumerate(messages) if "snmf SNP imputation:" in msg
    )
    prepared_idx = next(
        idx for idx, msg in enumerate(messages) if "snmf prepared SNP summary:" in msg
    )
    selected_idx = next(
        idx for idx, msg in enumerate(messages) if "selected sNMF K=" in msg
    )

    assert extraction_idx < imputation_idx < prepared_idx < selected_idx
    assert not any("subsampled " in msg for msg in messages)
    assert any(
        msg.strip() == f"wrote sNMF outputs to {(tmp_path / 'OUT').resolve()}"
        for msg in messages
    )


def test_run_dapc_method_writes_expected_outputs(tmp_path: Path) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")

    run_dapc_method(
        data=h5,
        name="dapc",
        outdir=tmp_path / "OUT",
        k=2,
        k_range=None,
        n_pcs=None,
        min_sample_coverage=2,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        impute_method="sample",
        subsample=True,
        random_seed=9,
        cores=1,
        force=True,
        log_level="INFO",
    )

    coords = pd.read_csv(tmp_path / "OUT" / "dapc.coords.tsv", sep="\t")
    membership = pd.read_csv(tmp_path / "OUT" / "dapc.membership.tsv", sep="\t")
    assignments = pd.read_csv(tmp_path / "OUT" / "dapc.assignments.tsv", sep="\t")
    k_scan = pd.read_csv(tmp_path / "OUT" / "dapc.k_scan.tsv", sep="\t")
    sample_summary = pd.read_csv(
        tmp_path / "OUT" / "dapc.sample_data_summary.tsv",
        sep="\t",
    )
    stats = (tmp_path / "OUT" / "dapc.stats.txt").read_text(encoding="utf-8")

    assert "sample" in coords.columns
    assert "axis1" in coords.columns
    assert {"sample", "cluster1", "cluster2"} == set(membership.columns)
    assert {"sample", "assigned_cluster", "assignment_score"} == set(assignments.columns)
    assert {"k", "bic", "selected"} == set(k_scan.columns)
    assert {
        "sample",
        "missing_fraction",
        "post_imputation_missing_fraction",
        "imputation_algorithm",
        "imputed_genotype_fraction",
    }.issubset(sample_summary.columns)
    assert sample_summary["imputation_algorithm"].eq("sample").all()
    assert "tool: dapc" in stats
    assert "n_pcs:" in stats


def test_run_dapc_method_supports_k_scan(tmp_path: Path) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")

    run_dapc_method(
        data=h5,
        name="scan",
        outdir=tmp_path / "OUT",
        k=None,
        k_range="2:3",
        n_pcs=3,
        min_sample_coverage=2,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        impute_method="none",
        subsample=False,
        random_seed=4,
        cores=1,
        force=True,
        log_level="INFO",
    )

    k_scan = pd.read_csv(tmp_path / "OUT" / "scan.k_scan.tsv", sep="\t")
    stats = (tmp_path / "OUT" / "scan.stats.txt").read_text(encoding="utf-8")
    assert k_scan["selected"].sum() == 1
    assert "impute_method: None" in stats
