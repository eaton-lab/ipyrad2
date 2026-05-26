from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pandas as pd
import pytest
from loguru import logger

from ipyrad2.analysis.methods import pca as pca_methods
from ipyrad2.analysis.methods.common import NumericalInput
from ipyrad2.analysis.methods.pca import run_pca_analysis, run_pca_method
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


def _write_assembly_style_phase2_snps_h5(path: Path) -> Path:
    """Write an SNP HDF5 that mimics assemble outputs with mismatched top-level names."""
    string_dtype = h5py.string_dtype(encoding="utf-8")
    names = ["a1", "a2", "b1", "b2"]
    doses = np.array(
        [
            [0, 0, 0, 1],
            [0, 0, 1, 0],
            [2, 2, 2, 1],
            [2, 2, 1, 2],
        ],
        dtype=np.uint8,
    )
    genos = np.zeros((doses.shape[0], doses.shape[1], 3), dtype=np.uint8)
    for row in range(doses.shape[0]):
        for col in range(doses.shape[1]):
            dose = int(doses[row, col])
            if dose == 0:
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
        ],
        dtype=np.uint32,
    )

    with h5py.File(path, "w") as io5:
        io5.attrs["version"] = 2.0
        io5.attrs["names"] = np.array(
            ["assembly_reference_sequence", *names],
            dtype=string_dtype,
        )
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


def _build_mock_umap_inputs() -> tuple[SimpleNamespace, dict[int, int], dict[int, NumericalInput]]:
    extracter = SimpleNamespace(
        snames=["a1", "a2"],
        sample_missing=pd.Series([0.0, 0.0], index=["a1", "a2"], dtype=float),
    )
    prepared = NumericalInput(
        extracter=extracter,
        view=None,
        matrix=np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.float64),
        imputation=None,
    )
    return extracter, {0: 12345}, {0: prepared}


def test_run_pca_method_writes_expected_outputs(tmp_path: Path) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")

    run_pca_method(
        data=h5,
        name="phase2",
        outdir=tmp_path / "OUT",
        method="pca",
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
        replicates=2,
        perplexity=5.0,
        max_iter=1000,
        n_neighbors=15,
        plot=False,
        cores=1,
        force=True,
        log_level="INFO",
    )

    coords = pd.read_csv(tmp_path / "OUT" / "phase2.coords.tsv", sep="\t")
    variance = pd.read_csv(tmp_path / "OUT" / "phase2.variance.tsv", sep="\t")
    sample_summary_path = tmp_path / "OUT" / "phase2.sample_data_summary.tsv"
    sample_summary = pd.read_csv(sample_summary_path, sep="\t")
    stats = (tmp_path / "OUT" / "phase2.stats.txt").read_text(encoding="utf-8")

    assert coords["method"].unique().tolist() == ["pca"]
    assert coords["replicate"].tolist().count(0) == 6
    assert coords["replicate"].tolist().count(1) == 6
    assert {"axis1", "axis2", "axis3"}.issubset(coords.columns)
    assert set(variance.columns) == {"replicate", "axis", "explained_variance_ratio"}
    assert {
        "sample",
        "population",
        "missing_fraction",
        "post_imputation_missing_fraction",
        "imputation_algorithm",
        "imputed_genotype_fraction",
    }.issubset(sample_summary.columns)
    assert sample_summary["population"].eq("all").all()
    assert sample_summary["imputation_algorithm"].eq("sample").all()
    assert np.allclose(sample_summary["post_imputation_missing_fraction"], 0.0)
    assert np.allclose(
        sample_summary["imputed_genotype_fraction"],
        sample_summary["missing_fraction"],
    )
    assert not (tmp_path / "OUT" / "phase2.plot.svg").exists()
    assert not (tmp_path / "OUT" / "phase2.sample_missing.tsv").exists()
    assert "tool: pca" in stats
    assert "method: pca" in stats
    assert "subsample: True" in stats
    assert "impute_method: sample" in stats
    assert "imputation_algorithm: sample" in stats
    assert "imputed_snp_fraction:" in stats
    assert "imputed_genotype_fraction:" in stats
    assert "samples_selected_initial_count: 6" in stats
    assert "samples_dropped_by_max_missing_count: 0" in stats
    assert "samples_final_count: 6" in stats
    assert "population_count: 1" in stats
    assert "Population sample counts" in stats
    assert "all: 6" in stats
    assert "samples_selected_initial: ['a1'" not in stats
    assert "samples_final: ['a1'" not in stats
    assert "imap: {'all'" not in stats
    assert "exported_snps: 3" in stats
    assert sample_summary_path.read_text(encoding="utf-8").count("0.000") > 0


def test_run_pca_method_writes_population_assignments_and_three_decimal_sample_summary(
    tmp_path: Path,
) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")
    imap = {
        "alpha": ["a1", "a2", "a3"],
        "beta": ["b1", "b2", "b3"],
    }

    run_pca_method(
        data=h5,
        name="phase2",
        outdir=tmp_path / "OUT",
        method="pca",
        min_sample_coverage=2,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=None,
        exclude=None,
        include_reference=False,
        impute_method="sample",
        subsample=False,
        random_seed=7,
        replicates=1,
        perplexity=5.0,
        max_iter=1000,
        n_neighbors=15,
        plot=False,
        cores=1,
        force=True,
        log_level="INFO",
    )

    sample_summary_path = tmp_path / "OUT" / "phase2.sample_data_summary.tsv"
    sample_summary = pd.read_csv(sample_summary_path, sep="\t")
    text = sample_summary_path.read_text(encoding="utf-8")

    assert sample_summary["population"].tolist() == [
        "alpha",
        "alpha",
        "alpha",
        "beta",
        "beta",
        "beta",
    ]
    assert "a3\talpha\t0.167\t0.000\tsample\t0.167" in text
    assert "b3\tbeta\t0.167\t0.000\tsample\t0.167" in text


def test_run_pca_method_stats_report_sample_counts_after_missing_filter(tmp_path: Path) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")
    imap = {
        "alpha": ["a1", "a2", "a3"],
        "beta": ["b1", "b2", "b3"],
    }

    run_pca_method(
        data=h5,
        name="phase2",
        outdir=tmp_path / "OUT",
        method="pca",
        min_sample_coverage=2,
        max_sample_missing=0.15,
        min_minor_allele_frequency=0.0,
        imap=imap,
        minmap=None,
        exclude=None,
        include_reference=False,
        impute_method="sample",
        subsample=False,
        random_seed=7,
        replicates=1,
        perplexity=5.0,
        max_iter=1000,
        n_neighbors=15,
        plot=False,
        cores=1,
        force=True,
        log_level="INFO",
    )

    sample_summary = pd.read_csv(
        tmp_path / "OUT" / "phase2.sample_data_summary.tsv",
        sep="\t",
    )
    stats = (tmp_path / "OUT" / "phase2.stats.txt").read_text(encoding="utf-8")

    assert sample_summary["sample"].tolist() == ["a1", "a2", "b1", "b2"]
    assert sample_summary["population"].tolist() == ["alpha", "alpha", "beta", "beta"]
    assert "samples_selected_initial_count: 6" in stats
    assert "samples_dropped_by_max_missing_count: 2" in stats
    assert "samples_final_count: 4" in stats
    assert "population_count: 2" in stats
    assert "alpha: 2" in stats
    assert "beta: 2" in stats
    assert "samples_selected_initial: ['a1'" not in stats
    assert "samples_final: ['a1'" not in stats


def test_run_pca_method_supports_tsne_and_umap_without_variance_file(tmp_path: Path) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")

    run_pca_method(
        data=h5,
        name="tsne",
        outdir=tmp_path / "TSNE",
        method="tsne",
        min_sample_coverage=2,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        impute_method="none",
        subsample=False,
        random_seed=5,
        replicates=1,
        perplexity=2.0,
        max_iter=250,
        n_neighbors=3,
        plot=False,
        cores=1,
        force=True,
        log_level="INFO",
    )
    run_pca_method(
        data=h5,
        name="umap",
        outdir=tmp_path / "UMAP",
        method="umap",
        min_sample_coverage=2,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        impute_method="sample",
        subsample=True,
        random_seed=5,
        replicates=1,
        perplexity=2.0,
        max_iter=250,
        n_neighbors=3,
        plot=False,
        cores=1,
        force=True,
        log_level="INFO",
    )

    tsne_coords = pd.read_csv(tmp_path / "TSNE" / "tsne.coords.tsv", sep="\t")
    umap_coords = pd.read_csv(tmp_path / "UMAP" / "umap.coords.tsv", sep="\t")
    tsne_stats = (tmp_path / "TSNE" / "tsne.stats.txt").read_text(encoding="utf-8")
    assert {"axis1", "axis2"}.issubset(tsne_coords.columns)
    assert {"axis1", "axis2"}.issubset(umap_coords.columns)
    assert not (tmp_path / "TSNE" / "tsne.variance.tsv").exists()
    assert not (tmp_path / "UMAP" / "umap.variance.tsv").exists()
    assert "impute_method: zero-fill" in tsne_stats
    assert "imputation_algorithm: zero-fill" in tsne_stats


def test_run_umap_analysis_uses_parallel_jobs_without_seed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    extracter, seeds, prepared_inputs = _build_mock_umap_inputs()
    calls: dict[str, object] = {}
    messages: list[str] = []

    monkeypatch.setattr(pca_methods, "_build_extracter", lambda **kwargs: extracter)

    def _fake_prepare_inputs(**kwargs):
        calls["prepare_random_seed"] = kwargs["random_seed"]
        return seeds, prepared_inputs

    def _fake_run_umap_once(
        matrix: np.ndarray,
        *,
        n_neighbors: int,
        embedding_random_state: int | None,
        n_jobs: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        calls["n_neighbors"] = n_neighbors
        calls["embedding_random_state"] = embedding_random_state
        calls["n_jobs"] = n_jobs
        return np.zeros((matrix.shape[0], 2), dtype=np.float64), np.array([], dtype=np.float64)

    monkeypatch.setattr(pca_methods, "_prepare_inputs", _fake_prepare_inputs)
    monkeypatch.setattr(pca_methods, "_run_umap_once", _fake_run_umap_once)
    sink_id = logger.add(messages.append, format="{message}", level="WARNING")
    try:
        result = pca_methods.run_umap_analysis(
            data=Path("dummy.hdf5"),
            random_seed=None,
            n_neighbors=3,
            cores=4,
        )
    finally:
        logger.remove(sink_id)

    assert result.method == "umap"
    assert calls["prepare_random_seed"] is None
    assert calls["n_neighbors"] == 3
    assert calls["embedding_random_state"] is None
    assert calls["n_jobs"] == 4
    assert not messages


def test_run_umap_analysis_honors_embedding_seed_when_serial(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    extracter, seeds, prepared_inputs = _build_mock_umap_inputs()
    calls: dict[str, object] = {}
    messages: list[str] = []

    monkeypatch.setattr(pca_methods, "_build_extracter", lambda **kwargs: extracter)

    def _fake_prepare_inputs(**kwargs):
        calls["prepare_random_seed"] = kwargs["random_seed"]
        return seeds, prepared_inputs

    def _fake_run_umap_once(
        matrix: np.ndarray,
        *,
        n_neighbors: int,
        embedding_random_state: int | None,
        n_jobs: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        calls["embedding_random_state"] = embedding_random_state
        calls["n_jobs"] = n_jobs
        return np.zeros((matrix.shape[0], 2), dtype=np.float64), np.array([], dtype=np.float64)

    monkeypatch.setattr(pca_methods, "_prepare_inputs", _fake_prepare_inputs)
    monkeypatch.setattr(pca_methods, "_run_umap_once", _fake_run_umap_once)
    sink_id = logger.add(messages.append, format="{message}", level="WARNING")
    try:
        pca_methods.run_umap_analysis(
            data=Path("dummy.hdf5"),
            random_seed=7,
            n_neighbors=3,
            cores=1,
        )
    finally:
        logger.remove(sink_id)

    assert calls["prepare_random_seed"] == 7
    assert calls["embedding_random_state"] == seeds[0]
    assert calls["n_jobs"] == 1
    assert not messages


def test_run_umap_analysis_warns_and_ignores_embedding_seed_when_parallel(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    extracter, seeds, prepared_inputs = _build_mock_umap_inputs()
    calls: dict[str, object] = {}
    messages: list[str] = []

    monkeypatch.setattr(pca_methods, "_build_extracter", lambda **kwargs: extracter)

    def _fake_prepare_inputs(**kwargs):
        calls["prepare_random_seed"] = kwargs["random_seed"]
        return seeds, prepared_inputs

    def _fake_run_umap_once(
        matrix: np.ndarray,
        *,
        n_neighbors: int,
        embedding_random_state: int | None,
        n_jobs: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        calls["embedding_random_state"] = embedding_random_state
        calls["n_jobs"] = n_jobs
        return np.zeros((matrix.shape[0], 2), dtype=np.float64), np.array([], dtype=np.float64)

    monkeypatch.setattr(pca_methods, "_prepare_inputs", _fake_prepare_inputs)
    monkeypatch.setattr(pca_methods, "_run_umap_once", _fake_run_umap_once)
    sink_id = logger.add(messages.append, format="{message}", level="WARNING")
    try:
        pca_methods.run_umap_analysis(
            data=Path("dummy.hdf5"),
            random_seed=7,
            n_neighbors=3,
            cores=4,
        )
    finally:
        logger.remove(sink_id)

    assert calls["prepare_random_seed"] == 7
    assert calls["embedding_random_state"] is None
    assert calls["n_jobs"] == 4
    assert any(
        "parallel UMAP does not support exact reproducibility" in msg
        and "Use --cores 1 for reproducible UMAP embeddings" in msg
        for msg in messages
    )


def test_run_pca_method_writes_svg_plot_when_requested(tmp_path: Path) -> None:
    pytest.importorskip("toyplot")
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")

    run_pca_method(
        data=h5,
        name="phase2",
        outdir=tmp_path / "OUT",
        method="pca",
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
        replicates=2,
        perplexity=5.0,
        max_iter=1000,
        n_neighbors=15,
        plot=True,
        cores=1,
        force=True,
        log_level="INFO",
    )

    plot = tmp_path / "OUT" / "phase2.plot.svg"
    assert plot.exists()
    plot_text = plot.read_text(encoding="utf-8")
    assert plot_text.lstrip().startswith("<svg")
    assert "<svg" in plot_text
    assert "<rect" in plot_text
    assert "stroke-width:2" in plot_text
    assert (tmp_path / "OUT" / "phase2.coords.tsv").exists()
    assert (tmp_path / "OUT" / "phase2.variance.tsv").exists()


def test_run_pca_method_accepts_assembly_style_hdf5_names(tmp_path: Path) -> None:
    h5 = _write_assembly_style_phase2_snps_h5(tmp_path / "assembly_snps.hdf5")

    run_pca_method(
        data=h5,
        name="assembly_phase2",
        outdir=tmp_path / "OUT",
        method="pca",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        impute_method="sample",
        subsample=False,
        random_seed=3,
        replicates=1,
        perplexity=5.0,
        max_iter=1000,
        n_neighbors=15,
        plot=False,
        cores=1,
        force=True,
        log_level="INFO",
    )

    coords = pd.read_csv(tmp_path / "OUT" / "assembly_phase2.coords.tsv", sep="\t")
    assert coords["sample"].tolist() == ["a1", "a2", "b1", "b2"]


def test_run_pca_method_include_reference_synthesizes_reference_row_for_assembly_style_hdf5(
    tmp_path: Path,
) -> None:
    h5 = _write_assembly_style_phase2_snps_h5(tmp_path / "assembly_snps.hdf5")

    run_pca_method(
        data=h5,
        name="assembly_with_ref",
        outdir=tmp_path / "OUT",
        method="pca",
        min_sample_coverage=1,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=True,
        impute_method="sample",
        subsample=False,
        random_seed=3,
        replicates=1,
        perplexity=5.0,
        max_iter=1000,
        n_neighbors=15,
        plot=False,
        cores=1,
        force=True,
        log_level="INFO",
    )

    coords = pd.read_csv(tmp_path / "OUT" / "assembly_with_ref.coords.tsv", sep="\t")
    assert coords["sample"].tolist() == [
        "assembly_reference_sequence",
        "a1",
        "a2",
        "b1",
        "b2",
    ]


def test_run_pca_analysis_public_api_returns_prepared_result(tmp_path: Path) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")

    result = run_pca_analysis(
        data=h5,
        min_sample_coverage=2,
        max_sample_missing=1.0,
        min_minor_allele_frequency=0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference=False,
        impute_method="none",
        subsample=True,
        random_seed=7,
        replicates=2,
        cores=1,
        log_level="INFO",
    )

    assert result.method == "pca"
    assert result.samples == ["a1", "a2", "a3", "b1", "b2", "b3"]
    assert sorted(result.coords_by_replicate) == [0, 1]
    assert sorted(result.variance_by_replicate) == [0, 1]
    assert result.primary_input.imputation.algorithm == "zero-fill"
    assert result.primary_input.imputation.imputed_snp_count == 1
    assert result.primary_input.imputation.total_snps == 3
    assert result.primary_input.matrix.shape == (6, 3)


def test_run_pca_method_logs_linewise_filter_stats_and_imputation_summary(tmp_path: Path) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")
    messages: list[str] = []
    sink_id = logger.add(messages.append, format="{message}", level="INFO")
    try:
        run_pca_method(
            data=h5,
            name="phase2",
            outdir=tmp_path / "OUT",
            method="pca",
            min_sample_coverage=2,
            max_sample_missing=1.0,
            min_minor_allele_frequency=0.0,
            imap=None,
            minmap=None,
            exclude=None,
            include_reference=False,
            impute_method="none",
            subsample=True,
            random_seed=7,
            replicates=1,
            perplexity=5.0,
            max_iter=1000,
            n_neighbors=15,
            plot=False,
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
        idx
        for idx, msg in enumerate(messages)
        if "pca SNP imputation: algorithm=zero-fill" in msg
    )
    prepared_idx = next(
        idx
        for idx, msg in enumerate(messages)
        if "pca prepared SNP summary:" in msg
    )

    assert extraction_idx < imputation_idx < prepared_idx
    assert any("filter statistic pre_filter_snps: 6" in msg for msg in messages)
    assert any("filter statistic post_filter_snps: 6" in msg for msg in messages)
    assert any(
        "filter statistic post_filter_percent_missing:" in msg
        and "linked post-filter genotype cells missing before optional subsampling" in msg
        for msg in messages
    )
    assert any(
        "pca SNP imputation: algorithm=zero-fill" in msg
        and "prepared_matrix_scope=subsampled_unlinked" in msg
        and "snp_columns_with_missing=1/3" in msg
        and "missing_genotype_cells=2/18" in msg
        for msg in messages
    )
    assert any(
        "pca prepared SNP summary:" in msg
        and "prepared_matrix_scope=subsampled_unlinked" in msg
        and "linked_post_filter_snps=6" in msg
        and "prepared_snps=3" in msg
        for msg in messages
    )
    assert not any("subsampled " in msg for msg in messages)


def test_run_pca_method_aggregates_multi_replicate_logging_and_keeps_details_at_debug(
    tmp_path: Path,
) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")
    messages: list[str] = []
    sink_id = logger.add(messages.append, format="{message}", level="DEBUG")
    try:
        run_pca_method(
            data=h5,
            name="phase2",
            outdir=tmp_path / "OUT",
            method="pca",
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
            replicates=2,
            perplexity=5.0,
            max_iter=1000,
            n_neighbors=15,
            plot=False,
            cores=1,
            force=True,
            log_level="INFO",
        )
    finally:
        logger.remove(sink_id)

    assert any(
        "pca SNP imputation across 2 replicates: algorithm=sample" in msg
        and "prepared_matrix_scope=subsampled_unlinked" in msg
        and "snp_columns_with_missing=" in msg
        and "missing_genotype_cells=" in msg
        for msg in messages
    )
    assert any("pca replicate 0 prepared SNP matrix:" in msg for msg in messages)
    assert any("pca replicate 1 prepared SNP matrix:" in msg for msg in messages)
    assert not any("replicate 0 matrix imputation" in msg for msg in messages)


def test_run_pca_method_rejects_vcf_input_and_multiple_tsne_replicates(tmp_path: Path) -> None:
    with pytest.raises(IPyradError, match="requires an SNP-capable HDF5 input"):
        run_pca_method(
            data=tmp_path / "variants.vcf.gz",
            name="bad",
            outdir=tmp_path / "OUT",
            method="pca",
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
            replicates=1,
            perplexity=2.0,
            max_iter=250,
            n_neighbors=3,
            plot=False,
            cores=1,
            force=True,
            log_level="INFO",
        )

    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")
    with pytest.raises(IPyradError, match="t-SNE supports exactly one run"):
        run_pca_method(
            data=h5,
            name="badtsne",
            outdir=tmp_path / "OUT2",
            method="tsne",
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
            replicates=2,
            perplexity=2.0,
            max_iter=250,
            n_neighbors=3,
            plot=False,
            cores=1,
            force=True,
            log_level="INFO",
        )


def test_run_pca_method_plot_is_supported_only_for_pca(tmp_path: Path) -> None:
    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")

    with pytest.raises(IPyradError, match="supported only with `-M pca`"):
        run_pca_method(
            data=h5,
            name="badtsneplot",
            outdir=tmp_path / "TSNE",
            method="tsne",
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
            replicates=1,
            perplexity=2.0,
            max_iter=250,
            n_neighbors=3,
            plot=True,
            cores=1,
            force=True,
            log_level="INFO",
        )

    with pytest.raises(IPyradError, match="supported only with `-M pca`"):
        run_pca_method(
            data=h5,
            name="badumapplot",
            outdir=tmp_path / "UMAP",
            method="umap",
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
            replicates=1,
            perplexity=2.0,
            max_iter=250,
            n_neighbors=3,
            plot=True,
            cores=1,
            force=True,
            log_level="INFO",
        )


def test_run_pca_method_without_plot_does_not_require_toyplot(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from ipyrad2.analysis.methods import common as common_methods

    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")

    def _fail_require_toyplot():
        raise AssertionError("toyplot should not be imported without --plot")

    monkeypatch.setattr(common_methods, "require_toyplot", _fail_require_toyplot)

    run_pca_method(
        data=h5,
        name="phase2",
        outdir=tmp_path / "OUT",
        method="pca",
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
        replicates=1,
        perplexity=5.0,
        max_iter=1000,
        n_neighbors=15,
        plot=False,
        cores=1,
        force=True,
        log_level="INFO",
    )


def test_run_pca_method_plot_fails_fast_when_toyplot_is_missing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from ipyrad2.analysis.methods import pca_drawing

    h5 = _write_phase2_snps_h5(tmp_path / "snps.hdf5")

    def _raise_missing():
        raise IPyradError("PCA plotting requires toyplot.")

    monkeypatch.setattr(pca_drawing, "require_toyplot", _raise_missing)

    with pytest.raises(IPyradError, match="requires toyplot"):
        run_pca_method(
            data=h5,
            name="phase2",
            outdir=tmp_path / "OUT",
            method="pca",
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
            replicates=1,
            perplexity=5.0,
            max_iter=1000,
            n_neighbors=15,
            plot=True,
            cores=1,
            force=True,
            log_level="INFO",
        )
