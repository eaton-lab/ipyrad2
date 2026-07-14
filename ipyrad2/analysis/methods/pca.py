#!/usr/bin/env python

"""Numerical PCA-family methods on SNP HDF5 data."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from loguru import logger
import numpy as np
import pandas as pd

from ...utils.exceptions import IPyradError
from ..extracters.snps_extracter import SNPsExtracter
from .common import (
    NumericalInput,
    aggregate_sample_data_summaries,
    build_imputed_sample_data_summary,
    count_linkage_blocks,
    ensure_output_paths,
    get_numerical_input,
    log_snp_imputation_summary,
    log_snp_replicate_details,
    log_snp_view_summary,
    normalize_impute_method,
    replicate_seeds,
    require_hdf5_input,
    require_sklearn,
    require_umap,
    run_snps_extracter_for_method,
    summarize_prepared_snp_view,
    write_sample_data_summary,
    write_stats_file,
)


@dataclass
class PCAFamilyResult:
    """Analysis-ready results for one PCA-family run."""

    method: str
    samples: list[str]
    coords_by_replicate: dict[int, np.ndarray]
    variance_by_replicate: dict[int, np.ndarray]
    sample_missing: pd.Series
    extracter: SNPsExtracter
    prepared_inputs_by_replicate: dict[int, NumericalInput]

    @property
    def primary_input(self) -> NumericalInput:
        """Return the first prepared input, used for stats and output summaries."""
        if not self.prepared_inputs_by_replicate:
            raise IPyradError("PCA-family analysis did not prepare any genotype matrix.")
        return self.prepared_inputs_by_replicate[min(self.prepared_inputs_by_replicate)]


@dataclass
class PCA:
    """Notebook-friendly wrapper around the PCA analysis runner."""

    result: PCAFamilyResult

    @classmethod
    def run(
        cls,
        *,
        data: Path | str,
        min_sample_coverage: float = 4,
        max_sample_missing: float = 1.0,
        min_minor_allele_frequency: float = 0.0,
        min_genotype_depth: int = 0,
        min_site_qual: float = 0.0,
        imap=None,
        minmap=None,
        exclude=None,
        include_reference: bool = False,
        impute_method: str | None = "sample",
        subsample: bool = True,
        random_seed: int | None = None,
        replicates: int = 1,
        cores: int = 1,
        log_level: str = "SUCCESS",
    ) -> "PCA":
        """Run PCA and keep the in-memory result for interactive use."""
        from ...utils.logger import set_api_log_level
        set_api_log_level(log_level) 
        
        result = run_pca_analysis(
            data=data,
            min_sample_coverage=min_sample_coverage,
            max_sample_missing=max_sample_missing,
            min_minor_allele_frequency=min_minor_allele_frequency,
            min_genotype_depth=min_genotype_depth,
            min_site_qual=min_site_qual,
            imap=imap,
            minmap=minmap,
            exclude=exclude,
            include_reference=include_reference,
            impute_method=impute_method,
            subsample=subsample,
            random_seed=random_seed,
            replicates=replicates,
            cores=cores,
            log_level=log_level,
        )
        return cls(result=result)

    def draw(
        self,
        *,
        width: int = 400,
        height: int = 300,
        marker_size: int = 10,
        colors: Path | str | None = None,
    ):
        """Return a Toyplot canvas for display in a notebook."""
        _validate_plot_args(
            plot_width=width,
            plot_height=height,
            plot_marker_size=marker_size,
        )
        from .pca_drawing import draw_pca_plot

        return draw_pca_plot(
            self.result,
            width=width,
            height=height,
            marker_size=marker_size,
            colors=colors,
        )


def _require_matrix(matrix: np.ndarray, label: str) -> tuple[int, int]:
    """Validate a 2D genotype matrix before running a numerical method."""
    if matrix.ndim != 2:
        raise IPyradError(f"{label} requires a 2D genotype matrix.")
    nsamples, nsites = matrix.shape
    if nsamples < 2:
        raise IPyradError(f"{label} requires at least two samples after filtering.")
    if nsites < 1:
        raise IPyradError(f"{label} requires at least one SNP after filtering.")
    return nsamples, nsites


def _normalize_pca_impute_method(impute_method: str | None) -> str:
    """Normalize PCA-family imputation names to explicit fill-all algorithms."""
    normalized = normalize_impute_method(impute_method)
    if normalized is None:
        return "zero-fill"
    return normalized


def _validate_method_args(
    *,
    method: str,
    replicates: int,
    perplexity: float | None = None,
    max_iter: int | None = None,
    n_neighbors: int | None = None,
) -> None:
    """Validate shared PCA-family method arguments before running extraction."""
    if method not in {"pca", "tsne", "umap"}:
        raise IPyradError(f"Unsupported PCA method: {method}")
    if replicates < 1:
        raise IPyradError("PCA replicate count must be at least 1.")
    if method == "tsne" and replicates != 1:
        raise IPyradError("t-SNE supports exactly one run; use --replicates 1.")
    if method == "umap" and replicates != 1:
        raise IPyradError("UMAP supports exactly one run; use --replicates 1.")
    if method == "tsne":
        if perplexity is None or perplexity <= 0:
            raise IPyradError("t-SNE perplexity must be greater than zero.")
        if max_iter is None or max_iter < 250:
            raise IPyradError("t-SNE max_iter must be at least 250.")
    if method == "umap" and (n_neighbors is None or n_neighbors < 2):
        raise IPyradError("UMAP n_neighbors must be at least 2.")


def _validate_plot_args(
    *,
    plot_width: int,
    plot_height: int,
    plot_marker_size: int,
) -> None:
    """Validate basic PCA plot settings."""
    if plot_width < 200:
        raise IPyradError("PCA plot width must be at least 200 pixels.")
    if plot_height < 200:
        raise IPyradError("PCA plot height must be at least 200 pixels.")
    if plot_marker_size < 1:
        raise IPyradError("PCA plot marker size must be at least 1.")


def _compute_pca(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return PCA coordinates and explained-variance ratios from a matrix."""
    nsamples, _nsites = _require_matrix(matrix, "PCA")
    centered = matrix.astype(np.float64, copy=False) - matrix.mean(axis=0, keepdims=True)
    u, singular_values, _vt = np.linalg.svd(centered, full_matrices=False)
    coords = u * singular_values
    if singular_values.size == 0:
        variance_ratio = np.zeros(0, dtype=np.float64)
    else:
        explained = (singular_values ** 2) / max(nsamples - 1, 1)
        total = explained.sum()
        variance_ratio = explained / total if total else np.zeros_like(explained)
    return coords, variance_ratio


def _run_tsne_once(
    matrix: np.ndarray,
    *,
    perplexity: float,
    max_iter: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Run one t-SNE embedding."""
    _decomposition, _KMeans, _LDA, TSNE = require_sklearn()
    nsamples, _nsites = _require_matrix(matrix, "t-SNE")
    if perplexity >= nsamples:
        raise IPyradError(
            "t-SNE perplexity must be smaller than the number of samples after filtering."
        )
    coords = TSNE(
        perplexity=perplexity,
        init="pca",
        max_iter=int(max_iter),
        random_state=seed,
    ).fit_transform(matrix)
    return coords, np.array([], dtype=np.float64)


def _run_umap_once(
    matrix: np.ndarray,
    *,
    n_neighbors: int,
    embedding_random_state: int | None,
    n_jobs: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Run one UMAP embedding."""
    umap = require_umap()
    _require_matrix(matrix, "UMAP")
    coords = umap.UMAP(
        n_neighbors=n_neighbors,
        init="spectral",
        random_state=embedding_random_state,
        n_jobs=n_jobs,
    ).fit_transform(matrix)
    return coords, np.array([], dtype=np.float64)


def _resolve_umap_embedding_config(
    *,
    random_seed: int | None,
    embedding_seed: int,
    cores: int,
) -> tuple[int | None, int]:
    """Return UMAP embedding random_state and thread count."""
    if random_seed is None:
        return None, cores
    if cores == 1:
        return embedding_seed, 1
    logger.warning(
        "parallel UMAP does not support exact reproducibility; ignoring the UMAP "
        "embedding seed while keeping the provided random seed for preprocessing "
        "(subsampling/imputation). Use --cores 1 for reproducible UMAP embeddings."
    )
    return None, cores


def _write_coords(
    path: Path,
    names: list[str],
    coords_by_rep: dict[int, np.ndarray],
    method: str,
) -> None:
    """Write one coordinate table covering all samples and replicates."""
    records = []
    for rep, coords in coords_by_rep.items():
        for idx, sample in enumerate(names):
            row = {
                "sample": sample,
                "replicate": rep,
                "method": method,
            }
            for axis_idx in range(coords.shape[1]):
                row[f"axis{axis_idx + 1}"] = float(coords[idx, axis_idx])
            records.append(row)
    pd.DataFrame.from_records(records).to_csv(path, sep="\t", index=False)


def _write_variance(path: Path, variances: dict[int, np.ndarray]) -> None:
    """Write explained-variance ratios for PCA replicates."""
    records = []
    for rep, values in variances.items():
        for axis_idx, value in enumerate(values, start=1):
            records.append(
                {
                    "replicate": rep,
                    "axis": axis_idx,
                    "explained_variance_ratio": float(value),
                }
            )
    pd.DataFrame.from_records(records).to_csv(path, sep="\t", index=False)


def _build_extracter(
    *,
    data: Path | str,
    min_sample_coverage: float,
    max_sample_missing: float,
    min_minor_allele_frequency: float,
    imap,
    minmap,
    exclude,
    include_reference: bool,
    cores: int,
    log_level: str,
    min_genotype_depth: int = 0,
    min_site_qual: float = 0.0,
) -> SNPsExtracter:
    """Create and run the canonical SNP extracter once for a PCA-family method."""
    return run_snps_extracter_for_method(
        data=require_hdf5_input(data, "pca"),
        min_sample_coverage=min_sample_coverage,
        max_sample_missing=max_sample_missing,
        min_minor_allele_frequency=min_minor_allele_frequency,
        imap=imap,
        minmap=minmap,
        min_genotype_depth=min_genotype_depth,
        min_site_qual=min_site_qual,
        exclude=exclude,
        include_reference=include_reference,
        cores=cores,
        log_level=log_level,
    )


def _prepare_inputs(
    *,
    method: str,
    extracter: SNPsExtracter,
    subsample: bool,
    random_seed: int | None,
    replicates: int,
    impute_method: str,
    log_level: str,
) -> tuple[dict[int, int], dict[int, NumericalInput]]:
    """Prepare imputed genotype matrices for one or more replicates."""
    seeds = {
        rep: rep_seed for rep, rep_seed in enumerate(replicate_seeds(random_seed, replicates))
    }
    prepared_inputs = {}
    for rep, rep_seed in seeds.items():
        prepared = get_numerical_input(
            extracter,
            subsample=subsample,
            random_seed=rep_seed,
            impute_method=impute_method,
            log_level=log_level,
        )
        prepared_inputs[rep] = prepared
    log_snp_imputation_summary(
        method,
        [prepared.imputation for prepared in prepared_inputs.values()],
        subsample=subsample,
    )
    log_snp_view_summary(
        method,
        summarize_prepared_snp_view(
            extracter,
            prepared_inputs[min(prepared_inputs)].view,
            subsample=subsample,
        ),
        view_label="prepared",
    )
    if len(prepared_inputs) > 1:
        log_snp_replicate_details(method, prepared_inputs, seeds=seeds)
    return seeds, prepared_inputs


def _build_result(
    *,
    method: str,
    extracter: SNPsExtracter,
    prepared_inputs_by_replicate: dict[int, NumericalInput],
    coords_by_replicate: dict[int, np.ndarray],
    variance_by_replicate: dict[int, np.ndarray],
) -> PCAFamilyResult:
    """Construct the normalized PCA-family result object."""
    return PCAFamilyResult(
        method=method,
        samples=list(extracter.snames),
        coords_by_replicate=coords_by_replicate,
        variance_by_replicate=variance_by_replicate,
        sample_missing=extracter.sample_missing.copy(),
        extracter=extracter,
        prepared_inputs_by_replicate=prepared_inputs_by_replicate,
    )


def _build_summary(
    result: PCAFamilyResult,
    *,
    replicates: int,
    perplexity: float,
    max_iter: int,
    n_neighbors: int,
) -> dict[str, object]:
    """Build the stats-file summary block for one PCA-family run."""
    prepared = result.primary_input
    exported_snps = int(prepared.view.snpsmap.shape[0])
    exported_blocks = count_linkage_blocks(prepared.view)
    axes_written = result.coords_by_replicate[min(result.coords_by_replicate)].shape[1]
    return {
        "method": result.method,
        "replicates": replicates,
        "linked_post_filter_snps": int(result.extracter.stats["post_filter_snps"]),
        "linked_post_filter_snp_containing_linkage_blocks": int(
            result.extracter.stats["post_filter_snp_containing_linkage_blocks"]
        ),
        "exported_snps": exported_snps,
        "exported_snp_containing_linkage_blocks": exported_blocks,
        "samples_retained": len(result.samples),
        "method_axes_written": int(axes_written),
        "imputation_algorithm": prepared.imputation.algorithm,
        "imputed_snp_count": prepared.imputation.imputed_snp_count,
        "imputed_snp_fraction": prepared.imputation.imputed_snp_fraction,
        "imputed_genotype_count": prepared.imputation.imputed_genotype_count,
        "imputed_genotype_fraction": prepared.imputation.imputed_genotype_fraction,
        "perplexity": perplexity if result.method == "tsne" else "NA",
        "max_iter": max_iter if result.method == "tsne" else "NA",
        "n_neighbors": n_neighbors if result.method == "umap" else "NA",
    }


def _add_population_assignments(
    sample_summary: pd.DataFrame,
    *,
    extracter: SNPsExtracter,
) -> pd.DataFrame:
    """Insert one final-population label column into the retained-sample summary."""
    sample_to_population: dict[str, str] = {}
    for population, names in extracter.imap.items():
        for name in names:
            if name in sample_to_population:
                raise IPyradError(f"Sample {name!r} was assigned to multiple populations.")
            sample_to_population[name] = population
    missing = [
        str(name)
        for name in sample_summary["sample"].tolist()
        if str(name) not in sample_to_population
    ]
    if missing:
        raise IPyradError(
            "Sample summary is missing population assignments for: " + ", ".join(missing)
        )
    result = sample_summary.copy()
    result.insert(
        1,
        "population",
        [sample_to_population[str(name)] for name in result["sample"].tolist()],
    )
    return result


def run_pca_analysis(
    *,
    data: Path | str,
    min_sample_coverage: float = 4,
    max_sample_missing: float = 1.0,
    min_minor_allele_frequency: float = 0.0,
    min_genotype_depth: int = 0,
    min_site_qual: float = 0.0,
    imap=None,
    minmap=None,
    exclude=None,
    include_reference: bool = False,
    impute_method: str | None = "sample",
    subsample: bool = True,
    random_seed: int | None = None,
    replicates: int = 1,
    cores: int = 1,
    log_level: str = "INFO",
) -> PCAFamilyResult:
    """Run one or more PCA replicates and return coordinates plus variance ratios."""
    _validate_method_args(method="pca", replicates=replicates)
    canonical_impute = _normalize_pca_impute_method(impute_method)
    if canonical_impute != impute_method and impute_method in {None, False, "none"}:
        logger.debug(
            "normalized PCA-family imputation mode {!r} to '{}'",
            impute_method,
            canonical_impute,
        )
    extracter = _build_extracter(
        data=data,
        min_sample_coverage=min_sample_coverage,
        max_sample_missing=max_sample_missing,
        min_minor_allele_frequency=min_minor_allele_frequency,
        imap=imap,
        minmap=minmap,
        min_genotype_depth=min_genotype_depth,
        min_site_qual=min_site_qual,
        exclude=exclude,
        include_reference=include_reference,
        cores=cores,
        log_level=log_level,
    )
    _seeds, prepared_inputs = _prepare_inputs(
        method="pca",
        extracter=extracter,
        subsample=subsample,
        random_seed=random_seed,
        replicates=replicates,
        impute_method=canonical_impute,
        log_level=log_level,
    )

    coords_by_replicate = {}
    variance_by_replicate = {}
    for rep, prepared in prepared_inputs.items():
        coords, variance = _compute_pca(prepared.matrix)
        coords_by_replicate[rep] = coords
        variance_by_replicate[rep] = variance
    return _build_result(
        method="pca",
        extracter=extracter,
        prepared_inputs_by_replicate=prepared_inputs,
        coords_by_replicate=coords_by_replicate,
        variance_by_replicate=variance_by_replicate,
    )


def run_tsne_analysis(
    *,
    data: Path | str,
    min_sample_coverage: float = 4,
    max_sample_missing: float = 1.0,
    min_minor_allele_frequency: float = 0.0,
    min_genotype_depth: int = 0,
    min_site_qual: float = 0.0,
    imap=None,
    minmap=None,
    exclude=None,
    include_reference: bool = False,
    impute_method: str | None = "sample",
    subsample: bool = True,
    random_seed: int | None = None,
    perplexity: float = 5.0,
    max_iter: int = 1000,
    cores: int = 1,
    log_level: str = "INFO",
) -> PCAFamilyResult:
    """Run one t-SNE embedding on a fully imputed genotype matrix."""
    _validate_method_args(
        method="tsne",
        replicates=1,
        perplexity=perplexity,
        max_iter=max_iter,
    )
    canonical_impute = _normalize_pca_impute_method(impute_method)
    if canonical_impute != impute_method and impute_method in {None, False, "none"}:
        logger.debug(
            "normalized PCA-family imputation mode {!r} to '{}'",
            impute_method,
            canonical_impute,
        )
    extracter = _build_extracter(
        data=data,
        min_sample_coverage=min_sample_coverage,
        max_sample_missing=max_sample_missing,
        min_minor_allele_frequency=min_minor_allele_frequency,
        imap=imap,
        minmap=minmap,
        min_genotype_depth=min_genotype_depth,
        min_site_qual=min_site_qual,
        exclude=exclude,
        include_reference=include_reference,
        cores=cores,
        log_level=log_level,
    )
    seeds, prepared_inputs = _prepare_inputs(
        method="tsne",
        extracter=extracter,
        subsample=subsample,
        random_seed=random_seed,
        replicates=1,
        impute_method=canonical_impute,
        log_level=log_level,
    )
    seed = seeds[0]
    prepared = prepared_inputs[0]
    coords, variance = _run_tsne_once(
        prepared.matrix,
        perplexity=perplexity,
        max_iter=max_iter,
        seed=seed,
    )
    return _build_result(
        method="tsne",
        extracter=extracter,
        prepared_inputs_by_replicate=prepared_inputs,
        coords_by_replicate={0: coords},
        variance_by_replicate={0: variance},
    )


def run_umap_analysis(
    *,
    data: Path | str,
    min_sample_coverage: float = 4,
    max_sample_missing: float = 1.0,
    min_minor_allele_frequency: float = 0.0,
    min_genotype_depth: int = 0,
    min_site_qual: float = 0.0,
    imap=None,
    minmap=None,
    exclude=None,
    include_reference: bool = False,
    impute_method: str | None = "sample",
    subsample: bool = True,
    random_seed: int | None = None,
    n_neighbors: int = 15,
    cores: int = 1,
    log_level: str = "INFO",
) -> PCAFamilyResult:
    """Run one UMAP embedding on a fully imputed genotype matrix."""
    _validate_method_args(method="umap", replicates=1, n_neighbors=n_neighbors)
    canonical_impute = _normalize_pca_impute_method(impute_method)
    if canonical_impute != impute_method and impute_method in {None, False, "none"}:
        logger.debug(
            "normalized PCA-family imputation mode {!r} to '{}'",
            impute_method,
            canonical_impute,
        )
    extracter = _build_extracter(
        data=data,
        min_sample_coverage=min_sample_coverage,
        max_sample_missing=max_sample_missing,
        min_minor_allele_frequency=min_minor_allele_frequency,
        imap=imap,
        minmap=minmap,
        min_genotype_depth=min_genotype_depth,
        min_site_qual=min_site_qual,
        exclude=exclude,
        include_reference=include_reference,
        cores=cores,
        log_level=log_level,
    )
    seeds, prepared_inputs = _prepare_inputs(
        method="umap",
        extracter=extracter,
        subsample=subsample,
        random_seed=random_seed,
        replicates=1,
        impute_method=canonical_impute,
        log_level=log_level,
    )
    seed = seeds[0]
    prepared = prepared_inputs[0]
    embedding_random_state, umap_n_jobs = _resolve_umap_embedding_config(
        random_seed=random_seed,
        embedding_seed=seed,
        cores=cores,
    )
    coords, variance = _run_umap_once(
        prepared.matrix,
        n_neighbors=n_neighbors,
        embedding_random_state=embedding_random_state,
        n_jobs=umap_n_jobs,
    )
    return _build_result(
        method="umap",
        extracter=extracter,
        prepared_inputs_by_replicate=prepared_inputs,
        coords_by_replicate={0: coords},
        variance_by_replicate={0: variance},
    )


def run_pca_method(
    *,
    data: Path | str,
    name: str,
    outdir: Path | str,
    method: str,
    min_sample_coverage: float,
    max_sample_missing: float,
    min_minor_allele_frequency: float,
    imap,
    minmap,
    exclude,
    include_reference: bool,
    impute_method: str | None,
    subsample: bool,
    random_seed: int | None,
    replicates: int,
    perplexity: float,
    max_iter: int,
    n_neighbors: int,
    plot: bool = False,
    plot_width: int = 400,
    plot_height: int = 300,
    plot_marker_size: int = 10,
    colors: Path | str | None = None,
    cores: int = 1,
    force: bool = False,
    log_level: str = "INFO",
    min_genotype_depth: int = 0,
    min_site_qual: float = 0.0,
) -> None:
    """CLI entrypoint for numerical PCA-family methods."""
    _validate_method_args(
        method=method,
        replicates=replicates,
        perplexity=perplexity,
        max_iter=max_iter,
        n_neighbors=n_neighbors,
    )
    if plot and method != "pca":
        raise IPyradError("PCA plotting is currently supported only with `-M pca`.")
    if colors is not None and not plot:
        raise IPyradError("PCA --plot-colors can only be used with --plot.")
    if colors is not None and imap is None:
        raise IPyradError("PCA --plot-colors requires --imap population assignments.")
    if plot:
        _validate_plot_args(
            plot_width=plot_width,
            plot_height=plot_height,
            plot_marker_size=plot_marker_size,
        )
    canonical_impute = _normalize_pca_impute_method(impute_method)
    if canonical_impute != impute_method and impute_method in {None, False, "none"}:
        logger.debug(
            "normalized PCA-family imputation mode {!r} to '{}'",
            impute_method,
            canonical_impute,
        )

    outdir = Path(outdir).expanduser().absolute()
    paths = {
        "coords": outdir / f"{name}.coords.tsv",
        "sample_data_summary": outdir / f"{name}.sample_data_summary.tsv",
        "stats": outdir / f"{name}.stats.txt",
    }
    if method == "pca":
        paths["variance"] = outdir / f"{name}.variance.tsv"
    if plot:
        from .pca_drawing import ensure_pca_plotting_available

        ensure_pca_plotting_available()
        paths["plot"] = outdir / f"{name}.plot.svg"

    ensure_output_paths(paths.values(), force=force)
    outdir.mkdir(parents=True, exist_ok=True)

    runners: dict[str, Callable[[], PCAFamilyResult]] = {
        "pca": lambda: run_pca_analysis(
            data=data,
            min_sample_coverage=min_sample_coverage,
            max_sample_missing=max_sample_missing,
            min_minor_allele_frequency=min_minor_allele_frequency,
            imap=imap,
            minmap=minmap,
            min_genotype_depth=min_genotype_depth,
            min_site_qual=min_site_qual,
            exclude=exclude,
            include_reference=include_reference,
            impute_method=canonical_impute,
            subsample=subsample,
            random_seed=random_seed,
            replicates=replicates,
            cores=cores,
            log_level=log_level,
        ),
        "tsne": lambda: run_tsne_analysis(
            data=data,
            min_sample_coverage=min_sample_coverage,
            max_sample_missing=max_sample_missing,
            min_minor_allele_frequency=min_minor_allele_frequency,
            imap=imap,
            minmap=minmap,
            min_genotype_depth=min_genotype_depth,
            min_site_qual=min_site_qual,
            exclude=exclude,
            include_reference=include_reference,
            impute_method=canonical_impute,
            subsample=subsample,
            random_seed=random_seed,
            perplexity=perplexity,
            max_iter=max_iter,
            cores=cores,
            log_level=log_level,
        ),
        "umap": lambda: run_umap_analysis(
            data=data,
            min_sample_coverage=min_sample_coverage,
            max_sample_missing=max_sample_missing,
            min_minor_allele_frequency=min_minor_allele_frequency,
            imap=imap,
            minmap=minmap,
            min_genotype_depth=min_genotype_depth,
            min_site_qual=min_site_qual,
            exclude=exclude,
            include_reference=include_reference,
            impute_method=canonical_impute,
            subsample=subsample,
            random_seed=random_seed,
            n_neighbors=n_neighbors,
            cores=cores,
            log_level=log_level,
        ),
    }
    result = runners[method]()

    sample_summary = aggregate_sample_data_summaries(
        build_imputed_sample_data_summary(
            samples=result.samples,
            matrix=prepared.view.genos,
            impute_method=canonical_impute,
        )
        for prepared in result.prepared_inputs_by_replicate.values()
    )
    sample_summary = _add_population_assignments(
        sample_summary,
        extracter=result.extracter,
    )
    _write_coords(paths["coords"], result.samples, result.coords_by_replicate, method)
    write_sample_data_summary(
        paths["sample_data_summary"],
        sample_summary,
        float_format="%.3f",
    )
    if method == "pca":
        _write_variance(paths["variance"], result.variance_by_replicate)
    if plot:
        from .pca_drawing import write_pca_svg_plot

        write_pca_svg_plot(
            result,
            paths["plot"],
            width=plot_width,
            height=plot_height,
            marker_size=plot_marker_size,
            colors=colors,
        )
    write_stats_file(
        paths["stats"],
        tool="pca",
        extracter=result.extracter,
        subsample=subsample,
        random_seed=random_seed,
        impute_method=canonical_impute,
        summary=_build_summary(
            result,
            replicates=replicates,
            perplexity=perplexity,
            max_iter=max_iter,
            n_neighbors=n_neighbors,
        ),
        sample_reporting="counts",
    )
    logger.info("wrote PCA-family coordinates to {}", paths["coords"])
    logger.info(
        "wrote PCA-family sample data summary to {}",
        paths["sample_data_summary"],
    )
    if method == "pca":
        logger.info("wrote PCA explained variance to {}", paths["variance"])
    if plot:
        logger.info("wrote PCA SVG plot to {}", paths["plot"])
    logger.info("wrote PCA-family stats to {}", paths["stats"])


__all__ = [
    "PCA",
    "PCAFamilyResult",
    "run_pca_analysis",
    "run_pca_method",
    "run_tsne_analysis",
    "run_umap_analysis",
]
