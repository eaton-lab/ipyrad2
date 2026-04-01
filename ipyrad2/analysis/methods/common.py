#!/usr/bin/env python

"""Shared helpers for phase-2 numerical analysis methods."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

from loguru import logger
import numpy as np
import pandas as pd

from ...utils.exceptions import IPyradError
from ..extractors.snp_extractor import SNPExportView, SNPExtractor
from .snp_imputer import _MISSING_GENO, SNPImputer


@dataclass
class ImputationSummary:
    """Summary statistics for matrix imputation on one prepared SNP view."""

    algorithm: str
    imputed_snp_count: int
    total_snps: int
    imputed_snp_fraction: float
    imputed_genotype_count: int
    total_genotypes: int
    imputed_genotype_fraction: float


@dataclass
class NumericalInput:
    """Filtered SNP data prepared for a numerical analysis method."""

    extractor: SNPExtractor
    view: SNPExportView
    matrix: np.ndarray
    imputation: ImputationSummary


@dataclass
class PreparedSNPViewSummary:
    """Shared counts for one prepared or exported SNP view."""

    samples_retained: int
    linked_post_filter_snps: int
    linked_post_filter_snp_containing_linkage_blocks: int
    selected_snps: int
    selected_snp_containing_linkage_blocks: int
    subsample: bool


def require_hdf5_input(data: Path | str, tool_name: str) -> Path:
    """Reject VCF input for phase-2 numerical tools."""
    path = Path(data).expanduser()
    text = str(path)
    if text.endswith((".vcf", ".vcf.gz")):
        raise IPyradError(
            f"`ipyrad2 analysis {tool_name}` requires an SNP-capable HDF5 input. "
            "Convert VCF first with `ipyrad2 analysis vcf-to-hdf5`."
        )
    return path


def normalize_impute_method(impute_method: str | None) -> str | None:
    """Normalize shared CLI imputation modes."""
    if impute_method == "sample":
        return "sample"
    if isinstance(impute_method, str) and impute_method.lower() in {"zero", "zero-fill"}:
        return "zero-fill"
    if impute_method in {None, False, "none"}:
        return None
    raise IPyradError("Unsupported imputation method. Use 'sample', 'zero-fill', or 'none'.")


def ensure_output_paths(paths: Iterable[Path], force: bool) -> None:
    """Fail early if any output already exists."""
    existing = next((path for path in paths if path.exists()), None)
    if existing is not None and not force:
        raise IPyradError(
            f"Output file already exists: {existing}. Use --force to overwrite."
        )


def run_snp_extractor_for_method(
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
) -> SNPExtractor:
    """Run the canonical SNP extractor once for a phase-2 method."""
    tool = SNPExtractor(
        data=Path(data),
        min_sample_coverage=min_sample_coverage,
        max_sample_missing=max_sample_missing,
        min_minor_allele_frequency=min_minor_allele_frequency,
        imap=imap,
        minmap=minmap,
        exclude=exclude,
        include_reference=include_reference,
        cores=cores,
    )
    tool.run(log_level=log_level)
    return tool


def get_numerical_input(
    extractor: SNPExtractor,
    *,
    subsample: bool,
    random_seed: int | None,
    impute_method: str | None,
    log_level: str,
) -> NumericalInput:
    """Return linked/unlinked SNP views and an analysis-ready matrix."""
    method = normalize_impute_method(impute_method)
    view = extractor.get_view(
        subsample=subsample,
        random_seed=random_seed,
        log_level=log_level,
    )
    imputation = summarize_imputation(view.genos, method)
    matrix = impute_genotype_matrix(
        view.genos,
        extractor,
        impute_method=method,
        random_seed=random_seed,
    )
    return NumericalInput(
        extractor=extractor,
        view=view,
        matrix=matrix,
        imputation=imputation,
    )


def summarize_prepared_snp_view(
    extractor: SNPExtractor,
    view: SNPExportView,
    *,
    subsample: bool,
) -> PreparedSNPViewSummary:
    """Return stable linked and selected SNP counts for one prepared view."""
    return PreparedSNPViewSummary(
        samples_retained=len(extractor.snames),
        linked_post_filter_snps=int(extractor.stats["post_filter_snps"]),
        linked_post_filter_snp_containing_linkage_blocks=int(
            extractor.stats["post_filter_snp_containing_linkage_blocks"]
        ),
        selected_snps=int(view.snpsmap.shape[0]),
        selected_snp_containing_linkage_blocks=count_linkage_blocks(view),
        subsample=bool(subsample),
    )


def _format_int_range(values: list[int]) -> str:
    """Render one integer or integer range compactly."""
    low = min(values)
    high = max(values)
    return str(low) if low == high else f"{low}-{high}"


def _format_fraction_range(values: list[float]) -> str:
    """Render one fraction or fraction range compactly."""
    low = min(values)
    high = max(values)
    return f"{low:.1%}" if np.isclose(low, high) else f"{low:.1%}-{high:.1%}"


def log_snp_imputation_summary(
    tool: str,
    summaries: ImputationSummary | Iterable[ImputationSummary] | None,
) -> None:
    """Log one shared SNP-imputation summary for one run."""
    if summaries is None:
        logger.info("{} SNP imputation: no imputation performed", tool)
        return

    if isinstance(summaries, ImputationSummary):
        items = [summaries]
    else:
        items = list(summaries)
    if not items:
        logger.info("{} SNP imputation: no imputation performed", tool)
        return

    algorithms = sorted({summary.algorithm for summary in items})
    if len(algorithms) != 1:
        raise IPyradError("SNP imputation logging requires one algorithm per run.")
    algorithm = algorithms[0]

    if len(items) == 1:
        summary = items[0]
        logger.info(
            "{} SNP imputation: algorithm={} snps={}/{} ({:.1%}) genotypes={}/{} ({:.1%})",
            tool,
            algorithm,
            summary.imputed_snp_count,
            summary.total_snps,
            summary.imputed_snp_fraction,
            summary.imputed_genotype_count,
            summary.total_genotypes,
            summary.imputed_genotype_fraction,
        )
        return

    logger.info(
        "{} SNP imputation across {} replicates: algorithm={} snps={}/{} ({}) genotypes={}/{} ({})",
        tool,
        len(items),
        algorithm,
        _format_int_range([summary.imputed_snp_count for summary in items]),
        _format_int_range([summary.total_snps for summary in items]),
        _format_fraction_range([summary.imputed_snp_fraction for summary in items]),
        _format_int_range([summary.imputed_genotype_count for summary in items]),
        _format_int_range([summary.total_genotypes for summary in items]),
        _format_fraction_range([summary.imputed_genotype_fraction for summary in items]),
    )


def log_snp_view_summary(
    tool: str,
    summary: PreparedSNPViewSummary,
    *,
    view_label: str,
) -> None:
    """Log one prepared or exported SNP-view count summary."""
    logger.info(
        "{} {} SNP summary: samples={} linked_post_filter_snps={} linked_post_filter_snp_containing_linkage_blocks={} {}_snps={} {}_snp_containing_linkage_blocks={} subsample={}",
        tool,
        view_label,
        summary.samples_retained,
        summary.linked_post_filter_snps,
        summary.linked_post_filter_snp_containing_linkage_blocks,
        view_label,
        summary.selected_snps,
        view_label,
        summary.selected_snp_containing_linkage_blocks,
        summary.subsample,
    )


def log_snp_replicate_details(
    tool: str,
    prepared_inputs: dict[int, NumericalInput],
    *,
    seeds: dict[int, int],
) -> None:
    """Log per-replicate SNP preparation details at DEBUG."""
    for rep, prepared in sorted(prepared_inputs.items()):
        summary = prepared.imputation
        logger.debug(
            "{} replicate {} prepared SNP matrix: seed={} shape={}x{} prepared_snps={} imputed_snps={}/{} imputed_genotypes={}/{}",
            tool,
            rep,
            seeds[rep],
            prepared.matrix.shape[0],
            prepared.matrix.shape[1],
            prepared.view.snpsmap.shape[0],
            summary.imputed_snp_count,
            summary.total_snps,
            summary.imputed_genotype_count,
            summary.total_genotypes,
        )


def resolve_imputation_algorithm_label(impute_method: str | None) -> str:
    """Return the user-facing label for one imputation mode."""
    if impute_method == "sample":
        return "sample"
    if impute_method in {None, "zero-fill"}:
        return "zero-fill"
    raise IPyradError("Unsupported imputation method label.")


def summarize_imputation(matrix: np.ndarray, impute_method: str | None) -> ImputationSummary:
    """Summarize missing-data imputation needs for one genotype matrix."""
    missing_mask = matrix == _MISSING_GENO
    imputed_genotype_count = int(np.sum(missing_mask))
    total_genotypes = int(matrix.size)
    total_snps = int(matrix.shape[1]) if matrix.ndim == 2 else 0
    imputed_snp_count = (
        int(np.count_nonzero(np.any(missing_mask, axis=0)))
        if total_snps
        else 0
    )
    return ImputationSummary(
        algorithm=resolve_imputation_algorithm_label(impute_method),
        imputed_snp_count=imputed_snp_count,
        total_snps=total_snps,
        imputed_snp_fraction=(imputed_snp_count / total_snps) if total_snps else 0.0,
        imputed_genotype_count=imputed_genotype_count,
        total_genotypes=total_genotypes,
        imputed_genotype_fraction=(
            imputed_genotype_count / total_genotypes if total_genotypes else 0.0
        ),
    )


def impute_genotype_matrix(
    matrix: np.ndarray,
    extractor: SNPExtractor,
    *,
    impute_method: str | None,
    random_seed: int | None,
) -> np.ndarray:
    """Impute a genotype matrix using the shared SNP imputer contract."""
    method = normalize_impute_method(impute_method)
    state = np.random.get_state()
    if random_seed is not None:
        np.random.seed(random_seed)
    try:
        result = SNPImputer(
            matrix.astype(np.uint8, copy=True),
            extractor.snames,
            imap=extractor.imap,
            impute_method=method,
            quiet=True,
        ).run()
        if np.any(result == _MISSING_GENO):
            raise IPyradError("Imputation failed to resolve all missing genotypes.")
        return result
    finally:
        np.random.set_state(state)


def calculate_sample_missing_fraction(
    matrix: np.ndarray,
    samples: list[str],
) -> pd.Series:
    """Return per-sample missingness on one selected genotype matrix."""
    values = np.asarray(matrix)
    if values.ndim != 2:
        raise IPyradError("Per-sample missingness requires a 2D genotype matrix.")
    if values.shape[0] != len(samples):
        raise IPyradError("Sample names must align to genotype-matrix rows.")
    if values.size == 0 or values.shape[1] == 0:
        missing = np.zeros(len(samples), dtype=float)
    else:
        missing = np.mean(values == _MISSING_GENO, axis=1)
    return pd.Series(missing, index=list(samples), dtype=float)


def _coerce_sample_fraction_series(
    values: pd.Series | np.ndarray | list[float],
    samples: list[str],
) -> pd.Series:
    """Coerce one per-sample numeric vector into the canonical sample order."""
    if isinstance(values, pd.Series):
        ordered = values.reindex(samples)
        if ordered.isna().any():
            raise IPyradError("Per-sample summaries must cover every retained sample.")
        return ordered.astype(float)
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.shape[0] != len(samples):
        raise IPyradError("Per-sample summaries must have one value per retained sample.")
    return pd.Series(array, index=list(samples), dtype=float)


def build_sample_data_summary(
    *,
    samples: list[str],
    missing_fraction: pd.Series | np.ndarray | list[float],
    post_imputation_missing_fraction: pd.Series | np.ndarray | list[float],
    imputation_algorithm: str,
) -> pd.DataFrame:
    """Build the shared sample-level data summary table."""
    pre = _coerce_sample_fraction_series(missing_fraction, samples)
    post = _coerce_sample_fraction_series(post_imputation_missing_fraction, samples)
    imputed = np.clip(pre.to_numpy() - post.to_numpy(), a_min=0.0, a_max=1.0)
    return pd.DataFrame(
        {
            "sample": samples,
            "missing_fraction": pre.to_numpy(),
            "post_imputation_missing_fraction": post.to_numpy(),
            "imputation_algorithm": [str(imputation_algorithm)] * len(samples),
            "imputed_genotype_fraction": imputed,
        }
    )


def build_imputed_sample_data_summary(
    *,
    samples: list[str],
    matrix: np.ndarray,
    impute_method: str | None,
) -> pd.DataFrame:
    """Build a sample summary for methods that fully impute missing data."""
    missing = calculate_sample_missing_fraction(matrix, samples)
    post = pd.Series(np.zeros(len(samples), dtype=float), index=samples, dtype=float)
    return build_sample_data_summary(
        samples=samples,
        missing_fraction=missing,
        post_imputation_missing_fraction=post,
        imputation_algorithm=resolve_imputation_algorithm_label(impute_method),
    )


def aggregate_sample_data_summaries(summaries: Iterable[pd.DataFrame]) -> pd.DataFrame:
    """Average one or more sample summary tables in shared sample order."""
    tables = list(summaries)
    if not tables:
        raise IPyradError("Cannot aggregate zero sample summary tables.")

    first = tables[0].reset_index(drop=True)
    sample_order = first["sample"].tolist()
    algorithm = first["imputation_algorithm"].tolist()
    numeric_cols = [
        "missing_fraction",
        "post_imputation_missing_fraction",
        "imputed_genotype_fraction",
    ]

    stacked = []
    for table in tables:
        current = table.reset_index(drop=True)
        if current["sample"].tolist() != sample_order:
            raise IPyradError("Sample summary aggregation requires identical sample order.")
        if current["imputation_algorithm"].tolist() != algorithm:
            raise IPyradError("Sample summary aggregation requires one imputation algorithm.")
        stacked.append(current[numeric_cols].to_numpy(dtype=float))

    mean_values = np.mean(np.stack(stacked, axis=0), axis=0)
    return pd.DataFrame(
        {
            "sample": sample_order,
            "missing_fraction": mean_values[:, 0],
            "post_imputation_missing_fraction": mean_values[:, 1],
            "imputation_algorithm": algorithm,
            "imputed_genotype_fraction": mean_values[:, 2],
        }
    )


def write_sample_data_summary(path: Path, data: pd.DataFrame) -> None:
    """Write the shared per-sample missingness/imputation summary table."""
    data.to_csv(path, sep="\t", index=False)


def write_membership(path: Path, samples: list[str], membership: np.ndarray) -> None:
    """Write sample-by-cluster membership coefficients."""
    data = {"sample": samples}
    for idx in range(membership.shape[1]):
        data[f"cluster{idx + 1}"] = membership[:, idx]
    pd.DataFrame(data).to_csv(path, sep="\t", index=False)


def write_assignments(path: Path, samples: list[str], membership: np.ndarray) -> None:
    """Write max-membership assignments from one membership matrix."""
    assignments = np.argmax(membership, axis=1) + 1
    scores = membership[np.arange(membership.shape[0]), assignments - 1]
    pd.DataFrame(
        {
            "sample": samples,
            "assigned_cluster": assignments,
            "assignment_score": scores,
        }
    ).to_csv(path, sep="\t", index=False)


def marker_ids_from_view(view: SNPExportView) -> list[str]:
    """Return stable SNP marker IDs from one filtered SNP view."""
    ids = []
    for row in view.snpsmap:
        loc = int(row[0])
        pos = int(row[4])
        ids.append(f"loc{loc}_pos{pos + 1}")
    return ids


def write_marker_cluster_matrix(
    path: Path,
    marker_ids: list[str],
    matrix: np.ndarray,
    *,
    marker_column: str = "marker_id",
) -> None:
    """Write a marker-by-cluster matrix in the shared clustering table shape."""
    data = {marker_column: marker_ids}
    for idx in range(matrix.shape[0]):
        data[f"cluster{idx + 1}"] = matrix[idx, :]
    pd.DataFrame(data).to_csv(path, sep="\t", index=False)


def write_stats_file(
    path: Path,
    *,
    tool: str,
    extractor: SNPExtractor,
    subsample: bool,
    random_seed: int | None,
    impute_method: str | None,
    summary: Dict[str, object],
    filter_stats: bool = True,
) -> None:
    """Write a human-readable stats file with shared SNP-extractor context."""
    with open(path, "w", encoding="utf-8") as out:
        out.write("Summary\n")
        out.write("-------\n")
        out.write(f"tool: {tool}\n")
        out.write(f"infile: {extractor.data}\n")
        out.write(f"samples_selected_initial: {extractor.initial_snames}\n")
        out.write(f"samples_dropped_by_max_missing: {extractor.dropped_samples_by_missing}\n")
        out.write(f"samples_final: {extractor.snames}\n")
        out.write(f"imap: {extractor.imap}\n")
        out.write(f"minmap: {extractor.minmap}\n")
        out.write(f"include_reference: {extractor.include_reference}\n")
        out.write(f"subsample: {subsample}\n")
        out.write(f"random_seed: {random_seed}\n")
        out.write(f"impute_method: {normalize_impute_method(impute_method)}\n")
        for key, value in summary.items():
            out.write(f"{key}: {value}\n")
        if filter_stats:
            out.write("\n")
            out.write("Filter statistics\n")
            out.write("-----------------\n")
            for key in extractor.stats.index:
                out.write(f"{key}: {extractor.stats[key]}\n")


def count_linkage_blocks(view: SNPExportView) -> int:
    """Return the number of linkage blocks represented in a SNP view."""
    if view.snpsmap.size == 0:
        return 0
    return int(np.unique(view.snpsmap[:, 0]).size)


def parse_k_range(value: str) -> tuple[int, int]:
    """Parse an inclusive K range of the form MIN:MAX."""
    try:
        lower_text, upper_text = value.split(":", maxsplit=1)
        lower = int(lower_text)
        upper = int(upper_text)
    except ValueError as exc:
        raise IPyradError("K ranges must use the form MIN:MAX.") from exc
    if lower < 1 or upper < lower:
        raise IPyradError("K ranges must satisfy 1 <= MIN <= MAX.")
    return lower, upper


def replicate_seeds(seed: int | None, nreplicates: int) -> List[int]:
    """Return deterministic per-replicate seeds."""
    base = int(seed) if seed is not None else int(np.random.default_rng().integers(2**31))
    rng = np.random.default_rng(base)
    return [int(rng.integers(2**31)) for _ in range(nreplicates)]


def require_sklearn():
    """Import scikit-learn lazily for methods that require it."""
    try:
        from sklearn import decomposition
        from sklearn.cluster import KMeans
        from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
        from sklearn.manifold import TSNE
    except ImportError as exc:
        raise IPyradError(
            "This analysis method requires scikit-learn. "
            "Install it with `pip install ipyrad2[analysis]` or "
            "`conda install scikit-learn -c conda-forge`."
        ) from exc
    return decomposition, KMeans, LinearDiscriminantAnalysis, TSNE


def require_umap():
    """Import umap-learn lazily for the UMAP method."""
    try:
        import umap
    except ImportError as exc:
        raise IPyradError(
            "UMAP requires umap-learn. Install it with "
            "`pip install ipyrad2[analysis]` or "
            "`conda install umap-learn -c conda-forge`."
        ) from exc
    return umap


def require_toyplot():
    """Import toyplot lazily for optional plotting helpers."""
    try:
        import toyplot
        import toyplot.svg
    except ImportError as exc:
        raise IPyradError(
            "PCA plotting requires toyplot. Install it with "
            "`pip install ipyrad2[analysis]` or "
            "`conda install toyplot -c conda-forge`."
        ) from exc
    return toyplot, toyplot.svg
