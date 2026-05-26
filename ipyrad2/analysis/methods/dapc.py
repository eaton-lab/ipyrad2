#!/usr/bin/env python

"""Sklearn-backed DAPC-style clustering on SNP HDF5 data."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ...utils.exceptions import IPyradError
from .common import (
    build_imputed_sample_data_summary,
    count_linkage_blocks,
    ensure_output_paths,
    get_numerical_input,
    log_snp_imputation_summary,
    log_snp_view_summary,
    normalize_impute_method,
    parse_k_range,
    require_hdf5_input,
    require_sklearn,
    run_snps_extracter_for_method,
    summarize_prepared_snp_view,
    write_assignments,
    write_membership,
    write_sample_data_summary,
    write_stats_file,
)


def _select_n_pcs(
    *,
    requested: int | None,
    n_samples: int,
    n_features: int,
    kmax: int,
) -> int:
    available = min(n_samples - 1, n_features)
    if available < 1:
        raise IPyradError("DAPC requires at least two samples and one SNP after filtering.")
    if requested is not None:
        if requested < kmax - 1:
            raise IPyradError("--n-pcs must be at least Kmax - 1.")
        if requested > available:
            raise IPyradError("--n-pcs cannot exceed the available principal components.")
        return requested
    base = min(20, available)
    return min(available, max(base, kmax - 1))


def _kmeans_bic(model, data: np.ndarray) -> float:
    """Return a lower-is-better BIC approximation for k-means."""
    labels = model.labels_
    centers = model.cluster_centers_
    n_clusters = model.n_clusters
    n_samples, n_features = data.shape
    cluster_sizes = np.bincount(labels, minlength=n_clusters)
    if np.any(cluster_sizes == 0):
        return float("inf")

    variance = np.sum((data - centers[labels]) ** 2) / max(n_samples - n_clusters, 1)
    variance /= max(n_features, 1)
    variance = max(float(variance), np.finfo(float).eps)

    log_likelihood = 0.0
    for size in cluster_sizes:
        log_likelihood += (
            size * np.log(size)
            - size * np.log(n_samples)
            - (size * n_features / 2.0) * np.log(2.0 * np.pi * variance)
            - ((size - 1.0) * n_features / 2.0)
        )
    n_params = n_clusters * (n_features + 1)
    return float(-2.0 * log_likelihood + n_params * np.log(n_samples))


def _fit_dapc(matrix: np.ndarray, *, k: int, n_pcs: int, seed: int | None):
    decomposition, KMeans, LinearDiscriminantAnalysis, _TSNE = require_sklearn()
    pca = decomposition.PCA(n_components=n_pcs, random_state=seed)
    pcs = pca.fit_transform(matrix.astype(np.float64, copy=False))
    kmeans = KMeans(n_clusters=k, random_state=seed, n_init=20)
    labels = kmeans.fit_predict(pcs)
    lda = LinearDiscriminantAnalysis()
    coords = lda.fit_transform(pcs, labels)
    membership = lda.predict_proba(pcs)
    return pcs, kmeans, coords, membership


def _write_coords(path: Path, samples: list[str], coords: np.ndarray) -> None:
    data = {"sample": samples}
    for idx in range(coords.shape[1]):
        data[f"axis{idx + 1}"] = coords[:, idx]
    pd.DataFrame(data).to_csv(path, sep="\t", index=False)


def run_dapc_method(
    *,
    data: Path | str,
    name: str,
    outdir: Path | str,
    k: int | None,
    k_range: str | None,
    n_pcs: int | None,
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
    cores: int,
    force: bool,
    log_level: str = "INFO",
    min_genotype_depth: int = 0,
    min_site_qual: float = 0.0,
) -> None:
    """CLI entrypoint for DAPC-style clustering."""
    require_hdf5_input(data, "dapc")
    normalized_impute = normalize_impute_method(impute_method)
    if (k is None) == (k_range is None):
        raise IPyradError("Specify exactly one of -k or --k-range.")

    outdir = Path(outdir).expanduser().absolute()
    paths = {
        "coords": outdir / f"{name}.coords.tsv",
        "membership": outdir / f"{name}.membership.tsv",
        "assignments": outdir / f"{name}.assignments.tsv",
        "k_scan": outdir / f"{name}.k_scan.tsv",
        "sample_data_summary": outdir / f"{name}.sample_data_summary.tsv",
        "stats": outdir / f"{name}.stats.txt",
    }
    ensure_output_paths(paths.values(), force=force)
    outdir.mkdir(parents=True, exist_ok=True)

    extracter = run_snps_extracter_for_method(
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
    prepared = get_numerical_input(
        extracter,
        subsample=subsample,
        random_seed=random_seed,
        impute_method=impute_method,
        log_level=log_level,
    )
    log_snp_imputation_summary("dapc", prepared.imputation, subsample=subsample)
    log_snp_view_summary(
        "dapc",
        summarize_prepared_snp_view(extracter, prepared.view, subsample=subsample),
        view_label="prepared",
    )

    if k_range is not None:
        lower, upper = parse_k_range(k_range)
        if lower < 2:
            raise IPyradError("K ranges for DAPC must start at 2 or greater.")
        if upper >= len(extracter.snames):
            raise IPyradError("Maximum K must be smaller than the number of retained samples.")
        kmax = upper
    else:
        if k is None or k < 2:
            raise IPyradError("K for DAPC must be 2 or greater.")
        if k >= len(extracter.snames):
            raise IPyradError("K must be smaller than the number of retained samples.")
        kmax = k

    retained_pcs = _select_n_pcs(
        requested=n_pcs,
        n_samples=prepared.matrix.shape[0],
        n_features=prepared.matrix.shape[1],
        kmax=kmax,
    )

    decomposition, KMeans, _LDA, _TSNE = require_sklearn()
    pca = decomposition.PCA(n_components=retained_pcs, random_state=random_seed)
    pcs = pca.fit_transform(prepared.matrix.astype(np.float64, copy=False))

    k_rows = []
    if k_range is not None:
        selected_row = None
        for current_k in range(lower, upper + 1):
            kmeans = KMeans(n_clusters=current_k, random_state=random_seed, n_init=20)
            kmeans.fit(pcs)
            bic = _kmeans_bic(kmeans, pcs)
            row = {"k": current_k, "bic": bic, "selected": False}
            k_rows.append(row)
            if selected_row is None or bic < selected_row["bic"]:
                selected_row = row
        selected_row["selected"] = True
        selected_k = int(selected_row["k"])
    else:
        selected_k = k
        kmeans = KMeans(n_clusters=selected_k, random_state=random_seed, n_init=20)
        kmeans.fit(pcs)
        k_rows.append(
            {"k": selected_k, "bic": _kmeans_bic(kmeans, pcs), "selected": True}
        )

    _pcs, _kmeans, coords, membership = _fit_dapc(
        prepared.matrix,
        k=selected_k,
        n_pcs=retained_pcs,
        seed=random_seed,
    )
    _write_coords(paths["coords"], extracter.snames, coords)
    write_membership(paths["membership"], extracter.snames, membership)
    write_assignments(paths["assignments"], extracter.snames, membership)
    pd.DataFrame.from_records(k_rows).to_csv(paths["k_scan"], sep="\t", index=False)
    sample_summary = build_imputed_sample_data_summary(
        samples=extracter.snames,
        matrix=prepared.view.genos,
        impute_method=normalized_impute,
    )
    write_sample_data_summary(paths["sample_data_summary"], sample_summary)
    write_stats_file(
        paths["stats"],
        tool="dapc",
        extracter=extracter,
        subsample=subsample,
        random_seed=random_seed,
        impute_method=impute_method,
        summary={
            "k_selected": selected_k,
            "k_range": k_range if k_range is not None else "NA",
            "n_pcs": retained_pcs,
            "linked_post_filter_snps": int(extracter.stats["post_filter_snps"]),
            "linked_post_filter_snp_containing_linkage_blocks": int(
                extracter.stats["post_filter_snp_containing_linkage_blocks"]
            ),
            "exported_snps": int(
                count_linkage_blocks(prepared.view)
                if subsample
                else prepared.view.snpsmap.shape[0]
            ),
            "exported_snp_containing_linkage_blocks": count_linkage_blocks(prepared.view),
            "samples_retained": len(extracter.snames),
        },
    )


__all__ = ["run_dapc_method"]
