#!/usr/bin/env python

"""External ADMIXTURE wrapper for SNP HDF5 inputs."""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from ...utils.exceptions import IPyradError
from .common import (
    build_imputed_sample_data_summary,
    count_linkage_blocks,
    ensure_output_paths,
    log_snp_imputation_summary,
    log_snp_view_summary,
    normalize_impute_method,
    parse_k_range,
    require_hdf5_input,
    run_snps_extracter_for_method,
    summarize_imputation,
    summarize_prepared_snp_view,
    write_assignments,
    write_marker_cluster_matrix,
    write_membership,
    write_sample_data_summary,
    write_stats_file,
)


_CV_ERROR_PATTERN = re.compile(r"CV error \(K=\d+\):\s*([0-9.eE+-]+)")


def _resolve_binary(binary: str | None) -> str:
    """Resolve the external ADMIXTURE binary path."""
    resolved = binary or shutil.which("admixture")
    if not resolved:
        raise IPyradError(
            "Could not find the `admixture` binary. Install it and ensure it is on PATH, "
            "or pass an explicit path with `--binary`."
        )
    return resolved


def _parse_cv_error(text: str) -> float | None:
    """Parse a CV error from combined admixture stdout/stderr."""
    match = _CV_ERROR_PATTERN.search(text)
    if not match:
        return None
    return float(match.group(1))


def _admixture_output_paths(stage_dir: Path, stem: str, k: int) -> dict[str, Path]:
    """Return the expected raw output paths for one ADMIXTURE run."""
    return {
        "p": stage_dir / f"{stem}.{k}.P",
        "q": stage_dir / f"{stem}.{k}.Q",
        "log": stage_dir / f"{stem}.{k}.log",
    }


def _run_admixture_once(
    *,
    binary: str,
    bed_path: Path,
    stage_dir: Path,
    k: int,
    cores: int,
    with_cv: bool,
) -> dict[str, object]:
    """Run one ADMIXTURE job and return parsed metadata."""
    cmd = [
        binary,
        f"-j{cores}",
    ]
    if with_cv:
        cmd.append("--cv=5")
    cmd.extend(
        [
            str(bed_path.name),
            str(k),
        ]
    )

    result = subprocess.run(
        cmd,
        cwd=stage_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    combined = "\n".join(part for part in (result.stdout, result.stderr) if part)
    paths = _admixture_output_paths(stage_dir, bed_path.stem, k)
    with open(paths["log"], "w", encoding="utf-8") as out:
        out.write(combined)
        if combined and not combined.endswith("\n"):
            out.write("\n")

    if result.returncode != 0:
        raise IPyradError(
            f"`admixture` failed for K={k} with exit code {result.returncode}.\n{combined.strip()}"
        )
    if not paths["p"].exists() or not paths["q"].exists():
        raise IPyradError(
            f"`admixture` completed for K={k} but did not write expected .P/.Q outputs."
        )

    return {
        "k": k,
        "cv_error": _parse_cv_error(combined) if with_cv else None,
        "paths": paths,
    }


def _read_p_matrix(path: Path, expected_k: int) -> np.ndarray:
    """Read the selected .P frequency matrix with stable 2D shape."""
    matrix = np.loadtxt(path)
    if matrix.ndim == 0:
        matrix = matrix.reshape(1, 1)
    elif matrix.ndim == 1:
        matrix = matrix.reshape(expected_k, -1)
    return np.asarray(matrix, dtype=float)


def _read_membership(path: Path, nsamples: int) -> np.ndarray:
    """Read the selected .Q membership matrix with stable 2D shape."""
    matrix = np.loadtxt(path)
    if matrix.ndim == 1:
        if nsamples == 1:
            matrix = matrix.reshape(1, -1)
        else:
            matrix = matrix.reshape(nsamples, -1)
    return np.asarray(matrix, dtype=float)


def _read_marker_ids(bim_path: Path) -> list[str]:
    """Read PLINK marker IDs in order."""
    ids = []
    with open(bim_path, encoding="utf-8") as infile:
        for line in infile:
            parts = line.rstrip().split("\t")
            if len(parts) < 2:
                raise IPyradError(f"Malformed BIM row in {bim_path}.")
            ids.append(parts[1])
    return ids


def run_admixture_method(
    *,
    data: Path | str,
    name: str,
    outdir: Path | str,
    k: int | None,
    k_range: str | None,
    binary: str | None,
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
    keep_intermediates: bool,
    cores: int,
    force: bool,
    log_level: str = "INFO",
    min_genotype_depth: int = 0,
    min_site_qual: float = 0.0,
) -> None:
    """CLI entrypoint for external ADMIXTURE runs."""
    require_hdf5_input(data, "admixture")
    normalized_impute = normalize_impute_method(impute_method)
    if (k is None) == (k_range is None):
        raise IPyradError("Specify exactly one of -k or --k-range.")

    outdir = Path(outdir).expanduser().absolute()
    paths = {
        "membership": outdir / f"{name}.membership.tsv",
        "allele_frequencies": outdir / f"{name}.allele_frequencies.tsv",
        "assignments": outdir / f"{name}.assignments.tsv",
        "k_scan": outdir / f"{name}.k_scan.tsv",
        "sample_data_summary": outdir / f"{name}.sample_data_summary.tsv",
        "stats": outdir / f"{name}.stats.txt",
    }
    ensure_output_paths(paths.values(), force=force)
    outdir.mkdir(parents=True, exist_ok=True)

    binary_path = _resolve_binary(binary)
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
    view = extracter.get_view(
        subsample=subsample,
        random_seed=random_seed,
        log_level="DEBUG",
    )
    imputation_summary = (
        summarize_imputation(view.genos, normalized_impute)
        if normalized_impute is not None
        else None
    )
    log_snp_imputation_summary("admixture", imputation_summary, subsample=subsample)
    log_snp_view_summary(
        "admixture",
        summarize_prepared_snp_view(extracter, view, subsample=subsample),
        view_label="prepared",
    )
    if view.reference is None:
        raise IPyradError(
            "ADMIXTURE export requires the HDF5 `reference` dataset. "
            "Rebuild the SNP HDF5 with a current assemble or `ipyrad2 vcf2hdf5` run."
        )

    nsamples = len(extracter.snames)
    if k_range is not None:
        lower, upper = parse_k_range(k_range)
        if lower < 2:
            raise IPyradError("ADMIXTURE K ranges must start at 2 or greater.")
        if upper >= nsamples:
            raise IPyradError("Maximum K must be smaller than the number of retained samples.")
        k_values = list(range(lower, upper + 1))
    else:
        if k is None or k < 2:
            raise IPyradError("K for ADMIXTURE must be 2 or greater.")
        if k >= nsamples:
            raise IPyradError("K must be smaller than the number of retained samples.")
        k_values = [k]

    stage_context = None
    if keep_intermediates:
        stage_dir = outdir / f"{name}.intermediates"
        if stage_dir.exists():
            if not force:
                raise IPyradError(
                    f"Output directory already exists: {stage_dir}. Use --force to overwrite."
                )
            shutil.rmtree(stage_dir)
        stage_dir.mkdir(parents=True, exist_ok=True)
    else:
        stage_context = tempfile.TemporaryDirectory(
            dir=outdir,
            prefix=f".{name}.admixture.",
        )
        stage_dir = Path(stage_context.name)

    try:
        bed_prefix = stage_dir / name
        plink_paths = extracter.write_plink(
            bed_prefix,
            view,
            impute_method=normalized_impute,
        )

        run_rows = []
        selected = None
        for current_k in k_values:
            run_info = _run_admixture_once(
                binary=binary_path,
                bed_path=plink_paths["bed"],
                stage_dir=stage_dir,
                k=current_k,
                cores=cores,
                with_cv=(k_range is not None),
            )
            row = {
                "k": current_k,
                "cv_error": run_info["cv_error"],
                "selected": False,
            }
            run_rows.append(row)
            if selected is None:
                selected = (row, run_info)
            elif k_range is not None and run_info["cv_error"] is not None:
                if selected[1]["cv_error"] is None or run_info["cv_error"] < selected[1]["cv_error"]:
                    selected = (row, run_info)

        selected[0]["selected"] = True
        selected_info = selected[1]

        membership = _read_membership(selected_info["paths"]["q"], nsamples)
        marker_ids = _read_marker_ids(plink_paths["bim"])
        allele_freqs = _read_p_matrix(selected_info["paths"]["p"], selected_info["k"])
        if allele_freqs.shape[1] != len(marker_ids):
            raise IPyradError(
                "Selected ADMIXTURE .P output does not match the staged BIM marker count."
            )

        write_membership(paths["membership"], extracter.snames, membership)
        write_assignments(paths["assignments"], extracter.snames, membership)
        write_marker_cluster_matrix(paths["allele_frequencies"], marker_ids, allele_freqs)
        pd.DataFrame.from_records(run_rows).to_csv(paths["k_scan"], sep="\t", index=False)
        sample_summary = build_imputed_sample_data_summary(
            samples=extracter.snames,
            matrix=view.genos,
            impute_method=normalized_impute,
        )
        write_sample_data_summary(paths["sample_data_summary"], sample_summary)
        write_stats_file(
            paths["stats"],
            tool="admixture",
            extracter=extracter,
            subsample=subsample,
            random_seed=random_seed,
            impute_method=normalized_impute,
            summary={
                "binary": binary_path,
                "k_selected": selected_info["k"],
                "requested_k": k if k is not None else "NA",
                "k_range": k_range if k_range is not None else "NA",
                "selected_cv_error": selected_info["cv_error"],
                "keep_intermediates": keep_intermediates,
                "linked_post_filter_snps": int(extracter.stats["post_filter_snps"]),
                "linked_post_filter_snp_containing_linkage_blocks": int(
                    extracter.stats["post_filter_snp_containing_linkage_blocks"]
                ),
                "exported_snps": int(view.snpsmap.shape[0]),
                "exported_snp_containing_linkage_blocks": count_linkage_blocks(view),
                "samples_retained": nsamples,
            },
        )
    finally:
        if stage_context is not None:
            stage_context.cleanup()


__all__ = ["run_admixture_method"]
