#!/usr/bin/env python

"""CLI runner for genome-wide population-genetic statistics."""

from __future__ import annotations

from io import StringIO
from pathlib import Path
from typing import Any

import h5py
from loguru import logger
from ....utils.exceptions import IPyradError
from ..common import ensure_output_paths
from .models import OUTPUT_TABLE_SPECS
from .models import ORDERED_STATS
from .models import PopgenRequest
from .models import PopgenResult
from .models import SEQUENCE_ONLY_STATS
from .models import SEQUENCE_STATS
from .models import SNP_STATS
from .models import WindowingConfig
from .models import normalize_stats
from .seq_backend import run_sequence_popgen
from .snp_backend import run_snp_popgen


POPGEN_FLOAT_FORMAT = "%.8f"


def _detect_hdf5_capabilities(data: Path) -> tuple[bool, bool]:
    """Return `(has_sequence, has_snp)` for one analysis HDF5 file."""
    if str(data).endswith((".vcf", ".vcf.gz")):
        raise IPyradError(
            "`ipyrad2 analysis popgen` requires an analysis HDF5 input. "
            "Convert VCF first with `ipyrad2 analysis vcf-to-hdf5`."
        )
    with h5py.File(data, "r") as io5:
        has_sequence = "phy" in io5 and "phymap" in io5
        has_snp = "genos" in io5 and "snpsmap" in io5
    return has_sequence, has_snp


def _resolve_requested_stats(
    stats,
    *,
    has_sequence: bool,
    has_snp: bool,
) -> list[str]:
    """Resolve `all` into the stats supported by the available input datasets."""
    requested = normalize_stats(stats)
    if requested == ["all"]:
        if has_sequence:
            return [name for name in ORDERED_STATS if name in SEQUENCE_STATS]
        if has_snp:
            return [name for name in ORDERED_STATS if name in SNP_STATS]
    return requested


def _choose_backend(
    requested_stats: list[str],
    *,
    has_sequence: bool,
    has_snp: bool,
    window_mode: str | None,
) -> str:
    """Choose the sequence or SNP backend from input capability and stat needs."""
    if window_mode is not None:
        if not has_sequence:
            raise IPyradError(
                "Windowed popgen statistics currently require sequence HDF5 with `phy`/`phymap`."
            )
        return "sequence"
    if any(stat in SEQUENCE_ONLY_STATS for stat in requested_stats):
        if not has_sequence:
            unsupported = ", ".join(
                stat for stat in requested_stats if stat in SEQUENCE_ONLY_STATS
            )
            raise IPyradError(
                "The requested popgen statistics require sequence HDF5 with `phy`/`phymap`: "
                f"{unsupported}"
            )
        return "sequence"
    if has_snp:
        return "snp"
    if has_sequence:
        return "sequence"
    raise IPyradError(
        "Input HDF5 does not contain a supported sequence or SNP dataset for popgen."
    )


def _resolve_windowing(
    *,
    window_size: int | None,
    step_size: int | None,
    loci_per_window: int | None,
    locus_step: int | None,
) -> WindowingConfig:
    """Normalize and validate windowing options."""
    if step_size is not None and window_size is None:
        raise IPyradError("--step-size requires --window-size.")
    if locus_step is not None and loci_per_window is None:
        raise IPyradError("--locus-step requires --loci-per-window.")

    if window_size is not None:
        if window_size <= 0:
            raise IPyradError("--window-size must be greater than zero.")
        resolved_step = window_size if step_size is None else step_size
        if resolved_step <= 0:
            raise IPyradError("--step-size must be greater than zero.")
        return WindowingConfig(
            mode="genomic",
            window_size=window_size,
            step_size=resolved_step,
        )

    if loci_per_window is not None:
        if loci_per_window <= 0:
            raise IPyradError("--loci-per-window must be greater than zero.")
        resolved_step = loci_per_window if locus_step is None else locus_step
        if resolved_step <= 0:
            raise IPyradError("--locus-step must be greater than zero.")
        return WindowingConfig(
            mode="locus",
            loci_per_window=loci_per_window,
            locus_step=resolved_step,
        )

    return WindowingConfig()


def _write_manifest_file(
    path: Path,
    sections: list[tuple[str, dict[str, Any] | str]],
) -> None:
    """Write a human-readable popgen manifest file."""
    with open(path, "w", encoding="utf-8") as out:
        for section_idx, (title, content) in enumerate(sections):
            out.write(f"{title}\n")
            out.write(f"{'-' * len(title)}\n")
            if isinstance(content, dict):
                for key, value in content.items():
                    out.write(f"{key}: {value}\n")
            else:
                out.write(content)
                if not content.endswith("\n"):
                    out.write("\n")
            if section_idx != len(sections) - 1:
                out.write("\n")


def _table_row_count(table) -> int:
    """Return the number of rows in one optional tabular output."""
    if table is None:
        return 0
    return int(len(table.index))


def _write_popgen_table(path: Path, table) -> None:
    """Write one popgen TSV with a shared readable float format."""
    table.to_csv(path, sep="\t", index=False, float_format=POPGEN_FLOAT_FORMAT)


def _render_popgen_table(table) -> str:
    """Return one popgen table rendered as tab-delimited text."""
    buffer = StringIO()
    table.to_csv(buffer, sep="\t", index=False, float_format=POPGEN_FLOAT_FORMAT)
    return buffer.getvalue()


def _build_popgen_request(
    *,
    data: Path,
    name: str,
    outdir: Path,
    requested_stats: list[str],
    backend: str,
    has_sequence: bool,
    has_snp: bool,
    min_sample_coverage: float,
    max_sample_missing: float,
    min_minor_allele_frequency: float,
    imap,
    minmap,
    exclude,
    include_reference: bool,
    subsample_unlinked: bool,
    random_seed: int | None,
    cores: int,
    force: bool,
    log_level: str,
    windowing: WindowingConfig,
) -> PopgenRequest:
    """Return one resolved typed request for a popgen run."""
    return PopgenRequest(
        data=data,
        name=name,
        outdir=outdir,
        requested_stats=tuple(requested_stats),
        backend=backend,
        has_sequence=has_sequence,
        has_snp=has_snp,
        min_sample_coverage=min_sample_coverage,
        max_sample_missing=max_sample_missing,
        min_minor_allele_frequency=min_minor_allele_frequency,
        imap=imap,
        minmap=minmap,
        exclude=tuple() if exclude is None else tuple(exclude),
        include_reference=include_reference,
        subsample_unlinked=subsample_unlinked,
        random_seed=random_seed,
        cores=cores,
        force=force,
        log_level=log_level,
        windowing=windowing,
    )


def _build_output_paths(request: PopgenRequest, result: PopgenResult) -> dict[str, Path]:
    """Return ordered output paths for one typed popgen result."""
    paths = {"manifest": request.outdir / f"{request.name}.manifest.txt"}
    for spec in OUTPUT_TABLE_SPECS:
        table = result.get_output_table(spec.key)
        if table is None or table.empty:
            continue
        paths[spec.key] = request.outdir / f"{request.name}.{spec.suffix}"
    return paths


def run_popgen_method(
    *,
    data: Path | str,
    name: str,
    outdir: Path | str,
    stats,
    min_sample_coverage: float,
    max_sample_missing: float,
    min_minor_allele_frequency: float,
    imap,
    minmap,
    exclude,
    include_reference: bool,
    subsample_unlinked: bool,
    random_seed: int | None,
    cores: int,
    force: bool,
    log_level: str,
    window_size: int | None = None,
    step_size: int | None = None,
    loci_per_window: int | None = None,
    locus_step: int | None = None,
) -> None:
    """Run the genome-wide population-genetic analysis workflow."""
    data = Path(data).expanduser().absolute()
    outdir = Path(outdir).expanduser().absolute()

    has_sequence, has_snp = _detect_hdf5_capabilities(data)
    windowing = _resolve_windowing(
        window_size=window_size,
        step_size=step_size,
        loci_per_window=loci_per_window,
        locus_step=locus_step,
    )
    requested_stats = _resolve_requested_stats(
        stats,
        has_sequence=has_sequence,
        has_snp=has_snp,
    )
    if not requested_stats:
        raise IPyradError("No popgen statistics were requested.")
    backend = _choose_backend(
        requested_stats,
        has_sequence=has_sequence,
        has_snp=has_snp,
        window_mode=windowing.mode,
    )

    if backend == "sequence":
        if subsample_unlinked:
            raise IPyradError("--subsample-unlinked is only supported on SNP-backed popgen runs.")
        if random_seed is not None:
            raise IPyradError("--seed is only supported on SNP-backed popgen runs.")
        if min_minor_allele_frequency:
            raise IPyradError(
                "--min-minor-allele-frequency is only supported on SNP-backed popgen runs."
            )
        result = run_sequence_popgen(
            data=data,
            requested_stats=requested_stats,
            min_sample_coverage=min_sample_coverage,
            max_sample_missing=max_sample_missing,
            imap=imap,
            minmap=minmap,
            exclude=exclude,
            include_reference=include_reference,
            cores=cores,
            log_level=log_level,
            window_size=windowing.window_size,
            step_size=windowing.step_size,
            loci_per_window=windowing.loci_per_window,
            locus_step=windowing.locus_step,
        )
    else:
        result = run_snp_popgen(
            data=data,
            requested_stats=requested_stats,
            min_sample_coverage=min_sample_coverage,
            max_sample_missing=max_sample_missing,
            min_minor_allele_frequency=min_minor_allele_frequency,
            imap=imap,
            minmap=minmap,
            exclude=exclude,
            include_reference=include_reference,
            subsample_unlinked=subsample_unlinked,
            random_seed=random_seed,
            cores=cores,
            log_level=log_level,
        )

    request = _build_popgen_request(
        data=data,
        name=name,
        outdir=outdir,
        requested_stats=requested_stats,
        backend=backend,
        has_sequence=has_sequence,
        has_snp=has_snp,
        min_sample_coverage=min_sample_coverage,
        max_sample_missing=max_sample_missing,
        min_minor_allele_frequency=min_minor_allele_frequency,
        imap=imap,
        minmap=minmap,
        exclude=exclude,
        include_reference=include_reference,
        subsample_unlinked=subsample_unlinked,
        random_seed=random_seed,
        cores=cores,
        force=force,
        log_level=log_level,
        windowing=windowing,
    )
    if "fit" in request.requested_stats and (result.global_stats is None or result.global_stats.empty):
        raise IPyradError("Popgen backend did not return global Fit statistics.")

    paths = _build_output_paths(request, result)
    ensure_output_paths(paths.values(), request.force)
    request.outdir.mkdir(parents=True, exist_ok=True)

    sample_payload = result.sample_data_summary
    if sample_payload is None:
        raise IPyradError("Popgen backend did not return a sample data summary.")
    sample_stats = result.sample_stats
    if sample_stats is None:
        raise IPyradError("Popgen backend did not return per-sample statistics.")
    _write_popgen_table(paths["sample_stats"], sample_stats)
    if "global_stats" in paths:
        _write_popgen_table(paths["global_stats"], result.global_stats)

    for label, table in result.iter_output_tables():
        if label == "sample_stats" or label == "global_stats":
            continue
        _write_popgen_table(paths[label], table)

    outputs = {"manifest_file": paths["manifest"]}
    output_tables = {
        "sample_stats": sample_stats,
        "global_stats": result.global_stats,
        "population_stats": result.population_stats,
        "pairwise_stats": result.pairwise_stats,
        "sfs": result.sfs,
        "window_population_stats": result.window_population_stats,
        "window_pairwise_stats": result.window_pairwise_stats,
    }
    outputs["sample_data_summary_rows_in_manifest"] = _table_row_count(sample_payload)
    for label, path in paths.items():
        if label == "manifest":
            continue
        outputs[f"{label}_file"] = path
        outputs[f"{label}_rows"] = _table_row_count(output_tables.get(label))

    backend_summary = dict(result.summary)
    backend_name = backend_summary.pop("input_backend", request.backend)
    requested_section = {
        "requested_stats": request.requested_stat_list,
        "sfs_mode": (
            "folded genome-wide only"
            if "sfs" in request.requested_stats
            else "not requested"
        ),
    }
    for key, value in request.requested_stat_formulas().items():
        requested_section[key] = value
    if request.windowing.mode is not None and "sfs" in request.requested_stats:
        requested_section["window_sfs_note"] = "Windowed SFS is not written in this phase."

    sections = [
        (
            "Inputs",
            {
                "tool": "popgen",
                "infile": request.data,
                "input_has_sequence": request.has_sequence,
                "input_has_snp": request.has_snp,
                "backend_used": backend_name,
            },
        ),
        (
            "Filters",
            {
                "min_sample_coverage": request.min_sample_coverage,
                "max_sample_missing": request.max_sample_missing,
                "min_minor_allele_frequency": request.min_minor_allele_frequency,
                "include_reference": request.include_reference,
                "exclude": list(request.exclude),
                "imap": request.imap,
                "minmap": request.minmap,
            },
        ),
        ("Requested Stats", requested_section),
        ("Windowing", request.windowing.as_manifest_dict()),
        ("Backend Summary", backend_summary),
        ("Sample Data Summary", _render_popgen_table(sample_payload)),
        ("Outputs", outputs),
    ]
    _write_manifest_file(paths["manifest"], sections)
    for path in paths.values():
        logger.info("wrote popgen file {}", path)


__all__ = ["run_popgen_method"]
