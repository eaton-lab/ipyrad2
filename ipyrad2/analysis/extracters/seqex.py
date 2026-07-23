"""Extract filtered, delimited loci from an ipyrad2 sequence HDF5 file."""

from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
import sys

import h5py
from loguru import logger
import numpy as np

from ...utils.exceptions import IPyradError
from ...utils.parallel import run_with_pool_iter
from .sequence_common import MISSING_BASE
from .sequence_common import build_sequence_extraction_context
from .sequence_common import filter_block_by_minmap
from .sequence_windows import intersect_phymap_locus
from .sequence_windows import resolve_sequence_windows

SEQEX_BATCH_SITES = 100_000
SEQEX_MIN_BATCH_SITES = 5_000


def _format_text_table(
    headers: tuple[str, ...],
    rows: list[tuple[str, ...]],
    right_align: set[int],
) -> str:
    """Format rows as a plain-text table with aligned columns."""
    widths = [
        max(len(header), *(len(row[idx]) for row in rows))
        for idx, header in enumerate(headers)
    ]

    def format_row(row: tuple[str, ...]) -> str:
        values = []
        for idx, (value, width) in enumerate(zip(row, widths)):
            values.append(
                value.rjust(width) if idx in right_align else value.ljust(width)
            )
        return "  ".join(values).rstrip()

    separator = tuple("-" * width for width in widths)
    return (
        "\n".join(
            [
                format_row(headers),
                format_row(separator),
                *(format_row(row) for row in rows),
            ]
        )
        + "\n"
    )


@dataclass(frozen=True)
class LocusSpec:
    """Coordinates needed to load and identify one locus or clipped fragment."""

    index: int
    scaffold: str
    phy0: int
    phy1: int
    pos0: int
    pos1: int
    source_pos0: int
    source_pos1: int
    selected_window: str
    clipped: bool

    @property
    def label(self) -> str:
        return f"{self.scaffold}:{self.pos0}-{self.pos1}"

    @property
    def length(self) -> int:
        return self.phy1 - self.phy0


@dataclass
class FilteredLocus:
    """One locus after site-wise and per-locus sample filtering."""

    spec: LocusSpec
    names: list[str]
    sequences: np.ndarray | None
    raw_samples: int
    raw_sites: int
    filtered_sites: int
    nonmissing_by_sample: dict[str, int]
    dropped_names: tuple[str, ...]


@dataclass
class FilterBatchResult:
    """Filtered loci and rejection counters returned by one batch."""

    loci: list[FilteredLocus]
    counts: dict[str, int]
    processed_loci: int


def _filter_locus_array(
    *,
    spec: LocusSpec,
    block: np.ndarray,
    snames: list[str],
    imap_row_indices: dict[str, np.ndarray],
    minmap: dict[str, int],
    max_sample_missing: float,
    min_length: int | None,
) -> tuple[FilteredLocus | None, str | None]:
    """Apply locus occupancy, site coverage, row missingness, and length filters."""
    present = np.any(block != MISSING_BASE, axis=1)
    if not all(
        int(np.sum(present[indices])) >= int(minmap[population])
        for population, indices in imap_row_indices.items()
    ):
        return None, "rejected_locus_coverage"

    filtered = filter_block_by_minmap(block, imap_row_indices, minmap)
    if filtered.shape[1] == 0:
        return None, "rejected_site_coverage"
    row_missing = np.mean(filtered == MISSING_BASE, axis=1)
    keep = row_missing <= max_sample_missing
    if not np.any(keep):
        return None, "rejected_sample_missing"
    if min_length is not None and filtered.shape[1] < min_length:
        return None, "rejected_filtered_length"

    names = [name for name, retain in zip(snames, keep) if retain]
    dropped = tuple(name for name, retain in zip(snames, keep) if not retain)
    kept = filtered[keep]
    nonmissing = np.sum(kept != MISSING_BASE, axis=1)
    return (
        FilteredLocus(
            spec=spec,
            names=names,
            sequences=kept,
            raw_samples=len(snames),
            raw_sites=block.shape[1],
            filtered_sites=filtered.shape[1],
            nonmissing_by_sample={
                name: int(nonmissing[idx]) for idx, name in enumerate(names)
            },
            dropped_names=dropped,
        ),
        None,
    )


def _filter_batch_from_phy(
    *,
    phy,
    specs: tuple[LocusSpec, ...],
    sidxs: list[int],
    snames: list[str],
    imap_row_indices: dict[str, np.ndarray],
    minmap: dict[str, int],
    max_sample_missing: float,
    min_length: int | None,
) -> FilterBatchResult:
    """Read and filter one physically contiguous batch of loci."""
    block = np.asarray(
        phy[sidxs, specs[0].phy0 : specs[-1].phy1],
        dtype=np.uint8,
    )
    block[block == 45] = MISSING_BASE
    loci = []
    counts: dict[str, int] = {}
    offset = 0
    for spec in specs:
        width = spec.length
        locus, rejection = _filter_locus_array(
            spec=spec,
            block=block[:, offset : offset + width],
            snames=snames,
            imap_row_indices=imap_row_indices,
            minmap=minmap,
            max_sample_missing=max_sample_missing,
            min_length=min_length,
        )
        offset += width
        if rejection is not None:
            counts[rejection] = counts.get(rejection, 0) + 1
        else:
            loci.append(locus)
    return FilterBatchResult(loci=loci, counts=counts, processed_loci=len(specs))


def _filter_locus_batch(
    *,
    data: Path,
    specs: tuple[LocusSpec, ...],
    sidxs: list[int],
    snames: list[str],
    imap_row_indices: dict[str, np.ndarray],
    minmap: dict[str, int],
    max_sample_missing: float,
    min_length: int | None,
) -> FilterBatchResult:
    """Process-pool entry point for one locus batch."""
    with h5py.File(data, "r") as io5:
        return _filter_batch_from_phy(
            phy=io5["phy"],
            specs=specs,
            sidxs=sidxs,
            snames=snames,
            imap_row_indices=imap_row_indices,
            minmap=minmap,
            max_sample_missing=max_sample_missing,
            min_length=min_length,
        )


def _safe_label(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", label).strip("_")


def _format_phylip(names: list[str], sequences: np.ndarray) -> str:
    width = max(len(name) for name in names) + 5
    rows = [
        f"{name.ljust(width)} {sequences[idx].tobytes().decode('utf-8')}"
        for idx, name in enumerate(names)
    ]
    return f"{len(names)} {sequences.shape[1]}\n{chr(10).join(rows)}\n"


def _format_fasta(
    names: list[str],
    sequences: np.ndarray,
    *,
    locus_label: str | None = None,
) -> str:
    rows = []
    for idx, name in enumerate(names):
        identifier = f"{name}|{locus_label}" if locus_label else name
        rows.append(f">{identifier}\n{sequences[idx].tobytes().decode('utf-8')}\n")
    return "".join(rows)


def _format_nexus_block(
    names: list[str],
    sequences: np.ndarray,
    *,
    title: str | None = None,
) -> str:
    width = max(len(name) for name in names) + 5
    lines = ["begin data;\n"]
    if title:
        escaped = title.replace("'", "''")
        lines.append(f"  title '{escaped}';\n")
    lines.extend(
        [
            f"  dimensions ntax={len(names)} nchar={sequences.shape[1]};\n",
            "  format datatype=dna missing=N gap=- interleave=yes;\n",
            "  matrix\n",
        ]
    )
    for start in range(0, sequences.shape[1], 100):
        stop = min(start + 100, sequences.shape[1])
        for idx, name in enumerate(names):
            seq = sequences[idx, start:stop].tobytes().decode("utf-8")
            lines.append(f"  {name.ljust(width)}{seq}\n")
        lines.append("\n")
    lines.append("  ;\nend;\n")
    return "".join(lines)


class _MultiLocusWriter:
    """Lazily write independent locus records to one file or stdout."""

    def __init__(self, path: Path, out_format: str, stdout: bool):
        self.path = path
        self.out_format = out_format
        self.stdout = stdout
        self.handle = None
        self.count = 0

    def _open(self):
        if self.handle is None:
            self.handle = (
                sys.stdout if self.stdout else self.path.open("w", encoding="utf-8")
            )
            if self.out_format == "nex":
                self.handle.write("#nexus\n")

    def write(self, locus: FilteredLocus, names: list[str]) -> None:
        self._open()
        if self.count and self.out_format == "phy":
            self.handle.write("\n")
        if self.out_format == "phy":
            text = _format_phylip(names, locus.sequences)
        elif self.out_format == "nex":
            text = _format_nexus_block(names, locus.sequences, title=locus.spec.label)
        else:
            text = _format_fasta(
                names,
                locus.sequences,
                locus_label=locus.spec.label,
            )
        self.handle.write(text)
        self.count += 1

    def close(self) -> None:
        if self.handle is not None and not self.stdout:
            self.handle.close()


class SeqexEngine:
    """Select, filter, sample, and write complete or clipped delimited loci."""

    def __init__(
        self,
        *,
        data,
        name,
        outdir,
        out_format,
        windows,
        max_loci,
        random_seed,
        min_length,
        min_sample_coverage,
        max_sample_missing,
        exclude=None,
        include_reference=False,
        imap=None,
        minmap=None,
        concatenate=False,
        split=False,
        append_population=False,
        stdout=False,
        force=False,
        logged_command=None,
        cores=1,
        log_level="INFO",
        clip: bool | None = None,
    ):
        self.data = Path(data)
        self.name = str(name)
        self.outdir = Path(outdir).expanduser().absolute()
        self.out_format = out_format
        self.max_loci = max_loci
        self.random_seed = random_seed
        self.min_length = min_length
        self.concatenate = bool(concatenate)
        self.split = bool(split)
        self.append_population = bool(append_population)
        self.stdout = bool(stdout)
        self.force = bool(force)
        self.logged_command = logged_command
        self.cores = int(cores)
        self.log_level = str(log_level)
        if windows is None:
            self.requested_windows = None
        elif isinstance(windows, str):
            self.requested_windows = [windows]
        else:
            self.requested_windows = [str(window) for window in windows]
        if clip not in (None, True, False):
            raise IPyradError("clip must be one of None, True, or False.")
        self.clip = clip
        self.context = build_sequence_extraction_context(
            data=data,
            min_sample_coverage=min_sample_coverage,
            max_sample_missing=max_sample_missing,
            exclude=exclude,
            include_reference=include_reference,
            imap=imap,
            minmap=minmap,
        )
        self.selected_windows = resolve_sequence_windows(
            self.context.scaffold_table,
            windows,
        )
        with h5py.File(self.data, "r") as io5:
            self.phymap = np.asarray(io5["phymap"], dtype=np.int64)
        self.coordinate_clipping_applied = self.clip is not False and any(
            window.explicit_coordinates for window in self.selected_windows
        )
        self.sample_to_population = {
            sample: population
            for population, samples in self.context.imap.items()
            for sample in samples
        }
        if self.append_population:
            invalid_samples = sorted(
                name for name in self.context.sample_names if "^" in name
            )
            invalid_populations = sorted(
                population
                for population in set(self.sample_to_population.values())
                if "^" in population
            )
            if invalid_samples or invalid_populations:
                details = []
                if invalid_samples:
                    details.append("samples=" + ", ".join(invalid_samples))
                if invalid_populations:
                    details.append("populations=" + ", ".join(invalid_populations))
                raise IPyradError(
                    "--append-population requires sample and population names without "
                    "the ^ delimiter (" + "; ".join(details) + ")."
                )
        self.counts = {
            "candidate_loci": 0,
            "rejected_raw_length": 0,
            "rejected_locus_coverage": 0,
            "rejected_site_coverage": 0,
            "rejected_sample_missing": 0,
            "rejected_filtered_length": 0,
            "accepted_before_sampling": 0,
            "written_loci": 0,
        }

    @property
    def alignment_path(self) -> Path:
        suffix = {"phy": "phy", "nex": "nex", "fa": "fa"}[self.out_format]
        return self.outdir / f"{self.name}.{suffix}"

    @property
    def stats_path(self) -> Path:
        return self.outdir / f"{self.name}.stats.txt"

    @property
    def stats_json_path(self) -> Path:
        return self.outdir / f"{self.name}.stats.json"

    def _output_names(self, names: list[str]) -> list[str]:
        if not self.append_population:
            return names
        return [f"{self.sample_to_population[name]}^{name}" for name in names]

    def _candidate_loci(self) -> list[LocusSpec]:
        specs = []
        scaffold_names = (
            self.context.scaffold_table["scaffold_name"].astype(str).tolist()
        )
        windows_by_scaffold = {}
        for window in self.selected_windows:
            windows_by_scaffold.setdefault(window.scaffold_index, []).append(window)
        complete_rows = set()
        for row_index, row in enumerate(self.phymap):
            scaff_idx = int(row[0])
            scaffold = scaffold_names[scaff_idx]
            source_pos0 = int(row[3])
            source_pos1 = int(row[4])
            for window in windows_by_scaffold.get(scaff_idx, []):
                should_clip = self.clip is True or (
                    self.clip is None and window.explicit_coordinates
                )
                intersection = intersect_phymap_locus(
                    row,
                    window,
                    clip=should_clip,
                )
                if intersection is None:
                    continue
                if not should_clip:
                    if row_index in complete_rows:
                        continue
                    complete_rows.add(row_index)
                specs.append(
                    LocusSpec(
                        index=row_index,
                        scaffold=scaffold,
                        phy0=intersection.phy0,
                        phy1=intersection.phy1,
                        pos0=intersection.pos0,
                        pos1=intersection.pos1,
                        source_pos0=source_pos0,
                        source_pos1=source_pos1,
                        selected_window=window.label,
                        clipped=intersection.clipped,
                    )
                )
        self.counts["candidate_loci"] = len(specs)
        if not specs:
            raise IPyradError("No loci overlap the selected windows.")
        return specs

    def _plan_filter_batches(
        self,
        specs: list[LocusSpec],
    ) -> list[tuple[LocusSpec, ...]]:
        """Return bounded batches of physically contiguous loci."""
        readable = []
        for spec in specs:
            if self.min_length is not None and spec.length < self.min_length:
                self.counts["rejected_raw_length"] += 1
            else:
                readable.append(spec)
        total_sites = sum(spec.length for spec in readable)
        target_sites = SEQEX_BATCH_SITES
        if self.cores > 1 and total_sites:
            target_sites = min(
                SEQEX_BATCH_SITES,
                max(
                    SEQEX_MIN_BATCH_SITES,
                    total_sites // (4 * self.cores),
                ),
            )
        batches = []
        start = 0
        while start < len(readable):
            stop = start + 1
            sites = readable[start].length
            while stop < len(readable):
                previous = readable[stop - 1]
                current = readable[stop]
                if current.phy0 != previous.phy1:
                    break
                if sites + current.length > target_sites:
                    break
                sites += current.length
                stop += 1
            batches.append(tuple(readable[start:stop]))
            start = stop
        return batches

    def _batch_kwargs(self, specs: tuple[LocusSpec, ...]) -> dict:
        return {
            "specs": specs,
            "sidxs": self.context.sample_indices,
            "snames": self.context.sample_names,
            "imap_row_indices": self.context.imap_row_indices,
            "minmap": self.context.minmap,
            "max_sample_missing": self.context.max_sample_missing,
            "min_length": self.min_length,
        }

    def _record_batch_counts(self, result: FilterBatchResult) -> None:
        for key, value in result.counts.items():
            self.counts[key] += value

    def _iter_filtered_batches(
        self,
        batches: list[tuple[LocusSpec, ...]],
    ):
        """Yield filter results in genomic order for serial or parallel runs."""
        if self.cores == 1 or len(batches) < 2:
            with h5py.File(self.data, "r") as io5:
                phy = io5["phy"]
                for specs in batches:
                    yield _filter_batch_from_phy(
                        phy=phy,
                        **self._batch_kwargs(specs),
                    )
            return

        jobs = (
            (
                idx,
                (
                    _filter_locus_batch,
                    {"data": self.data, **self._batch_kwargs(specs)},
                ),
            )
            for idx, specs in enumerate(batches)
        )
        pending: dict[int, FilterBatchResult] = {}
        next_idx = 0
        for idx, result in run_with_pool_iter(
            jobs,
            self.log_level,
            self.cores,
            max_inflight=2 * self.cores,
            msg="Filtering seqex loci",
            njobs=sum(len(batch) for batch in batches),
            progress_increment=lambda _key, value: value.processed_loci,
        ):
            pending[idx] = result
            while next_idx in pending:
                yield pending.pop(next_idx)
                next_idx += 1

    def _split_path(self, spec: LocusSpec) -> Path:
        suffix = {"phy": "phy", "nex": "nex", "fa": "fa"}[self.out_format]
        return self.outdir / f"{self.name}.{_safe_label(spec.label)}.{suffix}"

    def _check_outputs(self, specs: list[LocusSpec]) -> None:
        paths = [self.stats_path, self.stats_json_path]
        if not self.stdout:
            if self.split:
                paths.extend(self._split_path(spec) for spec in specs)
            else:
                paths.append(self.alignment_path)
        if not self.force:
            existing = next((path for path in paths if path.exists()), None)
            if existing is not None:
                raise IPyradError(
                    f"Output file already exists: {existing}. Use --force to overwrite."
                )

    def _write_one_file(self, path: Path, locus: FilteredLocus) -> None:
        names = self._output_names(locus.names)
        if self.out_format == "phy":
            text = _format_phylip(names, locus.sequences)
        elif self.out_format == "nex":
            text = "#nexus\n" + _format_nexus_block(names, locus.sequences)
        else:
            text = _format_fasta(names, locus.sequences)
        path.write_text(text, encoding="utf-8")

    def _build_concatenation(
        self,
        loci: list[FilteredLocus],
    ) -> tuple[list[str], np.ndarray, list[dict], dict[str, float]]:
        buffers = {name: bytearray() for name in self.context.sample_names}
        present = set()
        partitions = []
        position = 1
        for locus in loci:
            width = locus.filtered_sites
            by_name = {
                name: locus.sequences[idx].tobytes()
                for idx, name in enumerate(locus.names)
            }
            missing = b"N" * width
            for name in self.context.sample_names:
                buffers[name].extend(by_name.get(name, missing))
            present.update(locus.names)
            partitions.append(
                {
                    "locus": locus.spec.label,
                    "concat_start": position,
                    "concat_end": position + width - 1,
                }
            )
            position += width
        names = [name for name in self.context.sample_names if name in present]
        matrix = np.vstack(
            [np.frombuffer(bytes(buffers[name]), dtype=np.uint8) for name in names]
        )
        global_missing = np.mean(matrix == MISSING_BASE, axis=1)
        keep = global_missing <= self.context.max_sample_missing
        if not np.any(keep):
            raise IPyradError(
                "No samples passed max_sample_missing across the "
                "concatenated alignment."
            )
        dropped = {
            name: float(missing)
            for name, missing, retain in zip(names, global_missing, keep)
            if not retain
        }
        return (
            [name for name, retain in zip(names, keep) if retain],
            matrix[keep],
            partitions,
            dropped,
        )

    def _write_concatenated(self, loci: list[FilteredLocus]):
        names, matrix, partitions, dropped = self._build_concatenation(loci)
        output_names = self._output_names(names)
        if self.out_format == "phy":
            text = _format_phylip(output_names, matrix)
        elif self.out_format == "nex":
            text = "#nexus\n" + _format_nexus_block(output_names, matrix)
        else:
            text = _format_fasta(output_names, matrix)
        if self.stdout:
            sys.stdout.write(text)
        else:
            self.alignment_path.write_text(text, encoding="utf-8")
        return partitions, names, dropped

    def _summarize_output(
        self,
        loci: list[FilteredLocus],
        final_names: list[str],
    ) -> tuple[dict[str, int | float], list[dict[str, int | float | str]]]:
        """Return output-wide and per-sample occupancy statistics."""
        total_sites = sum(locus.filtered_sites for locus in loci)
        sample_rows = {
            name: {
                "sample": name,
                "population": self.sample_to_population.get(name, "unassigned"),
                "loci_written": 0,
                "loci_dropped_by_r": 0,
                "non_missing_bases": 0,
            }
            for name in self.context.sample_names
        }
        total_bases_written = 0
        for locus in loci:
            total_bases_written += len(locus.names) * locus.filtered_sites
            retained = set(locus.names)
            dropped = set(locus.dropped_names)
            for name, row in sample_rows.items():
                if name in retained:
                    row["loci_written"] += 1
                    row["non_missing_bases"] += locus.nonmissing_by_sample[name]
                if name in dropped:
                    row["loci_dropped_by_r"] += 1

        if self.concatenate:
            total_bases_written = len(final_names) * total_sites
        final_set = set(final_names)
        full_matrix_bases = len(final_names) * total_sites
        non_missing_bases = sum(
            int(sample_rows[name]["non_missing_bases"]) for name in final_names
        )
        occupancy = non_missing_bases / full_matrix_bases if full_matrix_bases else 0.0
        sample_counts = (
            [len(final_names)]
            if self.concatenate
            else [len(locus.names) for locus in loci]
        )
        output_summary = {
            "total_sites_written": total_sites,
            "total_bases_written": total_bases_written,
            "full_matrix_bases": full_matrix_bases,
            "non_missing_bases": non_missing_bases,
            "non_missing_occupancy": occupancy,
            "max_samples": max(sample_counts, default=0),
            "mean_samples": (
                sum(sample_counts) / len(sample_counts) if sample_counts else 0.0
            ),
        }
        rows = []
        for name in self.context.sample_names:
            row = sample_rows[name]
            matrix_bases = total_sites
            sample_occupancy = (
                int(row["non_missing_bases"]) / matrix_bases if matrix_bases else 0.0
            )
            rows.append(
                {
                    **row,
                    "written_final": name in final_set,
                    "matrix_bases": matrix_bases,
                    "non_missing_occupancy": sample_occupancy,
                }
            )
        return output_summary, rows

    def _write_stats(
        self,
        loci: list[FilteredLocus],
        partitions: list[dict],
        output_summary: dict[str, int | float],
        sample_rows: list[dict[str, bool | int | float | str]],
    ) -> None:
        output_layout = (
            "concatenated"
            if self.concatenate
            else "split"
            if self.split
            else "multi-locus"
        )
        seqex_summary = {
            "command": self.logged_command,
            "data": str(self.data),
            "output_layout": output_layout,
            "out_format": self.out_format,
            "cores": self.cores,
            "max_loci": self.max_loci,
            "random_seed": self.random_seed,
            "min_length": self.min_length,
            "clipping_mode": (
                "automatic"
                if self.clip is None
                else "always"
                if self.clip
                else "never"
            ),
            "coordinate_clipping_applied": self.coordinate_clipping_applied,
            "selected_windows": self.requested_windows,
        }
        partition_map = {row["locus"]: row for row in partitions}
        locus_rows = []
        for idx, locus in enumerate(loci, 1):
            partition = partition_map.get(locus.spec.label, {})
            locus_rows.append(
                {
                    "locus_index": idx,
                    "locus": locus.spec.label,
                    "source_locus": (
                        f"{locus.spec.scaffold}:"
                        f"{locus.spec.source_pos0}-{locus.spec.source_pos1}"
                    ),
                    "selected_window": locus.spec.selected_window,
                    "clipped": locus.spec.clipped,
                    "raw_samples": locus.raw_samples,
                    "raw_sites": locus.raw_sites,
                    "filtered_samples": len(locus.names),
                    "filtered_sites": locus.filtered_sites,
                    "concat_start": partition.get("concat_start"),
                    "concat_end": partition.get("concat_end"),
                }
            )

        if self.requested_windows is None:
            selected_windows_text = "none"
        elif self.requested_windows:
            selected_windows_text = ", ".join(self.requested_windows)
        else:
            selected_windows_text = "[]"

        lines = ["# Seqex Summary\n"]
        if self.logged_command:
            lines.append(f"command: {self.logged_command}\n")
        lines.extend(
            [
                f"data: {self.data}\n",
                f"output_layout: {output_layout}\n",
                f"out_format: {self.out_format}\n",
                f"cores: {self.cores}\n",
                f"max_loci: {self.max_loci if self.max_loci is not None else 'all'}\n",
                f"random_seed: {self.random_seed if self.random_seed is not None else 'none'}\n",
                f"min_length: {self.min_length if self.min_length is not None else 'none'}\n",
                f"clipping_mode: {seqex_summary['clipping_mode']}\n",
                "coordinate_clipping_applied: "
                f"{str(self.coordinate_clipping_applied).lower()}\n",
                f"windows_selected: {len(self.selected_windows)}\n",
                f"selected_windows: {selected_windows_text}\n",
            ]
        )
        lines.extend(f"{key}: {value}\n" for key, value in self.counts.items())
        lines.append("\n# Output Summary\n")
        lines.extend(
            [
                f"total_sites_written: {output_summary['total_sites_written']}\n",
                f"total_bases_written: {output_summary['total_bases_written']}\n",
                f"full_matrix_bases: {output_summary['full_matrix_bases']}\n",
                f"non_missing_bases: {output_summary['non_missing_bases']}\n",
                "non_missing_occupancy: "
                f"{float(output_summary['non_missing_occupancy']):.6f}\n",
                f"max_samples: {output_summary['max_samples']}\n",
                f"mean_samples: {float(output_summary['mean_samples']):.6f}\n",
            ]
        )
        lines.append("\n# Sample Occupancy\n")
        sample_headers = (
            "sample",
            "population",
            "written_final",
            "loci_written",
            "loci_dropped_by_r",
            "matrix_bases",
            "non_missing_bases",
            "non_missing_occupancy",
        )
        sample_table_rows = [
            (
                str(row["sample"]),
                str(row["population"]),
                "yes" if row["written_final"] else "no",
                str(row["loci_written"]),
                str(row["loci_dropped_by_r"]),
                str(row["matrix_bases"]),
                str(row["non_missing_bases"]),
                f"{float(row['non_missing_occupancy']):.6f}",
            )
            for row in sample_rows
        ]
        lines.append(
            _format_text_table(
                sample_headers,
                sample_table_rows,
                right_align={3, 4, 5, 6, 7},
            )
        )
        lines.append("\n# Written Loci\n")
        locus_headers = (
            "locus_index",
            "locus",
            "source_locus",
            "selected_window",
            "clipped",
            "raw_samples",
            "raw_sites",
            "filtered_samples",
            "filtered_sites",
            "concat_start",
            "concat_end",
        )
        locus_table_rows = [
            tuple(
                ""
                if row[key] is None
                else "yes"
                if key == "clipped" and row[key]
                else "no"
                if key == "clipped"
                else str(row[key])
                for key in locus_headers
            )
            for row in locus_rows
        ]
        lines.append(
            _format_text_table(
                locus_headers,
                locus_table_rows,
                right_align={0, 5, 6, 7, 8, 9, 10},
            )
        )
        self.stats_path.write_text("".join(lines), encoding="utf-8")
        stats_data = {
            "seqex_summary": seqex_summary,
            "filter_counts": self.counts,
            "output_summary": output_summary,
            "sample_occupancy": sample_rows,
            "written_loci": locus_rows,
        }
        self.stats_json_path.write_text(
            json.dumps(stats_data, indent=2) + "\n",
            encoding="utf-8",
        )

    def _log_results(
        self,
        loci: list[FilteredLocus],
        global_dropped: dict[str, float],
    ) -> None:
        """Log sample drops and describe the output destinations."""
        drop_counts = {name: 0 for name in self.context.sample_names}
        for locus in loci:
            for name in locus.dropped_names:
                drop_counts[name] += 1
        for name, count in drop_counts.items():
            if count:
                logger.debug(
                    "-r dropped {} from {} written locus/loci",
                    name,
                    count,
                )
        for name, missing in global_dropped.items():
            logger.info(
                "-r dropped {} from the concatenated alignment "
                "(missing={:.6f}, maximum={:.6f})",
                name,
                missing,
                self.context.max_sample_missing,
            )

        fmt = {"phy": "PHYLIP", "nex": "NEXUS", "fa": "FASTA"}[self.out_format]
        if self.stdout:
            destination = "stdout"
        elif self.split:
            destination = str(self.outdir)
        else:
            destination = str(self.alignment_path)
        if self.concatenate:
            description = f"concatenated into one {fmt} alignment"
        elif self.split:
            description = f"as separate {fmt} files"
        else:
            description = f"as independent records in one {fmt} file"
        logger.info(
            "wrote {} filtered loci {} to: {}",
            len(loci),
            description,
            destination,
        )
        logger.info("wrote stats report to: {}", self.stats_path)

    def run(self) -> list[FilteredLocus]:
        specs = self._candidate_loci()
        self._check_outputs(specs)
        self.outdir.mkdir(parents=True, exist_ok=True)
        batches = self._plan_filter_batches(specs)
        rng = np.random.default_rng(self.random_seed)
        reservoir: list[FilteredLocus] = []
        all_loci: list[FilteredLocus] = []
        stream_output = self.max_loci is None and not self.concatenate
        stream_writer = (
            _MultiLocusWriter(self.alignment_path, self.out_format, self.stdout)
            if stream_output and not self.split
            else None
        )

        try:
            for result in self._iter_filtered_batches(batches):
                self._record_batch_counts(result)
                for locus in result.loci:
                    self.counts["accepted_before_sampling"] += 1
                    if stream_output:
                        if self.split:
                            self._write_one_file(self._split_path(locus.spec), locus)
                        else:
                            stream_writer.write(
                                locus,
                                self._output_names(locus.names),
                            )
                        locus.sequences = None
                        all_loci.append(locus)
                        continue
                    if self.max_loci is None:
                        all_loci.append(locus)
                        continue
                    accepted = self.counts["accepted_before_sampling"]
                    if len(reservoir) < self.max_loci:
                        reservoir.append(locus)
                    else:
                        target = int(rng.integers(0, accepted))
                        if target < self.max_loci:
                            reservoir[target] = locus
        finally:
            if stream_writer is not None:
                stream_writer.close()

        selected = all_loci if self.max_loci is None else reservoir
        selected.sort(key=lambda locus: locus.spec.index)
        if not selected:
            raise IPyradError(
                "No loci passed coverage, missingness, and minimum-length filters."
            )
        if self.max_loci is not None and len(selected) < self.max_loci:
            logger.warning(
                "Requested {} loci, but only {} passed filtering.",
                self.max_loci,
                len(selected),
            )

        partitions = []
        global_dropped: dict[str, float] = {}
        final_names = [
            name
            for name in self.context.sample_names
            if any(name in locus.names for locus in selected)
        ]
        if stream_output:
            pass
        elif self.concatenate:
            partitions, final_names, global_dropped = self._write_concatenated(selected)
        elif self.split:
            for locus in selected:
                self._write_one_file(self._split_path(locus.spec), locus)
        else:
            writer = _MultiLocusWriter(
                self.alignment_path,
                self.out_format,
                self.stdout,
            )
            try:
                for locus in selected:
                    writer.write(locus, self._output_names(locus.names))
            finally:
                writer.close()
        self.counts["written_loci"] = len(selected)
        output_summary, sample_rows = self._summarize_output(
            selected,
            final_names,
        )
        self._write_stats(
            selected,
            partitions,
            output_summary,
            sample_rows,
        )
        self._log_results(selected, global_dropped)
        return selected


def run_seqex(**kwargs):
    """Validate arguments and run the locus-oriented sequence extractor."""
    print_scaffold_table = bool(kwargs.pop("print_scaffold_table", False))
    max_loci = kwargs.get("max_loci")
    min_length = kwargs.get("min_length")
    random_seed = kwargs.get("random_seed")
    cores = kwargs.get("cores", 1)
    if max_loci is not None and max_loci < 1:
        raise IPyradError("--max-loci must be at least 1.")
    if min_length is not None and min_length < 1:
        raise IPyradError("--min-length must be at least 1.")
    if random_seed is not None and random_seed < 0:
        raise IPyradError("--random-seed must be a non-negative integer.")
    if random_seed is not None and max_loci is None:
        raise IPyradError("--random-seed requires --max-loci.")
    if cores < 1:
        raise IPyradError("--cores must be at least 1.")
    if kwargs.get("concatenate") and kwargs.get("split"):
        raise IPyradError("--concatenate and --split are mutually exclusive.")
    if kwargs.get("split") and kwargs.get("stdout"):
        raise IPyradError("--split cannot be combined with --stdout.")
    if kwargs.get("append_population") and kwargs.get("imap") is None:
        raise IPyradError("--append-population requires --imap.")

    engine = SeqexEngine(**kwargs)
    if print_scaffold_table:
        try:
            engine.context.scaffold_table.to_csv(sys.stdout, sep="\t")
            sys.stdout.flush()
        except BrokenPipeError:
            sys.stdout = open(os.devnull, "w", encoding="utf-8")
        return None
    return engine.run()
