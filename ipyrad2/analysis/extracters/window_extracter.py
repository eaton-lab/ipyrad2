#!/usr/bin/env python

"""Extract/subset sequences from HDF5 database and write to a supermatrix.

Command
-------
$ ipyrad2 wex -d ... -w ... -o ... -O phy

Output file example
-------------------
10 100
A    NNNNNATCCGAGA...
B    NNNNNNNNNNNNN...
C    CCAGGATCCGAAA...
D    CCAGGATCCGAAA...
...

Stats file example
------------------
CMD: ipyrad2 wex -d ... -o ... ...
windows: Chr1:X-Y Chr1:A-B ...
populations: A B C
min_sample_coverage: A=1 B=2 C=3
max_sample_missing: 1.0
nsamples_before_filtering: 29
nsites_in_windows_before_filtering: 1000
nvariant_sites_in_windows_before_filtering: 100
nsamples_after_filtering: 29
nsites_in_windows_after_filtering: 300
nvariant_sites_in_windows_after_filtering: 20
outfile: alignment.phy
"""

from typing import Dict, List, Tuple
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
from loguru import logger

from ...utils.exceptions import IPyradError
from .sequence_common import build_sequence_imap_minmap
from .sequence_common import load_sequence_scaffold_table
from .sequence_common import load_sequence_chunk_from_phy
from .sequence_common import normalize_sequence_population_inputs
from .sequence_common import plan_sequence_chunk_spans
from .sequence_common import resolve_sequence_sample_subset


NEXHEADER = """#nexus
begin data;
  dimensions ntax={} nchar={};
  format datatype=dna missing=N gap=- interleave=yes;
  matrix
"""
REFERENCE_SAMPLE_NAME = "assembly_reference_sequence"
MISSING_BASE = 78


def _format_count(value: int) -> str:
    """Format integer counts consistently for text reports."""
    return f"{int(value):,}"


def _format_float(value: int | float, digits: int = 3) -> str:
    """Format floating-point values for text reports."""
    return f"{float(value):.{digits}f}"


def _format_fraction(value: int | float) -> str:
    """Format fraction-like values consistently for text reports."""
    return f"{float(value):.6f}"


def _format_percent(value: int | float, digits: int = 3) -> str:
    """Format one fraction as a percent string without a percent sign."""
    return f"{100 * float(value):.{digits}f}"


def _append_key_value_section(lines: list[str], title: str, rows: list[tuple[str, str]]) -> None:
    """Append one key/value report section to a list of output lines."""
    lines.append(f"# {title}")
    if rows:
        width = max(len(key) for key, _ in rows)
        for key, value in rows:
            lines.append(f"{key.ljust(width)}  {value}")
    lines.append("")


def _append_table_section(lines: list[str], title: str, headers: list[str], rows: list[list[str]]) -> None:
    """Append one simple whitespace-aligned table section."""
    lines.append(f"# {title}")
    widths = [len(header) for header in headers]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))
    lines.append("  ".join(header.ljust(widths[idx]) for idx, header in enumerate(headers)))
    for row in rows:
        lines.append("  ".join(value.ljust(widths[idx]) for idx, value in enumerate(row)))
    lines.append("")


def filter_block_by_minmap(
    block: np.ndarray,
    imap_row_indices: dict[str, np.ndarray],
    minmap: dict[str, int | float],
) -> np.ndarray:
    """Return one block after applying the per-population coverage filter."""
    if not block.size:
        return block[:, 0:0]

    mask = np.zeros(block.shape[1], dtype=np.bool_)
    for pop, pop_sidxs in imap_row_indices.items():
        pop_mincov = minmap[pop]
        mask |= np.sum(block[pop_sidxs, :] != MISSING_BASE, axis=0) < pop_mincov
    return block[:, np.invert(mask)]


class WindowExtracter:
    """Tool to extract sequences from one or more loci and write to a
    concatenated sequence file in phylip or nexus format.
    """
    def __init__(
        self,
        data: str,
        name: str,
        outdir: Path | str,
        out_format: str,
        windows: str | List[str],
        min_sample_coverage: int | float,
        max_sample_missing: float,
        exclude: List[str] | None = None,
        include_reference: bool = False,
        imap: Dict[str, List[str]] | None = None,
        minmap: Dict[str, int | float] | None = None,
        stdout: bool = False,
        force: bool = False,
        logged_command: str | None = None,
    ):
        # store params
        imap, minmap = normalize_sequence_population_inputs(imap, minmap)
        self.data = data
        self.name = name
        self.outdir = Path(outdir).expanduser().absolute()
        self.out_format = out_format
        self.windows = [] if windows is None else [windows] if isinstance(windows, str) else list(windows)
        self.exclude = set(exclude if exclude else [])
        self.include_reference = include_reference
        self.min_sample_coverage = min_sample_coverage
        self.max_sample_missing = min(1.0, max(0, max_sample_missing))
        self.stdout = stdout
        self.force = force
        self.logged_command = logged_command

        # data parsed from h5
        self.scaffold_table: pd.DataFrame = None
        self.snames: List[str] = None
        self.sidxs: List[str] = None
        self.pnames: Dict[str, str] = None
        self.phymap: pd.DataFrame = None
        self.phymap_windows: Dict[int, List[Tuple[int, int]]] = None
        self.imap: Dict[str, List[str]] = {}
        self.minmap: Dict[str, int] = {}
        self.selected_windows: List[str] = []
        self.seqarr: np.ndarray | None = None
        self._selected_phy_spans: tuple[tuple[int, int], ...] = tuple()
        self._chunk_spans: list[tuple[tuple[int, int], ...]] = []
        self._imap_row_indices: Dict[str, np.ndarray] = {}

        # fills: snames, sidxs, scaffold_table
        self.scaffold_table = load_sequence_scaffold_table(self.data)
        self.snames, self.sidxs, self.exclude, imap = resolve_sequence_sample_subset(
            self.data,
            exclude=self.exclude,
            include_reference=self.include_reference,
            imap=imap,
        )
        self.imap, self.minmap = build_sequence_imap_minmap(
            self.snames,
            min_sample_coverage=self.min_sample_coverage,
            imap=imap,
            minmap=minmap,
        )
        self._imap_row_indices = {
            pop: np.array([self.snames.index(name) for name in names], dtype=np.int64)
            for pop, names in self.imap.items()
        }

        # run commands
    def _run(self):
        # First two are fast
        self._get_phymap_windows()
        self._get_phymap()
        _, fnames, fseqarr = self._extract_filtered_alignment()
        return fnames, fseqarr

    def _get_scaffold_table(self) -> None:
        """Store table with scaffold names and lengths in the order they are stored in H5."""
        self.scaffold_table = load_sequence_scaffold_table(self.data)

    def _get_snames_and_sidxs_subset(self, imap) -> None:
        self.snames, self.sidxs, self.exclude, self.imap = resolve_sequence_sample_subset(
            self.data,
            exclude=self.exclude,
            include_reference=self.include_reference,
            imap=imap,
        )
        self._imap_row_indices = {
            pop: np.array([self.snames.index(name) for name in names], dtype=np.int64)
            for pop, names in self.imap.items()
        }

    def _parse_imap_minmap_inputs(self, imap, minmap):
        """Normalize imap/minmap inputs from dicts or files."""
        return normalize_sequence_population_inputs(imap, minmap)

    def _get_imap_minmap(self, imap, minmap):
        """Set _imap and _minmap for seqarr filtering."""
        self.imap, self.minmap = build_sequence_imap_minmap(
            self.snames,
            min_sample_coverage=self.min_sample_coverage,
            imap=imap,
            minmap=minmap,
        )
        self._imap_row_indices = {
            pop: np.array([self.snames.index(name) for name in names], dtype=np.int64)
            for pop, names in self.imap.items()
        }

    @staticmethod
    def _windows_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
        """Return True if 1-based inclusive windows overlap."""
        return not (end1 < start2 or end2 < start1)

    def _add_window(
        self,
        windows: Dict[str, List[Tuple[int, int]]],
        selected_windows: List[str],
        scaff: str,
        start: int,
        end: int,
        source: str,
    ) -> None:
        """Validate and store a 1-based inclusive window."""
        if start < 1 or end < start:
            raise IPyradError(
                f"Malformed window '{source}'. Windows must use valid positive coordinates."
            )

        existing = windows.setdefault(scaff, [])
        for existing_start, existing_end in existing:
            if self._windows_overlap(start, end, existing_start, existing_end):
                raise IPyradError(
                    f"windows cannot overlap. {source} overlaps "
                    f"{scaff}:{existing_start}-{existing_end}"
                )
        existing.append((start, end))
        selected_windows.append(f"{scaff}:{start}-{end}")

    def _get_phymap_windows(self) -> None:
        """Check each window for a matching scaffold name, and position within its bounds."""
        self._reset_selection_cache()
        windows: Dict[str, List[Tuple[int, int]]] = {}
        selected_windows: List[str] = []

        # set names in index for easy fetching
        t = self.scaffold_table.set_index("scaffold_name")

        if not self.windows:
            logger.info(
                "No windows specified; selecting the full length of all scaffolds. "
                "Use -w to subset scaffold windows and -P to view scaffold names."
            )
            self.windows = [r".*"]

        # Load windows from bed file if they are passed in this way
        if len(self.windows) == 1:
            bedfile = Path(self.windows[0])
            if bedfile.exists():
                logger.info(f"Loading windows from bed file: '{bedfile}'")
                for scaff, start, end in self._get_windows_from_bed(bedfile):
                    if scaff not in t.index:
                        raise IPyradError(
                            f"'{scaff}' from {bedfile} does not match to any scaffold names. Check with '-P'."
                        )
                    self._add_window(
                        windows,
                        selected_windows,
                        scaff,
                        start,
                        end,
                        f"{scaff}:{start}-{end}",
                    )
                self.selected_windows = selected_windows
                logger.debug(f"windows: {windows}")
                nwindows = sum(len(i) for i in windows.values())
                ws = 's' if nwindows > 1 else ''
                ss = 's' if len(windows) > 1 else ''
                logger.info(f"selected {nwindows} window{ws} from {len(windows)} scaffold{ss}")
                scaff_names = t.index.tolist()
                scaff_to_idx = {name: idx for idx, name in enumerate(scaff_names)}
                self.phymap_windows = {scaff_to_idx[i]: j for i, j in windows.items()}
                return

            logger.debug("Loading windows from command line arguments")

        # iterate over user-entered windows
        for window in self.windows:

            # sub-scaffold window
            if ":" in window:
                scaff, region = window.split(":", 1)
                mask = t.index.str.fullmatch(pat=scaff, na=False)
                scaffs = t.index[mask].tolist()
                if not scaffs:
                    raise IPyradError(
                        f"No scaffold names match '{window.split(':')[0]}'. Use -P to view scaffold names."
                    )
                if len(scaffs) > 1:
                    raise IPyradError("Cannot use regex with ':'. List windows separately: -w Chr1:1-1000 Chr2:1-1000")
                if region.count("-") != 1:
                    raise IPyradError(f"malformatted window '{window}'. Must be {{scaff}} or {{scaff}}:{{start}}-{{end}}")
                start, end = [int(i) for i in region.split("-")]
                self._add_window(windows, selected_windows, scaffs[0], start, end, window)

            # full scaffold window
            else:
                mask = t.index.str.fullmatch(pat=window, na=False)
                scaffs = t.index[mask].values
                if not scaffs.size:
                    raise IPyradError(f"'{window}' does not match to any scaffold names. Check with '-P'.")
                for scaff in scaffs:
                    if scaff not in windows:
                        length = int(t.loc[scaff, "scaffold_length"])
                        self._add_window(windows, selected_windows, scaff, 1, length, window)
                    else:
                        raise IPyradError(f"windows cannot overlap. {window} & {windows}")

        # log to INFO and DEBUG
        self.selected_windows = selected_windows
        logger.debug(f"windows: {windows}")
        nwindows = sum(len(i) for i in windows.values())
        ws = 's' if nwindows > 1 else ''
        ss = 's' if len(windows) > 1 else ''
        logger.info(f"selected {nwindows} window{ws} from {len(windows)} scaffold{ss}")

        # store as dict mapping {scaff_index: window, ...}
        scaff_names = t.index.tolist()
        scaff_to_idx = {name: idx for idx, name in enumerate(scaff_names)}
        self.phymap_windows = {
            scaff_to_idx[i]: j
            for i, j in windows.items()
        }


    def _get_windows_from_bed(self, bedfile: Path) -> List[Tuple[str, int, int]]:
        """Read windows from a BED file.

        BED uses 0-based, half-open coordinates. These are converted to the
        extracter's 1-based inclusive region semantics.
        """
        windows: List[Tuple[str, int, int]] = []
        with open(bedfile, "r", encoding="utf-8") as infile:
            for lineno, line in enumerate(infile, start=1):
                # Ignore comments and blank lines
                if line.startswith("#") or line.strip() == "":
                    continue

                chrom, start, end, *rest = line.rstrip("\t\n").split()
                start0 = int(start)
                end0 = int(end)
                if start0 < 0 or end0 <= start0:
                    raise IPyradError(
                        f"Malformed BED window at line {lineno}: {bedfile}"
                    )
                windows.append((chrom, start0 + 1, end0))
        return windows


    def _get_phymap(self) -> None:
        """Load the phymap for selecting windows from the seqs array."""
        self._reset_selection_cache()
        with h5py.File(self.data, 'r') as io5:
            phymap = io5["phymap"]
            colnames = phymap.attrs["columns"]
            mask = np.isin(phymap[:, 0], list(self.phymap_windows))
            phymap = pd.DataFrame(phymap[mask], columns=colnames)
        self.phymap = phymap

    def _reset_selection_cache(self) -> None:
        """Drop cached selection spans after window or phymap changes."""
        self.seqarr = None
        self._selected_phy_spans = tuple()
        self._chunk_spans = []

    def _get_selected_phy_spans(self) -> tuple[tuple[int, int], ...]:
        """Return ordered phy spans for the current window selection."""
        if self._selected_phy_spans:
            return self._selected_phy_spans

        rows_by_scaffold: Dict[int, List[np.ndarray]] = {}
        for row in self.phymap.to_numpy(copy=False):
            rows_by_scaffold.setdefault(int(row[0]), []).append(row)

        phy_spans: list[tuple[int, int]] = []
        for scaff_idx, windows in self.phymap_windows.items():
            rows = rows_by_scaffold.get(scaff_idx)
            if rows is None:
                continue
            scaff_rows = np.asarray(rows, dtype=np.int64)
            pos0 = scaff_rows[:, 3]
            pos1 = scaff_rows[:, 4]
            for start, end in windows:
                overlap = (pos1 >= start) & (pos0 <= end)
                overlap_rows = np.flatnonzero(overlap)
                if overlap_rows.size == 0:
                    continue

                first = scaff_rows[overlap_rows[0]]
                last = scaff_rows[overlap_rows[-1]]
                wmin_offset = max(0, start - int(first[3]))
                wmin = int(first[1]) + wmin_offset
                wmax_offset = max(0, int(last[4]) - end)
                wmax = int(last[2]) - wmax_offset
                if wmax > wmin:
                    phy_spans.append((wmin, wmax))

        if not phy_spans:
            raise IPyradError(
                "Selected windows contain zero data in the assembly. Try larger/different windows."
            )

        self._selected_phy_spans = tuple(phy_spans)
        self._chunk_spans = plan_sequence_chunk_spans(self._selected_phy_spans)
        logger.debug(
            "planned {} selected spans into {} HDF5 chunks",
            len(self._selected_phy_spans),
            len(self._chunk_spans),
        )
        return self._selected_phy_spans

    def _filter_block_sites(self, block: np.ndarray) -> np.ndarray:
        """Apply the minmap filter to one loaded sequence block."""
        return filter_block_by_minmap(block, self._imap_row_indices, self.minmap)

    def _summarize_filtered_selection(self) -> dict:
        """Return filtering decisions and stats using bounded-memory chunk scans."""
        self._get_selected_phy_spans()
        missing_counts = np.zeros(len(self.snames), dtype=np.int64)
        nsites_before = 0
        nvariants_before = 0
        nsites_after = 0
        nvariants_after = 0

        with h5py.File(self.data, "r") as io5:
            phy = io5["phy"]
            for spans in self._chunk_spans:
                block = load_sequence_chunk_from_phy(phy, self.sidxs, spans)
                nsites_before += int(block.shape[1])
                nvariants_before += count_snps(block)
                filtered = self._filter_block_sites(block)
                if not filtered.size:
                    continue
                nsites_after += int(filtered.shape[1])
                nvariants_after += count_snps(filtered)
                missing_counts += np.sum(filtered == MISSING_BASE, axis=1).astype(np.int64, copy=False)

        if nsites_after == 0:
            raise IPyradError("Selected windows contain zero data after filtering for coverage.")

        row_missing = missing_counts / nsites_after
        keep_mask = row_missing <= self.max_sample_missing
        keep_indices = np.flatnonzero(keep_mask).astype(np.int64, copy=False)
        fnames = [self.snames[idx] for idx in keep_indices]
        if not fnames:
            raise IPyradError("No samples passed max_sample_missing filter.")

        return {
            "fnames": fnames,
            "keep_indices": keep_indices,
            "row_missing": row_missing,
            "samples_selected_initial": list(self.snames),
            "samples_dropped_by_max_missing": [
                self.snames[idx] for idx in np.flatnonzero(~keep_mask)
            ],
            "nsites_before": nsites_before,
            "nvariants_before": nvariants_before,
            "nsites_after": nsites_after,
            "nvariants_after": nvariants_after,
        }

    def _iter_filtered_blocks(self, keep_indices: np.ndarray):
        """Yield filtered blocks for the retained samples."""
        self._get_selected_phy_spans()
        with h5py.File(self.data, "r") as io5:
            phy = io5["phy"]
            for spans in self._chunk_spans:
                block = load_sequence_chunk_from_phy(phy, self.sidxs, spans)
                filtered = self._filter_block_sites(block)
                if filtered.size:
                    yield filtered[keep_indices, :]

    def _collect_filtered_buffers(self, keep_indices: np.ndarray) -> list[bytearray]:
        """Return per-sample sequence buffers without materializing the full matrix."""
        buffers = [bytearray() for _ in range(len(keep_indices))]
        for block in self._iter_filtered_blocks(keep_indices):
            for idx in range(block.shape[0]):
                buffers[idx].extend(block[idx].tobytes())
        return buffers

    def _extract_filtered_alignment(self) -> tuple[dict, list[str], np.ndarray]:
        """Return summary stats and the fully filtered alignment matrix."""
        summary = self._summarize_filtered_selection()
        blocks = list(self._iter_filtered_blocks(summary["keep_indices"]))
        if len(blocks) == 1:
            fseqarr = blocks[0]
        else:
            fseqarr = np.concatenate(blocks, axis=1)
        self.seqarr = fseqarr
        return summary, summary["fnames"], fseqarr

    def _get_filtered_alignment_data(
        self,
        outfile: Path | str,
    ) -> tuple[list[str], np.ndarray, dict]:
        """Return filtered names, sequences, and stats without writing output."""
        self._get_phymap_windows()
        self._get_phymap()
        summary, fnames, fseqarr = self._extract_filtered_alignment()
        return fnames, fseqarr, self._build_stats_dict(summary, outfile)

    def _get_output_path(self, suffix: str) -> Path:
        return self.outdir / f"{self.name}.{suffix}"

    def _get_stats_path(self) -> Path:
        return self.outdir / f"{self.name}.stats.txt"

    def _prepare_output_paths(
        self,
        suffix: str | None,
        *,
        write_stats: bool,
        return_locus: bool = False,
    ) -> Path | str:
        paths: List[Path] = []

        if suffix and not self.stdout and not return_locus:
            paths.append(self._get_output_path(suffix))
        if write_stats:
            paths.append(self._get_stats_path())

        if paths:
            self.outdir.mkdir(exist_ok=True)
            if not self.force:
                existing = next((path for path in paths if path.exists()), None)
                if existing is not None:
                    raise IPyradError(
                        f"Output file already exists: {existing}. Use --force to overwrite."
                    )

        if self.stdout:
            return "STDOUT"
        if return_locus:
            return "RETURN"
        if suffix is None:
            raise IPyradError("Internal error: missing output suffix.")
        return self._get_output_path(suffix)

    def _build_stats_dict(self, summary: dict, outfile):
        """Build stats for the extracted windows without writing them."""
        return {
            "nsamples_before_filtering": len(self.snames),
            "nsites_in_windows_before_filtering": summary["nsites_before"],
            "nvariants_in_windows_before_filtering": summary["nvariants_before"],
            "nsamples_after_filtering": len(summary["fnames"]),
            "nsites_in_windows_after_filtering": summary["nsites_after"],
            "nvariants_in_windows_after_filtering": summary["nvariants_after"],
            "infile": self.data,
            "outfile": outfile,
            "windows": self.selected_windows,
            "imap": self.imap,
            "min_sample_coverage_filter": self.minmap,
            "max_sample_missing_filter": self.max_sample_missing,
        }

    def _format_windows_preview(self, limit: int = 10) -> str:
        """Return one readable preview of the selected windows."""
        if not self.selected_windows:
            return "(none)"
        if len(self.selected_windows) <= limit:
            return ", ".join(self.selected_windows)
        preview = ", ".join(self.selected_windows[:limit])
        return f"{preview}, ... ({len(self.selected_windows)} total)"

    def _format_minmap(self) -> str:
        """Return one readable `group=value` summary of the minmap filter."""
        rows = []
        for pop in sorted(self.minmap):
            value = self.minmap[pop]
            if isinstance(value, float) and not float(value).is_integer():
                formatted = _format_fraction(value)
            else:
                formatted = _format_count(int(value))
            rows.append(f"{pop}={formatted}")
        return ", ".join(rows) if rows else "(none)"

    def _build_sample_population_lookup(self) -> dict[str, str]:
        """Return one inverse sample-to-population lookup."""
        lookup = {}
        for pop, names in self.imap.items():
            for name in names:
                lookup[name] = pop
        return lookup

    def _write_stats_report(
        self,
        summary: dict,
        stats_dict: dict,
        stats_file: Path | None = None,
    ) -> None:
        """Write the human-readable wex stats report to disk."""
        stats_path = self._get_stats_path() if stats_file is None else Path(stats_file)
        self.outdir.mkdir(exist_ok=True)
        population_lookup = self._build_sample_population_lookup()
        dropped_by_missing = set(summary["samples_dropped_by_max_missing"])
        sample_rows: list[list[str]] = []
        for name, missing in sorted(
            zip(summary["samples_selected_initial"], summary["row_missing"], strict=True),
            key=lambda item: item[0],
        ):
            sample_rows.append([
                name,
                population_lookup.get(name, "unassigned"),
                _format_percent(missing),
                "yes" if name in dropped_by_missing else "no",
            ])

        extract_rows = [
            ("infile", str(stats_dict["infile"])),
            ("outfile", str(stats_dict["outfile"])),
            ("out_format", self.out_format),
            ("windows_selected", _format_count(len(self.selected_windows))),
            ("selected_windows_preview", self._format_windows_preview()),
        ]
        filtering_rows = [
            ("populations", ", ".join(sorted(self.imap)) if self.imap else "(none)"),
            ("min_sample_coverage_filter", self._format_minmap()),
            ("max_sample_missing", _format_fraction(self.max_sample_missing)),
            (
                "samples_selected_initial",
                _format_count(len(summary["samples_selected_initial"])),
            ),
            (
                "samples_dropped_by_max_missing",
                _format_count(len(summary["samples_dropped_by_max_missing"])),
            ),
            ("samples_final", _format_count(len(summary["fnames"]))),
        ]
        alignment_rows = [
            (
                "nsamples_before_filtering",
                _format_count(stats_dict["nsamples_before_filtering"]),
            ),
            (
                "nsites_in_windows_before_filtering",
                _format_count(stats_dict["nsites_in_windows_before_filtering"]),
            ),
            (
                "nvariants_in_windows_before_filtering",
                _format_count(stats_dict["nvariants_in_windows_before_filtering"]),
            ),
            (
                "nsamples_after_filtering",
                _format_count(stats_dict["nsamples_after_filtering"]),
            ),
            (
                "nsites_in_windows_after_filtering",
                _format_count(stats_dict["nsites_in_windows_after_filtering"]),
            ),
            (
                "nvariants_in_windows_after_filtering",
                _format_count(stats_dict["nvariants_in_windows_after_filtering"]),
            ),
        ]

        lines: list[str] = []
        _append_key_value_section(lines, "Extract Summary", extract_rows)
        _append_key_value_section(lines, "Filtering Summary", filtering_rows)
        _append_key_value_section(lines, "Alignment Summary", alignment_rows)
        _append_table_section(
            lines,
            "Sample Summary",
            ["sample", "population", "percent_missing", "dropped_by_max_missing"],
            sample_rows,
        )
        report_text = "\n".join(lines).rstrip() + "\n"
        if self.logged_command:
            report_text = f"CMD: {self.logged_command}\n\n{report_text}"
        with open(stats_path, "w", encoding="utf-8") as out:
            out.write(report_text)
        logger.info(f"wrote stats/log to: {stats_path}")

    def _write_to_phy(self, 
                      write_stats: bool = True,
                      prefix: str = None,
                      bpp_format: bool = False,
                      return_locus: bool = False,
                      return_alignment: bool = False,
                      return_stats: bool = False):
        """Writes the .seqarr matrix as a string to .outfile."""
        outfile = self._prepare_output_paths(
            "phy",
            write_stats=write_stats,
            return_locus=return_locus,
        )
        prefix = prefix if prefix else ""
        self._get_phymap_windows()
        self._get_phymap()

        if return_alignment or return_locus:
            summary, fnames, fseqarr = self._extract_filtered_alignment()

            # get padded names
            longname = max(len(i) for i in fnames)
            pnames = [i.ljust(longname + 5) for i in fnames]

            phy = []
            for idx, _ in enumerate(fnames):
                seq = fseqarr[idx].tobytes().decode("utf-8")
                phy.append(f"{prefix}{pnames[idx]} {seq}")

            ntaxa = len(fnames)
            nsites = fseqarr.shape[1]
            bpp_sep = "\n" if bpp_format else ""
            phy_text = "\n".join(phy)
            alignment = f"{ntaxa} {nsites}\n{bpp_sep}{phy_text}\n"
            stats_dict = self._build_stats_dict(summary, outfile)

            if return_alignment:
                if return_stats:
                    return alignment, stats_dict
                return alignment

            if write_stats:
                self._write_stats_report(summary, stats_dict)
            if return_locus and return_stats:
                return alignment, stats_dict
            if return_locus:
                return alignment
            if return_stats:
                return stats_dict
            return None

        summary = self._summarize_filtered_selection()
        fnames = summary["fnames"]
        seq_buffers = self._collect_filtered_buffers(summary["keep_indices"])
        stats_dict = self._build_stats_dict(summary, outfile)

        # get padded names
        longname = max(len(i) for i in fnames)
        pnames = [i.ljust(longname + 5) for i in fnames]
        ntaxa = len(fnames)
        nsites = summary["nsites_after"]

        # write to stdout
        if self.stdout:
            logger.debug("wrote alignment to stdout")
            out = sys.stdout
            out.write(f"{ntaxa} {nsites}\n")
            if bpp_format:
                out.write("\n")
            for idx, seq in enumerate(seq_buffers):
                out.write(f"{prefix}{pnames[idx]} {seq.decode('utf-8')}\n")
        else:
            with open(outfile, 'w', encoding="utf-8") as out:
                out.write(f"{ntaxa} {nsites}\n")
                if bpp_format:
                    out.write("\n")
                for idx, seq in enumerate(seq_buffers):
                    out.write(f"{prefix}{pnames[idx]} {seq.decode('utf-8')}\n")
            logger.info(f"wrote alignment ({ntaxa}, {nsites}) to: {outfile}")

        if write_stats:
            self._write_stats_report(summary, stats_dict)

        if return_stats:
            return stats_dict

    def _write_to_nex(
        self,
        write_stats: bool = True,
        return_alignment: bool = False,
        return_stats: bool = False,
    ):
        """Writes concatenated alignment to nex format..."""
        outfile = self._prepare_output_paths("nex", write_stats=write_stats)
        self._get_phymap_windows()
        self._get_phymap()
        if return_alignment:
            summary, fnames, fseqarr = self._extract_filtered_alignment()

            longname = max(len(i) for i in fnames)
            pnames = [i.ljust(longname + 5) for i in fnames]
            ntaxa = len(fnames)
            nsites = fseqarr.shape[1]

            lines = [NEXHEADER.format(ntaxa, nsites)]
            for block in range(0, fseqarr.shape[1], 100):
                stop = min(block + 100, fseqarr.shape[1])
                for idx, name in enumerate(pnames):
                    seq = fseqarr[idx, block:stop].tobytes().decode()
                    lines.append(f"  {name}{seq}\n")
                lines.append("\n")
            lines.append("  ;\nend;")
            alignment = "".join(lines)
            stats_dict = self._build_stats_dict(summary, outfile)

            if return_alignment:
                if return_stats:
                    return alignment, stats_dict
                return alignment

        summary = self._summarize_filtered_selection()
        fnames = summary["fnames"]
        seq_buffers = self._collect_filtered_buffers(summary["keep_indices"])
        longname = max(len(i) for i in fnames)
        pnames = [i.ljust(longname + 5) for i in fnames]
        ntaxa = len(fnames)
        nsites = summary["nsites_after"]
        stats_dict = self._build_stats_dict(summary, outfile)

        def _write_nexus(out) -> None:
            out.write(NEXHEADER.format(ntaxa, nsites))
            for block in range(0, nsites, 100):
                stop = min(block + 100, nsites)
                for idx, name in enumerate(pnames):
                    seq = bytes(seq_buffers[idx][block:stop]).decode("utf-8")
                    out.write(f"  {name}{seq}\n")
                out.write("\n")
            out.write("  ;\nend;")

        # write to stdout
        if self.stdout:
            logger.debug("wrote alignment to stdout")
            _write_nexus(sys.stdout)
        else:
            with open(outfile, 'w', encoding="utf-8") as out:
                _write_nexus(out)
            logger.info(f"wrote alignment ({ntaxa}, {nsites}) to: {outfile}")
        if write_stats:
            self._write_stats_report(summary, stats_dict)
        if return_stats:
            return stats_dict

    def _write_to_fa(self, write_stats: bool = True, return_stats: bool = False):
        """Write the extracted alignment as FASTA."""
        outfile = self._prepare_output_paths("fa", write_stats=write_stats)
        self._get_phymap_windows()
        self._get_phymap()
        summary = self._summarize_filtered_selection()
        fnames = summary["fnames"]
        seq_buffers = self._collect_filtered_buffers(summary["keep_indices"])
        stats_dict = self._build_stats_dict(summary, outfile)
        if self.stdout:
            logger.debug("wrote alignment to stdout")
            for idx, name in enumerate(fnames):
                sys.stdout.write(f">{name}\n{seq_buffers[idx].decode('utf-8')}\n")
        else:
            with open(outfile, "w", encoding="utf-8") as out:
                for idx, name in enumerate(fnames):
                    out.write(f">{name}\n{seq_buffers[idx].decode('utf-8')}\n")
            logger.info(
                f"wrote alignment ({len(fnames)}, {summary['nsites_after']}) to: {outfile}"
            )

        if write_stats:
            self._write_stats_report(summary, stats_dict)
        if return_stats:
            return stats_dict

    def _write_stats(self, fnames, fseqarr, outfile):
        """Write stats for the extracted windows."""
        summary = {
            "fnames": fnames,
            "row_missing": np.sum(fseqarr == MISSING_BASE, axis=1) / fseqarr.shape[1],
            "samples_selected_initial": list(fnames),
            "samples_dropped_by_max_missing": [],
            "nsites_before": self.seqarr.shape[1],
            "nvariants_before": count_snps(self.seqarr),
            "nsites_after": fseqarr.shape[1],
            "nvariants_after": count_snps(fseqarr),
        }
        self._write_stats_report(summary, self._build_stats_dict(summary, outfile))


def count_snps(arr):
    """Count variants to report in the stats for an alignment."""
    m = np.ma.masked_equal(arr, 78)
    multi_cols = (np.ma.ptp(m, axis=0) > 0).filled(False)
    return int(np.sum(multi_cols))


def run_window_extracter(**kwargs):
    """command line wrapper for window-extracter.

    Parameters:
    -----------
    data: Path | str
        A 'seqs.hdf5' database file from ipyrad2.
    name: str
        Prefix name used for outfiles. If None it is automatically set.
    outdir: Path | str
        Dir for output files. Created if it doesn't exist.
    out_format: str
        Format to write the alignments phy (default), nex, or fa.
    windows: str | List[str]:
        Subsample scaffold(s) by index number. If unsure, leave this
        empty when loading a file and then check the .scaffold_table
        to view the indices of scaffolds. Scaffolds are ordered by
        their order in the reference genome file.
    min_sample_coverage: int | float:
        Min number of individuals that must have data at a site
        for it to be included in the alignment (def=4).
    max_sample_missing: float
        Max proportion of sites that can be missing (N) in a sample.
        (def=1.0)
    exclude: List[str]
        A list of sample names to exclude from the data set. Samples
        can also be excluded by using an imap dictionary and not
        including them.
    imap: Dict
        A dictionary mapping group names (keys) to lists of sample
        names (values) to be included in the analysis. This can be
        used for 3 things: (1) to select samples to extract data for;
        (2) to filter based on sample coverage in groups (minmap);
        or (3) to use consensus_reduce=True to reduce the dataset to a
        consensus sequence for each group.
    minmap: Dict
        A dictionary mapping group names (keys) to integers or floats
        to act as a filter requiring that at least N (or N%) of samples
        in this group have data for a locus to be retained in the
        dataset. When using consensus_reduce=True the minmap applies to
        the reduced data set, i.e., it applies to the groups (keys) so
        that all values must be <= 1.
    stdout: bool
        ...
    force: bool
        ...
    """
    request_table = kwargs.pop("print_scaffold_table")

    tool = WindowExtracter(**kwargs)

    if request_table:
        try:
            tool.scaffold_table.to_csv(sys.stdout, sep="\t")
            sys.stdout.flush()
        except BrokenPipeError:
            sys.stdout = open(os.devnull, "w", encoding="utf-8")
        sys.exit(0)

    if tool.out_format == "phy":
        tool._write_to_phy()
    elif tool.out_format == "nex":
        tool._write_to_nex()
    elif tool.out_format == "fa":
        tool._write_to_fa()
    else:
        raise IPyradError(f"Unrecognized output format: {tool.out_format}")

    sys.exit(0)


if __name__ == "__main__":

    h5 = Path("/tmp/OUT_klmnop/assembly.seqs.hdf5")
    h5 = Path("/home/deren/Documents/ipyrad-tests/OUT/assembly.seqs.hdf5")
    h5 = Path("/home/deren/Documents/ipyrad-tests/Ped2_OUT/assembly.hdf5")
    assert h5.exists(), "h5 doesn't exist"

    with h5py.File(h5, 'r') as io5:
        print(len(io5.attrs['names']))
        print(io5["phymap"][:])
        print(io5["phy"].shape)

        # help(io5.create_dataset)


    # tool = WindowExtracter(
    #     data=h5,
    #     name='test',
    #     outdir=Path("/tmp/WEX"),
    #     windows=r"MT",
    #     min_sample_coverage=4,
    #     max_sample_missing=1.0,
    #     exclude=[],
    #     imap=None,
    #     minmap=None,
    #     stdout=True,
    #     force=True,
    # )
    # tool._write_to_phy()



    # print(tool.scaffold_table)
    # arr, stats = tool.run(return_data=True)
    # print(stats.T)
    # print(arr)
