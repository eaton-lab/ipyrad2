"""Resolve sequence-HDF5 windows and map them onto delimited loci."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

from loguru import logger
import pandas as pd

from ...utils.exceptions import IPyradError


@dataclass(frozen=True)
class SelectedWindow:
    """One resolved 1-based inclusive scaffold window."""

    scaffold_index: int
    scaffold: str
    start: int
    end: int
    explicit_coordinates: bool

    @property
    def label(self) -> str:
        return f"{self.scaffold}:{self.start}-{self.end}"


@dataclass(frozen=True)
class LocusIntersection:
    """One complete or window-clipped view into a phymap locus."""

    phy0: int
    phy1: int
    pos0: int
    pos1: int
    clipped: bool


def _overlaps(start1: int, end1: int, start2: int, end2: int) -> bool:
    """Return whether two 1-based inclusive coordinate intervals overlap."""
    return not (end1 < start2 or end2 < start1)


def _append_window(
    selected: list[SelectedWindow],
    *,
    scaffold_index: int,
    scaffold: str,
    start: int,
    end: int,
    explicit_coordinates: bool,
    source: str,
) -> None:
    """Validate and append one selected window."""
    if start < 1 or end < start:
        raise IPyradError(
            f"Malformed window '{source}'. Windows must use valid positive coordinates."
        )
    selected.append(
        SelectedWindow(
            scaffold_index=scaffold_index,
            scaffold=scaffold,
            start=start,
            end=end,
            explicit_coordinates=explicit_coordinates,
        )
    )


def _validate_nonoverlapping_windows(selected: list[SelectedWindow]) -> None:
    """Reject overlapping windows in one sorted pass."""
    ordered = sorted(
        selected,
        key=lambda window: (window.scaffold_index, window.start, window.end),
    )
    for previous, current in zip(ordered, ordered[1:]):
        if (
            previous.scaffold_index == current.scaffold_index
            and current.start <= previous.end
        ):
            raise IPyradError(
                f"windows cannot overlap. {current.label} overlaps {previous.label}"
            )


def _read_bed_windows(path: Path) -> list[tuple[str, int, int]]:
    """Return BED intervals converted to 1-based inclusive coordinates."""
    windows = []
    with path.open("r", encoding="utf-8") as infile:
        for lineno, line in enumerate(infile, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip().split()
            if len(fields) < 3:
                raise IPyradError(f"Malformed BED window at line {lineno}: {path}")
            scaffold, start, end = fields[:3]
            try:
                start0 = int(start)
                end0 = int(end)
            except ValueError as exc:
                raise IPyradError(
                    f"Malformed BED window at line {lineno}: {path}"
                ) from exc
            if start0 < 0 or end0 <= start0:
                raise IPyradError(f"Malformed BED window at line {lineno}: {path}")
            windows.append((scaffold, start0 + 1, end0))
    return windows


def resolve_sequence_windows(
    scaffold_table: pd.DataFrame,
    windows: str | list[str] | None,
) -> list[SelectedWindow]:
    """Resolve regex, region, or BED selectors in scaffold-table order."""
    requested = (
        []
        if windows is None
        else [windows]
        if isinstance(windows, str)
        else list(windows)
    )
    names = scaffold_table["scaffold_name"].astype(str).tolist()
    lengths = scaffold_table["scaffold_length"].astype(int).tolist()
    name_to_index = {name: idx for idx, name in enumerate(names)}
    selected: list[SelectedWindow] = []

    if not requested:
        logger.info(
            "No windows specified; selecting the full length of all scaffolds. "
            "Use -w to subset scaffold windows and -P to view scaffold names."
        )
        requested = [r".*"]

    if len(requested) == 1:
        bed_path = Path(requested[0])
        if bed_path.exists():
            logger.info("Loading windows from bed file: '{}'", bed_path)
            for scaffold, start, end in _read_bed_windows(bed_path):
                if scaffold not in name_to_index:
                    raise IPyradError(
                        f"'{scaffold}' from {bed_path} does not match any scaffold names. "
                        "Check with '-P'."
                    )
                _append_window(
                    selected,
                    scaffold_index=name_to_index[scaffold],
                    scaffold=scaffold,
                    start=start,
                    end=end,
                    explicit_coordinates=True,
                    source=f"{scaffold}:{start}-{end}",
                )
            _validate_nonoverlapping_windows(selected)
            _log_selection(selected)
            return sorted(selected, key=lambda item: (item.scaffold_index, item.start))

    for selector in requested:
        if ":" in selector:
            scaffold_pattern, region = selector.split(":", 1)
            matches = [name for name in names if re.fullmatch(scaffold_pattern, name)]
            if not matches:
                raise IPyradError(
                    f"No scaffold names match '{scaffold_pattern}'. Use -P to view scaffold names."
                )
            if len(matches) > 1:
                raise IPyradError(
                    "Cannot use regex with ':'. List windows separately: "
                    "-w Chr1:1-1000 Chr2:1-1000"
                )
            if region.count("-") != 1:
                raise IPyradError(
                    f"malformatted window '{selector}'. Must be "
                    "{scaff} or {scaff}:{start}-{end}"
                )
            try:
                start, end = (int(value) for value in region.split("-"))
            except ValueError as exc:
                raise IPyradError(
                    f"malformatted window '{selector}'. Must be "
                    "{scaff} or {scaff}:{start}-{end}"
                ) from exc
            scaffold = matches[0]
            _append_window(
                selected,
                scaffold_index=name_to_index[scaffold],
                scaffold=scaffold,
                start=start,
                end=end,
                explicit_coordinates=True,
                source=selector,
            )
            continue

        matches = [name for name in names if re.fullmatch(selector, name)]
        if not matches:
            raise IPyradError(
                f"'{selector}' does not match any scaffold names. Check with '-P'."
            )
        for scaffold in matches:
            idx = name_to_index[scaffold]
            _append_window(
                selected,
                scaffold_index=idx,
                scaffold=scaffold,
                start=1,
                end=lengths[idx],
                explicit_coordinates=False,
                source=selector,
            )

    _validate_nonoverlapping_windows(selected)
    selected.sort(key=lambda item: (item.scaffold_index, item.start))
    _log_selection(selected)
    return selected


def _log_selection(selected: list[SelectedWindow]) -> None:
    """Log a concise resolved-window summary."""
    scaffolds = {window.scaffold_index for window in selected}
    logger.debug("windows: {}", [window.label for window in selected])
    logger.info(
        "selected {} window{} from {} scaffold{}",
        len(selected),
        "s" if len(selected) != 1 else "",
        len(scaffolds),
        "s" if len(scaffolds) != 1 else "",
    )


def intersect_phymap_locus(
    row,
    window: SelectedWindow,
    *,
    clip: bool,
) -> LocusIntersection | None:
    """Return a complete or clipped locus intersection for one phymap row."""
    locus_start = int(row[3])
    locus_end = int(row[4])
    if not _overlaps(locus_start, locus_end, window.start, window.end):
        return None
    if not clip:
        return LocusIntersection(
            phy0=int(row[1]),
            phy1=int(row[2]),
            pos0=locus_start,
            pos1=locus_end,
            clipped=False,
        )

    start = max(locus_start, window.start)
    end = min(locus_end, window.end)
    phy0 = int(row[1]) + (start - locus_start)
    phy1 = phy0 + (end - start + 1)
    return LocusIntersection(
        phy0=phy0,
        phy1=phy1,
        pos0=start,
        pos1=end,
        clipped=start != locus_start or end != locus_end,
    )
