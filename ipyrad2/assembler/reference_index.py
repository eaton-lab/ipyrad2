"""Reuse FASTA indexes without rewriting shared reference inputs."""

from __future__ import annotations

import filecmp
import os
from pathlib import Path
import shlex
import sys
import tempfile

from loguru import logger

from ..utils.exceptions import IPyradError
from ..utils.parallel import run_pipeline

BIN_SAM = str(Path(sys.prefix) / "bin" / "samtools")


def _signature(path: Path) -> tuple[int, int, int, int]:
    stat = path.stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns


def _check_index(fai: Path, reference_size: int) -> None:
    """Check index structure and uncompressed-file bounds, without a genome scan."""
    seen: set[str] = set()
    previous_offset = -1
    with fai.open(encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 5 or not fields[0] or any(c.isspace() for c in fields[0]):
                raise ValueError(f"malformed FASTA index row {lineno}")
            name = fields[0]
            length, offset, line_bases, line_width = map(int, fields[1:])
            if name in seen:
                raise ValueError(f"duplicate contig {name!r}")
            seen.add(name)
            if length < 0 or not 0 < offset <= reference_size or offset <= previous_offset:
                raise ValueError(f"invalid length or offset for {name!r}")
            previous_offset = offset
            if length == 0:
                if line_bases != 0 or line_width != 0:
                    raise ValueError(f"invalid empty-contig line layout for {name!r}")
                continue
            if not (0 < line_bases <= length and line_bases <= line_width <= line_bases + 2):
                raise ValueError(f"invalid line layout for {name!r}")
            last_base = offset + ((length - 1) // line_bases) * line_width + (length - 1) % line_bases
            if last_base >= reference_size:
                raise ValueError(f"indexed sequence exceeds the FASTA size for {name!r}")
    if not seen:
        raise ValueError("empty FASTA index")


def ensure_reference_index(reference: Path, *, workspace: Path | None = None) -> Path:
    """Reuse a current .fai; verify an older one privately; create only if absent.

    Current-index checks use metadata, row structure and byte bounds, not a full
    sequence checksum. An older timestamp alone does not invalidate an index:
    compare it to a temporary samtools rebuild before deciding. Inconsistent
    existing indexes are never replaced automatically. Assembly's FASTA reader
    requires an uncompressed reference.
    """
    reference = Path(reference)
    fai = Path(str(reference) + ".fai")
    remedy = (
        "Check that the FASTA is the reference used for mapping, then explicitly "
        f"regenerate its index with {shlex.join(['samtools', 'faidx', str(reference)])}. "
        "The existing index has not been overwritten."
    )
    try:
        before = _signature(reference)
        if not reference.is_file():
            raise IPyradError(f"Reference FASTA is not a file: {reference}")
        with reference.open('rb') as handle:
            if handle.read(2) == b'\x1f\x8b':
                raise IPyradError("assemble requires an uncompressed reference FASTA.")
        exists = fai.exists() or fai.is_symlink()
        index_before = None
        if exists:
            if not fai.is_file():
                raise IPyradError(f"FASTA index is not a readable file: {fai}. {remedy}")
            index_before = _signature(fai)
            try:
                _check_index(fai, before[2])
            except (ValueError, UnicodeError) as exc:
                raise IPyradError(f"Invalid FASTA index {fai}: {exc}. {remedy}") from exc
            if index_before != _signature(fai) or before != _signature(reference):
                raise IPyradError("Reference FASTA or index changed during index validation.")
            if index_before[3] >= before[3]:
                logger.debug("reusing reference FASTA index {}", fai)
                return fai

        # Missing indexes are staged on the destination filesystem for atomic
        # no-clobber publication. Stale indexes can be checked in the assembly
        # workspace even if the shared reference directory is read-only.
        directory = reference.parent if not exists else workspace
        with tempfile.TemporaryDirectory(prefix=".ipyrad-faidx-", dir=directory) as temp:
            candidate = Path(temp) / "reference.fai"
            logger.info("{} reference FASTA index {}", "checking older" if exists else "creating missing", fai)
            run_pipeline([[BIN_SAM, "faidx", "--fai-idx", str(candidate), str(reference)]])
            if before != _signature(reference):
                raise IPyradError("Reference FASTA changed while generating its index.")
            try:
                _check_index(candidate, before[2])
            except (ValueError, UnicodeError) as exc:
                raise IPyradError(f"samtools generated an invalid FASTA index: {exc}") from exc
            if exists:
                if index_before != _signature(fai):
                    raise IPyradError("Reference index changed during index validation.")
                if not filecmp.cmp(candidate, fai, shallow=False):
                    raise IPyradError(f"FASTA index does not match the current reference: {fai}. {remedy}")
                logger.debug("older reference index matches the FASTA; preserving {}", fai)
            else:
                try:
                    os.link(candidate, fai)
                except FileExistsError:
                    # Another run may have just published the same index.
                    # Never replace its file, including a conflicting result.
                    if not fai.is_file() or not filecmp.cmp(candidate, fai, shallow=False):
                        raise IPyradError(f"A conflicting FASTA index appeared during creation: {fai}. {remedy}") from None
        return fai
    except OSError as exc:
        raise IPyradError(f"Cannot prepare reference FASTA index {fai}: {exc}") from exc
