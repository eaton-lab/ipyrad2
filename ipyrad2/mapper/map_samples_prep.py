#!/usr/bin/env python

"""Sample parsing and materialization helpers for ipyrad2 map."""

from __future__ import annotations

from collections import defaultdict
import gzip
import shutil
from pathlib import Path
from typing import Dict, Tuple

from loguru import logger

from ..utils.exceptions import IPyradError
from ..utils.names import get_name_to_fastq_dict
from ..utils.pops import expand_imap_patterns


def _require(condition: bool, message: str) -> None:
    """Raise a user-facing sample-prep error when a precondition is not met."""
    if not condition:
        raise IPyradError(message)


def _open_mapper_fastq(path: Path):
    """Open one mapper FASTQ input in binary mode."""
    if path.suffix == ".bz2":
        raise IPyradError(
            f"ipyrad2 map supports only plain FASTQ or .gz-compressed FASTQ inputs: {path}"
        )
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        return opener(path, "rb")
    except OSError as err:
        raise IPyradError(f"Failed to read FASTQ input: {path}") from err


def _concat_fastqs(inputs: list[Path], outfile: Path) -> None:
    """Concatenate plain or gzipped FASTQs into one gzipped output file."""
    outfile.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(outfile, "wb") as out:
        for path in inputs:
            with _open_mapper_fastq(path) as infile:
                shutil.copyfileobj(infile, out, length=1024 * 1024)


def validate_fastq_inputs(
    fastq_dict: Dict[str, Tuple[Path, Path | None]],
) -> None:
    """Validate mapper FASTQ inputs after parsing or imap materialization."""
    for sname, fastq_tuple in fastq_dict.items():
        for path in fastq_tuple:
            if path is None:
                continue
            _require(path.exists(), f"FASTQ path does not exist for sample '{sname}': {path}")
            _require(
                path.suffix != ".bz2",
                f"ipyrad2 map supports only plain FASTQ or .gz-compressed FASTQ inputs: {path}",
            )


def detect_is_paired(
    fastq_dict: Dict[str, Tuple[Path, Path | None]],
) -> bool:
    """Validate pairedness consistency across parsed mapper samples."""
    paired_states = {fastq_tuple[1] is not None for fastq_tuple in fastq_dict.values()}
    if not paired_states:
        return False
    if len(paired_states) != 1:
        raise IPyradError("some but not all files have R1 and R2 pairs. Check inputs.")
    return next(iter(paired_states))


def apply_imap_to_samples(
    imap: Path,
    tmpdir: Path,
    fastq_dict: Dict[str, Tuple[Path, Path | None]],
) -> Dict[str, Tuple[Path, Path | None]]:
    """Return a new fastq_dict after applying subset, rename, and merge rules."""
    source_fastqs = dict(fastq_dict)
    original_total = len(source_fastqs)
    pname_to_path_tuples: dict[str, list[Tuple[Path, Path | None]]] = defaultdict(list)
    pname_to_snamelist: dict[str, list[str]] = defaultdict(list)
    raw_imap: dict[str, list[str]] = defaultdict(list)
    identity_prefix = "__imap_identity__"

    with open(imap, "r", encoding="utf-8") as infile:
        for lineno, line in enumerate(infile, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            sname, *data = line.split()
            pname = data[0] if data else None
            if pname is None:
                raw_imap[f"{identity_prefix}{lineno}"].append(sname)
            else:
                raw_imap[pname].append(sname)

    expanded_imap, warn_list = expand_imap_patterns(
        raw_imap,
        source_fastqs,
        mapping_name="imap file",
        available_name="parsed sample names",
        strict_unmatched=False,
    )
    for pname, snames in expanded_imap.items():
        if pname.startswith(identity_prefix):
            for sname in snames:
                pname_to_snamelist[sname].append(sname)
                pname_to_path_tuples[sname].append(source_fastqs[sname])
            continue
        for sname in snames:
            pname_to_snamelist[pname].append(sname)
            pname_to_path_tuples[pname].append(source_fastqs[sname])

    if warn_list:
        logger.warning(
            "One or more names or glob patterns in imap file did not match sample names and will be skipped: {}",
            " ".join(warn_list),
        )
    if not pname_to_snamelist:
        raise IPyradError("No samples in imap file match parsed sample names. Revise imap file or name parsing args.")

    selected_total = sum(len(i) for i in pname_to_snamelist.values())
    logger.info(
        "subselecting, renaming, or merging {}/{} samples into {} samples according to: {}",
        selected_total,
        original_total,
        len(pname_to_snamelist),
        imap.name,
    )

    result: Dict[str, Tuple[Path, Path | None]] = {}
    maxlen = max(len(i) for i in pname_to_snamelist)
    for pname, path_tuples in pname_to_path_tuples.items():
        snames = pname_to_snamelist[pname]
        logger.info(f"{pname}{' ' * (maxlen - len(pname))} <- {' + '.join(snames)}")

        if len(path_tuples) == 1:
            result[pname] = path_tuples[0]
            continue

        paired_states = {fastq_tuple[1] is not None for fastq_tuple in path_tuples}
        if len(paired_states) != 1:
            raise IPyradError(
                f"Cannot merge technical replicates with mixed SE and PE inputs for sample '{pname}'."
            )

        out1 = tmpdir / f"{pname}.tmp.R1.fastq.gz"
        _concat_fastqs([fastq_tuple[0] for fastq_tuple in path_tuples], out1)
        if path_tuples[0][1] is not None:
            out2 = tmpdir / f"{pname}.tmp.R2.fastq.gz"
            _concat_fastqs(
                [fastq_tuple[1] for fastq_tuple in path_tuples if fastq_tuple[1] is not None],
                out2,
            )
            result[pname] = (out1, out2)
        else:
            result[pname] = (out1, None)
    return result


def unmate_paired_samples(
    fastq_dict: Dict[str, Tuple[Path, Path | None]],
    tmpdir: Path,
) -> Dict[str, Tuple[Path, Path | None]]:
    """Return one effective SE FASTQ per paired sample by concatenating R1 then R2."""
    result: Dict[str, Tuple[Path, Path | None]] = {}
    for sname, (r1, r2) in fastq_dict.items():
        if r2 is None:
            raise IPyradError(
                f"--unmate can only be used with paired-end FASTQ inputs; sample '{sname}' is single-end."
            )
        out = tmpdir / f"{sname}.tmp.unmated.fastq.gz"
        _concat_fastqs([r1, r2], out)
        result[sname] = (out, None)
    return result


def prepare_map_samples(
    fastqs,
    delim_str: str | None,
    delim_idx: int,
    imap: Path | None,
    tmpdir: Path,
    unmate: bool = False,
) -> Tuple[Dict[str, Tuple[Path, Path | None]], bool]:
    """Parse, optionally materialize, and validate mapper samples."""
    fastq_dict = get_name_to_fastq_dict(fastqs, delim_str, delim_idx)
    if imap is not None:
        fastq_dict = apply_imap_to_samples(imap, tmpdir, fastq_dict)
    validate_fastq_inputs(fastq_dict)
    is_paired = detect_is_paired(fastq_dict)
    if unmate:
        if not is_paired:
            raise IPyradError("--unmate can only be used with paired-end FASTQ inputs.")
        fastq_dict = unmate_paired_samples(fastq_dict, tmpdir)
        is_paired = False
    return fastq_dict, is_paired
