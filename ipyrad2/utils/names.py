#!/usr/bin/env python

"""Parse names from PE or SE filenames.

PE
-----
- custom delim
- split on all '.' and keep longest left substring w/ proper pairing.
- split on all '_' and keep longest left substring w/ proper pairing.
- cut from right and longest with proper pairing.
- else assume SE

SE
----
- custom delim
- cut to find longest unique left substrings
"""


import os
import re
from typing import Dict, List, Tuple, Union
from pathlib import Path
from collections import defaultdict
from loguru import logger
from .exceptions import IPyradError


KNOWN_FILE_SUFFIXES = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq.bz2",
    ".fq.bz2",
    ".fastq",
    ".fq",
    ".bam",
)

MATE_PATTERNS = (
    re.compile(
        r"^(?P<sample>.+?)[._-](?P<mate>r[12]|read[12])"
        r"(?:[._-]001)?(?P<trailing>(?:[._-].+)?)$",
        re.IGNORECASE,
    ),
    re.compile(
        r"^(?P<sample>.+?)[._-](?P<mate>[12])"
        r"(?:[._-]001)?(?P<trailing>(?:[._-].+)?)$",
        re.IGNORECASE,
    ),
    re.compile(
        r"^(?P<sample>.+?)(?P<mate>r[12]|read[12])"
        r"(?:[._-]001)?(?P<trailing>(?:[._-].+)?)$",
        re.IGNORECASE,
    ),
)


def expand_path(p: str | Path) -> Path:
    """Returns an absolute path after expanding ~ and env variables"""
    return Path(os.path.expandvars(str(p))).expanduser().absolute()


def _strip_known_suffix(name: str) -> str:
    """Remove common sequencing/compression suffixes from a filename."""
    lower = name.lower()
    for suffix in KNOWN_FILE_SUFFIXES:
        if lower.endswith(suffix):
            return name[:-len(suffix)]
    return name


def _validate_delim_args(delim: str | None, delim_index: int | None) -> None:
    """Validate user-supplied delimiter arguments before parsing names."""
    if delim is None:
        return
    if delim == "":
        raise IPyradError("delim cannot be an empty string.")
    if not isinstance(delim_index, int):
        raise IPyradError("delim_index must be an integer when delim is set.")
    if delim_index == 0:
        raise IPyradError("delim_index cannot be 0 when delim is set.")


def _normalize_paths(paths: List[Path | str]) -> List[Path]:
    """Normalize mixed string/Path inputs and require at least one path."""
    normalized = [Path(path) if isinstance(path, str) else path for path in paths]
    if not normalized:
        raise IPyradError("No fastq data were provided.")
    return normalized


def _group_paths_by_delim(
    fastqs: List[Path],
    delim: str,
    delim_index: int,
) -> Dict[str, List[Path]]:
    """Group paths by the substring defined by delim and delim_index."""
    groups = defaultdict(list)
    for path in fastqs:
        parts = path.name.split(delim)
        sample_name = delim.join(parts[:delim_index])
        if not sample_name:
            raise IPyradError(
                f"Delimiter parsing produced an empty sample name for '{path.name}'. "
                "Adjust -dx/-di."
            )
        groups[sample_name].append(path)
    return groups


def _parse_mate_token(path: Path) -> tuple[str, int, str] | None:
    """Extract a sample name and mate orientation from common PE filenames."""
    stem = _strip_known_suffix(path.name)
    for pattern in MATE_PATTERNS:
        match = pattern.match(stem)
        if not match:
            continue
        sample_name = match.group("sample").rstrip("._-")
        if not sample_name:
            return None
        mate = 1 if match.group("mate").lower().endswith("1") else 2
        trailing = match.group("trailing") or ""
        return sample_name, mate, trailing
    return None


def _pair_group(paths: List[Path], sample_name: str) -> tuple[Path, Path]:
    """Return a deterministic (R1, R2) tuple for a paired-end sample."""
    ordered = {}
    trailing = None
    for path in sorted(paths, key=lambda item: str(item)):
        parsed = _parse_mate_token(path)
        if parsed is None:
            names = ", ".join(sorted(item.name for item in paths))
            raise IPyradError(
                f"Cannot determine read orientation for paired files in sample "
                f"'{sample_name}': {names}"
            )
        mate = parsed[1]
        current_trailing = parsed[2]
        if trailing is None:
            trailing = current_trailing
        elif trailing != current_trailing:
            names = ", ".join(sorted(item.name for item in paths))
            raise IPyradError(
                f"Paired files for sample '{sample_name}' do not share the same "
                f"trailing suffix after the mate token: {names}"
            )
        if mate in ordered:
            names = ", ".join(sorted(item.name for item in paths))
            raise IPyradError(
                f"Found duplicate read {mate} files for sample '{sample_name}': {names}"
            )
        ordered[mate] = path
    if set(ordered) != {1, 2}:
        names = ", ".join(sorted(item.name for item in paths))
        raise IPyradError(f"Missing read pair for sample '{sample_name}': {names}")
    return ordered[1], ordered[2]


def _build_paired_result(
    groups: Dict[str, List[Path]],
    fastqs: List[Path],
) -> Dict[str, Tuple[Path, Path | None]] | None:
    """Build a paired-end result dict if all groups are valid pairs."""
    if not perfect_pairs(groups, fastqs):
        return None
    return {name: _pair_group(groups[name], name) for name in sorted(groups)}


def _build_single_end_result(
    groups: Dict[str, List[Path]],
    fastqs: List[Path],
) -> Dict[str, Tuple[Path, Path | None]] | None:
    """Build a single-end result dict if all groups are unique."""
    if not all_unique(groups, fastqs):
        return None
    return {name: (groups[name][0], None) for name in sorted(groups)}


def _get_default_single_end_groups(fastqs: List[Path]) -> Dict[str, List[Path]]:
    """Group files by stripped basename for conservative SE parsing."""
    groups = defaultdict(list)
    for path in fastqs:
        sample_name = _strip_known_suffix(path.name)
        if not sample_name:
            raise IPyradError(f"Cannot parse a sample name from '{path.name}'.")
        groups[sample_name].append(path)
    return groups


def _group_paths_by_detected_mates(
    fastqs: List[Path],
) -> tuple[Dict[str, List[Path]], int]:
    """Group paths by detected mate token and count how many parsed as PE."""
    groups = defaultdict(list)
    parsed_count = 0
    for path in fastqs:
        parsed = _parse_mate_token(path)
        if parsed is None:
            continue
        parsed_count += 1
        groups[parsed[0]].append(path)
    return groups, parsed_count


def get_paths_list_from_fastq_str(fastq_paths: Union[Path, List[Path]]) -> List[Path]:
    """Expand fastq_paths str (e.g., 'data/*.gz') into List[Path].
    """
    expanded = []
    # ensure paths is a List[Path] but where the Path elements may be
    # regex path names that have not yet been expanded.

    # ensure it is a list
    if isinstance(fastq_paths, (str, Path)):
        fastq_paths = [fastq_paths]
    else:
        fastq_paths = list(fastq_paths)

    # ensure each is a Path object
    paths = _normalize_paths(fastq_paths)

    # for each Path in paths list expand into a list of Paths
    seen = set()
    for path in paths:
        # expand to a full path
        path = expand_path(path)

        # raise if path is a dir.
        if path.is_dir():
            raise IPyradError(f"{path} is a dir. Use regex to select files in the dir (e.g., './path/*.fastq.gz')")

        # expand a regex operator to possibly match multiple files
        # such as paired-end files.
        try:
            fastqs = sorted(path.parent.glob(path.name), key=lambda item: str(item))
        except ValueError as err:
            raise IPyradError(f"No fastq data match input: {path}") from err
        if not fastqs:
            raise IPyradError(f"No fastq data match input: {path}")
        for fastq in fastqs:
            fastq = expand_path(fastq)
            if fastq.is_dir():
                raise IPyradError(
                    f"{fastq} is a dir. Use regex to select files in the dir "
                    "(e.g., './path/*.fastq.gz')"
                )
            if fastq in seen:
                raise IPyradError(f"Input fastq matched more than once: {fastq}")
            seen.add(fastq)
            expanded.append(fastq)
    return expanded


def get_name_to_fastq_dict(
    fastqs: List[Path | str],
    delim: str | None,
    delim_index: int | None,
    suffix: str | None = None,
) -> Dict[str, Tuple[Path, Path | None]]:
    """
    """
    # expand str to List[Path]
    paths_list = get_paths_list_from_fastq_str(fastqs)

    # parse List[Path] to {name: (Path, Path)} or {name: (Path, None)}
    fastq_dict = get_pairs_or_single_by_trim(paths_list, delim, delim_index)

    # add optional suffix
    if suffix is not None:
        fastq_dict = {f"{i}{suffix}": j for (i, j) in fastq_dict.items()}

    # report to logger
    total = len(fastq_dict)
    fmax = min(10, total)
    logger.info(f"showing first {fmax}/{total} names parsed from file paths")
    max_len = max(len(i) for i in fastq_dict)
    for fidx, name in enumerate(sorted(fastq_dict)):
        paths = fastq_dict[name]
        key_padded = name + " " * (max_len - len(name))
        if paths[1]:
            if fidx < fmax:
                logger.info(f"{key_padded} <- {(paths[0].name, paths[1].name)}")
        else:
            if fidx < fmax:
                logger.info(f"{key_padded} <- {paths[0].name}")
    return fastq_dict


def perfect_pairs(ndict: Dict[str, List[Path]], paths: List[Path]) -> bool:
    """valid PE name"""
    if not paths or not ndict:
        return False
    if sum(len(v) for v in ndict.values()) != len(paths):
        return False
    return all(name and len(v) == 2 for name, v in ndict.items())


def all_unique(ndict: Dict[str, List[Path]], paths: List[Path]) -> bool:
    """valid SE name"""
    if not paths or not ndict:
        return False
    if sum(len(v) for v in ndict.values()) != len(paths):
        return False
    return all(name and len(v) == 1 for name, v in ndict.items())


def get_pairs_or_single_by_trim(
    fastqs: List[Path | str],
    delim: str | None,
    delim_index: int | None,
) -> Dict[str, Tuple[Path, Path | None]]:
    """..."""
    fastqs = _normalize_paths(fastqs)
    _validate_delim_args(delim, delim_index)
    delim_groups = None

    # try to pair sample by delim args
    if delim is not None:
        delim_groups = _group_paths_by_delim(fastqs, delim, delim_index)
        paired = _build_paired_result(delim_groups, fastqs)
        if paired is not None:
            logger.info(f"paired files by user args: -dx={delim} -di={delim_index}")
            return paired
        logger.info("pairing files by name-delim failed. Falling back to auto-detection.")

    # conservatively pair only filenames with recognizable mate tokens.
    names_to_paths, parsed_count = _group_paths_by_detected_mates(fastqs)
    paired = _build_paired_result(names_to_paths, fastqs)
    if paired is not None:
        logger.info("paired files by auto-detecting mate tokens in filenames")
        return paired
    if parsed_count:
        raise IPyradError(
            "Cannot safely pair files by name. Some filenames look paired-end "
            "but do not form complete R1/R2 pairs. Try setting the delim args "
            "explicitly, or enter paths with different naming conventions as "
            "separate inputs."
        )

    # --------------------------------------------------------------
    logger.info("failed to pair files, assuming data is single-end")
    # --------------------------------------------------------------

    # try to get unique SE names using delim args
    if delim is not None:
        if all_unique(delim_groups, fastqs):
            logger.info(f"parsed names by user args: -dx={delim} -di={delim_index}")
            return {i: (j[0], None) for i, j in sorted(delim_groups.items())}
        logger.info("parsing names by user args failed. Falling back to auto-detection.")

    names_to_paths = _get_default_single_end_groups(fastqs)
    single_end = _build_single_end_result(names_to_paths, fastqs)
    if single_end is not None:
        logger.info("parsed names by stripping known file suffixes")
        return single_end

    raise IPyradError(
        "Cannot parse names from file names, likely because they do not share "
        "unique sample names after removing common sequencing suffixes. Try "
        "setting the delim args explicitly, or entering paths with different "
        "naming conventions as separate inputs, e.g., '-d *.fq.gz *.fastq.gz'."
    )


if __name__ == "__main__":

    pass
    # from ipyrad3.utils.logger import set_log_level
    # set_log_level("DEBUG")

    # path = Path("/home/deren/Documents/tools/ipyrad2/examples/Pedic-PE-ddRAD/*_R*")
    # get_fastq_tuples_dict_from_paths_list(path)

    a = Path("/home/deren/Documents/ipyrad-tests/examples/Pedic-PE-ddRAD/c*_R*")
    b = Path("/home/deren/Documents/ipyrad-tests/examples/Pedic-PE-ddRAD/m*_R*")
    fastqs = [a, b]
    # fastqs = Path("/tmp/TRIM_SE/*.gz")
    # fastqs = Path("/tmp/MAP_PE/*.bam")
    pdict = get_name_to_fastq_dict(fastqs, None, None)
    # print(pdict)
