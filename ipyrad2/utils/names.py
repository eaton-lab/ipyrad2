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
from collections import Counter, defaultdict
from loguru import logger
from .exceptions import IPyradError


KNOWN_FILE_SUFFIXES = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq",
    ".fq",
)

WORKFLOW_SAMPLE_SUFFIXES = (
    ".trimmed",
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

LITERAL_MATE_PATTERNS = (
    re.compile(
        r"^(?P<prefix>.+?)[._-](?P<mate>r[12]|read[12])"
        r"(?P<suffix>(?:[._-].+)?)$",
        re.IGNORECASE,
    ),
    re.compile(
        r"^(?P<prefix>.+?)[._-](?P<mate>[12])"
        r"(?P<suffix>(?:[._-].+)?)$",
        re.IGNORECASE,
    ),
    re.compile(
        r"^(?P<prefix>.+?)(?P<mate>r[12]|read[12])"
        r"(?P<suffix>(?:[._-].+)?)$",
        re.IGNORECASE,
    ),
)

PAIR_WARNING_MIN_COMPLETE_FRACTION = 0.8


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


def _strip_workflow_sample_suffix(name: str) -> str:
    """Remove one recognized workflow suffix from a parsed FASTQ sample name."""
    for suffix in WORKFLOW_SAMPLE_SUFFIXES:
        if not name.endswith(suffix):
            continue
        trimmed = name[:-len(suffix)]
        if trimmed:
            return trimmed
    return name


def normalize_workflow_sample_name(name: str) -> str:
    """Return one canonical sample name after internal workflow normalization."""
    return _strip_workflow_sample_suffix(name)


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
    parsed = _parse_mate_token_from_stem(stem)
    if parsed is not None:
        return parsed
    trimmed_stem = stem.rstrip("._-")
    if trimmed_stem == stem:
        return None
    return _parse_mate_token_from_stem(trimmed_stem)


def _parse_mate_token_from_stem(stem: str) -> tuple[str, int, str] | None:
    """Extract a sample name and mate orientation from a filename stem."""
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


def _parse_mate_token_for_delim_pairing(path: Path) -> tuple[str, int, str] | None:
    """Parse mate tokens for user-delimited grouping, tolerating trailing separators."""
    parsed = _parse_mate_token(path)
    if parsed is not None:
        return parsed
    stem = _strip_known_suffix(path.name).rstrip("._-")
    if not stem:
        return None
    return _parse_mate_token_from_stem(stem)


def _parse_literal_mate_token(path: Path) -> tuple[str, int] | None:
    """Extract a mate orientation and stem with only the mate token removed."""
    stem = _strip_known_suffix(path.name)
    for pattern in LITERAL_MATE_PATTERNS:
        match = pattern.match(stem)
        if not match:
            continue
        prefix = match.group("prefix").rstrip("._-")
        if not prefix:
            return None
        suffix = match.group("suffix") or ""
        sample_name = f"{prefix}{suffix}"
        if not sample_name:
            return None
        mate = 1 if match.group("mate").lower().endswith("1") else 2
        return sample_name, mate
    return None


def _pair_group(
    paths: List[Path],
    sample_name: str,
    parser=_parse_mate_token,
) -> tuple[Path, Path]:
    """Return a deterministic (R1, R2) tuple for a paired-end sample."""
    ordered = {}
    trailing = None
    for path in sorted(paths, key=lambda item: str(item)):
        parsed = parser(path)
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
    try:
        return {name: _pair_group(groups[name], name) for name in sorted(groups)}
    except IPyradError:
        return None


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


def _group_paths_by_literal_mate_names(
    fastqs: List[Path],
) -> tuple[Dict[str, List[Path]], int]:
    """Group paths by literal filename stem after removing only the mate token."""
    groups = defaultdict(list)
    parsed_count = 0
    for path in fastqs:
        parsed = _parse_literal_mate_token(path)
        if parsed is None:
            continue
        parsed_count += 1
        groups[parsed[0]].append(path)
    return groups, parsed_count


def _format_group_examples(groups: Dict[str, List[Path]], limit: int = 5) -> str:
    """Render a short sample-to-files summary for user-facing parse errors."""
    examples = []
    for sample_name in sorted(groups)[:limit]:
        files = ", ".join(
            path.name for path in sorted(groups[sample_name], key=lambda item: str(item))
        )
        examples.append(f"{sample_name}: {files}")
    if len(groups) > limit:
        examples.append(f"... and {len(groups) - limit} more")
    return "; ".join(examples)


def _paired_name_error_message(
    fastqs: List[Path],
    groups: Dict[str, List[Path]],
    parsed_count: int,
    parser,
    source: str,
) -> str:
    """Return a detailed error for incomplete or inconsistent paired-name evidence."""
    details = []
    for sample_name in sorted(groups):
        paths = sorted(groups[sample_name], key=lambda item: str(item))
        mate_counts = Counter()
        parsed_names = []
        for path in paths:
            parsed = parser(path)
            if parsed is None:
                continue
            mate_counts[parsed[1]] += 1
            parsed_names.append(path.name)
        if not mate_counts:
            continue
        if mate_counts.get(1, 0) == 0 or mate_counts.get(2, 0) == 0:
            missing = "R1" if mate_counts.get(1, 0) == 0 else "R2"
            details.append(
                f"{sample_name}: missing {missing} ({', '.join(parsed_names)})"
            )
            continue
        duplicates = [
            f"R{mate} x{count}"
            for mate, count in sorted(mate_counts.items())
            if count > 1
        ]
        if duplicates:
            details.append(f"{sample_name}: duplicate mates ({', '.join(duplicates)})")

    if parsed_count < len(fastqs):
        unparsed = ", ".join(
            sorted(path.name for path in fastqs if parser(path) is None)
        )
        details.append(f"unrecognized alongside paired-end-looking files: {unparsed}")

    detail_text = "; ".join(details[:5])
    if len(details) > 5:
        detail_text += f"; and {len(details) - 5} more"

    return (
        f"Cannot safely pair files by {source}. Some filenames look paired-end "
        f"but do not form complete consistent R1/R2 pairs. {detail_text}"
    )


def _is_complete_pair_group(paths: List[Path], parser) -> bool:
    """Return whether paths contain exactly one matching R1/R2 pair."""
    if len(paths) != 2:
        return False
    mates = set()
    trailing = None
    for path in paths:
        parsed = parser(path)
        if parsed is None:
            return False
        mate = parsed[1]
        if mate in mates:
            return False
        mates.add(mate)
        if len(parsed) > 2:
            current_trailing = parsed[2]
            if trailing is None:
                trailing = current_trailing
            elif trailing != current_trailing:
                return False
    return mates == {1, 2}


def _count_complete_pair_files(groups: Dict[str, List[Path]], parser) -> int:
    """Return the number of files in groups that form valid complete pairs."""
    complete = 0
    for paths in groups.values():
        if _is_complete_pair_group(paths, parser):
            complete += len(paths)
    return complete


def _warn_incomplete_pair_evidence(
    fastqs: List[Path],
    evidences,
) -> None:
    """Warn for incomplete pairing only when complete-pair coverage is strong."""
    scored = []
    for groups, parsed_count, parser, source in evidences:
        complete_pair_files = _count_complete_pair_files(groups, parser)
        if not complete_pair_files:
            continue
        scored.append((complete_pair_files, parsed_count, groups, parser, source))
    if not scored:
        return

    complete_pair_files, parsed_count, groups, parser, source = max(
        scored,
        key=lambda item: (item[0], item[1]),
    )
    if (complete_pair_files / len(fastqs)) < PAIR_WARNING_MIN_COMPLETE_FRACTION:
        return
    message = _paired_name_error_message(
        fastqs=fastqs,
        groups=groups,
        parsed_count=parsed_count,
        parser=parser,
        source=source,
    )
    logger.warning(
        "{} Proceeding as single-end; complete auto-detected pairs cover {}/{} "
        "input files.",
        message,
        complete_pair_files,
        len(fastqs),
    )


def get_paths_list_from_fastq_str(
    fastq_paths: Union[str, Path, List[str | Path]],
) -> List[Path]:
    """Expand FASTQ paths or glob patterns into a concrete ordered list of paths."""
    expanded = []
    # ensure paths is a List[Path] but where the Path elements may still be
    # glob patterns that have not yet been expanded.

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
            raise IPyradError(
                f"{path} is a dir. Use a glob pattern to select files in the dir "
                "(e.g., './path/*.fastq.gz')"
            )

        # expand a glob pattern to possibly match multiple files
        try:
            fastqs = sorted(path.parent.glob(path.name), key=lambda item: str(item))
        except ValueError as err:
            raise IPyradError(f"No fastq data match input: {path}") from err
        if not fastqs:
            raise IPyradError(f"No fastq data match input: {path}")
        for fastq in fastqs:
            fastq = expand_path(fastq)
            if fastq.is_symlink():
                try:
                    fastq.stat()
                except OSError as err:
                    try:
                        target = os.readlink(fastq)
                    except OSError:
                        target = "<unreadable>"
                    raise IPyradError(
                        "FASTQ input symlink cannot be resolved: "
                        f"{fastq} -> {target}. Filesystem error: {err}"
                    ) from err
            if fastq.is_dir():
                raise IPyradError(
                    f"{fastq} is a dir. Use a glob pattern to select files in the dir "
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
    """Return parsed sample names mapped to single-end or paired-end FASTQ tuples."""
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


def normalize_parsed_fastq_sample_names(
    fastq_dict: Dict[str, Tuple[Path, Path | None]],
) -> Dict[str, Tuple[Path, Path | None]]:
    """Normalize parsed FASTQ names to canonical workflow sample names."""
    normalized: Dict[str, Tuple[Path, Path | None]] = {}
    sources: dict[str, list[str]] = defaultdict(list)

    for parsed_name, fastq_tuple in fastq_dict.items():
        canonical_name = normalize_workflow_sample_name(parsed_name)
        sources[canonical_name].append(parsed_name)
        if canonical_name not in normalized:
            normalized[canonical_name] = fastq_tuple

    collisions = {
        canonical_name: sorted(set(parsed_names))
        for canonical_name, parsed_names in sources.items()
        if len(set(parsed_names)) > 1
    }
    if collisions:
        detail = "; ".join(
            f"{canonical_name} <- {', '.join(parsed_names)}"
            for canonical_name, parsed_names in sorted(collisions.items())
        )
        raise IPyradError(
            "Parsed FASTQ sample names collide after internal workflow-suffix normalization. "
            f"Revise input filenames or parsing args so canonical sample names stay unique: {detail}"
        )

    renamed = {
        parsed_name: canonical_name
        for canonical_name, parsed_names in sources.items()
        for parsed_name in parsed_names
        if parsed_name != canonical_name
    }
    if renamed:
        logger.info(
            "normalized {} parsed FASTQ sample name(s) by stripping recognized workflow suffixes",
            len(renamed),
        )
        for parsed_name, canonical_name in sorted(renamed.items()):
            logger.info("{} -> {}", parsed_name, canonical_name)

    return normalized


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
    """Parse FASTQ filenames into deterministic SE or PE sample groupings."""
    fastqs = _normalize_paths(fastqs)
    _validate_delim_args(delim, delim_index)
    delim_groups = None

    # try to pair samples by explicit delimiter args
    if delim is not None:
        delim_groups = _group_paths_by_delim(fastqs, delim, delim_index)
        if perfect_pairs(delim_groups, fastqs):
            try:
                paired = {
                    name: _pair_group(
                        delim_groups[name],
                        name,
                        parser=_parse_mate_token_for_delim_pairing,
                    )
                    for name in sorted(delim_groups)
                }
            except IPyradError as err:
                raise IPyradError(
                    f"Delimiter parsing with -dx={delim} -di={delim_index} grouped files into "
                    f"ambiguous non-paired sample names. Adjust -dx/-di. "
                    f"{_format_group_examples(delim_groups)}"
                ) from err
            logger.info(f"paired files by user args: -dx={delim} -di={delim_index}")
            return paired
        if not all_unique(delim_groups, fastqs):
            raise IPyradError(
                f"Delimiter parsing with -dx={delim} -di={delim_index} grouped files into "
                f"ambiguous non-paired sample names. Adjust -dx/-di. "
                f"{_format_group_examples(delim_groups)}"
            )
        logger.info("pairing files by name-delim failed. Falling back to auto-detection.")

    # conservatively pair only filenames with recognizable mate tokens.
    names_to_paths, parsed_count = _group_paths_by_detected_mates(fastqs)
    paired = _build_paired_result(names_to_paths, fastqs)
    if paired is not None:
        logger.info("paired files by auto-detecting mate tokens in filenames")
        return paired

    literal_names_to_paths, literal_parsed_count = _group_paths_by_literal_mate_names(fastqs)
    paired = _build_paired_result(literal_names_to_paths, fastqs)
    if paired is not None:
        logger.info("paired files by secondary mate-token fallback")
        return paired
    evidences = []
    if parsed_count:
        evidences.append(
            (
                names_to_paths,
                parsed_count,
                _parse_mate_token,
                "auto-detected mate tokens",
            )
        )
    if literal_parsed_count:
        evidences.append(
            (
                literal_names_to_paths,
                literal_parsed_count,
                _parse_literal_mate_token,
                "secondary mate-token fallback",
            )
        )
    if evidences:
        _warn_incomplete_pair_evidence(fastqs, evidences)

    # --------------------------------------------------------------
    logger.info("failed to pair files, assuming data in single-end")
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
