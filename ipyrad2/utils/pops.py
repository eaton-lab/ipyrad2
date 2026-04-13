#!/usr/bin/env python

import pandas as pd
from fnmatch import fnmatchcase
from pathlib import Path
from typing import Dict, List, Sequence, Tuple
from .exceptions import IPyradError


def parse_pops_file(popfile: Path) -> Tuple[Dict[str, List[str]], Dict[str, int]]:
    """Parse an ipyrad1-style pop_assign_file and return a dictionary
    mapping pop names to tuples of minsamples and a list of sample
    names. The value of minsamples is used by grouped-calling and
    analysis workflows that enforce per-population coverage rules.

    popfile format:

            ind1 pop1
            ind2 pop2
            ind3 pop3
            ...
            indN popN
            # pop1:3 pop2:3 pop3:3

    Returns 2 dictionaries keying populations to minsamples (minmap)
    and populations to lists of sample names or glob patterns (imap).
    """
    # Offload some of the parsing to the `parse_imap()` function which is used
    # by the new CLI mode and expects a file with _only_ the sample/pop mapping
    popdict = parse_imap(popfile)

    ## parse minsamples per population if present (line with #)
    mindat = [
        i.lstrip("#").lstrip().rstrip() for i in
        open(popfile, 'r').readlines() if i.startswith("#")]

    if mindat:
        popmins = {}
        for i in range(len(mindat)):
            minlist = mindat[i].replace(",", "").split()
            popmins.update({i.split(':')[0]: int(i.split(':')[1])
                                    for i in minlist})
        if set(popdict.keys()) != set(popmins.keys()):
            raise IPyradError(IMAP_MINMAP_DISAGREE.format(popdict.keys(),
                                                          popmins.keys()))
    else:
        raise IPyradError(MIN_SAMPLES_PER_POP_MALFORMED)

    imap = {i:popdict[i] for i in popdict}
    minmap = {i:popmins[i] for i in popdict}

    return imap, minmap


def parse_imap(popfile: Path) -> Dict[str, List[str]]:
    """Return an imap dict mapping pop names to list of sample names.

    Each line is a sample name or glob pattern mapped to a str group.
    Many samples can be mapped to the same group name.

    Format
    ------
    sample\tgroup
    sample\tgroup
    """

    try:
        ## parse populations file
        popdat = pd.read_csv(
            popfile, header=None,
            sep=r"\s+",
            names=["inds", "pops"],
            comment="#",
            dtype=str)

        popdict = {
            key: group.inds.values.tolist() for key, group in
            popdat.groupby("pops")}


    except FileNotFoundError as exc:
        raise IPyradError(f"Populations file not found - {popfile}") from exc
    except OSError as exc:
        raise IPyradError(f"Failed to read populations file - {popfile}") from exc
    except ValueError as exc:
        raise IPyradError(f"  Populations file malformed - {popfile}") from exc

    return popdict


def _has_glob_magic(value: str) -> bool:
    """Return True when a sample token should be treated as a glob pattern."""
    return any(char in value for char in "*?[")


def expand_imap_patterns(
    imap: Dict[str, List[str]],
    available_names: Sequence[str],
    *,
    mapping_name: str = "IMAP",
    available_name: str = "the available samples",
    strict_unmatched: bool = True,
) -> tuple[Dict[str, List[str]], List[str]]:
    """Expand exact names or globs in an IMAP against available sample names."""
    available_names = [str(name) for name in available_names]
    available_set = set(available_names)
    sample_to_groups: dict[str, set[str]] = {}
    expanded: Dict[str, List[str]] = {}
    unmatched: list[str] = []

    for group, entries in imap.items():
        resolved: list[str] = []
        seen: set[str] = set()
        for entry in entries:
            if _has_glob_magic(entry):
                matches = [
                    name for name in available_names
                    if fnmatchcase(name, entry)
                ]
            elif entry in available_set:
                matches = [entry]
            else:
                matches = []

            if not matches:
                unmatched.append(entry)
                continue

            for name in matches:
                if name in seen:
                    continue
                seen.add(name)
                resolved.append(name)
                sample_to_groups.setdefault(name, set()).add(group)

        if resolved:
            expanded[group] = resolved

    duplicate_samples = sorted(
        name for name, groups in sample_to_groups.items()
        if len(groups) > 1
    )
    if duplicate_samples:
        raise IPyradError(
            f"{mapping_name} assigns sample(s) multiple times: "
            + ", ".join(duplicate_samples)
        )

    unmatched = list(dict.fromkeys(unmatched))
    if strict_unmatched and unmatched:
        raise IPyradError(
            f"{mapping_name} contains sample names or glob patterns that were not found in "
            f"{available_name}: {', '.join(unmatched[:10])}"
        )

    return expanded, unmatched


def parse_minmap(path: Path) -> Dict[str, int]:
    """Return a minmap dict mapping groups to int values.

    Each line is a unique group mapped to an int value.

    Format
    ------
    group\tsize
    group\tsize
    """
    minmap: Dict[str, int] = {}

    try:
        with open(path, "r", encoding="utf-8") as infile:
            for lineno, line in enumerate(infile, start=1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split()
                if len(parts) != 2:
                    raise IPyradError(
                        f"Minmap file malformed at line {lineno}: {path}"
                    )

                group, size = parts
                if group in minmap:
                    raise IPyradError(
                        f"Minmap file contains duplicate group '{group}': {path}"
                    )

                minmap[group] = int(size)

    except ValueError as exc:
        raise IPyradError(f"Minmap file malformed - {path}") from exc
    except OSError as exc:
        raise IPyradError(f"Failed to read minmap file - {path}") from exc

    if not minmap:
        raise IPyradError(f"Minmap file is empty - {path}")

    return minmap

IMAP_MINMAP_DISAGREE = """\n
    The populations specified in the main body of the populations file do
    not agree with the populations specified on the minmap line (final line
    of the imap file.

    populations for all samples: {}
    populations in minmap: {}
"""


MIN_SAMPLES_PER_POP_MALFORMED = """\n\
    Population assignment file must include a line indicating the minimum
    number of samples per population. This line should come at the end
    of the file and should be preceded by a hash sign (#), e.g.:

    # pop1:3 pop2:3 pop3:3
    """


if __name__ == "__main__":
    pass

    test1 = Path("../tests/test_map1_imap.tsv")
    test2 = Path("../tests/test_map2_imap.tsv")

    imap1 = parse_imap(test1)
    print(imap1)
