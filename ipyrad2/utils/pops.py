#!/usr/bin/env python

import pandas as pd
from typing import Dict, List, Tuple
from pathlib import Path
from .exceptions import IPyradError


def parse_pops_file(popfile: Path) -> Tuple[Dict[str, List[str]], Dict[str, int]]:
    """Parse an ipyrad1-style pop_assign_file and return a dictionary
    mapping pop names to tuples of minsamples and a list of sample
    names. The value of minsamples is used by denovo for selecting
    the number of samples per population for constructing the
    pseudo-reference and also by assemble for specifying minsamples
    per locus per population for filtering.

    popfile format:

            ind1 pop1
            ind2 pop2
            ind3 pop3
            ...
            indN popN
            # pop1:3 pop2:3 pop3:3

    Returns 2 dictionaries keying populations to minsamples (minmap) and
    populations to lists of sample names (imap).
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
    else:
        raise IPyradError(MIN_SAMPLES_PER_POP_MALFORMED)

    imap = {i:popdict[i] for i in popdict}
    minmap = {i:popmins[i] for i in popdict}

    return imap, minmap


def parse_imap(popfile: Path) -> Dict[str, List[str]]:
    """Return an imap dict mapping pop names to list of sample names.

    Each line is a unique sample mapped to a str group. Many samples
    can be mapped to the same group name.

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


    except (ValueError, IOError):
        raise IPyradError("  Populations file malformed - {}".format(popfile))

    return popdict


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
