#!/usr/bin/env python

import pandas as pd
from typing import Dict, List
from pathlib import Path
from .exceptions import IPyradError


def parse_pops_file(popfile: Path) -> Dict[str, (int, List[str])]:
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

    return  {i: (popmins[i], popdict[i]) for i in popdict}


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
    pass


MIN_SAMPLES_PER_POP_MALFORMED = """\n\
    Population assignment file must include a line indicating the minimum
    number of samples per population. This line should come at the end
    of the file and should be preceded by a hash sign (#), e.g.:

    # pop1:3 pop2:3 pop3:3
    """
