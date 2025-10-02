#!/usr/bin/env python

import string

COMP = str.maketrans("ACGTacgtnN", "TGCAtgcanN")

def comp(seq: str) -> str:
    return seq.translate(COMP)

def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]

# used in demux to fix sample names.
BADCHARS = (
    string.punctuation
    .replace("_", "")
    .replace("-", "")
    .replace(".", "") + " "
)

# used in demux to resolve ambiguous cutters
AMBIGS = {
    "R": ("G", "A"),
    "K": ("G", "T"),
    "S": ("G", "C"),
    "Y": ("T", "C"),
    "W": ("T", "A"),
    "M": ("C", "A"),
}

IUPAC = {
    frozenset(("A","G")) : "R",
    frozenset(("C","T")) : "Y",
    frozenset(("G","C")) : "S",
    frozenset(("A","T")) : "W",
    frozenset(("G","T")) : "K",
    frozenset(("A","C")) : "M",
}


if __name__ == "__main__":

    S = "ACTGGnnnnAAATTTCCCGGG"
    print(S)
    print(comp(S))
    print(revcomp(S))
