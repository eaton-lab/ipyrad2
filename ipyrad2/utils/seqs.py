#!/usr/bin/env python

import string

def comp(seq: str) -> str:
    """Returns a seq complemented. Preserves little n's for denovo inserts."""
    return seq.replace("A", 't')\
              .replace('T', 'a')\
              .replace('C', 'g')\
              .replace('G', 'c')\
              .replace('n', 'Z')\
              .upper()\
              .replace("Z", "n")

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
