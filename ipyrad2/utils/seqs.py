#!/usr/bin/env python


def comp(seq: str) -> str:
    """Returns a seq complemented. Preserves little n's for denovo inserts."""
    return seq.replace("A", 't')\
              .replace('T', 'a')\
              .replace('C', 'g')\
              .replace('G', 'c')\
              .replace('n', 'Z')\
              .upper()\
              .replace("Z", "n")
