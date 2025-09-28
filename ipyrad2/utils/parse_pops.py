#!/usr/bin/env python

from typing import Dict, List
from pathlib import Path


def parse_imap(path: Path) -> Dict[str, List[str]]:
    """Return an imap dict mapping pop names to list of sample names.

    Each line is a unique sample mapped to a str group. Many samples
    can be mapped to the same group name.

    Format
    ------
    sample\tgroup
    sample\tgroup
    """
    pass


def parse_minmap(path: Path) -> Dict[str, int]:
    """Return a minmap dict mapping groups to int values.

    Each line is a unique group mapped to an int value.

    Format
    ------
    group\tsize
    group\tsize
    """
    pass


