#!/usr/bin/env python

from argparse import ArgumentTypeError


def make_wide(formatter, width: int = 120, max_help_position: int = 36):
    """Return a wider HelpFormatter, if possible."""
    try:
        # https://stackoverflow.com/a/5464440
        # beware: "Only the name of this class is considered a public API."
        kwargs = {'width': width, 'max_help_position': max_help_position}
        formatter(None, **kwargs)
        return lambda prog: formatter(prog, **kwargs)
    except TypeError:
        return formatter


def intlike(s: str) -> int:
    """Allows int or float to be entered but converts to an int"""
    try:
        return int(round(float(s)))                # 3.5 -> 4; 1e5 -> 100_000
    except ValueError:
        raise ArgumentTypeError(f"{s!r} is not a number")


def floatlike(s: str) -> float:
    """Allows int or float to be entered but converts to a float."""
    try:
        return float(s)
    except ValueError:
        raise ArgumentTypeError(f"{s!r} is not a number")
