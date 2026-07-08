#!/usr/bin/env python

"""Logger for ipyrad to stderr.

logging to stderr
-----------------
DEBUG: used by developers to examine extra details.
INFO: info reported to users, including progress bars. (DEFAULT)
WARNING: warnings to users, if set to default then progress bars are not shown.
ERROR: sometimes printed along with raised errors.

Examples
--------
>>> import ipyrad2 as ip2
>>> ip2.set_log_level("DEBUG")
"""

from __future__ import annotations
import sys
from loguru import logger


_CURRENT_LOG_LEVEL = "INFO"


def formatter(record):
    """Custom formatter that allows for progress bar."""
    end = record["extra"].get("end", "\n")

    if record["extra"].get("progress", False):
        return "{message}" + end

    fmessage = (
        "{time:YYYY-MM-DD HH:mm:ss} | "
        "<level>{level:<8}</level> <white>|</white> "
        "<magenta>{file:<20}</magenta> <white>|</white> "
        "{message}"
    ) + end
    return fmessage


def color_support():
    """Check for color support in stderr as a notebook or terminal/tty."""
    return sys.stderr.isatty()
    # check if we're in IPython/jupyter
    # try:
    #     import IPython
    #     tty1 = bool(IPython.get_ipython())
    # except ImportError:
    #     tty1 = False
    # # check if we're in a terminal
    # tty2 = sys.stderr.isatty()
    # return tty1 or tty2


def normalize_log_level(level: str) -> str:
    """Convert level input to full string if it is a substring, else return user value"""
    OPTIONS = ("TRACE", "DEBUG", "INFO", "WARNING", "ERROR")
    level = level.upper()
    for lvl in OPTIONS:
        if lvl.startswith(level):
            return lvl
    return level


def is_log_level_enabled(level: str) -> bool:
    """Return True when the configured stderr sink includes this level."""
    current_no = logger.level(_CURRENT_LOG_LEVEL).no
    requested_no = logger.level(normalize_log_level(level)).no
    return requested_no >= current_no


def set_log_level(log_level: str = "DEBUG"):
    """Add the shared ipyrad logger to stderr.

    logger.info("...")
    """
    global _CURRENT_LOG_LEVEL
    _CURRENT_LOG_LEVEL = normalize_log_level(log_level)
    logger.remove()

    logger.add(
        sink=sys.stderr,
        level=_CURRENT_LOG_LEVEL,
        colorize=color_support(),
        format=formatter,
        enqueue=False,
    )
    return logger


def setup_loguru_worker(log_level: str) -> None:
    """initialized on parallel Worker processes."""
    global _CURRENT_LOG_LEVEL
    _CURRENT_LOG_LEVEL = normalize_log_level(log_level)

    logger.remove()
    logger.add(
        sys.stderr,
        level=_CURRENT_LOG_LEVEL,
        colorize=color_support(),
        format=formatter,
        enqueue=False,
    )


if __name__ == "__main__":
    pass
