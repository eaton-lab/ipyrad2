#!/usr/bin/env python

"""Logger for ipyrad to STDERR and optionally also to a LOGFILE.

logging to STDERR
-----------------
DEBUG: used by developers to examine extra details.
INFO: info reported to users, including progress bars. (DEFAULT)
WARNING: warnings to users, if set to default then progress bars are not shown.
ERROR: sometimes printed along with raised errors.

logging to LOGFILE
------------------
DEBUG: developer stuff
INFO: same as above, w/ some extra info, but not progress bars. (DEFAULT)
same
same

Examples
--------
>>> import ipyrad as ip
>>> ip.set_log_level("DEBUG")
>>> ip.set_log_level("DEBUG", log_file="/tmp/ip-log.txt")

Note
----
Exceptions written to the logfile have color support, which
can be viewed using `less -R logfile.txt`
"""

from __future__ import annotations
from typing import Optional
import sys
from pathlib import Path
from loguru import logger


_CURRENT_LOG_LEVEL = "INFO"


def formatter(record):
    """Custom formatter that allows for progress bar."""
    end = record["extra"].get("end", "\n")
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


def set_log_level(log_level: str = "DEBUG", log_file: Optional[Path] = None):
    """Add logger for ipyrad to stderr and optionally to file.

    logger.info("...")
    logger.bind(to_file=True).info("...")
    """
    global _CURRENT_LOG_LEVEL
    _CURRENT_LOG_LEVEL = normalize_log_level(log_level)
    logger.remove()

    # always log to stderr
    logger.add(
        sink=sys.stderr,
        level=_CURRENT_LOG_LEVEL,
        colorize=color_support(),
        format=formatter,
        enqueue=False,
        # traceback=True,
    )
    # optionally log to file
    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(exist_ok=True)
        log_file.touch(exist_ok=True)
        logger.add(
            sink=str(log_file),
            level=_CURRENT_LOG_LEVEL,
            colorize=False,
            format=formatter,
            enqueue=True,
            rotation="50 MB",
        )
    return logger


def setup_loguru_worker(log_level: str) -> None:
    """initialized on parallel Worker processes."""
    from loguru import logger
    import sys

    global _CURRENT_LOG_LEVEL
    _CURRENT_LOG_LEVEL = normalize_log_level(log_level)
    logger.remove()
    logger.add(
        sys.stderr,
        level=_CURRENT_LOG_LEVEL,
        colorize=color_support(),
        format=formatter,
        enqueue=True,
    )


if __name__ == "__main__":
    pass
