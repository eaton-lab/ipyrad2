#!/usr/bin/env python

from __future__ import annotations
from typing import Final

try: # Python 3.8+
    from importlib.metadata import PackageNotFoundError, version
except Exception: # pragma: no cover - fallback for very old Pythons
    from importlib_metadata import PackageNotFoundError, version # type: ignore

# Keep this in sync with [project].name in pyproject.toml
DIST_NAME: Final[str] = "ipyrad2"

def get_version() -> str:
    """Return installed distribution version or a safe fallback.

    The fallback avoids runtime crashes when metadata is unavailable (e.g., running
    from a source tree without installation).
    """
    try:
        return version(DIST_NAME)
    except PackageNotFoundError: # e.g., not installed, or running from a checkout
        return "0.0.0+unknown"

__version__: Final[str] = get_version()
