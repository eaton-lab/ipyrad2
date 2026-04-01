#!/usr/bin/env python

from __future__ import annotations

try:
    from ._version import __version__
except ImportError:
    try:
        from importlib.metadata import PackageNotFoundError, version
    except Exception:  # pragma: no cover - fallback for very old Pythons
        from importlib_metadata import PackageNotFoundError, version  # type: ignore

    try:
        __version__ = version("ipyrad2")
    except PackageNotFoundError:
        __version__ = "0.0.0+unknown"
