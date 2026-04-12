#!/usr/bin/env python

"""Internal runtime profiling helpers."""

from __future__ import annotations

from contextlib import contextmanager
import resource
import sys
import time

from loguru import logger


def _format_peak_rss(value: int) -> str:
    """Format ru_maxrss values across Linux/macOS into a compact string."""
    factor = 1 if sys.platform == "darwin" else 1024
    nbytes = max(0, int(value)) * factor
    mib = nbytes / (1024 * 1024)
    return f"{mib:.1f} MiB"


def _current_self_rss_bytes() -> int | None:
    """Return the current RSS for the main process when the platform exposes it."""
    if sys.platform.startswith("linux"):
        try:
            with open("/proc/self/status", "rt", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("VmRSS:"):
                        parts = line.split()
                        if len(parts) >= 2:
                            return int(parts[1]) * 1024
        except OSError:
            return None
    return None


def _format_bytes(nbytes: int | None) -> str:
    """Format byte counts consistently for debug profiling logs."""
    if nbytes is None:
        return "n/a"
    mib = nbytes / (1024 * 1024)
    return f"{mib:.1f} MiB"


@contextmanager
def profile_stage(stage: str):
    """Log elapsed time plus RSS deltas for one stage."""
    start_time = time.perf_counter()
    start_self_rss = _current_self_rss_bytes()
    start_self_peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    start_child_peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start_time
        end_self_rss = _current_self_rss_bytes()
        end_self_peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        end_child_peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss

        rss_delta = None
        if start_self_rss is not None and end_self_rss is not None:
            rss_delta = end_self_rss - start_self_rss

        peak_self_delta = end_self_peak - start_self_peak
        peak_child_delta = end_child_peak - start_child_peak
        logger.debug(
            (
                "stage profile {}: elapsed={:.2f}s, "
                "self_rss_start={}, self_rss_end={}, self_rss_delta={}, "
                "self_peak_delta={}, child_peak_delta={}, "
                "self_peak_total={}, child_peak_total={}"
            ),
            stage,
            elapsed,
            _format_bytes(start_self_rss),
            _format_bytes(end_self_rss),
            _format_bytes(rss_delta),
            _format_peak_rss(peak_self_delta),
            _format_peak_rss(peak_child_delta),
            _format_peak_rss(end_self_peak),
            _format_peak_rss(end_child_peak),
        )
