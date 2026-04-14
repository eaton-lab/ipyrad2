#!/usr/bin/env python

"""Benchmark denovo split-stage memory and runtime on an existing _denovo_work."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess as sp
import sys
import tempfile
import time
from pathlib import Path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Copy an existing denovo _denovo_work into a temporary run directory, "
            "run ipyrad2.denovo.graph.make_global_tables on it, and report wall time "
            "plus sampled peak RSS across the process tree."
        )
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        required=True,
        help="Path to an existing _denovo_work directory.",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="Cores passed to make_global_tables. [default=%(default)s]",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Log level passed to make_global_tables. [default=%(default)s]",
    )
    parser.add_argument(
        "--sample-interval",
        type=float,
        default=0.1,
        help="Seconds between RSS samples while the child is running. [default=%(default)s]",
    )
    parser.add_argument(
        "--keep-copy",
        action="store_true",
        help="Keep the copied benchmark run directory instead of deleting it.",
    )
    parser.add_argument(
        "--copy-root",
        type=Path,
        default=None,
        help="Optional directory under which the temporary benchmark copy is created.",
    )
    parser.add_argument(
        "--json-out",
        type=Path,
        default=None,
        help="Optional path to write the benchmark result JSON.",
    )
    return parser.parse_args()


def _read_children(pid: int) -> list[int]:
    path = Path("/proc") / str(pid) / "task" / str(pid) / "children"
    try:
        text = path.read_text(encoding="utf-8").strip()
    except OSError:
        return []
    if not text:
        return []
    return [int(part) for part in text.split()]


def _read_rss_kb(pid: int) -> int:
    path = Path("/proc") / str(pid) / "status"
    try:
        with open(path, "rt", encoding="utf-8") as infile:
            for line in infile:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except OSError:
        return 0
    return 0


def _collect_tree_rss_kb(root_pid: int) -> tuple[int, list[int]]:
    seen: set[int] = set()
    stack = [root_pid]
    total = 0
    ordered: list[int] = []
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        rss_kb = _read_rss_kb(pid)
        if rss_kb <= 0:
            continue
        total += rss_kb
        ordered.append(pid)
        stack.extend(_read_children(pid))
    return total, sorted(ordered)


def _copy_workdir(source: Path, copy_root: Path | None) -> Path:
    if copy_root is None:
        tmpdir = Path(tempfile.mkdtemp(prefix="denovo-split-bench-"))
    else:
        copy_root.mkdir(parents=True, exist_ok=True)
        tmpdir = Path(tempfile.mkdtemp(prefix="denovo-split-bench-", dir=copy_root))
    copied = tmpdir / source.name
    shutil.copytree(source, copied)
    return copied


def _build_worker_code(
    workdir: Path,
    cores: int,
    log_level: str,
) -> str:
    return (
        "from pathlib import Path\n"
        "from ipyrad2.denovo.graph import make_global_tables\n"
        f"make_global_tables(Path({workdir.as_posix()!r}), "
        f"cores={cores}, "
        f"log_level={log_level!r})\n"
    )


def main() -> int:
    args = _parse_args()
    source = args.workdir.expanduser().resolve()
    if source.name != "_denovo_work":
        raise SystemExit(f"--workdir must point to a _denovo_work directory: {source}")
    if not source.is_dir():
        raise SystemExit(f"--workdir does not exist: {source}")

    copied_workdir = _copy_workdir(source, args.copy_root.expanduser().resolve() if args.copy_root else None)
    run_root = copied_workdir.parent

    code = _build_worker_code(
        copied_workdir,
        cores=args.cores,
        log_level=args.log_level,
    )
    cmd = [sys.executable, "-c", code]

    started = time.time()
    proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    peak_rss_kb = 0
    peak_pids: list[int] = [proc.pid]
    sample_count = 0

    while proc.poll() is None:
        rss_kb, pids = _collect_tree_rss_kb(proc.pid)
        if rss_kb > peak_rss_kb:
            peak_rss_kb = rss_kb
            peak_pids = pids
        sample_count += 1
        time.sleep(args.sample_interval)

    stdout, stderr = proc.communicate()
    elapsed = time.time() - started
    rss_kb, pids = _collect_tree_rss_kb(proc.pid)
    if rss_kb > peak_rss_kb:
        peak_rss_kb = rss_kb
        peak_pids = pids

    result = {
        "source_workdir": str(source),
        "copied_workdir": str(copied_workdir),
        "run_root": str(run_root),
        "cores": int(args.cores),
        "log_level": args.log_level,
        "sample_interval_seconds": float(args.sample_interval),
        "wall_seconds": round(elapsed, 3),
        "peak_rss_kb": int(peak_rss_kb),
        "peak_rss_mb": round(peak_rss_kb / 1024.0, 2),
        "peak_pids": peak_pids,
        "rss_samples_taken": int(sample_count),
        "returncode": int(proc.returncode),
        "stdout": stdout,
        "stderr": stderr,
    }

    if args.json_out is not None:
        args.json_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(result, indent=2))

    if not args.keep_copy:
        shutil.rmtree(run_root, ignore_errors=True)

    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
