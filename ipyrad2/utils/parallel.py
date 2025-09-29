#!/usr/bin/env python

# POSIX-only robust ProcessPool with hard Ctrl-C and subprocess cleanup.
# Formatted with Black-style line breaks/indentation.

from typing import Any, List, Optional, Sequence, Tuple, Callable, Dict
import atexit
import multiprocessing as mp
import queue
import os
import signal
import subprocess as sp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from .logger import setup_loguru_worker

# ---------- Worker-side (signal-safe; kills subprocess groups) ----------

_CHILD_PROCS: set[sp.Popen] = set()


def safe_popen(argv: Sequence[str], **kwargs) -> sp.Popen:
    """start Popen with process group and store in global set."""
    kwargs.setdefault("preexec_fn", os.setsid)  # start new process group
    p = sp.Popen(argv, **kwargs)
    _CHILD_PROCS.add(p)
    return p


def _kill_all_children(sig: int = signal.SIGTERM) -> None:
    """kill all Popens in global set."""
    for p in list(_CHILD_PROCS):
        try:
            if p.poll() is None:
                os.killpg(p.pid, sig)
        except Exception:
            pass
    for p in list(_CHILD_PROCS):
        try:
            if p.poll() is None:
                os.killpg(p.pid, signal.SIGKILL)
        except Exception:
            pass
    _CHILD_PROCS.clear()


def _worker_signal_handler(signum, _frame) -> None:
    """Callback function to kill on signal."""
    _kill_all_children()
    raise SystemExit(128 + signum)


def _init_worker_with_pid(pid_queue: "mp.queues.Queue") -> None:
    """initialization function to store pids and register killer."""
    pid_queue.put(os.getpid())
    signal.signal(signal.SIGINT, _worker_signal_handler)
    signal.signal(signal.SIGTERM, _worker_signal_handler)
    atexit.register(_kill_all_children)

    setup_loguru_worker()


# ---------- Pipeline runner (supports optional outfile) ----------


def run_pipeline(cmds: List[Sequence[str]], outfile: Optional[Path] = None) -> Tuple[int, bytes, bytes]:
    """Run a stdout->stdin pipeline (list of argv lists).

    If 'outfile' is provided, final stage stdout is streamed to that file,
    and the returned 'out' is b''.
    Otherwise, 'out' captures the final stage's stdout in memory.

    Returns: (rc_of_last_stage, out_bytes, err_bytes_of_last_stage)
    Raises RuntimeError on failure (last stage rc!=0).
    """
    procs: List[sp.Popen] = []
    fout = None
    try:
        prev = None
        for i, argv in enumerate(cmds):
            is_last = i == len(cmds) - 1
            if is_last and outfile is not None:
                outfile.parent.mkdir(parents=True, exist_ok=True)
                fout = outfile.open("wb")
                p = safe_popen(
                    argv,
                    stdin=None if prev is None else prev.stdout,
                    stdout=fout,
                    stderr=sp.PIPE,
                    text=False,
                )
            else:
                p = safe_popen(
                    argv,
                    stdin=None if prev is None else prev.stdout,
                    stdout=sp.PIPE,
                    stderr=sp.PIPE,
                    text=False,
                )
            if prev is not None and prev.stdout is not None:
                prev.stdout.close()
            procs.append(p)
            prev = p

        last = procs[-1]
        if outfile is not None:
            _, err = last.communicate()
            rc = last.returncode
            for p in procs[:-1]:
                try:
                    if p.stderr:
                        p.stderr.read()
                except Exception:
                    pass
            if fout is not None:
                fout.flush()
                fout.close()
                fout = None
            if rc != 0:
                raise RuntimeError(
                    f"pipeline failed (rc={rc}): {cmds[-1]}\n"
                    f"{err.decode(errors='replace')}"
                )
            return rc, b"", err

        out, err = last.communicate()
        rc = last.returncode
        for p in procs[:-1]:
            try:
                if p.stderr:
                    p.stderr.read()
            except Exception:
                pass
        if rc != 0:
            raise RuntimeError(
                f"pipeline failed (rc={rc}): {cmds[-1]}\n"
                f"{err.decode(errors='replace')}"
            )
        return rc, out, err

    except Exception:
        _kill_all_children()
        raise
    finally:
        try:
            if fout is not None:
                fout.close()
        except Exception:
            pass
        for p in procs:
            try:
                if p.poll() is None:
                    p.wait(timeout=0.1)
            except Exception:
                pass


# ---------- Parent-side wrapper (always uses 'spawn'; no option exposed) ----------


def run_with_pool(func: Callable[[Any], Any], jobs: Dict[Any, Any], log_level: str, max_workers: int | None = None) -> List[Any]:
    """
    Execute func(**job) over 'jobs' on a ProcessPoolExecutor.
    - Ctrl-C in parent: SIGTERM then SIGKILL to workers; workers kill their child procs.
    - Exceptions in workers propagate via fut.result().
    Returns results in submission order.
    """
    # jobs_dictst = list(jobs)
    results = {}
    ctx = mp.get_context("spawn")  # single, safe choice for Linux/macOS
    pid_queue: mp.queues.Queue = ctx.Queue()
    worker_pids: set[int] = set()

    def _collect_pids_nonblock() -> None:
        while True:
            try:
                pid = pid_queue.get_nowait()
                worker_pids.add(pid)
            except queue.Empty:
                break
            except Exception:
                break

    ex = None
    try:
        with ProcessPoolExecutor(
            max_workers=max_workers,
            mp_context=ctx,
            initializer=_init_worker_with_pid,
            initargs=(pid_queue,),
        ) as ex:
            fut2idx = {ex.submit(func, **kwargs): key for key, kwargs in jobs.items()}
            for fut in as_completed(fut2idx):
                _collect_pids_nonblock()
                results[fut2idx[fut]] = fut.result()
    except KeyboardInterrupt:
        # Grab any PIDs already reported
        _collect_pids_nonblock()

        # Best effort: also ask the executor for its child processes (private API)
        try:
            if ex is not None and hasattr(ex, "_processes"):
                for proc in ex._processes.values():
                    if proc is not None and proc.pid is not None:
                        worker_pids.add(proc.pid)
        except Exception:
            pass

        # Send SIGTERM, brief grace, then SIGKILL
        for pid in list(worker_pids):
            try:
                os.kill(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            except Exception:
                pass
        try:
            import time

            time.sleep(0.4)
        except Exception:
            pass
        for pid in list(worker_pids):
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            except Exception:
                pass
        raise
    return results


if __name__ == "__main__":
    pass


