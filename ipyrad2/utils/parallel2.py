#!/usr/bin/env python
from __future__ import annotations

# POSIX-only robust ProcessPool with hard Ctrl-C and subprocess cleanup.

from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Callable
import atexit
import multiprocessing as mp
import os
import queue
import signal
import subprocess as sp
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, as_completed, wait
from pathlib import Path

from .logger import setup_loguru_worker


# ---------- Worker-side (signal-safe; kills subprocess groups) ----------

_CHILD_PROCS: set[sp.Popen] = set()
_CHILD_PG_QUEUE: mp.queues.Queue[int] | None = None  # parent learns child PGIDs


def safe_popen(argv: Sequence[str], **kwargs) -> sp.Popen:
    """Start Popen in its own process group and register it for cleanup."""
    kwargs.setdefault("preexec_fn", os.setsid)  # start new PGID
    p = sp.Popen(argv, **kwargs)
    _CHILD_PROCS.add(p)
    # Inform parent of this subprocess PGID so parent can nuke it on Ctrl-C.
    try:
        if _CHILD_PG_QUEUE is not None:
            _CHILD_PG_QUEUE.put(os.getpgid(p.pid))
    except Exception:
        pass
    return p


def _kill_all_children(sig: int = signal.SIGTERM) -> None:
    """Kill all registered child processes; escalate to SIGKILL."""
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
    """Worker signal handler: kill children then exit."""
    _kill_all_children()
    raise SystemExit(128 + signum)


def _init_worker_with_pid(
    pid_queue: mp.queues.Queue[int],
    child_pg_queue: mp.queues.Queue[int],
    log_level: str,
) -> None:
    """Initializer: register handlers, send worker PID, and set logger."""
    global _CHILD_PG_QUEUE
    _CHILD_PG_QUEUE = child_pg_queue
    pid_queue.put(os.getpid())
    signal.signal(signal.SIGINT, _worker_signal_handler)
    signal.signal(signal.SIGTERM, _worker_signal_handler)
    atexit.register(_kill_all_children)
    setup_loguru_worker(log_level)


# ---------- Pipeline runner (supports optional outfile) ----------

def run_pipeline(
    cmds: List[Sequence[str]],
    outfile: Optional[Path] = None,
    stdin_text: Optional[str] = None,
    stdin_encoding: str = "utf-8",
) -> Tuple[int, bytes, bytes]:
    """Run a shell-like pipeline of cmds; kill safely on exceptions."""
    procs: List[sp.Popen] = []
    fout = None
    try:
        prev = None
        for i, argv in enumerate(cmds):
            is_first = i == 0
            is_last = i == len(cmds) - 1

            if is_last and outfile is not None:
                outfile.parent.mkdir(parents=True, exist_ok=True)
                fout = outfile.open("wb")
                stdout_target = fout
            else:
                stdout_target = sp.PIPE

            stdin_source = sp.PIPE if (is_first and stdin_text is not None) else (None if prev is None else prev.stdout)

            p = safe_popen(
                argv,
                stdin=stdin_source,
                stdout=stdout_target,
                stderr=sp.PIPE,
                text=False,
            )

            if prev is not None and prev.stdout is not None and not is_last:
                prev.stdout.close()

            procs.append(p)
            prev = p

        if stdin_text is not None and procs and procs[0].stdin is not None:
            data = stdin_text.encode(stdin_encoding, errors="strict")
            procs[0].stdin.write(data)
            procs[0].stdin.close()
            procs[0].stdin = None

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
                try:
                    fout.flush()
                finally:
                    fout.close()
                    fout = None
            if rc != 0:
                raise RuntimeError(
                    f"pipeline failed (rc={rc}): {cmds[-1]}\n{err.decode(errors='replace')}"
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
                f"pipeline failed (rc={rc}): {cmds[-1]}\n{err.decode(errors='replace')}"
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
            try:
                _CHILD_PROCS.discard(p)
            except Exception:
                pass


# ---------- Parent-side helpers ----------


def _collect_nonblock(
    pid_queue: mp.queues.Queue[int],
    child_pg_queue: mp.queues.Queue[int],
    worker_pids: set[int],
    child_pgids: set[int],
) -> None:
    """Drain PID/PGID queues without blocking."""
    while True:
        try:
            worker_pids.add(pid_queue.get_nowait())
        except queue.Empty:
            break
        except Exception:
            break
    while True:
        try:
            child_pgids.add(child_pg_queue.get_nowait())
        except queue.Empty:
            break
        except Exception:
            break


def _hard_shutdown(
    ex: Optional[ProcessPoolExecutor],
    worker_pids: set[int],
    child_pgids: set[int],
) -> None:
    """Kill tool subprocess groups first, then workers, then cancel futures."""
    # Try to pick up any still-unknown worker PIDs from executor internals
    try:
        if ex is not None and hasattr(ex, "_processes"):
            for proc in ex._processes.values():
                if proc and proc.pid:
                    worker_pids.add(proc.pid)
    except Exception:
        pass

    # 1) Kill subprocess groups (fastp/bwa/etc.)
    for pgid in list(child_pgids):
        try:
            os.killpg(pgid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        except Exception:
            pass
    for pgid in list(child_pgids):
        try:
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except Exception:
            pass

    # 2) Kill workers
    for pid in list(worker_pids):
        try:
            os.kill(pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        except Exception:
            pass
    for pid in list(worker_pids):
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except Exception:
            pass


# ---------- Pool wrappers ----------

def run_with_pool(
    jobs: Dict[Any, Tuple[Callable[[Any], Any], Dict[str, Any]]],
    log_level: str,
    max_workers: int | None = None,
    max_inflight: int | None = None,
) -> Dict[Any, Any]:
    """Run jobs in parallel with bounded submission; return {key: result}."""
    results: Dict[Any, Any] = {}
    ctx = mp.get_context("spawn")

    # Use SimpleQueue to avoid SemLock-backed Queue semaphores.
    pid_queue = ctx.SimpleQueue()
    child_pg_queue = ctx.SimpleQueue()

    worker_pids: set[int] = set()
    child_pgids: set[int] = set()

    if max_inflight is None:
        max_inflight = max(1, (max_workers or (os.cpu_count() or 2)))

    jobs_iter = iter(jobs.items())
    futures_to_key: Dict[Any, Any] = {}
    inflight: set = set()
    ex: ProcessPoolExecutor | None = None

    try:
        with ProcessPoolExecutor(
            max_workers=max_workers,
            mp_context=ctx,
            initializer=_init_worker_with_pid,
            initargs=(pid_queue, child_pg_queue, log_level),
        ) as ex:
            # Pre-fill bounded window
            while len(inflight) < max_inflight:
                try:
                    key, (func, kwargs) = next(jobs_iter)
                except StopIteration:
                    break
                fut = ex.submit(func, **kwargs)
                futures_to_key[fut] = key
                inflight.add(fut)

            # Consume and backfill one-by-one
            while inflight:
                _collect_nonblock(pid_queue, child_pg_queue, worker_pids, child_pgids)
                done, _ = wait(inflight, return_when=FIRST_COMPLETED)
                for fut in done:
                    inflight.remove(fut)
                    key = futures_to_key.pop(fut)
                    results[key] = fut.result()

                    try:
                        key2, (func2, kwargs2) = next(jobs_iter)
                        fut2 = ex.submit(func2, **kwargs2)
                        futures_to_key[fut2] = key2
                        inflight.add(fut2)
                    except StopIteration:
                        pass

    except KeyboardInterrupt:
        _collect_nonblock(pid_queue, child_pg_queue, worker_pids, child_pgids)
        _hard_shutdown(ex, worker_pids, child_pgids)
        try:
            if ex is not None:
                ex.shutdown(wait=False, cancel_futures=True)
        except Exception:
            pass
        raise
    finally:
        # Close queues so their FDs/threads are torn down; avoids leaked semaphores.
        try:
            pid_queue.close()
        except Exception:
            pass
        try:
            child_pg_queue.close()
        except Exception:
            pass
    return results



def run_with_pool_iter(
    jobs_iter: Iterable[Tuple[Any, Tuple[Callable[..., Any], Dict[str, Any]]]],
    log_level: str,
    max_workers: Optional[int] = None,
    max_inflight: Optional[int] = None,
) -> Iterator[Tuple[Any, Any]]:
    """Yield (key, result) as jobs complete; bounds inflight work."""
    ctx = mp.get_context("spawn")
    pid_queue: mp.queues.Queue[int] = ctx.Queue()
    child_pg_queue: mp.queues.Queue[int] = ctx.Queue()
    worker_pids: set[int] = set()
    child_pgids: set[int] = set()

    if max_inflight is None:
        max_inflight = max(1, 2 * (max_workers or (os.cpu_count() or 2)))

    job_it = iter(jobs_iter)
    futures_to_key: Dict[Any, Any] = {}
    inflight: set = set()
    ex: ProcessPoolExecutor | None = None

    try:
        with ProcessPoolExecutor(
            max_workers=max_workers,
            mp_context=ctx,
            initializer=_init_worker_with_pid,
            initargs=(pid_queue, child_pg_queue, log_level),
        ) as ex:
            # Pre-fill window
            while len(inflight) < max_inflight:
                try:
                    key, (func, kwargs) = next(job_it)
                except StopIteration:
                    break
                fut = ex.submit(func, **kwargs)
                futures_to_key[fut] = key
                inflight.add(fut)

            # Consume as they finish
            while inflight:
                _collect_nonblock(pid_queue, child_pg_queue, worker_pids, child_pgids)
                done, _ = wait(inflight, return_when=FIRST_COMPLETED)
                for fut in done:
                    inflight.remove(fut)
                    key = futures_to_key.pop(fut)
                    yield key, fut.result()

                    # Backfill one job
                    try:
                        key2, (func2, kwargs2) = next(job_it)
                        fut2 = ex.submit(func2, **kwargs2)
                        futures_to_key[fut2] = key2
                        inflight.add(fut2)
                    except StopIteration:
                        pass

    except KeyboardInterrupt:
        _collect_nonblock(pid_queue, child_pg_queue, worker_pids, child_pgids)
        _hard_shutdown(ex, worker_pids, child_pgids)
        # Cancel any pending futures (best-effort)
        try:
            if ex is not None and hasattr(ex, "_futures"):
                for f in list(ex._futures):
                    f.cancel()
        except Exception:
            pass
        raise


if __name__ == "__main__":
    pass
