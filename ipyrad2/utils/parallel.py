#!/usr/bin/env python

# POSIX-only robust ProcessPool with hard Ctrl-C and subprocess cleanup.
# Formatted with Black-style line breaks/indentation.

from typing import Any, List, Optional, Sequence, Tuple, Callable, Dict, Iterable, Iterator
import atexit
import multiprocessing as mp
import queue
import os
import signal
import subprocess as sp
from concurrent.futures import ProcessPoolExecutor, as_completed, FIRST_COMPLETED, wait
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


def _init_worker_with_pid(pid_queue: "mp.queues.Queue", log_level: str) -> None:
    """initialization function to store pids and register killer."""
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
    procs: List[sp.Popen] = []
    fout = None
    try:
        prev = None
        for i, argv in enumerate(cmds):
            is_first = (i == 0)
            is_last  = (i == len(cmds) - 1)

            # stdout target
            if is_last and outfile is not None:
                outfile.parent.mkdir(parents=True, exist_ok=True)
                fout = outfile.open("wb")
                stdout_target = fout
            else:
                stdout_target = sp.PIPE

            # stdin source
            if is_first and stdin_text is not None:
                stdin_source = sp.PIPE
            else:
                stdin_source = None if prev is None else prev.stdout

            p = safe_popen(
                argv,
                stdin=stdin_source,
                stdout=stdout_target,
                stderr=sp.PIPE,
                text=False,
            )

            # IMPORTANT: only close prev.stdout if we are NOT wiring it into the last stage
            if prev is not None and prev.stdout is not None and not is_last:
                prev.stdout.close()

            procs.append(p)
            prev = p

        # feed stdin to first stage if requested
        if stdin_text is not None:
            first = procs[0]
            if first.stdin is not None:
                data = stdin_text.encode(stdin_encoding, errors="strict")
                first.stdin.write(data)
                first.stdin.close()
                # Prevent later flush attempts on a closed handle
                first.stdin = None

        last = procs[-1]

        if outfile is not None:
            _, err = last.communicate()
            rc = last.returncode
            # drain earlier stderrs
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


# ---------- Parent-side wrapper (always uses 'spawn'; no option exposed) ----------


def run_with_pool(
    jobs: Dict[Any, Tuple[Callable[[Any], Any], Dict[Any, Any]]],
    log_level: str,
    max_workers: int | None = None,
) -> List[Any]:
    """Distribute jobs in parallel and collect results.

    Submit jobs as a {key: (func, **kwargs), ...}.

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
        # start process pool
        with ProcessPoolExecutor(
            max_workers=max_workers,
            mp_context=ctx,
            initializer=_init_worker_with_pid,
            initargs=(pid_queue, log_level),
        ) as ex:

            # submit jobs to the pool
            futures_to_jnames = {}
            for jname, job in jobs.items():
                func, kwargs = job
                future = ex.submit(func, **kwargs)
                futures_to_jnames[future] = jname

            # wait for jobs to finish while storing pids in case we need to kill
            for fut in as_completed(futures_to_jnames):
                _collect_pids_nonblock()
                results[futures_to_jnames[fut]] = fut.result()

    # stop jobs on interrupt
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



def run_with_pool_iter(
    jobs_iter: Iterable[Tuple[Any, Tuple[Callable[..., Any], Dict[str, Any]]]],
    log_level: str,
    max_workers: Optional[int] = None,
    max_inflight: Optional[int] = None,
) -> Iterator[Tuple[Any, Any]]:
    """
    Yield (key, result) as jobs complete, with robust SIGINT cleanup.
    Accepts a streaming iterator of jobs instead of a pre-built dict.

    Parameters
    ----------
    jobs_iter
        An iterator yielding (key, (func, kwargs_dict)) for each job.
        Example item: ("locus0001", (worker_build_consensus, {"locus_id": "locus0001", ...}))
    log_level
        Passed to your worker initializer (loguru setup).
    max_workers
        Max worker processes; passed to ProcessPoolExecutor.
    max_inflight
        Max number of in-flight futures. Defaults to 2 * max_workers (or max_workers if None).

    Yields
    ------
    (key, result)
        Each job's key and the value returned by the worker function.
        Exceptions raised by a worker will propagate when that job completes.

    Behavior
    --------
    - Submits up to `max_inflight` jobs at a time to bound memory.
    - Yields results as soon as futures complete.
    - Ctrl-C: sends SIGTERM then SIGKILL to workers; workers kill their child procs.
    """
    # --- setup worker context and pid collection (same pattern as run_with_pool) ---
    ctx = mp.get_context("spawn")
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

    if max_inflight is None:
        if max_workers is None:
            max_inflight = os.cpu_count() or 2
        else:
            max_inflight = max(1, 2 * max_workers)

    # --- submission/consumption state ---
    job_it = iter(jobs_iter)
    futures_to_key: Dict[Any, Any] = {}
    inflight: set = set()
    ex = None

    try:
        with ProcessPoolExecutor(
            max_workers=max_workers,
            mp_context=ctx,
            initializer=_init_worker_with_pid,
            initargs=(pid_queue, log_level),
        ) as ex:
            # Pre-fill the window
            while len(inflight) < max_inflight:
                try:
                    key, (func, kwargs) = next(job_it)
                except StopIteration:
                    break
                fut = ex.submit(func, **kwargs)
                futures_to_key[fut] = key
                inflight.add(fut)

            # As soon as any completes, yield it and submit the next job
            while inflight:
                _collect_pids_nonblock()
                done, _pending = wait(inflight, return_when=FIRST_COMPLETED)
                for fut in done:
                    inflight.remove(fut)
                    key = futures_to_key.pop(fut)
                    # propagate exceptions here
                    result = fut.result()
                    yield key, result

                    # backfill one new job per completion (bounded window)
                    try:
                        key2, (func2, kwargs2) = next(job_it)
                        fut2 = ex.submit(func2, **kwargs2)
                        futures_to_key[fut2] = key2
                        inflight.add(fut2)
                    except StopIteration:
                        pass

    except KeyboardInterrupt:
        # Collect any worker pids already registered
        _collect_pids_nonblock()
        # Best effort: ask executor for its child processes (private API)
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



if __name__ == "__main__":
    pass


