#!/usr/bin/env python

"""Shared process-pool helpers with hard Ctrl-C and subprocess cleanup."""

from __future__ import annotations

from concurrent.futures import FIRST_COMPLETED
from concurrent.futures import Future
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import wait
from dataclasses import dataclass
import multiprocessing as mp
import os
import queue
import signal
import time
from typing import Any, Callable, Dict, Iterable, Iterator, Tuple

from loguru import logger

from ..progress import ProgressBar
from .pipeline import _init_worker_with_pid


POOL_WAIT_POLL_SECONDS = 0.25


def _validate_positive_int(name: str, value: int | None) -> int | None:
    """Validate optional positive integer parameters."""
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _func_name(func: Callable[..., Any]) -> str:
    """Return the most useful callable name for logging."""
    return getattr(func, "__qualname__", getattr(func, "__name__", repr(func)))


def _summarize_value(value: Any) -> Any:
    """Summarize large values for exception messages and logs."""
    if isinstance(value, (bytes, bytearray)):
        return f"<{type(value).__name__} len={len(value)}>"
    if isinstance(value, (list, tuple, set, frozenset)):
        return f"<{type(value).__name__} len={len(value)}>"
    if isinstance(value, dict):
        return f"<dict len={len(value)}>"
    text = repr(value)
    if len(text) > 120:
        return text[:117] + "..."
    return value


def _summarize_kwargs(kwargs: Dict[str, Any]) -> Dict[str, Any]:
    """Return a log-safe summary of kwargs."""
    return {key: _summarize_value(value) for key, value in kwargs.items()}


class ParallelJobError(RuntimeError):
    """Raised when a job submitted to the shared process pool fails."""

    def __init__(
        self,
        key: Any,
        func_name: str,
        kwargs: Dict[str, Any],
        original: BaseException,
    ) -> None:
        self.key = key
        self.func_name = func_name
        self.kwargs = kwargs
        self.kwargs_repr = _summarize_kwargs(kwargs)
        self.original = original
        super().__init__(
            f"job failed (key={key!r}, func={func_name}, kwargs={self.kwargs_repr!r})"
        )


def _log_parallel_job_error(exc: ParallelJobError) -> None:
    """Emit a concise structured log line for a failed pool job."""
    logger.error(
        f"job failed key={exc.key!r} func={exc.func_name} kwargs={exc.kwargs_repr!r}: "
        f"{type(exc.original).__name__}: {exc.original}"
    )


def _collect_nonblock(
    pid_queue: mp.queues.Queue[int],
    child_pg_queue: mp.queues.Queue[int],
    worker_pids: set[int],
    child_pgids: set[int],
) -> None:
    """Drain PID and PGID queues without blocking."""
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


def _close_queue(qobj: Any) -> None:
    """Close a multiprocessing queue-like object if possible."""
    try:
        qobj.close()
    except Exception:
        pass


def _pid_is_alive(pid: int) -> bool:
    """Return True if a PID still exists."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except Exception:
        return False
    return True


def _kill_pid(pid: int, sig: int) -> None:
    """Send a signal to a PID if it still exists."""
    try:
        os.kill(pid, sig)
    except ProcessLookupError:
        pass
    except Exception:
        pass


def _kill_pgid(pgid: int, sig: int) -> None:
    """Send a signal to a process group if it still exists."""
    try:
        os.killpg(pgid, sig)
    except ProcessLookupError:
        pass
    except Exception:
        pass


def _executor_processes(executor: ProcessPoolExecutor | None) -> list[Any]:
    """Return worker process objects tracked by ProcessPoolExecutor."""
    if executor is None or not hasattr(executor, "_processes"):
        return []
    try:
        return [proc for proc in executor._processes.values() if proc is not None]
    except Exception:
        return []


def _staged_shutdown(
    executor: ProcessPoolExecutor | None,
    worker_pids: set[int],
    child_pgids: set[int],
    child_grace: float = 0.05,
    worker_grace: float = 0.05,
) -> None:
    """Kill tool subprocess groups first, then workers, with brief grace."""
    procs = _executor_processes(executor)
    for proc in procs:
        pid = getattr(proc, "pid", None)
        if pid:
            worker_pids.add(pid)

    def _alive_worker_pids() -> set[int]:
        alive: set[int] = set()
        for proc in procs:
            try:
                proc.join(timeout=0)
            except Exception:
                pass
            pid = getattr(proc, "pid", None)
            if not pid:
                continue
            try:
                if proc.is_alive():
                    alive.add(pid)
            except Exception:
                if pid and _pid_is_alive(pid):
                    alive.add(pid)
        for pid in list(worker_pids):
            if _pid_is_alive(pid):
                alive.add(pid)
        return alive

    for pgid in list(child_pgids):
        _kill_pgid(pgid, signal.SIGTERM)

    deadline = time.monotonic() + child_grace
    alive_pids = _alive_worker_pids()
    while alive_pids and time.monotonic() < deadline:
        time.sleep(0.01)
        alive_pids = _alive_worker_pids()

    for pid in list(alive_pids):
        _kill_pid(pid, signal.SIGTERM)

    deadline = time.monotonic() + worker_grace
    alive_pids = _alive_worker_pids()
    while alive_pids and time.monotonic() < deadline:
        time.sleep(0.01)
        alive_pids = _alive_worker_pids()

    for pgid in list(child_pgids):
        _kill_pgid(pgid, signal.SIGKILL)
    for pid in list(_alive_worker_pids()):
        _kill_pid(pid, signal.SIGKILL)

    for proc in procs:
        try:
            proc.join(timeout=0.05)
        except Exception:
            pass


def _hard_shutdown(
    executor: ProcessPoolExecutor | None,
    worker_pids: set[int],
    child_pgids: set[int],
) -> None:
    """Kill tool subprocess groups and workers immediately."""
    procs = _executor_processes(executor)
    for proc in procs:
        pid = getattr(proc, "pid", None)
        if pid:
            worker_pids.add(pid)

    for pgid in list(child_pgids):
        _kill_pgid(pgid, signal.SIGTERM)
    for pid in list(worker_pids):
        _kill_pid(pid, signal.SIGTERM)

    time.sleep(0.02)

    for pgid in list(child_pgids):
        _kill_pgid(pgid, signal.SIGKILL)
    for pid in list(worker_pids):
        _kill_pid(pid, signal.SIGKILL)


@dataclass(frozen=True)
class _FutureInfo:
    """Metadata stored for each submitted future."""

    key: Any
    func: Callable[..., Any]
    kwargs: Dict[str, Any]


class _ManagedProcessPool:
    """Internal spawn-based process pool with shared cleanup logic."""

    def __init__(self, log_level: str, max_workers: int | None = None) -> None:
        self.log_level = log_level
        self.max_workers = _validate_positive_int("max_workers", max_workers)
        self.ctx = mp.get_context("spawn")
        self.pid_queue = self.ctx.SimpleQueue()
        self.child_pg_queue = self.ctx.SimpleQueue()
        self.worker_pids: set[int] = set()
        self.child_pgids: set[int] = set()
        self.executor = ProcessPoolExecutor(
            max_workers=self.max_workers,
            mp_context=self.ctx,
            initializer=_init_worker_with_pid,
            initargs=(self.pid_queue, self.child_pg_queue, log_level),
        )
        self._aborted = False
        self._executor_shutdown = False
        self._queues_closed = False

    def collect(self) -> None:
        """Drain worker PID and subprocess PGID queues."""
        _collect_nonblock(
            self.pid_queue,
            self.child_pg_queue,
            self.worker_pids,
            self.child_pgids,
        )

    def submit(self, func: Callable[..., Any], kwargs: Dict[str, Any]) -> Future[Any]:
        """Submit a keyword-argument job to the underlying executor."""
        return self.executor.submit(func, **kwargs)

    def abort(self, inflight: Iterable[Future[Any]] = (), fast: bool = False) -> None:
        """Cancel pending futures, terminate workers, and stop the executor."""
        if self._aborted:
            return
        self._aborted = True

        for fut in tuple(inflight):
            try:
                fut.cancel()
            except Exception:
                pass

        self.collect()
        if fast:
            _hard_shutdown(self.executor, self.worker_pids, self.child_pgids)
        else:
            _staged_shutdown(self.executor, self.worker_pids, self.child_pgids)

    def close(self, wait: bool = True) -> None:
        """Close the executor and IPC queues."""
        self.collect()
        self._shutdown_executor(wait=wait, cancel_futures=self._aborted)
        self._close_queues()

    def _shutdown_executor(self, wait: bool, cancel_futures: bool) -> None:
        if self._executor_shutdown:
            return
        try:
            self.executor.shutdown(wait=wait, cancel_futures=cancel_futures)
        except Exception:
            pass
        finally:
            self._executor_shutdown = True

    def _close_queues(self) -> None:
        if self._queues_closed:
            return
        _close_queue(self.pid_queue)
        _close_queue(self.child_pg_queue)
        self._queues_closed = True

    def iter_results(
        self,
        jobs_iter: Iterable[Tuple[Any, Tuple[Callable[..., Any], Dict[str, Any]]]],
        max_inflight: int,
    ) -> Iterator[Tuple[Any, Any]]:
        """Yield `(key, result)` pairs as jobs finish, aborting on first error."""
        _validate_positive_int("max_inflight", max_inflight)

        inflight: dict[Future[Any], _FutureInfo] = {}
        job_it = iter(jobs_iter)

        def _submit_next() -> bool:
            try:
                key, (func, kwargs) = next(job_it)
            except StopIteration:
                return False

            info = _FutureInfo(key=key, func=func, kwargs=kwargs)
            try:
                future = self.submit(func, kwargs)
            except BaseException as exc:
                self.abort(inflight)
                raise ParallelJobError(key, _func_name(func), kwargs, exc) from exc

            inflight[future] = info
            return True

        while len(inflight) < max_inflight and _submit_next():
            pass

        while inflight:
            done: set[Future[Any]] = set()
            while not done:
                self.collect()
                done, _ = wait(
                    tuple(inflight),
                    timeout=POOL_WAIT_POLL_SECONDS,
                    return_when=FIRST_COMPLETED,
                )

            for fut in done:
                info = inflight.pop(fut)
                try:
                    result = fut.result()
                except BaseException as exc:
                    self.abort(inflight)
                    raise ParallelJobError(
                        info.key,
                        _func_name(info.func),
                        info.kwargs,
                        exc,
                    ) from exc

                yield info.key, result

                while len(inflight) < max_inflight and _submit_next():
                    pass


def _default_max_workers() -> int:
    """Return the default process count used for inflight heuristics."""
    return os.cpu_count() or 2


def _prime_jobs_iter(
    jobs_iter: Iterable[Tuple[Any, Tuple[Callable[..., Any], Dict[str, Any]]]],
) -> Iterator[Tuple[Any, Tuple[Callable[..., Any], Dict[str, Any]]]] | None:
    """Return an iterator with its first item restored, or None if empty."""
    job_it = iter(jobs_iter)
    try:
        first = next(job_it)
    except StopIteration:
        return None

    def _yield_first() -> Iterator[Tuple[Any, Tuple[Callable[..., Any], Dict[str, Any]]]]:
        yield first
        yield from job_it

    return _yield_first()


def run_with_pool(
    jobs: Dict[Any, Tuple[Callable[..., Any], Dict[str, Any]]],
    log_level: str,
    max_workers: int | None = None,
    max_inflight: int | None = None,
    msg: str = "Processing",
) -> Dict[Any, Any]:
    """Run jobs in parallel with bounded submission; fail fast on first error."""
    _validate_positive_int("max_workers", max_workers)
    if not jobs:
        return {}

    if max_inflight is None:
        max_inflight = max(1, max_workers or _default_max_workers())
    _validate_positive_int("max_inflight", max_inflight)

    results: Dict[Any, Any] = {}
    prog = ProgressBar(len(jobs), None, f"{msg} - total jobs: {len(jobs)}")
    prog.update()

    pool = _ManagedProcessPool(log_level=log_level, max_workers=max_workers)
    close_wait = True
    try:
        for key, result in pool.iter_results(jobs.items(), max_inflight=max_inflight):
            results[key] = result
            prog.finished += 1
            prog.update()
    except ParallelJobError as exc:
        prog.close()
        _log_parallel_job_error(exc)
        raise
    except KeyboardInterrupt:
        prog.close()
        logger.warning("interrupted by user. Cleaning up.")
        pool.abort(fast=True)
        close_wait = False
        raise SystemExit(130)
    finally:
        prog.close()
        pool.close(wait=close_wait)

    return results


def run_with_pool_iter(
    jobs_iter: Iterable[Tuple[Any, Tuple[Callable[..., Any], Dict[str, Any]]]],
    log_level: str,
    max_workers: int | None = None,
    max_inflight: int | None = None,
    msg: str | None = None,
    njobs: int | None = None,
    progress_increment: Callable[[Any, Any], int] | None = None,
) -> Iterator[Tuple[Any, Any]]:
    """Yield `(key, result)` as jobs complete; fail fast on first error."""
    _validate_positive_int("max_workers", max_workers)
    if (msg is None) != (njobs is None):
        raise ValueError("msg and njobs must be provided together for progress reporting")
    if njobs is not None:
        _validate_positive_int("njobs", njobs)
    if progress_increment is None:
        def progress_increment(_key, _result) -> int:
            return 1

    primed_iter = _prime_jobs_iter(jobs_iter)
    if primed_iter is None:
        return

    if max_inflight is None:
        max_inflight = max(1, 2 * (max_workers or _default_max_workers()))
    _validate_positive_int("max_inflight", max_inflight)

    prog = None
    if msg is not None and njobs is not None:
        prog = ProgressBar(njobs, None, f"{msg} - total jobs: {njobs}")
        prog.update()

    pool = _ManagedProcessPool(log_level=log_level, max_workers=max_workers)
    close_wait = True
    try:
        for key, result in pool.iter_results(primed_iter, max_inflight=max_inflight):
            if prog is not None:
                increment = progress_increment(key, result)
                if isinstance(increment, bool) or not isinstance(increment, int) or increment < 0:
                    raise ValueError("progress_increment must return a non-negative integer")
                prog.finished = min(njobs, prog.finished + increment)
                prog.update()
            yield key, result
    except ParallelJobError as exc:
        if prog is not None:
            prog.close()
        _log_parallel_job_error(exc)
        raise
    except KeyboardInterrupt:
        if prog is not None:
            prog.close()
        logger.warning("interrupted by user. Cleaning up.")
        pool.abort(fast=True)
        close_wait = False
        raise SystemExit(130)
    finally:
        if prog is not None:
            prog.close()
        pool.close(wait=close_wait)
