#!/usr/bin/env python

"""Pipeline subprocess helpers with shared cleanup and streaming support."""

from __future__ import annotations

import atexit
from contextlib import ExitStack
import multiprocessing as mp
import math
import os
import select
import selectors
import signal
import subprocess as sp
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterator, Optional, Sequence, Tuple

from ..logger import setup_loguru_worker


_CHILD_PROCS: set[sp.Popen] = set()
_CHILD_PROCS_LOCK = threading.RLock()
_CHILD_PG_QUEUE: mp.queues.Queue[int] | None = None
_WORKER_EXITING = False


class PipelineTimeoutError(RuntimeError):
    """Raised when a subprocess pipeline exceeds a requested timeout."""


@dataclass
class _LaunchedPipeline:
    """Tracked subprocesses and stderr tempfiles for one pipeline."""

    procs: list[sp.Popen]
    stderr_tmpfiles: list[tempfile.NamedTemporaryFile]


def safe_popen(argv: Sequence[str], **kwargs) -> sp.Popen:
    """Start Popen in its own process group and register it for cleanup."""
    # `start_new_session=True` gives us the same kill-by-process-group behavior
    # as `os.setsid`, but avoids `preexec_fn`, which is brittle under threads.
    kwargs.pop("preexec_fn", None)
    kwargs.setdefault("start_new_session", True)
    proc = sp.Popen(argv, **kwargs)
    with _CHILD_PROCS_LOCK:
        # A prior canceled child might have left kernel I/O only after its
        # cleanup deadline. Reap and forget it before registering new work.
        _CHILD_PROCS.difference_update(p for p in tuple(_CHILD_PROCS) if p.poll() is not None)
        _CHILD_PROCS.add(proc)
    try:
        if _CHILD_PG_QUEUE is not None:
            _CHILD_PG_QUEUE.put(os.getpgid(proc.pid))
    except Exception:
        pass
    return proc


def _kill_all_children(sig: int = signal.SIGTERM) -> None:
    """Kill all registered child processes and clear the registry."""
    with _CHILD_PROCS_LOCK:
        procs = list(_CHILD_PROCS)
        _CHILD_PROCS.clear()

    alive: list[sp.Popen] = []
    for proc in procs:
        try:
            if proc.poll() is None:
                alive.append(proc)
        except Exception:
            pass

    for proc in alive:
        try:
            os.killpg(proc.pid, sig)
        except Exception:
            pass

    if sig != signal.SIGKILL and alive:
        time.sleep(0.005)
        for proc in alive:
            try:
                if proc.poll() is None:
                    os.killpg(proc.pid, signal.SIGKILL)
            except Exception:
                pass


def _worker_signal_handler(signum, _frame) -> None:
    """Worker signal handler: kill children then exit immediately."""
    global _WORKER_EXITING
    if _WORKER_EXITING:
        os._exit(128 + signum)
    _WORKER_EXITING = True
    try:
        _kill_all_children()
    finally:
        try:
            signal.signal(signal.SIGINT, signal.SIG_DFL)
            signal.signal(signal.SIGTERM, signal.SIG_DFL)
        except Exception:
            pass
        os._exit(128 + signum)


def _init_worker_with_pid(
    pid_queue: mp.queues.Queue[int],
    child_pg_queue: mp.queues.Queue[int],
    log_level: str,
) -> None:
    """Initializer: register handlers, send worker PID, and set logger."""
    global _CHILD_PG_QUEUE, _WORKER_EXITING
    _CHILD_PG_QUEUE = child_pg_queue
    _WORKER_EXITING = False
    with _CHILD_PROCS_LOCK:
        _CHILD_PROCS.clear()
    pid_queue.put(os.getpid())
    signal.signal(signal.SIGINT, _worker_signal_handler)
    signal.signal(signal.SIGTERM, _worker_signal_handler)
    atexit.register(_kill_all_children)
    setup_loguru_worker(log_level)


def _is_sigpipe_rc(rc: int, sigpipe: int) -> bool:
    """Return True when a return code reflects SIGPIPE."""
    return rc == -sigpipe or rc == 128 + sigpipe


def _terminate_proc(proc: sp.Popen, sig: int) -> None:
    """Terminate one process group if it is still alive."""
    try:
        if proc.poll() is None:
            os.killpg(proc.pid, sig)
    except Exception:
        pass


def _collect_pipeline_stderr(
    stderr_tmpfiles: list[tempfile.NamedTemporaryFile],
) -> bytes:
    """Collect and return combined stderr from temporary files."""
    err_all = bytearray()
    for tmpf in stderr_tmpfiles:
        try:
            tmpf.flush()
            tmpf.seek(0)
            err_all.extend(tmpf.read())
        except Exception:
            pass
    return bytes(err_all)


def _cleanup_pipeline_tempfiles(
    stderr_tmpfiles: list[tempfile.NamedTemporaryFile],
) -> None:
    """Close and unlink pipeline stderr temporary files."""
    for tmpf in stderr_tmpfiles:
        try:
            name = getattr(tmpf, "name", None)
            tmpf.close()
            if name and os.path.exists(name):
                os.unlink(name)
        except Exception:
            pass


def _launch_pipeline(
    cmds: list[Sequence[str]],
    stdin_text: Optional[str],
    stdin_encoding: str,
    last_stdout_target,
    last_stderr_pipe: bool,
    env: dict[str, str] | None = None,
    cwd: str | Path | None = None,
) -> _LaunchedPipeline:
    """Launch a pipeline and optionally feed text into its first stage."""
    procs: list[sp.Popen] = []
    stderr_tmpfiles: list[tempfile.NamedTemporaryFile] = []
    prev: sp.Popen | None = None

    try:
        with ExitStack() as stack:
            input_file = None
            if stdin_text is not None:
                # A file lets the child read input while the parent drains output.
                input_file = stack.enter_context(tempfile.TemporaryFile())
                input_file.write(stdin_text.encode(stdin_encoding, errors="strict"))
                input_file.seek(0)

            for idx, argv in enumerate(cmds):
                is_first = idx == 0
                is_last = idx == len(cmds) - 1

                stdout_target = last_stdout_target if is_last else sp.PIPE
                stdin_source = input_file if is_first else prev.stdout

                if is_last and last_stderr_pipe:
                    stderr_target = sp.PIPE
                else:
                    tmpf = tempfile.NamedTemporaryFile(prefix="pipe-stderr-", delete=False)
                    stderr_target = tmpf
                    stderr_tmpfiles.append(tmpf)

                child_context = {}
                if env is not None:
                    child_context["env"] = env
                if cwd is not None:
                    child_context["cwd"] = cwd
                proc = safe_popen(
                    argv,
                    stdin=stdin_source,
                    stdout=stdout_target,
                    stderr=stderr_target,
                    text=False,
                    **child_context,
                )

                procs.append(proc)
                if prev is not None and prev.stdout is not None:
                    prev.stdout.close()
                prev = proc
    except BaseException as exc:
        # A later executable can fail to start after earlier children launched.
        _terminate_pipeline(procs)
        try:
            _wait_for_pipeline(procs, timeout_s=PIPELINE_FAILURE_REAP_TIMEOUT_S)
        except Exception as cleanup:
            if hasattr(exc, "add_note"):
                exc.add_note(str(cleanup))
        finally:
            _cleanup_pipeline_processes(procs, stderr_tmpfiles)
        raise

    return _LaunchedPipeline(procs=procs, stderr_tmpfiles=stderr_tmpfiles)


def _terminate_pipeline(procs: list[sp.Popen]) -> None:
    """Terminate one pipeline quickly and forcefully."""
    for proc in procs:
        _terminate_proc(proc, signal.SIGTERM)
    time.sleep(0.02)
    for proc in procs:
        _terminate_proc(proc, signal.SIGKILL)


# Keep a shared short teardown budget for actual failures/deadlines; a
# deliberately truncated successful scan can tolerate delayed filesystem I/O.
PIPELINE_FAILURE_REAP_TIMEOUT_S = 1.2
STREAM_EARLY_CLOSE_REAP_TIMEOUT_S = 30.0


def _wait_for_pipeline(
    procs: list[sp.Popen], *, timeout_s: float = PIPELINE_FAILURE_REAP_TIMEOUT_S,
) -> list[int]:
    """Reap a whole pipeline against one deadline, not one timeout per child.

    Callers send termination signals before waiting on canceled pipelines.
    SIGKILL cannot immediately release a process blocked in kernel I/O.
    """
    deadline = time.monotonic() + timeout_s
    while True:
        rcs = [proc.poll() for proc in procs]
        if all(rc is not None for rc in rcs):
            return [int(rc) for rc in rcs]
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            pending = [str(proc.args) for proc, rc in zip(procs, rcs) if rc is None]
            raise PipelineTimeoutError(
                f"pipeline cleanup exceeded {timeout_s:g}s after termination; "
                f"unreaped children: {', '.join(pending)}"
            )
        time.sleep(min(0.01, remaining))


def _cleanup_pipeline_processes(
    procs: list[sp.Popen],
    stderr_tmpfiles: list[tempfile.NamedTemporaryFile],
) -> None:
    """Close pipes and temporary files; unregister only reaped children."""
    for proc in procs:
        with _CHILD_PROCS_LOCK:
            # Retain a still-live child for later reaping/worker shutdown.
            if proc.poll() is not None:
                _CHILD_PROCS.discard(proc)
        try:
            if proc.stdin is not None:
                proc.stdin.close()
        except Exception:
            pass
        try:
            if proc.stdout is not None:
                proc.stdout.close()
        except Exception:
            pass
        try:
            if proc.stderr is not None:
                proc.stderr.close()
        except Exception:
            pass

    _cleanup_pipeline_tempfiles(stderr_tmpfiles)


def _bad_pipeline_rcs(
    rcs: list[int],
    ncmds: int,
    ignore_sigpipe: bool,
    sigpipe: int,
) -> list[tuple[int, int]]:
    """Return the non-allowed nonzero return codes from one pipeline."""
    bad: list[tuple[int, int]] = []
    for idx, rc in enumerate(rcs):
        if rc == 0:
            continue
        if ignore_sigpipe and idx < ncmds - 1 and _is_sigpipe_rc(rc, sigpipe):
            continue
        bad.append((idx, rc))
    return bad


def _pipeline_error(
    cmds: list[Sequence[str]],
    err_all: bytes,
    primary_fail: tuple[int, int] | None,
    bad: list[tuple[int, int]],
) -> RuntimeError:
    """Build the canonical pipeline failure error."""
    if primary_fail is not None:
        fail_idx, fail_rc = primary_fail
    else:
        fail_idx, fail_rc = bad[0]
    stderr_text = err_all.decode(errors="replace")
    detail = _summarize_disk_space_error(stderr_text) or stderr_text
    return RuntimeError(
        f"pipeline failed (rc={fail_rc}): {cmds[fail_idx]}\n"
        f"{detail}"
    )


def _pipeline_timeout_error(
    cmds: list[Sequence[str]],
    timeout_s: float,
    elapsed_s: float,
) -> PipelineTimeoutError:
    """Build the canonical pipeline timeout error."""
    rendered = " | ".join(" ".join(map(str, cmd)) for cmd in cmds)
    return PipelineTimeoutError(
        f"pipeline timed out after {elapsed_s:.1f}s (timeout={timeout_s:.1f}s): {rendered}"
    )


def _summarize_disk_space_error(stderr_text: str) -> str | None:
    """Return a concise message for disk-full pipeline failures."""
    disk_markers = (
        "No space left on device",
        "Disk quota exceeded",
        "ENOSPC",
    )
    context_markers = (
        "File write failed",
        "failed writing to",
    )
    lines = [line.rstrip() for line in stderr_text.splitlines() if line.strip()]
    if not any(any(marker in line for marker in disk_markers) for line in lines):
        return None

    selected: list[str] = []
    for line in lines:
        if any(marker in line for marker in disk_markers + context_markers):
            if line not in selected:
                selected.append(line)
            if len(selected) >= 6:
                break

    if not selected:
        selected = [
            line for line in lines if any(marker in line for marker in disk_markers)
        ][:3]

    detail = "\n".join(selected) if selected else "disk space error"
    return f"disk space error:\n{detail}"


def _make_atomic_outfile(outfile: Path) -> tuple[BinaryIO, Path]:
    """Create a temporary sibling file for atomic final output replacement."""
    outfile.parent.mkdir(parents=True, exist_ok=True)
    fout = tempfile.NamedTemporaryFile(
        mode="wb",
        prefix=f".{outfile.name}.",
        suffix=".tmp",
        dir=outfile.parent,
        delete=False,
    )
    return fout, Path(fout.name)


STREAM_FINALIZE_TIMEOUT_S = 60.0


class _PipelineLineStream:
    """Incremental stdout reader for a subprocess pipeline."""

    def __init__(
        self,
        cmds: list[Sequence[str]],
        stdin_text: Optional[str],
        stdin_encoding: str,
        ignore_sigpipe: bool,
        line_encoding: str,
    ) -> None:
        import signal as _signal

        self.cmds = cmds
        self.ignore_sigpipe = ignore_sigpipe
        self.line_encoding = line_encoding
        self.sigpipe = _signal.SIGPIPE
        self.buffer = bytearray()
        self.primary_fail: tuple[int, int] | None = None
        self.saw_eof = False
        self.closed = False
        self.finalized = False

        launched = _launch_pipeline(
            cmds=cmds,
            stdin_text=stdin_text,
            stdin_encoding=stdin_encoding,
            last_stdout_target=sp.PIPE,
            last_stderr_pipe=False,
        )
        self.procs = launched.procs
        self.stderr_tmpfiles = launched.stderr_tmpfiles
        self.last = self.procs[-1]
        if self.last.stdout is None:
            raise RuntimeError("pipeline was started without a stdout pipe")
        self.stdout_fd = self.last.stdout.fileno()

    def __iter__(self) -> "_PipelineLineStream":
        return self

    def __enter__(self) -> "_PipelineLineStream":
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        try:
            self.close(fast=exc_type is not None)
        except Exception as cleanup:
            if exc is None:
                raise
            if hasattr(exc, "add_note"):
                exc.add_note(str(cleanup))
        return False

    def __del__(self) -> None:
        try:
            self.close(fast=True)
        except Exception:
            pass

    def __next__(self) -> str:
        if self.closed:
            raise StopIteration

        while True:
            line = self._pop_line()
            if line is not None:
                return line

            if self.primary_fail is not None:
                self._terminate_live()
                if self.buffer:
                    return self._drain_buffer()
                self._finalize()
                raise StopIteration

            if self.saw_eof:
                if self.buffer:
                    return self._drain_buffer()
                self._finalize()
                raise StopIteration

            ready, _, _ = select.select([self.stdout_fd], [], [], 0.05)
            if ready:
                chunk = os.read(self.stdout_fd, 65536)
                if not chunk:
                    self.saw_eof = True
                else:
                    self.buffer.extend(chunk)

            self._poll_failures()

    def close(self, *, fast: bool = False) -> None:
        """Cancel a deliberately truncated stream and release resources."""
        if self.closed:
            return
        self.closed = True
        if self.finalized:
            return
        # If a producer has already failed before cancellation, do not
        # suppress that genuine failure as an intentional truncation.
        if self.primary_fail is None:
            self._poll_failures()
        already_failed = self.primary_fail is not None
        # Early callers intentionally stop after a bounded sample of lines.
        # Release the pipe before signaling, so producers blocked on output
        # can exit promptly rather than waiting for the pipe to drain.
        if self.last.stdout is not None:
            self.last.stdout.close()
        self._terminate_live()
        budget = (PIPELINE_FAILURE_REAP_TIMEOUT_S if fast or already_failed
                  else STREAM_EARLY_CLOSE_REAP_TIMEOUT_S)
        self._finalize(suppress_errors=fast or not already_failed, reap_timeout_s=budget)

    def _pop_line(self) -> str | None:
        newline = self.buffer.find(b"\n")
        if newline < 0:
            return None
        line = bytes(self.buffer[:newline])
        del self.buffer[: newline + 1]
        return line.decode(self.line_encoding, errors="replace")

    def _drain_buffer(self) -> str:
        line = bytes(self.buffer).decode(self.line_encoding, errors="replace")
        self.buffer.clear()
        return line

    def _poll_failures(self) -> None:
        for idx, proc in enumerate(self.procs):
            rc = proc.poll()
            if rc is None or rc == 0:
                continue
            if self.ignore_sigpipe and idx < len(self.procs) - 1 and _is_sigpipe_rc(rc, self.sigpipe):
                continue
            self.primary_fail = (idx, rc)
            break

    def _terminate_live(self) -> None:
        _terminate_pipeline(self.procs)

    def _finalize(self, suppress_errors: bool = False, *,
                  reap_timeout_s: float = PIPELINE_FAILURE_REAP_TIMEOUT_S) -> None:
        if self.finalized:
            return
        self.finalized = True
        start = time.monotonic()
        timed_out = False
        try:
            # EOF is not process exit. Allow bounded normal finalization rather
            # than killing a successful command one second after its last byte.
            if not suppress_errors and self.primary_fail is None:
                while any(proc.poll() is None for proc in self.procs):
                    self._poll_failures()
                    if self.primary_fail is not None:
                        self._terminate_live()
                        break
                    if time.monotonic() - start >= STREAM_FINALIZE_TIMEOUT_S:
                        timed_out = True
                        self._terminate_live()
                        break
                    time.sleep(0.01)
            try:
                rcs = _wait_for_pipeline(self.procs, timeout_s=reap_timeout_s)
            except PipelineTimeoutError as cleanup:
                if timed_out:
                    raise _pipeline_timeout_error(
                        self.cmds, STREAM_FINALIZE_TIMEOUT_S, time.monotonic() - start,
                    ) from cleanup
                if self.primary_fail is not None and not suppress_errors:
                    raise _pipeline_error(
                        self.cmds, _collect_pipeline_stderr(self.stderr_tmpfiles),
                        self.primary_fail, [self.primary_fail],
                    ) from cleanup
                raise
            if timed_out:
                raise _pipeline_timeout_error(
                    self.cmds, STREAM_FINALIZE_TIMEOUT_S, time.monotonic() - start,
                )
            err_all = _collect_pipeline_stderr(self.stderr_tmpfiles)
            bad = _bad_pipeline_rcs(
                rcs=rcs,
                ncmds=len(self.cmds),
                ignore_sigpipe=self.ignore_sigpipe,
                sigpipe=self.sigpipe,
            )
            if bad and not suppress_errors:
                raise _pipeline_error(
                    cmds=self.cmds,
                    err_all=err_all,
                    primary_fail=self.primary_fail,
                    bad=bad,
                )
        finally:
            _cleanup_pipeline_processes(self.procs, self.stderr_tmpfiles)


def run_pipeline(
    cmds: list[Sequence[str]],
    outfile: Optional[Path] = None,
    stdin_text: Optional[str] = None,
    stdin_encoding: str = "utf-8",
    ignore_sigpipe: bool = True,
    timeout_s: float | None = None,
    env: dict[str, str] | None = None,
    cwd: str | Path | None = None,
) -> Tuple[int, bytes, bytes]:
    """Run a shell-like pipeline with pipefail and an optional operation deadline.

    Once the deadline expires, all children are canceled and given one short,
    shared cleanup budget. The longer intentional stream-close grace does not
    apply to this API.
    """
    import signal as _signal

    if not cmds:
        raise ValueError("run_pipeline(): cmds must be a non-empty list of argv sequences")
    if timeout_s is not None:
        if (isinstance(timeout_s, bool) or not isinstance(timeout_s, (int, float))
                or not math.isfinite(timeout_s) or timeout_s <= 0):
            raise ValueError("run_pipeline(): timeout_s must be a positive number")

    sigpipe = _signal.SIGPIPE

    procs: list[sp.Popen] = []
    stderr_tmpfiles: list[tempfile.NamedTemporaryFile] = []
    fout: BinaryIO | None = None
    outfile_tmp: Path | None = None

    stdout_buffer = bytearray()
    stderr_buffer = bytearray()
    timed_out = False
    elapsed_s = 0.0
    primary_fail: tuple[int, int] | None = None
    cleanup_attempted = False

    def drain_ready(selector, timeout: float) -> bool:
        ready = selector.select(timeout)
        for key, _ in ready:
            block = os.read(key.fd, 65536)
            if block:
                key.data.extend(block)
            else:
                selector.unregister(key.fileobj)
        return bool(ready)

    try:
        if outfile is not None:
            fout, outfile_tmp = _make_atomic_outfile(outfile)

        launched = _launch_pipeline(
            cmds=cmds,
            stdin_text=stdin_text,
            stdin_encoding=stdin_encoding,
            last_stdout_target=fout if fout is not None else sp.PIPE,
            last_stderr_pipe=True,
            env=env,
            cwd=cwd,
        )
        procs = launched.procs
        stderr_tmpfiles = launched.stderr_tmpfiles

        last = procs[-1]
        start_time = time.monotonic()
        deadline = start_time + timeout_s if timeout_s is not None else None

        with selectors.DefaultSelector() as selector:
            if last.stdout is not None:
                selector.register(last.stdout, selectors.EVENT_READ, stdout_buffer)
            if last.stderr is not None:
                selector.register(last.stderr, selectors.EVENT_READ, stderr_buffer)

            # EOF is not proof of upstream exit. Every process and stdout/stderr
            # pipe remains subject to the operation deadline until completed.
            while True:
                rcs = [proc.poll() for proc in procs]
                bad = [
                    (i, rc) for i, rc in enumerate(rcs)
                    if rc not in (None, 0) and not (
                        ignore_sigpipe and i < len(procs) - 1 and _is_sigpipe_rc(rc, sigpipe)
                    )
                ]
                if bad:
                    primary_fail = bad[0]
                    break
                if all(rc is not None for rc in rcs) and not selector.get_map():
                    break
                remaining = deadline - time.monotonic() if deadline is not None else None
                if remaining is not None and remaining <= 0:
                    timed_out = True
                    elapsed_s = time.monotonic() - start_time
                    break
                wait = min(0.01, remaining) if remaining is not None else 0.01
                if selector.get_map():
                    drain_ready(selector, wait)
                else:
                    time.sleep(wait)

            if primary_fail is not None or timed_out:
                _terminate_pipeline(procs)

            # On failure, collect already-buffered stderr without waiting for
            # an uninterruptibly blocked child or a slow descendant.
            cleanup_attempted = True
            try:
                rcs = _wait_for_pipeline(procs, timeout_s=PIPELINE_FAILURE_REAP_TIMEOUT_S)
            except PipelineTimeoutError as cleanup:
                # Do not replace the operation timeout or primary exit error.
                if timed_out:
                    raise _pipeline_timeout_error(cmds, float(timeout_s), elapsed_s) from cleanup
                if primary_fail is not None:
                    drain_ready(selector, 0)
                    err = bytes(stderr_buffer) + _collect_pipeline_stderr(stderr_tmpfiles)
                    raise _pipeline_error(cmds, err, primary_fail, [primary_fail]) from cleanup
                raise
            drain_deadline = time.monotonic() + 0.05
            while selector.get_map() and time.monotonic() < drain_deadline:
                # Children have exited; drain buffered output, but never wait
                # for an unrelated descendant that inherited a pipe.
                if not drain_ready(selector, 0):
                    break

        if fout is not None:
            try:
                fout.flush()
            finally:
                fout.close()
                fout = None

        if timed_out:
            raise _pipeline_timeout_error(cmds, float(timeout_s), elapsed_s or (time.monotonic() - start_time))

        err_all = stderr_buffer
        err_all.extend(_collect_pipeline_stderr(stderr_tmpfiles))
        bad = _bad_pipeline_rcs(
            rcs=rcs,
            ncmds=len(cmds),
            ignore_sigpipe=ignore_sigpipe,
            sigpipe=sigpipe,
        )
        if bad:
            raise _pipeline_error(
                cmds=cmds,
                err_all=bytes(err_all),
                primary_fail=primary_fail,
                bad=bad,
            )

        if outfile is not None and outfile_tmp is not None:
            os.replace(outfile_tmp, outfile)
            outfile_tmp = None

        return 0, bytes(stdout_buffer), bytes(err_all)

    except BaseException as exc:
        if not cleanup_attempted:
            _terminate_pipeline(procs)
            try:
                _wait_for_pipeline(procs, timeout_s=PIPELINE_FAILURE_REAP_TIMEOUT_S)
            except Exception as cleanup:
                if hasattr(exc, "add_note"):
                    exc.add_note(str(cleanup))
        raise

    finally:
        try:
            if fout is not None:
                fout.close()
        except Exception:
            pass

        if outfile_tmp is not None:
            try:
                outfile_tmp.unlink()
            except FileNotFoundError:
                pass

        _cleanup_pipeline_processes(procs, stderr_tmpfiles)


def stream_pipeline_lines(
    cmds: list[Sequence[str]],
    stdin_text: Optional[str] = None,
    stdin_encoding: str = "utf-8",
    ignore_sigpipe: bool = True,
    line_encoding: str = "utf-8",
) -> Iterator[str]:
    """Yield decoded stdout lines incrementally from a pipeline."""
    if not cmds:
        raise ValueError("stream_pipeline_lines(): cmds must be a non-empty list of argv sequences")
    return _PipelineLineStream(
        cmds=cmds,
        stdin_text=stdin_text,
        stdin_encoding=stdin_encoding,
        ignore_sigpipe=ignore_sigpipe,
        line_encoding=line_encoding,
    )
