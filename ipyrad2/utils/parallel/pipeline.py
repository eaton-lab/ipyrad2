#!/usr/bin/env python

"""Pipeline subprocess helpers with shared cleanup and streaming support."""

from __future__ import annotations

import atexit
import multiprocessing as mp
import os
import select
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
) -> _LaunchedPipeline:
    """Launch a pipeline and optionally feed text into its first stage."""
    procs: list[sp.Popen] = []
    stderr_tmpfiles: list[tempfile.NamedTemporaryFile] = []
    prev: sp.Popen | None = None

    for idx, argv in enumerate(cmds):
        is_first = idx == 0
        is_last = idx == len(cmds) - 1

        stdout_target = last_stdout_target if is_last else sp.PIPE
        stdin_source = (
            sp.PIPE
            if is_first and stdin_text is not None
            else None if prev is None else prev.stdout
        )

        if is_last and last_stderr_pipe:
            stderr_target = sp.PIPE
        else:
            tmpf = tempfile.NamedTemporaryFile(prefix="pipe-stderr-", delete=False)
            stderr_target = tmpf
            stderr_tmpfiles.append(tmpf)

        proc = safe_popen(
            argv,
            stdin=stdin_source,
            stdout=stdout_target,
            stderr=stderr_target,
            text=False,
        )

        if prev is not None and prev.stdout is not None:
            prev.stdout.close()

        procs.append(proc)
        prev = proc

    if stdin_text is not None and procs and procs[0].stdin is not None:
        data = stdin_text.encode(stdin_encoding, errors="strict")
        procs[0].stdin.write(data)
        procs[0].stdin.close()
        procs[0].stdin = None

    return _LaunchedPipeline(procs=procs, stderr_tmpfiles=stderr_tmpfiles)


def _terminate_pipeline(procs: list[sp.Popen]) -> None:
    """Terminate one pipeline quickly and forcefully."""
    for proc in procs:
        _terminate_proc(proc, signal.SIGTERM)
    time.sleep(0.02)
    for proc in procs:
        _terminate_proc(proc, signal.SIGKILL)


def _wait_for_pipeline(procs: list[sp.Popen]) -> list[int]:
    """Wait for pipeline processes to exit and return their return codes."""
    rcs: list[int] = []
    for proc in procs:
        try:
            rc = proc.poll()
            if rc is None:
                rc = proc.wait(timeout=1.0)
        except sp.TimeoutExpired:
            _terminate_proc(proc, signal.SIGKILL)
            rc = proc.wait(timeout=0.2)
        except Exception:
            rc = proc.returncode if proc.returncode is not None else 1
        rcs.append(int(rc) if rc is not None else 0)
    return rcs


def _cleanup_pipeline_processes(
    procs: list[sp.Popen],
    stderr_tmpfiles: list[tempfile.NamedTemporaryFile],
) -> None:
    """Close process pipes, unregister children, and clean stderr tempfiles."""
    for proc in procs:
        with _CHILD_PROCS_LOCK:
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
        self.close()
        return False

    def __del__(self) -> None:
        try:
            self.close()
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

    def close(self) -> None:
        """Terminate any live subprocesses and release resources."""
        if self.closed:
            return
        self.closed = True
        self._terminate_live()
        self._finalize(suppress_errors=True)

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

    def _finalize(self, suppress_errors: bool = False) -> None:
        if self.finalized:
            return
        self.finalized = True

        rcs = _wait_for_pipeline(self.procs)
        err_all = _collect_pipeline_stderr(self.stderr_tmpfiles)
        bad = _bad_pipeline_rcs(
            rcs=rcs,
            ncmds=len(self.cmds),
            ignore_sigpipe=self.ignore_sigpipe,
            sigpipe=self.sigpipe,
        )
        _cleanup_pipeline_processes(self.procs, self.stderr_tmpfiles)
        if bad and not suppress_errors:
            raise _pipeline_error(
                cmds=self.cmds,
                err_all=err_all,
                primary_fail=self.primary_fail,
                bad=bad,
            )


def run_pipeline(
    cmds: list[Sequence[str]],
    outfile: Optional[Path] = None,
    stdin_text: Optional[str] = None,
    stdin_encoding: str = "utf-8",
    ignore_sigpipe: bool = True,
    timeout_s: float | None = None,
) -> Tuple[int, bytes, bytes]:
    """Run a shell-like pipeline of cmds with pipefail + fail-fast semantics."""
    import signal as _signal

    if not cmds:
        raise ValueError("run_pipeline(): cmds must be a non-empty list of argv sequences")
    if timeout_s is not None:
        if isinstance(timeout_s, bool) or not isinstance(timeout_s, (int, float)) or timeout_s <= 0:
            raise ValueError("run_pipeline(): timeout_s must be a positive number")

    sigpipe = _signal.SIGPIPE

    procs: list[sp.Popen] = []
    stderr_tmpfiles: list[tempfile.NamedTemporaryFile] = []
    fout: BinaryIO | None = None
    outfile_tmp: Path | None = None

    last_out = b""
    last_err = b""
    comm_exc: BaseException | None = None
    comm_done = threading.Event()
    timed_out = False
    elapsed_s = 0.0

    def _communicate_last(last: sp.Popen) -> None:
        nonlocal last_out, last_err, comm_exc
        try:
            if outfile is not None:
                _out, last_err = last.communicate()
                last_out = b""
            else:
                last_out, last_err = last.communicate()
        except BaseException as exc:
            comm_exc = exc
        finally:
            comm_done.set()

    try:
        if outfile is not None:
            fout, outfile_tmp = _make_atomic_outfile(outfile)

        launched = _launch_pipeline(
            cmds=cmds,
            stdin_text=stdin_text,
            stdin_encoding=stdin_encoding,
            last_stdout_target=fout if fout is not None else sp.PIPE,
            last_stderr_pipe=True,
        )
        procs = launched.procs
        stderr_tmpfiles = launched.stderr_tmpfiles

        last = procs[-1]
        start_time = time.monotonic()
        comm_thread = threading.Thread(target=_communicate_last, args=(last,), daemon=True)
        comm_thread.start()

        primary_fail: tuple[int, int] | None = None

        while not comm_done.is_set():
            for idx, proc in enumerate(procs):
                rc = proc.poll()
                if rc is None or rc == 0:
                    continue
                if ignore_sigpipe and idx < len(procs) - 1 and _is_sigpipe_rc(rc, sigpipe):
                    continue
                primary_fail = (idx, rc)
                break
            if primary_fail is not None:
                break
            if timeout_s is not None:
                elapsed_s = time.monotonic() - start_time
                if elapsed_s >= timeout_s:
                    timed_out = True
                    break
            time.sleep(0.005)

        if primary_fail is None:
            for idx, proc in enumerate(procs):
                rc = proc.poll()
                if rc is None or rc == 0:
                    continue
                if ignore_sigpipe and idx < len(procs) - 1 and _is_sigpipe_rc(rc, sigpipe):
                    continue
                primary_fail = (idx, rc)
                break

        if primary_fail is not None or timed_out:
            _terminate_pipeline(procs)

        comm_done.wait(timeout=1.0)
        comm_thread.join(timeout=1.0)

        if fout is not None:
            try:
                fout.flush()
            finally:
                fout.close()
                fout = None

        rcs = _wait_for_pipeline(procs)

        if timed_out:
            raise _pipeline_timeout_error(cmds, float(timeout_s), elapsed_s or (time.monotonic() - start_time))

        if comm_exc is not None:
            raise comm_exc

        err_all = bytearray()
        if last_err:
            err_all.extend(last_err)
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

        return 0, last_out, bytes(err_all)

    except Exception:
        _terminate_pipeline(procs)
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
