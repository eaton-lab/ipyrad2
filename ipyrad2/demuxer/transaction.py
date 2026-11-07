"""Recoverable publication of privately staged demux outputs.

Renames are atomic per file, not across the entire output set. A durable journal
lets a subsequent invocation roll back an interrupted publication.
"""

from __future__ import annotations

from contextlib import contextmanager
import json
import os
from pathlib import Path
import shutil
import signal
import uuid

from loguru import logger
from ..utils.exceptions import IPyradError

RUN_PREFIX = ".ipyrad_demux_run-"
SCRATCH_PREFIX = ".ipyrad_demux_scratch-"


def pending_runs(outdir):
    return sorted(outdir.glob(RUN_PREFIX + "*"))


def _sync_dir(path):
    if os.name == "posix":
        fd = os.open(path, os.O_RDONLY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)


def _save(root, state):
    path = root / "manifest.json"
    with (root / "manifest.tmp").open("w") as handle:
        json.dump(state, handle)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(root / "manifest.tmp", path)
    _sync_dir(root)


@contextmanager
def output_lock(outdir):
    """Serialize runs sharing an output directory; the lock file is persistent."""
    path = outdir / ".ipyrad_demux.lock"
    with path.open("a+b") as handle:
        if os.name == "posix":
            import fcntl

            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise IPyradError(f"Another demux run is using {outdir}") from exc
        else:
            import msvcrt

            handle.seek(0)
            handle.write(b"0")
            handle.flush()
            handle.seek(0)
            try:
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            except OSError as exc:
                raise IPyradError(f"Another demux run is using {outdir}") from exc
        try:
            yield
        finally:
            if os.name == "posix":
                fcntl.flock(handle, fcntl.LOCK_UN)
            else:
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)


class DemuxTransaction:
    def __init__(self, outdir, scratch_root=None):
        self.outdir = outdir
        token = uuid.uuid4().hex
        self.root = outdir / (RUN_PREFIX + token)
        self.root.mkdir(mode=0o700)
        self.state = {"phase": "processing", "entries": [], "scratch": None}
        self.outputs = self.root / "outputs"
        self.scratch = None
        scratch_created = False
        try:
            self.outputs.mkdir()
            (self.root / "backup").mkdir()
            _save(self.root, self.state)
            if scratch_root is not None:
                self.scratch = scratch_root / (SCRATCH_PREFIX + token)
                self.state["scratch"] = str(self.scratch)
                _save(self.root, self.state)
                self.scratch.mkdir(mode=0o700)
                scratch_created = True
                (self.scratch / ".owner").write_text(token)
        except BaseException:
            try:
                if scratch_created:
                    shutil.rmtree(self.scratch)
                shutil.rmtree(self.root)
            except OSError as exc:
                logger.warning(
                    "Could not clean failed demux staging setup at {}: {}",
                    self.root,
                    exc,
                )
            raise

    @classmethod
    def load(cls, root):
        self = cls.__new__(cls)
        self.root, self.outdir = root, root.parent
        self.outputs = root / "outputs"
        try:
            self.state = json.loads((root / "manifest.json").read_text())
            if self.state["phase"] not in ("processing", "publishing", "committed"):
                raise ValueError("unknown phase")
            for entry in self.state["entries"]:
                for value in (entry["target"], entry["source"]):
                    if value is not None and (
                        not value or Path(value).name != value or value in (".", "..")
                    ):
                        raise ValueError("unsafe artifact name")
            self.scratch = (
                Path(self.state["scratch"]) if self.state["scratch"] else None
            )
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise IPyradError(
                f"Cannot recover demux manifest {root}: {exc}; inspect it manually."
            ) from exc
        return self

    def publish(self, replacements, removals=()):
        for source, _ in replacements:
            with source.open("rb") as handle:
                os.fsync(handle.fileno())
        _sync_dir(self.outputs)
        targets = {target.name: source.name for source, target in replacements}
        for target in removals:
            targets.setdefault(target.name, None)
        self.state["entries"] = [
            dict(target=target, source=source)
            for target, source in sorted(targets.items())
        ]
        self.state["phase"] = "publishing"
        _save(self.root, self.state)
        for entry in self.state["entries"]:
            target = self.outdir / entry["target"]
            if target.exists():
                os.replace(target, self.root / "backup" / entry["target"])
                _sync_dir(self.outdir)
                _sync_dir(self.root / "backup")
            if entry["source"] is not None:
                os.replace(self.outputs / entry["source"], target)
                _sync_dir(self.outputs)
                _sync_dir(self.outdir)
        self.state["phase"] = "committed"
        try:
            _save(self.root, self.state)
        except BaseException:
            self.state["phase"] = "publishing"
            raise

    def rollback(self):
        if self.state["phase"] != "publishing":
            return
        for entry in reversed(self.state["entries"]):
            target = self.outdir / entry["target"]
            backup = self.root / "backup" / entry["target"]
            source = (
                self.outputs / entry["source"] if entry["source"] is not None else None
            )
            # Moving new outputs back makes rollback itself restartable.
            if source is not None and not source.exists() and target.exists():
                os.replace(target, source)
                _sync_dir(self.outputs)
                _sync_dir(self.outdir)
            if backup.exists():
                os.replace(backup, target)
                _sync_dir(self.root / "backup")
                _sync_dir(self.outdir)
        self.state["phase"] = "processing"
        _save(self.root, self.state)

    def cleanup(self):
        if self.scratch is not None and self.scratch.exists():
            token = self.root.name.removeprefix(RUN_PREFIX)
            marker = self.scratch / ".owner"
            if (
                self.scratch.name != SCRATCH_PREFIX + token
                or self.scratch.is_symlink()
                or not marker.is_file()
                or marker.read_text() != token
            ):
                raise IPyradError(
                    f"Refusing to remove scratch without matching ownership marker: {self.scratch}"
                )
            shutil.rmtree(self.scratch)
        shutil.rmtree(self.root)
        _sync_dir(self.outdir)

    def abort(self):
        previous = signal.getsignal(signal.SIGINT)
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        try:
            self.rollback()
            self.cleanup()
        except Exception as exc:
            logger.warning(
                "demux cleanup failed; recovery directory retained at {}: {}",
                self.root,
                exc,
            )
        finally:
            signal.signal(signal.SIGINT, previous)


def recover_runs(outdir):
    for root in pending_runs(outdir):
        if root.is_symlink():
            raise IPyradError(f"Refusing to recover symlink {root}")
        transaction = DemuxTransaction.load(root)
        transaction.rollback()
        transaction.cleanup()
