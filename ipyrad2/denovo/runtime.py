"""Stage measurements and recoverable publication for denovo runs."""
from contextlib import contextmanager
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import time
from uuid import uuid4

from ..utils.exceptions import IPyradError
from ..utils.reporting import atomic_text, jsonable, utc_now


MEMORY_DEFINITION = (
    'Sampled sum of resident bytes for the main process and its live descendants, '
    'including pool workers and native tools. Shared pages may be counted more than once; '
    'short-lived peaks between samples may be missed. Not a scheduler MaxRSS or an exact peak. '
    'Unavailable measurements are null, never zero.'
)


def process_tree_rss(pid=None):
    """Read Linux procfs without adding a required monitoring dependency."""
    root = os.getpid() if pid is None else pid
    if not Path('/proc/self/statm').exists():
        return None
    pending, seen, total = [root], set(), 0
    try:
        page_size = os.sysconf('SC_PAGE_SIZE')
        while pending:
            child = pending.pop()
            if child in seen:
                continue
            seen.add(child)
            base = Path('/proc') / str(child)
            try:
                total += int((base/'statm').read_text().split()[1]) * page_size
                # Native subprocesses may be children of non-main worker threads.
                for task in (base/'task').iterdir():
                    try:
                        pending.extend(int(x) for x in (task/'children').read_text().split())
                    except FileNotFoundError:
                        continue  # exited during sampling
            except FileNotFoundError:
                if child == root:
                    return None
                continue
        return total
    except (OSError, ValueError, IndexError):
        return None


class StageRecorder:
    def __init__(self, interval=0.25, sampler=process_tree_rss):
        self.stages = []
        self.interval = interval
        self.sampler = sampler

    @contextmanager
    def stage(self, name):
        row = dict(name=name, status='RUNNING', started_utc=utc_now(),
                   finished_utc=None, elapsed_seconds=None, sampled_peak_tree_rss_bytes=None,
                   rss_samples=0, unavailable_rss_samples=0, sample_interval_seconds=self.interval)
        self.stages.append(row)
        stop = threading.Event()

        def sample():
            try:
                value = self.sampler()
            except Exception:
                value = None
            if value is None:
                row['unavailable_rss_samples'] += 1
            else:
                row['rss_samples'] += 1
                row['sampled_peak_tree_rss_bytes'] = max(row['sampled_peak_tree_rss_bytes'] or 0, value)

        def loop():
            while not stop.wait(self.interval):
                sample()

        start = time.monotonic()
        sample()
        thread = threading.Thread(target=loop, daemon=True, name='denovo-memory')
        thread.start()
        try:
            yield row
        except BaseException:
            row['status'] = 'FAILED'
            raise
        else:
            row['status'] = 'COMPLETE'
        finally:
            stop.set()
            thread.join()
            sample()
            row.update(finished_utc=utc_now(), elapsed_seconds=time.monotonic()-start)


def tool_version(binary):
    """Capture versions only; a missing version string is not a tool failure."""
    try:
        result = subprocess.run([str(binary), '--version'], capture_output=True, text=True, timeout=5)
        lines = (result.stdout + '\n' + result.stderr).strip().splitlines()
        return dict(path=str(binary), version=lines[0] if lines else None,
                    version_exit_code=result.returncode)
    except (OSError, subprocess.TimeoutExpired) as error:
        return dict(path=str(binary), version=None, version_error=str(error))


class OutputTransaction:
    """Keep old outputs until success; roll back ordinary publication errors.

    Files are renamed atomically one at a time, with reports published last.
    This is not an atomic filesystem snapshot. A publication journal and backups
    remain on abrupt interruption; there is no automatic or blind recovery.
    """
    def __init__(self, outdir, names, force):
        self.outdir = Path(outdir)
        self.names = list(names)
        self.force = force
        self.id = uuid4().hex
        self.attempt = self.outdir / f'.denovo-attempt-{self.id}'
        self.lock = self.outdir / '.denovo.lock'
        self.publication_status = 'NOT_STARTED'

    def __enter__(self):
        self.outdir.mkdir(parents=True, exist_ok=True)
        try:
            fd = os.open(self.lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        except FileExistsError as error:
            raise IPyradError(f'Denovo output directory is locked: {self.lock}. Inspect the recorded run before removing a stale lock.') from error
        try:
            with os.fdopen(fd, 'w') as stream:
                json.dump(dict(pid=os.getpid(), run_id=self.id, attempt=str(self.attempt)), stream)
            existing = [self.outdir/name for name in self.names
                        if (self.outdir/name).exists() or (self.outdir/name).is_symlink()]
            if existing and not self.force:
                raise IPyradError(f'denovo outputs already exist in {self.outdir}. Use --force to overwrite: ' + ', '.join(p.name for p in existing))
            self.attempt.mkdir()
        except BaseException:
            self.lock.unlink(missing_ok=True)
            raise
        return self

    def publish(self):
        # Stats are completion markers and go last, after biological artifacts.
        names = sorted(self.names, key=lambda n: (n.startswith('denovo.stats'), n))
        backup = self.attempt / '_previous'
        backup.mkdir()
        journal = self.attempt / 'PUBLICATION.json'
        moved_old, moved_new = [], []
        state = dict(status='PUBLISHING', old=moved_old, new=moved_new, target=str(self.outdir))
        self.publication_status = 'PUBLISHING'
        atomic_text(journal, json.dumps(state, indent=2)+'\n')
        try:
            for name in names:
                old, new = self.outdir/name, self.attempt/name
                if old.exists() or old.is_symlink():
                    old.rename(backup/name)
                    moved_old.append(name)
                if new.exists():
                    new.rename(old)
                    moved_new.append(name)
                atomic_text(journal, json.dumps(state, indent=2)+'\n')
            state['status'] = 'PUBLISHED'
            atomic_text(journal, json.dumps(state, indent=2)+'\n')
        except BaseException as publication_error:
            try:
                for name in reversed(moved_new):
                    (self.outdir/name).rename(self.attempt/name)
                for name in reversed(moved_old):
                    (backup/name).rename(self.outdir/name)
            except BaseException as rollback_error:
                self.publication_status = state['status'] = 'RECOVERY_REQUIRED'
                state.update(publication_error=str(publication_error), rollback_error=str(rollback_error))
                try:
                    atomic_text(journal, json.dumps(state, indent=2)+'\n')
                except OSError:
                    pass
                raise IPyradError(
                    f'Denovo publication rollback failed: {rollback_error}. '
                    f'The output directory may contain mixed generations. The lock and backups '
                    f'are retained; inspect {journal} before recovery.'
                ) from publication_error
            self.publication_status = 'ROLLED_BACK'
            state['status'] = 'ROLLED_BACK'
            atomic_text(journal, json.dumps(state, indent=2)+'\n')
            raise
        self.publication_status = 'PUBLISHED'
        # Publication is complete. A cleanup error must not label it a failed run.
        try:
            # The retired experimental output has no replacement. Preserve it
            # outside the disposable transaction backup on a forced rerun.
            legacy = backup / 'denovo.read_audit'
            if legacy.exists() or legacy.is_symlink():
                archive = self.outdir / f'.denovo-previous-{self.id}'
                archive.mkdir()
                legacy.rename(archive / legacy.name)
            shutil.rmtree(self.attempt)
        except OSError:
            # In particular, a failed legacy archive keeps the attempt backup.
            pass

    def __exit__(self, *exc):
        if self.publication_status != 'RECOVERY_REQUIRED':
            self.lock.unlink(missing_ok=True)


def write_runtime_report(outputs, run, stages, body=None):
    """Augment existing scientific summaries; keep measured quantities numeric."""
    txt, js = outputs['run_stats'], outputs['run_stats_json']
    if body is None:
        data = json.loads(js.read_text()) if js.exists() else {}
        body = txt.read_text() if txt.exists() else ''
    else:
        data = {}
    if 'outputs' in data:
        data['outputs'] = {key: value.replace(run['attempt_directory'], run['output_directory'])
                           if isinstance(value, str) else value for key, value in data['outputs'].items()}
    body = body.replace(run['attempt_directory'], run['output_directory'])
    data.update(stats_schema_version='2.0', tool='denovo', run=run,
                stages=stages, memory_definition=MEMORY_DEFINITION)
    from .reporting import DEFINITIONS, render_report
    data['definitions'] = dict(DEFINITIONS)
    atomic_text(txt, render_report(data, failure_body=body if run.get('error') else ''))
    atomic_text(js, json.dumps(jsonable(data), indent=2, allow_nan=False)+'\n')
