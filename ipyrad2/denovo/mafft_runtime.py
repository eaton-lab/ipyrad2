"""Child-local MAFFT runtime selection, independent of executable PATH lookup.

The scientific MAFFT executable and arguments are unchanged. Only the child
working directory, temporary directory, and shell-helper search path vary.
"""
from __future__ import annotations

from contextlib import contextmanager
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import threading

from loguru import logger

_NETWORK_FS = {"nfs", "nfs4", "cifs", "smbfs", "lustre", "gpfs", "ceph", "glusterfs"}


def _unescape_mount(value: str) -> str:
    return re.sub(r"\\([0-7]{3})", lambda m: chr(int(m[1], 8)), value)


@lru_cache(maxsize=256)
def filesystem_type(path: str) -> str | None:
    """Find the longest matching Linux mount, resolving executable symlinks."""
    resolved = Path(path).resolve()
    try:
        lines = Path('/proc/self/mountinfo').read_text().splitlines()
    except OSError:
        return None
    best = (-1, None)
    for line in lines:
        fields = line.split()
        try:
            mount = Path(_unescape_mount(fields[4]))
            kind = fields[fields.index('-') + 1]
            if resolved == mount or mount in resolved.parents:
                if len(str(mount)) > best[0]:
                    best = (len(str(mount)), kind)
        except (ValueError, IndexError):
            continue
    return best[1]


def _network(path: str) -> bool:
    kind = filesystem_type(path)
    return kind in _NETWORK_FS or bool(kind and kind.startswith('fuse.'))


@lru_cache(maxsize=32)
def _helper_path(inherited_path: str) -> str:
    """Use local system helpers when the environment's awk is on shared storage.

    Do not infer a defective awk version or change helpers after a timeout.
    Unknown/non-Linux filesystems retain the environment's selection.
    """
    awk = shutil.which('awk', path=inherited_path)
    if not awk or not _network(awk):
        return inherited_path
    system_path = os.pathsep.join(p for p in ('/usr/bin', '/bin') if Path(p).is_dir())
    candidate = shutil.which('awk', path=system_path)
    if not candidate or _network(candidate):
        raise RuntimeError('MAFFT awk is on shared storage and no local system awk is available')
    # Exercise the arithmetic and field parsing used by MAFFT before selecting it.
    check = subprocess.run(
        [candidate, 'BEGIN {print -1.0 * "1.53"} {print $1, $3}'],
        input='7 x 143\n', capture_output=True, text=True, timeout=5,
    )
    if check.returncode or check.stdout.split() != ['-1.53', '7', '143']:
        raise RuntimeError(f'Local MAFFT awk compatibility check failed: {candidate}')
    logger.info('MAFFT shell helpers: {} (environment awk {} is on shared storage)', system_path, awk)
    return system_path + os.pathsep + inherited_path


def _scratch_parent() -> str:
    # Respect existing user/scheduler scratch choices when they are local.
    candidates = [os.environ.get(k) for k in ('MAFFT_TMPDIR', 'SLURM_TMPDIR', 'TMPDIR')]
    candidates += [tempfile.gettempdir(), '/tmp']
    for value in dict.fromkeys(v for v in candidates if v):
        path = Path(value).expanduser().absolute()
        if path.is_dir() and not _network(str(path)) and os.access(path, os.W_OK | os.X_OK):
            return str(path)
    raise RuntimeError('No writable local MAFFT scratch directory; set MAFFT_TMPDIR to node-local scratch')


_RUNTIME_LOCK = threading.Lock()
_LOCAL_RUNTIMES: dict[tuple[str, str, str], tuple[str, str]] = {}


def _local_runtime(binary: str, env: dict[str, str], parent: str) -> str:
    """Stage a shared Conda MAFFT runtime once, without changing its contents.

    Bash reads its launcher while running, and native helpers also live on NFS.
    Moving scratch and awk alone does not remove those shared-file accesses.
    Retain the small private runtime bundle for failed-call diagnostics.
    """
    supplied = env.get('MAFFT_BINARIES')
    prefix = Path(supplied).expanduser().absolute() if supplied else Path(binary).parent.parent/'libexec/mafft'
    if not _network(binary) and not (prefix.is_dir() and _network(str(prefix))):
        return binary
    if not prefix.is_dir() or not (prefix/'version').is_file():
        raise RuntimeError('Cannot stage shared MAFFT installation: set MAFFT_BINARIES to its native program directory')
    key = (binary, str(prefix), parent)
    with _RUNTIME_LOCK:
        cached = _LOCAL_RUNTIMES.get(key)
        if cached is None:
            root = Path(tempfile.mkdtemp(prefix='ipyrad2-mafft-runtime-', dir=parent))
            launcher = root/Path(binary).name
            shutil.copy2(binary, launcher)
            native = root/'native'
            shutil.copytree(prefix, native)
            files = [(Path(binary), launcher)] + [(p, native/p.relative_to(prefix)) for p in prefix.rglob('*') if p.is_file()]
            manifest = {}
            for source, copied in files:
                digest = hashlib.sha256(source.read_bytes()).hexdigest()
                if hashlib.sha256(copied.read_bytes()).hexdigest() != digest:
                    raise RuntimeError(f'MAFFT runtime copy checksum mismatch: {source}')
                manifest[str(source)] = digest
            (root/'source-sha256.json').write_text(json.dumps(manifest, indent=2)+'\n')
            cached = (str(launcher), str(native))
            _LOCAL_RUNTIMES[key] = cached
            logger.info('MAFFT byte-identical local runtime retained at {}', root)
        env['MAFFT_BINARIES'] = cached[1]
        return cached[0]


@contextmanager
def mafft_process_context(binary: str, input_text: str):
    """Prepare one isolated invocation without mutating process-wide environment.

    Preserve failed inputs and scratch for replay. Successful private scratch
    is removed; existing user inputs and directories are never removed.
    """
    env = dict(os.environ)
    env['PATH'] = _helper_path(env.get('PATH', os.defpath))
    parent = _scratch_parent()
    work = Path(tempfile.mkdtemp(prefix='ipyrad2-mafft-', dir=parent))
    env.update(TMPDIR=str(work), MAFFT_TMPDIR=str(work), PWD=str(work))
    # Resolve before changing cwd, including caller-supplied relative paths.
    executable = str(Path(shutil.which(binary) or binary).absolute())
    try:
        executable = _local_runtime(executable, env, parent)
        yield executable, env, work
    except Exception as exc:
        try:
            (work / 'input.fa').write_text(input_text)
            (work / 'failure.json').write_text(json.dumps({
                'executable': executable, 'helper_path': env['PATH'],
                'mafft_binaries': env.get('MAFFT_BINARIES'),
                'scratch': str(work), 'error': str(exc),
            }, indent=2) + '\n')
        except OSError as diagnostic_error:
            logger.warning('Could not write MAFFT replay metadata: {}', diagnostic_error)
        logger.error('MAFFT failed; input and scratch retained at {}', work)
        raise
    else:
        shutil.rmtree(work)
