from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import shutil

import pytest

from ipyrad2.denovo import align
from ipyrad2.denovo import mafft_runtime as runtime


@pytest.fixture(autouse=True)
def clear_policy_cache():
    runtime._helper_path.cache_clear()
    runtime.filesystem_type.cache_clear()
    yield
    runtime._helper_path.cache_clear()
    runtime.filesystem_type.cache_clear()


def test_absolute_conda_launcher_uses_local_helpers(monkeypatch, tmp_path):
    """PATH wrappers cannot intercept an absolute binary; child env must change."""
    conda = tmp_path / 'conda bin'; conda.mkdir()
    awk = conda / 'awk'; awk.write_text('#!/bin/sh\nexit 77\n'); awk.chmod(0o700)
    binary = conda / 'mafft'
    binary.write_text('#!/bin/sh\nexec awk "{print}"\n'); binary.chmod(0o700)
    scratch = tmp_path / 'local scratch'; scratch.mkdir()
    inherited = str(conda) + os.pathsep + os.environ['PATH']
    monkeypatch.setenv('PATH', inherited)
    monkeypatch.setenv('MAFFT_TMPDIR', str(scratch))
    monkeypatch.setattr(runtime, '_network', lambda p: Path(p) == awk)
    result = align.mafft_align_one([('a', 'ACGT'), ('b', 'ACGA')], str(binary))
    assert result == [('a', 'ACGT'), ('b', 'ACGA')]
    assert os.environ['PATH'] == inherited
    assert list(scratch.iterdir()) == []


def test_local_environment_helpers_unchanged(monkeypatch):
    monkeypatch.setattr(runtime, '_network', lambda p: False)
    value = '/some/environment/bin:' + os.environ['PATH']
    assert runtime._helper_path(value) == value


def test_unusable_local_awk_fails_before_alignment(monkeypatch):
    monkeypatch.setattr(runtime.shutil, 'which', lambda *a, **k: '/awk')
    monkeypatch.setattr(runtime, '_network', lambda p: True)
    with pytest.raises(RuntimeError, match='no local system awk'):
        runtime._helper_path('/shared/bin')


def test_shared_tmpdir_falls_back_to_local(monkeypatch, tmp_path):
    shared = tmp_path / 'shared'; shared.mkdir()
    local = tmp_path / 'local'; local.mkdir()
    monkeypatch.setenv('MAFFT_TMPDIR', str(shared))
    monkeypatch.setenv('SLURM_TMPDIR', str(local))
    monkeypatch.setattr(runtime, '_network', lambda p: p == str(shared))
    assert runtime._scratch_parent() == str(local)


def test_concurrent_contexts_do_not_mutate_parent(monkeypatch, tmp_path):
    monkeypatch.setenv('MAFFT_TMPDIR', str(tmp_path))
    env = dict(os.environ); original_cwd = Path.cwd()
    def one(i):
        with runtime.mafft_process_context('/bin/cat', 'data') as (_, child, work):
            assert child['MAFFT_TMPDIR'] == str(work)
            assert child['PWD'] == str(work)
            return str(work)
    with ThreadPoolExecutor(max_workers=8) as pool:
        paths = list(pool.map(one, range(24)))
    assert len(set(paths)) == 24
    assert all(not Path(p).exists() for p in paths)
    assert dict(os.environ) == env
    assert Path.cwd() == original_cwd


def test_failed_alignment_preserves_replay_input(monkeypatch, tmp_path):
    monkeypatch.setenv('MAFFT_TMPDIR', str(tmp_path))
    with pytest.raises(RuntimeError, match='deliberate'):
        with runtime.mafft_process_context('/bin/false', '>a\nACGT\n') as (_, _, work):
            raise RuntimeError('deliberate')
    assert (work/'input.fa').read_text() == '>a\nACGT\n'
    assert json.loads((work/'failure.json').read_text())['executable'] == '/bin/false'


def test_real_mafft_same_alignment_with_child_context(tmp_path, monkeypatch):
    binary = shutil.which('mafft')
    if binary is None:
        pytest.skip('MAFFT not installed')
    monkeypatch.setenv('MAFFT_TMPDIR', str(tmp_path))
    result = align.mafft_align_one([('a', 'ACGTACGT'), ('b', 'ACGTTCGT')], binary)
    assert [(name, seq.upper()) for name, seq in result] == [('a', 'ACGTACGT'), ('b', 'ACGTTCGT')]


def test_relative_executable_resolved_before_child_cwd(monkeypatch, tmp_path):
    binary = tmp_path/'mafft'; binary.write_text('#!/bin/sh\ncat\n'); binary.chmod(0o700)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('MAFFT_TMPDIR', str(tmp_path))
    result = align.mafft_align_one([('a','ACGT'), ('b','ACGA')], './mafft')
    assert result == [('a','ACGT'), ('b','ACGA')]
