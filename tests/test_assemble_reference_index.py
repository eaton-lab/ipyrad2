"""Native index reuse and failure/publication checks for shared references."""

import os
from pathlib import Path
import subprocess

import pytest

from ipyrad2.assembler import reference_index as indexing
from ipyrad2.assembler.beds import get_reference_sort_order
from ipyrad2.utils.exceptions import IPyradError


def reference_file(tmp_path, content=b'>chr2 description\nACGTAC\n>chr1\nTGCATG\n'):
    path = tmp_path / 'reference.fa'
    path.write_bytes(content)
    return path


def native_index(reference):
    subprocess.run([indexing.BIN_SAM, 'faidx', str(reference)], check=True, capture_output=True)
    return Path(str(reference) + '.fai')


def snapshot(path):
    return path.read_bytes(), path.stat().st_mtime_ns, path.stat().st_ino


@pytest.mark.parametrize('content', [
    b'>chr2 description\nACGTAC\n>chr1\nTGCATG\n',
    b'>chr2\r\nACGT\r\nAC\r\n>chr1\r\nTGCA\r\nTG\r\n',
    b'>chr2\nACGTAC\n>chr1\nTGCATG',
])
def test_existing_native_index_and_sort_order_unchanged(monkeypatch, tmp_path, content):
    reference = reference_file(tmp_path, content)
    fai = native_index(reference)
    before = snapshot(reference), snapshot(fai)

    def no_reindex(*args, **kwargs):
        pytest.fail('An existing current index must not be regenerated')

    monkeypatch.setattr(indexing, 'run_pipeline', no_reindex)
    for _ in range(2):
        result = get_reference_sort_order(reference, tmp_path)
        assert result.read_text() == 'chr2\t6\nchr1\t6\n'
    assert (snapshot(reference), snapshot(fai)) == before


def test_missing_index_created_once_and_native_fetch_matches(monkeypatch, tmp_path):
    reference = reference_file(tmp_path)
    before = snapshot(reference)
    original = indexing.run_pipeline
    calls = []

    def track(cmds):
        calls.append(cmds)
        return original(cmds)

    monkeypatch.setattr(indexing, 'run_pipeline', track)
    fai = indexing.ensure_reference_index(reference)
    first = snapshot(fai)
    assert indexing.ensure_reference_index(reference) == fai
    assert len(calls) == 1
    assert snapshot(fai) == first and snapshot(reference) == before
    fetched = subprocess.run([indexing.BIN_SAM, 'faidx', str(reference), 'chr2:2-5'],
                             check=True, capture_output=True, text=True).stdout
    assert fetched == '>chr2:2-5\nCGTA\n'
    assert not list(tmp_path.glob('.ipyrad-faidx-*'))


def test_older_but_correct_index_is_verified_without_rewriting(monkeypatch, tmp_path):
    reference = reference_file(tmp_path)
    fai = native_index(reference)
    now = fai.stat().st_mtime_ns + 10_000_000
    os.utime(reference, ns=(now, now))
    before = snapshot(reference), snapshot(fai)
    workspace = tmp_path / 'job'
    workspace.mkdir()
    original = indexing.run_pipeline
    calls = []

    def track(cmds):
        candidate = Path(cmds[0][cmds[0].index('--fai-idx') + 1])
        assert workspace in candidate.parents
        assert candidate != fai
        calls.append(cmds)
        return original(cmds)

    monkeypatch.setattr(indexing, 'run_pipeline', track)
    assert indexing.ensure_reference_index(reference, workspace=workspace) == fai
    assert len(calls) == 1
    assert (snapshot(reference), snapshot(fai)) == before
    assert list(workspace.iterdir()) == []


def test_stale_plausible_index_rejected_and_preserved(tmp_path):
    reference = reference_file(tmp_path, b'>chr1\nACGT\n')
    fai = native_index(reference)
    reference.write_bytes(b'>chr1\nACGTACGT\n')
    newer = fai.stat().st_mtime_ns + 10_000_000
    os.utime(reference, ns=(newer, newer))
    before = snapshot(reference), snapshot(fai)
    with pytest.raises(IPyradError, match='does not match the current reference'):
        indexing.ensure_reference_index(reference, workspace=tmp_path)
    assert (snapshot(reference), snapshot(fai)) == before
    assert not list(tmp_path.glob('.ipyrad-faidx-*'))


@pytest.mark.parametrize('contents', [
    '', 'broken\n', 'chr1\tNaN\t6\t4\t5\n', 'chr1\t100\t6\t4\t5\n',
    'chr1\t4\t6\t0\t1\n', 'chr1\t4\t6\t4\t9\n',
    'chr1\t4\t6\t4\t5\nchr1\t1\t7\t1\t2\n',
])
def test_invalid_index_never_overwritten(monkeypatch, tmp_path, contents):
    reference = reference_file(tmp_path, b'>chr1\nACGT\n')
    fai = Path(str(reference) + '.fai')
    fai.write_text(contents)
    before = snapshot(fai)

    def no_reindex(*args, **kwargs):
        pytest.fail('A malformed existing index must not be silently replaced')

    monkeypatch.setattr(indexing, 'run_pipeline', no_reindex)
    with pytest.raises(IPyradError, match='Invalid FASTA index'):
        indexing.ensure_reference_index(reference)
    assert snapshot(fai) == before


def test_native_failure_does_not_publish_partial_index(monkeypatch, tmp_path):
    reference = reference_file(tmp_path)
    before = snapshot(reference)

    def fail(cmds):
        Path(cmds[0][cmds[0].index('--fai-idx')+1]).write_text('partial')
        raise IPyradError('native indexing failed')

    monkeypatch.setattr(indexing, 'run_pipeline', fail)
    with pytest.raises(IPyradError, match='native indexing failed'):
        indexing.ensure_reference_index(reference)
    assert not Path(str(reference)+'.fai').exists()
    assert not list(tmp_path.glob('.ipyrad-faidx-*'))
    assert snapshot(reference) == before


@pytest.mark.parametrize('conflict', [False, True])
def test_concurrent_publication_never_clobbers_existing_index(monkeypatch, tmp_path, conflict):
    reference = reference_file(tmp_path)
    link = os.link
    published = []

    def another_run_wins(candidate, fai):
        fai.write_bytes(b'other index\n' if conflict else candidate.read_bytes())
        published.append(snapshot(fai))
        link(candidate, fai)

    monkeypatch.setattr(indexing.os, 'link', another_run_wins)
    if conflict:
        with pytest.raises(IPyradError, match='conflicting FASTA index'):
            indexing.ensure_reference_index(reference)
    else:
        indexing.ensure_reference_index(reference)
    assert snapshot(Path(str(reference)+'.fai')) == published[0]
    assert not list(tmp_path.glob('.ipyrad-faidx-*'))


def test_reference_change_during_indexing_prevents_publication(monkeypatch, tmp_path):
    reference = reference_file(tmp_path)
    original = indexing.run_pipeline

    def changed(cmds):
        result = original(cmds)
        stat = reference.stat()
        os.utime(reference, ns=(stat.st_atime_ns, stat.st_mtime_ns+10_000_000))
        return result

    monkeypatch.setattr(indexing, 'run_pipeline', changed)
    with pytest.raises(IPyradError, match='FASTA changed while generating'):
        indexing.ensure_reference_index(reference)
    assert not Path(str(reference)+'.fai').exists()


def test_existing_index_can_be_reused_through_symlinks(monkeypatch, tmp_path):
    source = reference_file(tmp_path)
    fai = native_index(source)
    linked = tmp_path/'linked.fa'
    linked.symlink_to(source)
    Path(str(linked)+'.fai').symlink_to(fai)
    before = snapshot(source), snapshot(fai)
    monkeypatch.setattr(indexing, 'run_pipeline', lambda *a, **k: pytest.fail('reindex'))
    assert indexing.ensure_reference_index(linked) == Path(str(linked)+'.fai')
    assert (snapshot(source), snapshot(fai)) == before
