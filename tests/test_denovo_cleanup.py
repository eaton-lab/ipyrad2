"""Removed experimental interfaces fail explicitly; native output semantics remain."""
import inspect

import pytest

from ipyrad2.cli.cli_main import setup_parsers
from ipyrad2.denovo import align, denovo
from ipyrad2.denovo.runtime import OutputTransaction


@pytest.mark.parametrize('flag,value', [
    ('--read-refinement', 'none'),
    ('--read-refinement', 'audit'),
    ('--audit-max-components', '64'),
    ('--audit-max-samples', '6'),
    ('--audit-max-fragments', '25000'),
    ('--audit-scan-fragments', '250000'),
])
def test_removed_audit_cli_rejected(flag, value):
    with pytest.raises(SystemExit) as error:
        setup_parsers().parse_args(['denovo', '-d', 'sample.fastq', flag, value])
    assert error.value.code == 2


def test_removed_python_arguments_fail_before_execution():
    for name, value in [('read_refinement', 'none'), ('audit_max_components', 64),
                        ('audit_max_samples', 6), ('audit_max_fragments', 25000),
                        ('audit_scan_fragments', 250000)]:
        with pytest.raises(TypeError, match=name):
            denovo.run_denovo(**{name: value})
    for function in (align.consensus_from_aligned, align.worker_build_consensus,
                     align.iter_alignment_jobs, align.write_ordered_consensus_stream_to_file):
        assert 'min_prop' not in inspect.signature(function).parameters
        with pytest.raises(TypeError, match='min_prop'):
            function(min_prop=.9)


@pytest.mark.parametrize('function,args', [
    (align.consensus_from_aligned, ([], .9)),
    (align.worker_build_consensus, (1, [], '', .9)),
    (align.iter_alignment_jobs, (None, None, '', .9)),
    (align.write_ordered_consensus_stream_to_file, (None, None, None, '', 1, .9)),
])
def test_removed_positional_consensus_argument_cannot_become_another_setting(function, args):
    with pytest.raises(TypeError, match='positional'):
        function(*args)


def test_consensus_ignores_gaps_and_unknowns_and_preserves_ties():
    assert align.consensus_from_aligned([
        ('a', 'a-N-T'), ('b', 'c-N-G'), ('c', 'c-N-C'), ('d', '--N-A'),
    ]) == 'CA'
    assert align.consensus_from_aligned([('empty', '-NN-')]) == ''


def test_legacy_audit_preserved_in_force_backup(tmp_path):
    legacy = tmp_path/'denovo.read_audit'
    legacy.mkdir()
    (legacy/'REPORT.txt').write_text('old experimental evidence')
    names = [p.name for p in denovo._iter_denovo_outputs(tmp_path)]
    with OutputTransaction(tmp_path, names, True) as transaction:
        (transaction.attempt/'denovo_reference.fa').write_text('>new\nACGT\n')
        transaction.publish()
        assert not legacy.exists()
        archive = tmp_path/f'.denovo-previous-{transaction.id}'
        assert (archive/'denovo.read_audit/REPORT.txt').read_text() == 'old experimental evidence'
    assert not legacy.exists()


def test_legacy_archive_failure_retains_transaction_backup(tmp_path, monkeypatch):
    legacy = tmp_path/'denovo.read_audit'
    legacy.mkdir()
    (legacy/'REPORT.txt').write_text('preserve this')
    names = [p.name for p in denovo._iter_denovo_outputs(tmp_path)]
    with OutputTransaction(tmp_path, names, True) as transaction:
        original = type(tmp_path).rename
        def fail_archive(path, target):
            if path == transaction.attempt/'_previous/denovo.read_audit':
                raise OSError('archive unavailable')
            return original(path, target)
        monkeypatch.setattr(type(tmp_path), 'rename', fail_archive)
        transaction.publish()
        assert transaction.publication_status == 'PUBLISHED'
        assert (transaction.attempt/'_previous/denovo.read_audit/REPORT.txt').read_text() == 'preserve this'
