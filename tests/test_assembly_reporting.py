"""Report provenance, bounded display, completion and native-evidence accounting."""
import json
from pathlib import Path
import shlex

import h5py
import numpy as np
import pytest

from ipyrad2.cli.command_log import format_stats_command
from ipyrad2.assembler.reporting import (
    record_assembly_run, activate_run_reporting, publish_assembly_report, record_inputs,
)
from ipyrad2.assembler.evidence_hdf5 import write_snp_evidence_hdf5
from ipyrad2.assembler.rad_evidence_hdf5 import write_rad_evidence_hdf5
from ipyrad2.analysis.rad_evidence import RadSiteEvidenceReader
from test_hdf5_evidence import fixture
from test_rad_site_evidence import setup


def test_command_summary_and_exact_quoted_invocation():
    argv = ['ipyrad2', 'assemble', '-d', 'space sample.bam', '-r', "ref's.fa"]
    assert shlex.split(format_stats_command(argv).removeprefix('Command: ')) == argv
    long = ['ipyrad2', 'assemble', '-d', *[f'RAD/{i}.bam' for i in range(240)],
            '-w', *[f'WGS/{i}.bam' for i in range(135)], '--evidence', 'loci']
    line = format_stats_command(long)
    assert line.startswith('Command (summary):')
    assert '<240 BAMs>' in line and '<135 BAMs>' in line
    assert len(line) <= 120 and len(line.splitlines()) == 1
    # A single huge argument cannot evade the bound; no token is cut in half.
    line = format_stats_command(['ipyrad2', 'assemble', '-r', 'x'*1000])
    assert len(line) <= 120 and line.endswith(' ...') and 'xxxx' not in line
    line = format_stats_command(['ipyrad2', 'assemble', '-n', 'new\nline\tname\r\v\f\x85\u2028\u2029'])
    assert len(line.splitlines()) == 1 and '\t' not in line and '\r' not in line
    assert line.startswith('Command (summary):')


def _minimal_stats(mode='none'):
    return {'evidence': dict(mode=mode, status='not_requested', snps=None, loci=None)}


def test_success_delayed_until_cleanup_and_embedded_provenance(tmp_path):
    argv = ['ipyrad2', 'assemble', '-d', *[f'folder /{i}.bam' for i in range(300)]]
    invocation = dict(kind='cli', argv=argv, display_argv=argv)
    groups = tmp_path/'groups.tsv'
    groups.write_text('s1\tp1\ns2\tp1\n')

    @record_assembly_run
    def run(outdir, name, evidence='none', invocation=None, keep_tmpdir=False, populations=groups):
        activate_run_reporting()
        record_inputs({'s1': tmp_path/'s1.bam'}, {'s2': tmp_path/'s2.bam'}, groups)
        with h5py.File(outdir/f'{name}.hdf5', 'w') as h:
            h.attrs['format_version'] = '3.0'
            h.create_dataset('genos', data=[1, 2, 3])
        publish_assembly_report(outdir, name, _minimal_stats(), '# Assemble Summary\nSamples 2\n')
        assert not (outdir/f'{name}.stats.json').exists()
        (outdir/'cleanup_done').touch()

    run(tmp_path, 'test', invocation=invocation)
    stats = json.loads((tmp_path/'test.stats.json').read_text())
    assert stats['run']['status'] == 'COMPLETE'
    assert stats['run']['invocation']['argv'] == argv
    assert shlex.split(stats['run']['invocation']['command']) == argv
    assert stats['command'] == shlex.join(argv)
    assert stats['run']['parameters']['evidence'] == 'none'
    assert stats['run']['working_directory'] == str(Path.cwd())
    assert stats['run']['calling']['sample_population_assignments'] == [['s1', 'p1'], ['s2', 'p1']]
    assert len(stats['run']['configuration_files']['populations']['sha256']) == 64
    assert stats['evidence']['snps'] is None and stats['storage']['evidence_dataset_bytes']['snps'] is None
    assert stats['storage']['hdf5_bytes'] == (tmp_path/'test.hdf5').stat().st_size
    assert stats['run']['finished_utc'] >= stats['run']['started_utc']
    assert stats['run']['elapsed_seconds'] > 0
    with h5py.File(tmp_path/'test.hdf5') as h:
        embedded = json.loads(h['assembly_report/summary_json'][:].tobytes())
        assert embedded['run'] == stats['run']
        assert embedded['evidence'] == stats['evidence']
        np.testing.assert_array_equal(h['genos'][:], [1, 2, 3])
    text = (tmp_path/'test.stats.txt').read_text()
    commands = [line for line in text.splitlines() if line.startswith('Command')]
    assert len(commands) == 1 and len(commands[0]) <= 120
    assert '<300 BAMs>' in commands[0] and 'not_requested' in text


def test_failure_never_publishes_complete_or_modifies_earlier_report(tmp_path):
    previous = tmp_path/'test.stats.json'
    previous.write_text('existing results')

    @record_assembly_run
    def refuse(outdir, name, evidence='none'):
        raise ValueError('overwrite refused')
    with pytest.raises(ValueError):
        refuse(tmp_path, 'test')
    assert previous.read_text() == 'existing results'

    @record_assembly_run
    def fail(outdir, name, evidence='loci'):
        activate_run_reporting()
        # Also handles unreadable output and failure after report preparation.
        (outdir/f'{name}.hdf5').write_bytes(b'incomplete')
        publish_assembly_report(outdir, name, _minimal_stats(), '')
        raise RuntimeError('evidence writing failed')
    with pytest.raises(RuntimeError, match='evidence writing failed'):
        fail(tmp_path, 'test')
    stats = json.loads(previous.read_text())
    assert stats['run']['status'] == 'FAILED'
    assert stats['evidence']['status'] == 'not_completed'
    assert 'summary' not in stats
    assert 'COMPLETE' not in (tmp_path/'test.stats.txt').read_text()


def test_python_api_does_not_invent_command(tmp_path):
    @record_assembly_run
    def run(outdir, name, evidence='none'):
        activate_run_reporting()
        publish_assembly_report(outdir, name, _minimal_stats(), '')
    run(tmp_path, 'test')
    stats = json.loads((tmp_path/'test.stats.json').read_text())
    assert stats['run']['invocation']['kind'] == 'python_api'
    assert stats['run']['invocation']['argv'] is None
    assert stats['command'] is None


def test_snp_status_and_field_denominators(tmp_path):
    db, raw, final = fixture(tmp_path)
    result = write_snp_evidence_hdf5(database=db, raw_vcf=raw, final_vcf=final, chunk_sites=2)
    assert result['sites'] == 4 and result['matched_native_records'] == 3
    assert result['cells'] == 8 and result['usable_pl_cells'] == 4
    assert result['status_counts'] == dict(PASS=3, NoHardCall=1, NoRead=1, HardMask=0, NoMatch=2, NoLikelihood=1)
    assert result['complete_field_cells'] == dict(DP=6, AD=5, GQ=4, PL=4)
    with h5py.File(db) as h:
        assert json.loads(h['evidence/summary_json'][:].tobytes()) == result
        statuses = h['evidence/status'][:]
        for i, row in enumerate(result['per_sample']):
            assert row['cells'] == 4
            assert list(row['status_counts'].values()) == np.bincount(statuses[i], minlength=6).tolist()
            assert row['usable_pl_cells'] == int(np.isin(statuses[i], [0, 1]).sum())


def test_locus_counts_distinguish_bases_records_masks_and_field_availability(tmp_path):
    db, ref, bed, raw, mask = setup(tmp_path)
    result = write_rad_evidence_hdf5(database=db, reference=ref, bed=bed, native_files=[raw],
                                   native_regions=[bed], hard_masks=({'RAD': mask},), chunk_sites=2)
    assert result['sites'] == 8 and result['records'] == 7 and result['cells'] == 14
    assert sum(result['site_status_counts'].values()) == 8
    assert result['site_status_counts'] == dict(NoRecord=1, HasRecord=6, NonACGT=1)
    assert result['symbolic_allele_records'] == 5
    assert result['complete_field_cells']['GQ'] == 0
    assert result['complete_field_cells']['PL'] > result['usable_pl_cells']
    with RadSiteEvidenceReader(db) as rd:
        values = rd.read(np.arange(rd.nsites))
        assert result['usable_pl_cells'] == int(values['usable_PL'].sum())
        for i, row in enumerate(result['per_sample']):
            assert row['cells'] == rd.nrecords
            assert row['usable_pl_cells'] == int(values['usable_PL'][i].sum())
            assert list(row['status_counts'].values()) == np.bincount(values['status'][i], minlength=5).tolist()
    with h5py.File(db) as h:
        assert json.loads(h['rad_evidence/summary_json'][:].tobytes()) == result
