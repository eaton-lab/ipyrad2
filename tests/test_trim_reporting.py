"""Statistical denominators, provenance, and failure semantics of trim reports."""
import copy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from ipyrad2.cli.command_log import format_stats_command
from ipyrad2.trimmer import reporting
from ipyrad2.trimmer import trim_fastqs as trim
from ipyrad2.utils.exceptions import IPyradError
from ipyrad2.utils.parallel.pool import ParallelJobError


def native(reads=10, bases=1000, q30=900, paired=False):
    stage = dict(total_reads=reads, total_bases=bases, q20_bases=bases, q30_bases=q30,
                 q30_rate=0.123, read1_mean_length=bases/reads if reads else 0)
    if paired:
        stage['read2_mean_length'] = stage['read1_mean_length']
    return dict(summary=dict(fastp_version='test', before_filtering=stage, after_filtering=copy.deepcopy(stage)),
                filtering_result=dict(passed_filter_reads=reads, **{key: 0 for key in reporting.FILTERS}))


def write_native(tmp, name, value):
    path = tmp / f'{name}.stats.json'
    path.write_text(json.dumps(value))
    return path


def api_params(tmp):
    return dict(fastqs=[tmp/'reads.fastq'], outdir=tmp/'out', cutsite_motifs=None,
                max_reads=1, min_trimmed_length=1, max_unqualified_percent=15,
                min_quality=20, min_mean_window_quality=30, cut_window_size=5,
                phred64=False, max_reads_kmer=100, max_ns=5, disable_infer_cutsite_motifs=True,
                disable_adapter_trimming=True, disable_quality_filtering=False,
                cores=1, threads=1, delim_str=None, delim_idx=1, suffix=None,
                umi_tag_in_i5=False, force=False, log_level='ERROR')


def test_weighted_quality_optional_fields_mixed_layout_and_native_preservation(tmp_path):
    a = native(reads=1, bases=10, q30=0)
    b = native(reads=20, bases=1000, q30=1000, paired=True)
    b['adapter_cutting'] = dict(adapter_trimmed_reads=3, adapter_trimmed_bases=15)
    paths = [write_native(tmp_path, 'se', a), write_native(tmp_path, 'pe', b)]
    original = [p.read_bytes() for p in paths]
    trim.write_stats_summary(['se', 'pe'], tmp_path)
    data = json.loads((tmp_path/'ipyrad_trim_stats_0.json').read_text())
    assert data['totals']['before']['q30_rate'] == pytest.approx(1000/1010)
    assert data['totals']['before']['q30_bases'] == 1000
    pe, se = data['sample_summary']
    assert pe['before']['read2_mean_length'] == 50
    assert se['before']['read2_mean_length'] is None
    assert pe['operations']['adapter_trimmed_reads'] == 3
    assert se['operations']['adapter_trimmed_reads'] is None
    assert data['totals']['operations']['adapter_trimmed_reads'] is None
    assert data['totals']['filtering']['too_short_reads'] == 0
    assert data['totals']['filtering_reconciled'] is True
    assert [p.read_bytes() for p in paths] == original
    for row in data['sample_summary']:
        p = Path(row['native_report']['path'])
        assert hashlib.sha256(p.read_bytes()).hexdigest() == row['native_report']['sha256']
    text = (tmp_path/'ipyrad_trim_stats_0.txt').read_text()
    assert 'R2 before' in text and 'N/A' in text
    assert '-->' not in text and '→' not in text
    assert data['run']['status'] == 'REPORT_ONLY'
    first = (tmp_path/'ipyrad_trim_stats_0.json').read_bytes()
    trim.write_stats_summary(['se'], tmp_path)
    assert (tmp_path/'ipyrad_trim_stats_1.json').exists()
    assert (tmp_path/'ipyrad_trim_stats_0.json').read_bytes() == first


def test_zero_denominator_and_all_reads_filtered(tmp_path):
    zero = native(0, 0, 0)
    removed = native()
    removed['summary']['after_filtering'] = zero['summary']['after_filtering']
    removed['filtering_result'].update(passed_filter_reads=0, too_short_reads=10)
    write_native(tmp_path, 'zero', zero)
    write_native(tmp_path, 'removed', removed)
    trim.write_stats_summary(['zero', 'removed'], tmp_path)
    data = json.loads((tmp_path/'ipyrad_trim_stats_0.json').read_text())
    removed, zero = data['sample_summary']
    assert removed['retention']['reads'] == 0
    assert zero['retention']['reads'] is None
    assert data['totals']['after']['q30_rate'] is None
    assert data['totals']['filtering_reconciled'] is True


def test_incomplete_filter_counts_are_not_assumed_zero(tmp_path):
    report = native()
    del report['filtering_result']['too_long_reads']
    write_native(tmp_path, 'a', report)
    trim.write_stats_summary(['a'], tmp_path)
    data = json.loads((tmp_path/'ipyrad_trim_stats_0.json').read_text())
    assert data['totals']['filtering']['too_long_reads'] is None
    assert data['totals']['filtering_reconciled'] is None


@pytest.mark.parametrize('value', [{}, {'summary': None}, native(bases=-1), native(q30=float('nan'))])
def test_invalid_native_report_does_not_publish_success(tmp_path, value):
    write_native(tmp_path, 'bad', value)
    with pytest.raises(IPyradError, match='Invalid fastp stats report'):
        trim.write_stats_summary(['bad'], tmp_path)
    assert not list(tmp_path.glob('ipyrad_trim_stats_*'))


def test_full_cli_invocation_fastq_summary_line_and_control_characters(monkeypatch):
    import ipyrad2.cli.cli_main as cli
    argv = ['trim', '--fastqs', *[f'data with spaces/sample{i}.fastq.gz' for i in range(375)], '-o', 'OUT']
    args = cli.setup_parsers().parse_args(argv)
    seen = {}
    monkeypatch.setattr(cli.importlib, 'import_module', lambda *a, **k: SimpleNamespace(run_trimmer=lambda **kw: seen.update(kw)))
    monkeypatch.setattr(cli.sys, 'argv', ['/env/bin/ipyrad2', *argv])
    cli.run_subcommand(args, _exit=False)
    assert seen['invocation']['argv'] == ['/env/bin/ipyrad2', *argv]
    line = format_stats_command(seen['invocation']['display_argv'])
    assert '<375 FASTQs>' in line and 'BAMs' not in line
    assert len(line) <= 120 and len(line.splitlines()) == 1
    line = format_stats_command(['ipyrad2', 'trim', '-d', 'embedded\nnewline.fastq'])
    assert '\\n' in line and len(line.splitlines()) == 1


def test_failed_pool_does_not_reuse_stale_native_reports(tmp_path, monkeypatch):
    params = api_params(tmp_path)
    params['force'] = True
    params['fastqs'][0].write_text('@read\nACGTACGT\n+\nIIIIIIII\n')
    params['outdir'].mkdir()
    stale = write_native(params['outdir'], 'reads', native())
    original = stale.read_bytes()
    (params['outdir']/'ipyrad_trim_stats_0.json').write_text('old run')
    def fail(*args, **kwargs):
        raise ParallelJobError('reads', 'fastp', {}, RuntimeError('native failure'))
    monkeypatch.setattr(trim, 'run_with_pool', fail)
    with pytest.raises(ParallelJobError):
        trim.run_trimmer(**params)
    data = json.loads((params['outdir']/'ipyrad_trim_stats_1.json').read_text())
    assert data['run']['status'] == 'FAILED'
    assert data['counts']['failed_samples'] == 1
    assert data['sample_summary'] == [] and data['totals'] is None
    assert stale.read_bytes() == original
    assert (params['outdir']/'ipyrad_trim_stats_0.json').read_text() == 'old run'
    assert data['run']['cutsite_junctions']['R1']['inference'] is None
    assert data['run']['native_commands']['reads']['argv'][0].endswith('/fastp')


def test_early_refusal_preserves_report_history(tmp_path):
    params = api_params(tmp_path)
    params['fastqs'][0].write_text('@read\nACGTACGT\n+\nIIIIIIII\n')
    params['outdir'].mkdir()
    write_native(params['outdir'], 'reads', native())
    with pytest.raises(IPyradError, match='exists in outdir'):
        trim.run_trimmer(**params)
    assert not list(params['outdir'].glob('ipyrad_trim_stats_*'))
