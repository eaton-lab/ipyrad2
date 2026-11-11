"""Trim statistics schema 2.0; native fastp reports remain unmodified."""
from contextvars import ContextVar
from dataclasses import asdict
from functools import wraps
import hashlib
import inspect
import json
import math
from pathlib import Path
import shlex
import time
from uuid import uuid4

from loguru import logger

from .. import __version__
from ..cli.command_log import format_stats_command
from ..utils.exceptions import IPyradError
from ..utils.reporting import atomic_text, jsonable, utc_now

STATS_SCHEMA_VERSION = '2.0'
_RUN = ContextVar('trim_report', default=None)
FILTERS = ('low_quality_reads', 'too_many_N_reads', 'low_complexity_reads',
           'too_short_reads', 'too_long_reads')
DEFINITIONS = {
    'reads': 'Individual FASTQ records; a paired-end pair contributes two reads.',
    'before': 'Reads processed by fastp, before trimming/filtering; with max_reads this is a subset, not an input-file census.',
    'retention': 'Retained reads / processed reads; retained bases / processed bases.',
    'quality': 'Q20/Q30 bases / total bases at that stage; aggregate rates use summed base counts, not mean sample rates.',
    'filtering': 'Native fastp rejection categories; paired-end rejection is pairwise. Categories describe filtering decisions, not independent defect prevalence.',
    'operations': 'Adapter trimming and overlap correction can overlap each other and filtering; do not sum them as rejected reads.',
    'missing': 'JSON null and text N/A mean unavailable, inapplicable, or an undefined zero denominator; zero is a measured zero.',
}


def paired(paths):
    return paths[1] is not None and Path(paths[1]).name != '-null-'


def next_paths(outdir):
    idx = 0
    while True:
        paths = tuple(Path(outdir) / f'ipyrad_trim_stats_{idx}.{ext}' for ext in ('txt', 'json'))
        if not any(p.exists() for p in paths):
            return paths
        idx += 1


def record_trim_run(function):
    """Delay COMPLETE publication until the run and native-report checks succeed."""
    signature = inspect.signature(function)

    @wraps(function)
    def wrapped(*args, **kwargs):
        bound = signature.bind(*args, **kwargs)
        bound.apply_defaults()
        params = dict(bound.arguments)
        invocation = params.pop('invocation', None)
        supplied = params.pop('logged_command', None)
        run = dict(id=uuid4().hex, status='RUNNING', started_utc=utc_now(),
                   finished_utc=None, elapsed_seconds=None, ipyrad2_version=__version__,
                   working_directory=str(Path.cwd()), parameters=jsonable(params),
                   invocation=invocation or dict(kind='python_api', argv=None, command=None),
                   supplied_command=supplied if invocation is None else None,
                   inputs={}, cutsite_junctions={}, native_commands={})
        if invocation is not None:
            run['invocation'] = dict(invocation, command=shlex.join(invocation['argv']))
        context = dict(run=run, active=False, pending=None)
        token = _RUN.set(context)
        start = time.monotonic()
        try:
            result = function(*args, **kwargs)
            run.update(status='COMPLETE', finished_utc=utc_now(), elapsed_seconds=time.monotonic()-start)
            if context['pending'] is not None:
                _publish(context['pending'], run, context['paths'])
            return result
        except BaseException as error:
            if context['active']:
                run.update(status='FAILED', finished_utc=utc_now(), elapsed_seconds=time.monotonic()-start,
                           error=dict(type=type(error).__name__, message=str(error)))
                failed = getattr(error, 'key', None)
                run['sample_outcomes'] = {name: ('failed' if name == failed else 'not_confirmed')
                                          for name in run['inputs'] if name not in run['skipped']}
                try:
                    _publish(None, run, context['paths'])
                except Exception:
                    logger.exception('Could not write failure statistics; original trim error follows')
            raise
        finally:
            _RUN.reset(token)
    return wrapped


def activate(fastqs, skipped, outdir):
    """Activate only after input and overwrite guards have passed."""
    context = _RUN.get()
    if context is None:
        return
    context.update(active=True, paths=next_paths(outdir))
    context['run'].update(inputs={name: dict(layout='PE' if paired(paths) else 'SE',
                                            fastqs=jsonable(paths)) for name, paths in fastqs.items()},
                          skipped=dict(skipped))


def record_junctions(junctions, requested, disabled, has_paired):
    context = _RUN.get()
    if context is None:
        return
    for idx, junction in enumerate(junctions):
        applicable = idx == 0 or has_paired
        manual = bool(requested and requested[idx])
        source = ('not_applicable' if not applicable else 'manual' if manual else
                  'disabled' if disabled else 'inferred')
        info = dict(source=source, motifs=list(junction.motifs) if applicable else [],
                    offset=junction.offset if applicable else None,
                    trim_front_bases=junction.trim_length if applicable else None,
                    inference=asdict(junction) if source == 'inferred' else None)
        context['run']['cutsite_junctions'][f'R{idx+1}'] = jsonable(info)


def record_command(name, argv):
    context = _RUN.get()
    if context is not None:
        context['run']['native_commands'][name] = dict(argv=list(argv), command=shlex.join(argv))
        if 'adapter_fasta' not in context['run']:
            path = Path(argv[argv.index('--adapter_fasta') + 1])
            context['run']['adapter_fasta'] = dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def _number(mapping, key, required=False, integer=False):
    value = mapping.get(key)
    if value is None and not required:
        return None
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or value < 0 or (integer and int(value) != value)):
        raise ValueError(f'invalid or missing {key}: {value!r}')
    return int(value) if integer else value


def _ratio(numerator, denominator):
    return numerator / denominator if numerator is not None and denominator else None


def _stage(raw, layout):
    result = {key: _number(raw, key, required=True, integer=True)
              for key in ('total_reads', 'total_bases')}
    for q in ('q20', 'q30'):
        result[q+'_bases'] = _number(raw, q+'_bases', integer=True)
        if result[q+'_bases'] is not None and result[q+'_bases'] > result['total_bases']:
            raise ValueError(f'{q}_bases exceeds total_bases')
        result[q+'_rate'] = _ratio(result[q+'_bases'], result['total_bases'])
    for end in ('read1', 'read2'):
        result[end+'_mean_length'] = (_number(raw, end+'_mean_length')
                                     if end == 'read1' or layout == 'PE' else None)
    return result


def summarize_native(snames, outdir, loader):
    """Check native summary fields without altering or duplicating native JSON."""
    context = _RUN.get()
    inputs = context['run']['inputs'] if context else {}
    rows = []
    for name in sorted(snames):
        path = Path(outdir) / f'{name}.stats.json'
        native = loader(path)
        try:
            layout = inputs.get(name, {}).get('layout')
            if layout is None:
                layout = 'PE' if 'read2_mean_length' in native['summary']['before_filtering'] else 'SE'
            before = _stage(native['summary']['before_filtering'], layout)
            after = _stage(native['summary']['after_filtering'], layout)
            if context is not None:
                suffixes = ('.R1.trimmed.fastq.gz', '.R2.trimmed.fastq.gz') if layout == 'PE' else ('.trimmed.fastq.gz',)
                for suffix in suffixes:
                    if not (Path(outdir) / (name + suffix)).is_file():
                        raise ValueError(f'missing trimmed FASTQ: {name + suffix}')
            if after['total_reads'] > before['total_reads'] or after['total_bases'] > before['total_bases']:
                raise ValueError('retained totals exceed processed totals')
            filtering = {key: _number(native.get('filtering_result', {}), key, integer=True)
                         for key in ('passed_filter_reads', *FILTERS)}
            known = all(filtering[key] is not None for key in FILTERS)
            rejected = sum(filtering[key] for key in FILTERS) if known else None
            reconciled = (filtering['passed_filter_reads'] == after['total_reads']
                          and rejected + after['total_reads'] == before['total_reads']) if known else None
            operations = {key: _number(native.get('adapter_cutting', {}), key, integer=True)
                          for key in ('adapter_trimmed_reads', 'adapter_trimmed_bases')}
            operations.update({key: _number(native.get('filtering_result', {}), key, integer=True)
                               for key in ('corrected_reads', 'corrected_bases')})
            rows.append(dict(sample=name, layout=layout, before=before, after=after,
                             retention=dict(reads=_ratio(after['total_reads'], before['total_reads']),
                                            bases=_ratio(after['total_bases'], before['total_bases'])),
                             filtering=filtering, filtering_reconciled=reconciled, operations=operations,
                             fastp_version=native['summary'].get('fastp_version'),
                             native_report=dict(path=str(path.absolute()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())))
        except (KeyError, TypeError, ValueError) as err:
            raise IPyradError(f'Invalid fastp stats report {path}: {err}') from err
    if context is not None:
        expected = set(inputs) - set(context['run']['skipped'])
        if set(snames) != expected:
            raise IPyradError('Trim result sample IDs do not match the scheduled samples.')
    return rows


def _sum_known(rows, key):
    values = [row.get(key) for row in rows]
    return sum(values) if values and all(value is not None for value in values) else None


def aggregate(rows):
    out = {}
    for stage in ('before', 'after'):
        records = [r[stage] for r in rows]
        out[stage] = {key: _sum_known(records, key)
                      for key in ('total_reads', 'total_bases', 'q20_bases', 'q30_bases')}
        for q in ('q20', 'q30'):
            out[stage][q+'_rate'] = _ratio(out[stage][q+'_bases'], out[stage]['total_bases'])
    out['retention'] = {key: _ratio(out['after']['total_'+key], out['before']['total_'+key])
                        for key in ('reads', 'bases')}
    out['filtering'] = {key: _sum_known([r['filtering'] for r in rows], key)
                        for key in ('passed_filter_reads', *FILTERS)}
    out['operations'] = {key: _sum_known([r['operations'] for r in rows], key)
                         for key in ('adapter_trimmed_reads', 'adapter_trimmed_bases', 'corrected_reads', 'corrected_bases')}
    out['filtering_reconciled'] = (all(r['filtering_reconciled'] for r in rows)
                                   if rows and all(r['filtering_reconciled'] is not None for r in rows) else None)
    return out


def publish_summary(rows, outdir, logged_command=None):
    context = _RUN.get()
    if context is not None:
        context['pending'] = rows
        return
    # Standalone summarization has no observed run timing or process exit status.
    run = dict(status='REPORT_ONLY', id=uuid4().hex, ipyrad2_version=__version__,
               started_utc=None, finished_utc=None, elapsed_seconds=None,
               invocation=dict(kind='report_only', argv=None, command=None),
               supplied_command=logged_command, parameters={}, inputs={}, skipped={},
               cutsite_junctions={}, native_commands={}, working_directory=str(Path.cwd()))
    _publish(rows, run, next_paths(outdir))


def _safe(value):
    return ''.join(c if c.isprintable() else c.encode('unicode_escape').decode('ascii') for c in str(value))


def _fmt(value, kind='count'):
    if value is None:
        return 'N/A'
    if kind == 'percent':
        return f'{value:.2%}'
    if kind == 'float':
        return f'{value:,.1f}'
    return f'{value:,}' if isinstance(value, (int, float)) else _safe(value)


def _table(headers, rows):
    cells = [list(headers), *[[str(c) for c in row] for row in rows]]
    widths = [max(len(row[i]) for row in cells) for i in range(len(headers))]
    return '\n'.join('  '.join(cell.ljust(widths[i]) if i == 0 else cell.rjust(widths[i])
                                for i, cell in enumerate(row)).rstrip() for row in cells)


def _render(data):
    run, rows, totals = data['run'], data['sample_summary'], data['totals']
    inv = run['invocation']
    if inv.get('argv'):
        command = format_stats_command(inv.get('display_argv', inv['argv']))
    else:
        command = 'Command: Python API (resolved parameters in JSON)'
        if run.get('supplied_command'):
            command = 'Command (supplied): ' + _safe(run['supplied_command'])
            if len(command) > 120:
                command = command[:116] + ' ...'
    lines = ['ipyrad2 trim statistics', '=' * 72, command,
             f"Status: {run['status']}    Stats schema: {STATS_SCHEMA_VERSION}    ipyrad2: {__version__}",
             f"Started (UTC): {_fmt(run['started_utc'])}", f"Finished (UTC): {_fmt(run['finished_utc'])}",
             f"Elapsed seconds: {_fmt(run['elapsed_seconds'], 'float')}",
             'fastp version(s): ' + (', '.join(data['fastp_versions']) or 'N/A'), '', 'Inputs and processing']
    c = data['counts']
    lines += [f"  Input files: {_fmt(c['input_files'])}    Identified samples: {_fmt(c['identified_samples'])}",
              f"  Completed: {c['completed_samples']}    Skipped: {c['skipped_samples']}    Failed: {c['failed_samples']}    Not confirmed: {c['unconfirmed_samples']}",
              f"  Completed layouts: {c['single_end_samples']} SE, {c['paired_end_samples']} PE",
              '  Reads are individual FASTQ records; each PE pair contributes two reads.']
    params = run['parameters']
    if params:
        lines += [f"  Parallel jobs: {max(1, params['cores']//params['threads'])}    Threads per job: {params['threads']}",
                  f"  Read limit per file: {_fmt(params['max_reads']) if params['max_reads'] else 'all'}",
                  '', 'Trimming settings',
                  f"  Minimum retained length: {params['min_trimmed_length']} bp    Maximum Ns: {params['max_ns']}",
                  f"  Quality filtering: {'disabled' if params['disable_quality_filtering'] else 'enabled'}    Input Phred: {64 if params['phred64'] else 33}",
                  f"  Adapter trimming: {'disabled' if params['disable_adapter_trimming'] else 'enabled'}    UMI index2: {'PE only' if params['umi_tag_in_i5'] else 'disabled'}",
                  '  Poly-G, poly-X and low-complexity filters enabled; overlap correction for PE.']
        if not params['disable_quality_filtering']:
            lines += [f"  Qualified base Q: {params['min_quality']}    Unqualified bases allowed: {params['max_unqualified_percent']}%",
                      f"  Tail window: 5 bp    Minimum mean window Q: {params['min_mean_window_quality']}"]
        lines += ['  Full resolved parameters and native fastp commands are in JSON.']
    if run['cutsite_junctions']:
        lines += ['', 'Restriction-junction trimming', _table(
            ['End', 'Source', 'Motifs', 'Offset', 'Front bp', 'Supporting reads'],
            [[end, j['source'], _safe(','.join(j['motifs'])) or 'N/A', _fmt(j['offset']),
              _fmt(j['trim_front_bases']), _fmt(j['inference']['total_support'] if j['inference'] else None)]
             for end, j in run['cutsite_junctions'].items()])]
    if run.get('error'):
        lines += ['', 'Failure', '  ' + _safe(run['error']['message']),
                  '  No aggregate success statistics published; unconfirmed sample outcomes remain unknown.']
    if totals is not None:
        b, a = totals['before'], totals['after']
        lines += ['', 'Overall retention and quality',
                  'Before counts cover reads processed by fastp, not a census of input files.',
                  _table(['Metric', 'Before', 'After', 'Retained'], [
                      ['Reads', _fmt(b['total_reads']), _fmt(a['total_reads']), _fmt(totals['retention']['reads'], 'percent')],
                      ['Bases', _fmt(b['total_bases']), _fmt(a['total_bases']), _fmt(totals['retention']['bases'], 'percent')],
                      ['Q20 bases (%)', _fmt(b['q20_rate'], 'percent'), _fmt(a['q20_rate'], 'percent'), 'N/A'],
                      ['Q30 bases (%)', _fmt(b['q30_rate'], 'percent'), _fmt(a['q30_rate'], 'percent'), 'N/A']]),
                  '', 'Per-sample retention and quality',
                  _table(['Sample', 'Layout', 'Reads before', 'Reads after', 'Retained', 'Q30 before', 'Q30 after'],
                         [[_safe(r['sample']), r['layout'], _fmt(r['before']['total_reads']), _fmt(r['after']['total_reads']),
                           _fmt(r['retention']['reads'], 'percent'), _fmt(r['before']['q30_rate'], 'percent'),
                           _fmt(r['after']['q30_rate'], 'percent')] for r in rows]),
                  '', 'Mean read length (bp)',
                  _table(['Sample', 'R1 before', 'R1 after', 'R2 before', 'R2 after'],
                         [[_safe(r['sample']), *[_fmt(r[stage][end+'_mean_length'], 'float')
                                                for end in ('read1', 'read2') for stage in ('before', 'after')]] for r in rows]),
                  '', 'Read filtering (native fastp categories)',
                  _table(['Category', 'Reads', '% processed'],
                         [[key, _fmt(value), _fmt(_ratio(value, b['total_reads']), 'percent')]
                          for key, value in totals['filtering'].items()]),
                  '  Passed plus rejection categories reconcile: ' + ('N/A' if totals['filtering_reconciled'] is None else str(totals['filtering_reconciled'])),
                  '  PE rejection is pairwise; category totals are not independent defect counts.',
                  '', 'Per-sample rejected reads',
                  _table(['Sample', 'Low quality', 'Too many N', 'Low complexity', 'Too short', 'Too long'],
                         [[_safe(r['sample']), *[_fmt(r['filtering'][key]) for key in FILTERS]] for r in rows]),
                  '', 'Trimming and correction operations',
                  _table(['Operation', 'Count'], [[key, _fmt(value)] for key, value in totals['operations'].items()]),
                  '  These operations can overlap filtering and each other; do not add them to rejected reads.']
    if run['skipped']:
        lines += ['', 'Skipped samples', *['  '+_safe(reason) for reason in run['skipped'].values()]]
    lines += ['', 'Outputs', '  Text report: ' + data['outputs']['text'], '  JSON report: ' + data['outputs']['json'],
              '  Per-sample native JSON and actual FASTQ paths are listed in the run JSON.',
              '  Native fastp JSON is preserved; temporary HTML reports are removed.',
              '', 'N/A means unavailable, inapplicable, or an undefined denominator. Zero is a measured zero.']
    return '\n'.join(lines) + '\n'


def _publish(rows, run, paths):
    txt, js = paths
    rows = rows or []
    completed = {row['sample'] for row in rows}
    inputs = run['inputs']
    skipped = run['skipped']
    outcomes = run.get('sample_outcomes', {})
    counts = dict(input_files=sum(sum(p is not None and Path(p).name != '-null-' for p in info['fastqs'])
                                  for info in inputs.values()) if inputs else None,
                  identified_samples=len(inputs) if inputs else None, completed_samples=len(rows),
                  skipped_samples=len(skipped), failed_samples=sum(v == 'failed' for v in outcomes.values()),
                  unconfirmed_samples=len(set(inputs)-set(skipped)-completed)-sum(v == 'failed' for v in outcomes.values()),
                  single_end_samples=sum(r['layout'] == 'SE' for r in rows), paired_end_samples=sum(r['layout'] == 'PE' for r in rows))
    counts['unconfirmed_samples'] = max(0, counts['unconfirmed_samples'])
    outputs = dict(text=str(txt.absolute()), json=str(js.absolute()), samples={})
    for row in rows:
        name = row['sample']
        suffixes = ('.R1.trimmed.fastq.gz', '.R2.trimmed.fastq.gz') if row['layout'] == 'PE' else ('.trimmed.fastq.gz',)
        outputs['samples'][name] = dict(native_json=row['native_report'],
                                       fastqs=[str((txt.parent / (name+s)).absolute()) for s in suffixes
                                               if (txt.parent / (name+s)).is_file()])
    data = dict(stats_schema_version=STATS_SCHEMA_VERSION, tool='trim', run=run, counts=counts,
                fastp_versions=sorted({str(r['fastp_version']) for r in rows if r['fastp_version'] is not None}),
                definitions=DEFINITIONS, sample_summary=rows,
                totals=aggregate(rows) if rows else None, outputs=outputs)
    # JSON is serialized first, including a strict finite-value check, before publishing either file.
    encoded = json.dumps(data, indent=2, allow_nan=False)+'\n'
    text = _render(data)
    atomic_text(txt, text)
    atomic_text(js, encoded)
    logger.info('trimming stats written to {} and {}', txt, js)
