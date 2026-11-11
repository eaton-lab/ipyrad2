"""Versioned assembly reports, independent of biological output formats."""
from contextvars import ContextVar
from functools import wraps
import hashlib
import inspect
import json
from pathlib import Path
import shlex
import time
from uuid import uuid4

import h5py
import numpy as np
from loguru import logger

from .. import __version__
from ..cli.command_log import format_stats_command
from ..utils.reporting import atomic_text, jsonable as _jsonable, utc_now as _utc

STATS_SCHEMA_VERSION = '2.0'
_RUN = ContextVar('assembly_report', default=None)


def record_assembly_run(function):
    """Publish successful stats only after all output work and cleanup finish.

    Failures after output-directory preparation get an explicit failure report;
    early validation/overwrite refusals never replace an earlier run's report.
    """
    signature = inspect.signature(function)

    @wraps(function)
    def wrapped(*args, **kwargs):
        bound = signature.bind(*args, **kwargs)
        bound.apply_defaults()
        params = dict(bound.arguments)
        invocation = params.pop('invocation', None)
        logged = params.pop('logged_command', None)
        run = dict(id=uuid4().hex, status='RUNNING', started_utc=_utc(), finished_utc=None,
                   elapsed_seconds=None, ipyrad2_version=__version__,
                   working_directory=str(Path.cwd()), parameters=_jsonable(params),
                   invocation=invocation or dict(kind='python_api', argv=None, command=None),
                   supplied_command=logged if invocation is None else None,
                   inputs=None, calling=None)
        if invocation is not None:
            run['invocation'] = dict(invocation, command=shlex.join(invocation['argv']))
        context = dict(run=run, active=False, pending=None)
        token = _RUN.set(context)
        start = time.monotonic()
        try:
            result = function(*args, **kwargs)
            run.update(status='COMPLETE', finished_utc=_utc(), elapsed_seconds=time.monotonic()-start)
            if context['pending'] is not None:
                _publish(*context['pending'], run=run)
            return result
        except BaseException as error:
            if context['active']:
                run.update(status='FAILED', finished_utc=_utc(), elapsed_seconds=time.monotonic()-start,
                           error=dict(type=type(error).__name__, message=str(error)))
                try:
                    evidence = dict(mode=params['evidence'], status='not_completed', snps=None, loci=None)
                    _publish(Path(params['outdir']).expanduser().absolute(), params['name'],
                             {'evidence': evidence}, '', run=run, embed=False)
                except Exception:
                    logger.exception('Could not write failure statistics; original assembly error follows')
            raise
        finally:
            _RUN.reset(token)
    return wrapped


def activate_run_reporting():
    context = _RUN.get()
    if context is not None:
        context['active'] = True


def record_inputs(rad, wgs, group_file=None):
    context = _RUN.get()
    if context is None:
        return
    run = context['run']
    run['inputs'] = {'RAD': _jsonable(rad), 'WGS': _jsonable(wgs)}
    assignments = []
    if group_file is not None:
        assignments = [line.split() for line in Path(group_file).read_text().splitlines() if line.strip()]
    run['calling'] = dict(mode='grouped' if group_file is not None else 'independent',
                          sample_population_assignments=assignments)
    configs = {}
    for key in ('populations', 'rename', 'subsample', 'loci_bed'):
        value = run['parameters'].get(key)
        if value and Path(value).is_file():
            configs[key] = dict(path=value, sha256=hashlib.sha256(Path(value).read_bytes()).hexdigest())
    run['configuration_files'] = configs


class EvidenceCounts:
    """Accumulate small per-sample summaries while native arrays are in memory."""
    def __init__(self, names, statuses):
        self.names, self.statuses = list(names), tuple(statuses)
        self.counts = np.zeros((len(names), len(statuses)), dtype='i8')
        self.fields = {k: np.zeros(len(names), dtype='i8') for k in ('AD', 'DP', 'GQ', 'PL')}

    def add(self, data, allele_counts=None):
        status = data['status']
        for i in range(len(self.statuses)):
            self.counts[:, i] += np.sum(status == i, axis=1)
        for key in self.fields:
            values = data[key]
            if values.ndim == 3:
                widths = allele_counts if key == 'AD' else allele_counts*(allele_counts+1)//2
                meaningful = np.arange(values.shape[2])[None, :] < widths[:, None]
                present = np.all((values >= 0) | ~meaningful[None, :, :], axis=2)
            else:
                present = values >= 0
            self.fields[key] += np.sum(present, axis=1)

    def add_record(self, data):
        for i in range(len(self.statuses)):
            self.counts[:, i] += data['status'][:, 0] == i
        for key in self.fields:
            self.fields[key] += np.all(data[key] >= 0, axis=1)

    def add_ragged(self, data, allele_starts, likelihood_starts):
        """Summarize a whole batch with variable allele/genotype widths."""
        for i in range(len(self.statuses)):
            self.counts[:, i] += np.sum(data['status'] == i, axis=1)
        for key in self.fields:
            present = data[key] >= 0
            if key in ('AD', 'PL'):
                starts = allele_starts if key == 'AD' else likelihood_starts
                present = np.logical_and.reduceat(present, starts, axis=1)
            self.fields[key] += np.sum(present, axis=1)

    def summary(self):
        rows = []
        for i, name in enumerate(self.names):
            statuses = dict(zip(self.statuses, map(int, self.counts[i])))
            rows.append(dict(sample=name, cells=sum(statuses.values()), status_counts=statuses,
                             usable_pl_cells=statuses['PASS']+statuses.get('NoHardCall', 0),
                             complete_field_cells={k: int(v[i]) for k, v in self.fields.items()}))
        return dict(samples=len(rows), cells=int(self.counts.sum()),
                    status_counts=dict(zip(self.statuses, map(int, self.counts.sum(axis=0)))),
                    usable_pl_cells=sum(r['usable_pl_cells'] for r in rows),
                    complete_field_cells={k: int(v.sum()) for k, v in self.fields.items()},
                    per_sample=rows)


def store_summary(group, summary):
    group.create_dataset('summary_json', data=np.frombuffer(json.dumps(summary).encode(), dtype='u1'),
                         compression='gzip')


def publish_assembly_report(outdir, name, stats, body):
    context = _RUN.get()
    if context is not None:
        context['pending'] = (outdir, name, stats, body)
    else:
        _publish(outdir, name, stats, body)


def _storage(database, embed=None):
    if not database.is_file():
        return dict(hdf5_bytes=None, evidence_dataset_bytes={'snps': None, 'loci': None}), None
    with h5py.File(database, 'a' if embed is not None else 'r') as h:
        version = str(h.attrs.get('format_version', h.attrs.get('version', 'unknown')))
        if embed is not None:
            if 'assembly_report' in h:
                del h['assembly_report']
            group = h.create_group('assembly_report')
            group.attrs['stats_schema_version'] = STATS_SCHEMA_VERSION
            store_summary(group, embed)
        sizes = {}
        for scope, key in [('snps', 'evidence'), ('loci', 'rad_evidence')]:
            values = []
            if key in h:
                h[key].visititems(lambda n, d: values.append(int(d.id.get_storage_size()))
                                 if isinstance(d, h5py.Dataset) else None)
            sizes[scope] = sum(values) if key in h else None
        h.flush()
    return dict(hdf5_bytes=database.stat().st_size, evidence_dataset_bytes=sizes,
                dataset_bytes_note='Physical dataset allocation; excludes group metadata and some variable-length payload overhead.'), version


def _n(value):
    return 'N/A' if value is None else f'{value:,}'


def _bytes(value):
    return 'N/A' if value is None else f'{value:,} bytes ({value / 1024**3:.3f} GiB)'


def _fraction(numerator, denominator):
    return f'{_n(numerator)} / {_n(denominator)}' + (f' ({numerator/denominator:.1%})' if denominator else ' (N/A)')


def _section(title, rows):
    return [f'# {title}', *[f'  {label:<36} {value}' for label, value in rows], '']


def _evidence_text(evidence):
    mode = evidence.get('mode')
    lines = _section('Evidence Storage', [('Requested mode', mode or 'N/A'), ('Storage status', evidence['status'])])
    if mode == 'none' or evidence['status'] != 'complete':
        return lines
    lines += ['  Native AD/DP/PL and available GQ; no posterior probabilities are stored.',
              '  Evidence storage does not change standard genotype calls.', '']
    for scope, title in [('snps', 'SNP Evidence'), ('loci', 'All-locus Evidence')]:
        data = evidence.get(scope)
        if data is None:
            continue
        rows = []
        if scope == 'snps':
            rows += [('Retained SNP positions', _n(data['sites'])),
                     ('Exact native matches', _fraction(data['matched_native_records'], data['sites'])),
                     ('Unmatched SNP records', _fraction(data['sites']-data['matched_native_records'], data['sites']))]
        else:
            rows += [('Scope', 'Every base in the final retained locus BED'), ('Locus bases', _n(data['sites'])),
                     *[(f'Bases: {k}', _fraction(v, data['sites'])) for k,v in data['site_status_counts'].items()],
                     ('Native records', _n(data['records'])),
                     ('Symbolic-allele records', _fraction(data['symbolic_allele_records'], data['records']))]
        rows += [('Sample x '+('SNP' if scope=='snps' else 'record')+' cells', _n(data['cells'])),
                 *[(k, _fraction(v, data['cells'])) for k,v in data['status_counts'].items()],
                 ('Usable native PL' if scope == 'loci' else 'Usable PL incl. NoHardCall',
                  _fraction(data['usable_pl_cells'], data['cells'])),
                 *[(f'Complete native {k}', _fraction(v, data['cells'])) for k,v in data['complete_field_cells'].items()]]
        lines += _section(title, rows)
        if scope == 'loci':
            lines += ['  A base can have multiple records. Symbolic alleles are literal native tokens.',
                      '  Raw masked values remain available for audit but are not usable evidence.', '']
        for kind, values in data.get('by_sample_type', {}).items():
            lines.append(f"  {kind} usable PL: {_fraction(values['usable_pl_cells'], values['cells'])}")
        lines += ['']
    lines += ['  Status categories are mutually exclusive; field availability is separate.',
              '  Usable PL requires positive depth, complete PL and no hard mask.',
              '  Uncovered cells remain missing, not homozygous reference.',
              '  Per-sample counts and exact denominators are in stats.json.', '']
    return lines


def _publish(outdir, name, stats, body, run=None, embed=True):
    stats = dict(stats)
    run = run or dict(status='COMPLETE', started_utc=None, finished_utc=None, elapsed_seconds=None,
                     ipyrad2_version=__version__, working_directory=None, parameters=None,
                     invocation=dict(kind='not_recorded', argv=None, command=None))
    evidence = stats.setdefault('evidence', dict(mode=None, status='not_recorded', snps=None, loci=None))
    if evidence.get('temporary_pileup_bytes') is not None and run.get('parameters') is not None:
        evidence['temporary_pileup_disposition'] = 'retained (--keep-tmpdir)' if run['parameters']['keep_tmpdir'] else 'removed after success'
    stats.update(stats_schema_version=STATS_SCHEMA_VERSION, run=run)
    database = Path(outdir)/f'{name}.hdf5'
    # Embedded provenance omits filesystem byte counts to avoid self-referential sizes.
    if embed:
        storage, version = _storage(database, dict(stats_schema_version=STATS_SCHEMA_VERSION, run=run, evidence=evidence))
    else:
        # Failure reporting must also work for an incomplete/unreadable HDF5.
        storage, version = dict(hdf5_bytes=None, evidence_dataset_bytes={'snps': None, 'loci': None}), None
    stats['storage'], stats['hdf5_format_version'] = storage, version
    stats['outputs'] = {kind: str(Path(outdir)/f'{name}.{suffix}') for kind, suffix in
                        [('database', 'hdf5'), ('locus_sequences', 'loci.gz'), ('variants', 'vcf.gz'),
                         ('locus_intervals', 'bed'), ('text_stats', 'stats.txt'), ('json_stats', 'stats.json')]}

    inv = run['invocation']
    command = inv.get('command') or stats.get('command')
    stats['command'] = command
    command_tokens = inv.get('display_argv') or inv.get('argv')
    if command_tokens is None and command:
        try:
            command_tokens = shlex.split(command)
        except ValueError:
            command_tokens = [command]
    command_line = format_stats_command(command_tokens) if command_tokens else 'Command: N/A (Python API; see resolved parameters in stats.json)'
    run_lines = _section('Run', [('Status', run['status']), ('Started UTC', run['started_utc'] or 'N/A'),
                              ('Finished UTC', run['finished_utc'] or 'N/A'),
                              ('Elapsed seconds', 'N/A' if run['elapsed_seconds'] is None else f"{run['elapsed_seconds']:.3f}"),
                              ('ipyrad2 version', run['ipyrad2_version']), ('HDF5 format', version or 'N/A'),
                              ('Stats schema', STATS_SCHEMA_VERSION)])
    if run.get('calling'):
        run_lines += _section('Calling', [('Mode', run['calling']['mode']),
                              ('Populations file', run['parameters']['populations'] or 'N/A')])
    timings = evidence.get('timing_seconds', {})
    storage_lines = _section('Storage and Time', [('HDF5 file size', _bytes(storage['hdf5_bytes'])),
             *[(f'{k} evidence datasets', _bytes(v)) for k,v in storage['evidence_dataset_bytes'].items()],
             *[(k.replace('_',' '), 'N/A' if v is None else f'{v:.3f} s') for k,v in timings.items()],
             ('Temporary pileup BCF size', _bytes(evidence.get('temporary_pileup_bytes'))),
             ('Temporary pileup disposition', evidence.get('temporary_pileup_disposition') or 'N/A')])
    storage_lines += ['  Dataset bytes exclude metadata/variable-length overhead; file size includes it.',
                      '  Timings are measured operations, not added cost versus a lean assembly.',
                      '  N/A means not measured or not applicable, not zero.', '']
    extra = '\n'.join(_evidence_text(evidence)+storage_lines)
    # Keep the main biological summary first, ahead of evidence and detailed tables.
    split = body.find('# Mixed RAD/WGS Diagnostics')
    if split < 0:
        split = body.find('# Locus Filtering')
    if split < 0:
        split = len(body)
    text = '\n'.join(run_lines)+command_line+'\n\n'+body[:split]+extra+'\n'+body[split:]
    if run['status'] == 'COMPLETE':
        text += '\n'+'\n'.join(_section('Outputs', [(k.replace('_', ' '), Path(v).name) for k,v in stats['outputs'].items()]))
    text += '\nFull invocation, resolved parameters and per-sample statistics: '+f'{name}.stats.json\n'
    if run['status'] == 'FAILED':
        text += '\nFailure: '+run['error']['message']+'\n'
    atomic_text(Path(outdir)/f'{name}.stats.json', json.dumps(stats, indent=2)+'\n')
    atomic_text(Path(outdir)/f'{name}.stats.txt', text)
