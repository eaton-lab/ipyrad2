"""Human-readable de novo reports; scientific JSON fields remain stable."""
from pathlib import Path
import shlex

from .. import __version__
from ..cli.command_log import format_stats_command

DEFINITIONS = {
    'builder': 'A selected sample used to construct the reference library. Representation here refers only to builders, before downstream mapping of all samples.',
    'input_sequence': 'A representative of a within-sample read cluster, called a seed or core in diagnostic tables. Counted in selected_sample_summary.consensus_records.',
    'similarity_group': 'A connected group of input sequences linked by sequence similarity; called a component in diagnostic tables. It can contain several loci.',
    'output_membership': 'A sequence-to-locus mapping row. denovo_summary.consensus_records counts these rows, not the input sequences or FASTA entries.',
    'same_sample_merging': 'Combining compatible sequences from the same builder before splitting a similarity group; called reconciliation in diagnostic tables. This does not certify orthology.',
    'audit_scope': 'component_qc and component_node_summary describe only audited groups (same-builder duplicates or oversized groups), not all similarity groups.',
    'placeholder': 'An oversized group retained as one unsplit reference entry using representative sequences. It has not been resolved into individual loci.',
    'fragment_support': 'Sum of native abundance counts for retained within-sample clusters; paired mates contribute once after merging/joining. Not a census of raw FASTQ reads or a retention rate.',
    'representation': 'Builder samples represented in each output locus; not downstream mapping coverage, depth or genotype completeness.',
    'output_spacer': 'N padding actually inserted between paired-end arms in a completed output sequence. Joined-only, mixed PE and no-spacer counts partition written loci. Natural ambiguity Ns and planned layout tags are not inserted spacers.',
    'missing': 'Text N/A means unavailable, inapplicable or an undefined denominator. Legacy scientific JSON values are preserved.',
}

STAGE_LABELS = {
    'preflight': 'Check inputs', 'within_samples': 'Cluster reads within builders',
    'combine_summaries': 'Combine builder sequences', 'across_samples': 'Find sequence similarities',
    'graph': 'Separate similarity groups', 'alignment': 'Build reference sequences',
    'final_qc': 'Summarize library', 'summary': 'Write reports', 'cleanup': 'Clean temporary files',
}


def _safe(value):
    return ''.join(c if c.isprintable() else c.encode('unicode_escape').decode('ascii') for c in str(value))


def _number(value):
    if value is None:
        return None
    return float(str(value).replace(',', ''))


def _count(value):
    number = _number(value)
    return 'N/A' if number is None else f'{number:,.0f}'


def _float(value, digits=2):
    number = _number(value)
    return 'N/A' if number is None else f'{number:,.{digits}f}'


def _percent(numerator, denominator=1):
    a, b = _number(numerator), _number(denominator)
    return 'N/A' if a is None or not b else f'{a / b:.2%}'


def _share(numerator, denominator):
    return f'{_count(numerator)} / {_count(denominator)} ({_percent(numerator, denominator)})'


def _section(title, rows):
    width = max((len(k) for k, _ in rows), default=0)
    return [f'# {title}', *[f'  {k:<{width}}  {_safe(v)}' for k, v in rows], '']


def _table(headers, rows):
    cells = [headers, *[[str(v) for v in row] for row in rows]]
    widths = [max(len(row[i]) for row in cells) for i in range(len(headers))]
    return ['  '.join(v.ljust(widths[i]) if i == 0 else v.rjust(widths[i])
                     for i, v in enumerate(row)).rstrip() for row in cells]


def render_report(data, *, failure_body=''):
    """Render existing counters without changing their definitions or values."""
    run = data.get('run', {})
    inv = run.get('invocation', {})
    argv = inv.get('display_argv') or inv.get('argv')
    supplied = run.get('supplied_command') or data.get('command')
    if argv:
        command = format_stats_command(argv)
    elif supplied:
        try:
            command = format_stats_command(shlex.split(supplied))
        except ValueError:
            command = 'Command (supplied): ' + _safe(supplied)[:100]
    else:
        command = 'Command: Python API (resolved parameters in JSON)'
    lines = ['ipyrad2 denovo statistics', '=' * 72, command, '']
    lines += _section('Run', [
        ('Status', run.get('status', 'REPORT_ONLY')),
        ('Started UTC', run.get('started_utc') or 'N/A'),
        ('Finished UTC', run.get('finished_utc') or 'N/A'),
        ('Elapsed seconds', _float(run.get('elapsed_seconds'), 3)),
        ('ipyrad2 version', __version__), ('Stats schema', data.get('stats_schema_version', '2.0')),
    ])
    failed = run.get('status') in {'FAILED', 'RECOVERY_REQUIRED'}
    if failed:
        error = run.get('error', {})
        lines += _section('Failure', [('Error', f"{error.get('type', 'Error')}: {error.get('message', 'N/A')}")])
        lines += ['  Reference library not finalized; available stage measurements follow.', '']
        if failure_body:
            lines += [_safe(failure_body.strip()), '']
    elif 'denovo_summary' in data:
        lines += _library_sections(data)
    lines += _execution_sections(data)
    lines += ['Full parameters, scientific counters, field definitions and paths: denovo.stats.json',
              'N/A means unavailable, inapplicable or an undefined denominator; zero is a count.', '']
    return '\n'.join(lines)


def _library_sections(data):
    inputs, summary = data['inputs'], data['denovo_summary']
    qc, audit = data['locus_qc'], data['component_qc']
    samples = data['selected_sample_summary']
    total = summary.get('loci_written')
    builders = inputs.get('selected_sample_count')
    graph = next((s.get('details', {}) for s in data.get('stages', []) if s['name'] == 'graph'), {})
    seq_count = sum(_number(s['consensus_records']) for s in samples)
    lines = _section('Reference Library', [
        ('Reference sequences (loci)', _count(total)),
        ('Represented by 2+ builders', _share(qc.get('loci_with_2plus_samples'), total)),
        ('Represented by one builder', _share(qc.get('singleton_loci'), total)),
        ('Represented by all builders', _share(qc.get('loci_with_all_selected_samples'), total)),
        ('Unresolved large-group placeholders', _count(audit.get('oversize_unsplit_components'))),
    ])
    if _number(audit.get('oversize_unsplit_components')):
        lines += ['  Placeholders are included in the reference count; they are unsplit groups,',
                  '  not resolved individual loci.', '']
    lines += _section('Inputs and Builder Samples', [
        ('Input samples available', _count(inputs.get('total_input_sample_count'))),
        ('Builder samples selected', _count(builders)),
        ('Builder FASTQ files', _count(inputs.get('fastq_count'))),
        ('Read layout', inputs.get('paired_mode', 'N/A')),
        ('Sample selection mode', inputs.get('sample_selection_mode', 'N/A')),
        ('Input representative sequences', _count(seq_count)),
        ('Sequence-to-locus output memberships', _count(summary.get('consensus_records'))),
    ])
    lines += ['  Builders are the selected samples used to construct this library.',
              '  Each input representative (seed/core) comes from a within-sample read cluster.',
              '  Output memberships count mapping-table rows, not FASTA entries.', '']
    lines += _table(['Sample', 'Input sequences', 'Fragment support', 'Joined PE', 'Merged PE', 'SE'],
                    [[_safe(s['sample']), *[_count(s.get(k)) for k in
                      ['consensus_records','n_reads_sum','joined_records','merged_records','single_records']]] for s in samples])
    lines += ['  Joined PE: non-overlapping mates; merged PE: overlapping mates combined.',
              '  Fragment support sums native abundance counts in retained clusters;',
              '  paired mates count once. It is not raw-read retention.', '']
    lines += _section('Builder Representation', [
        ('Builders per locus: mean / median / max', ' / '.join([
            _float(qc.get('mean_samples_per_locus') if _number(total) else None),
            _float(qc.get('median_samples_per_locus') if _number(total) else None),
            _count(qc.get('max_samples_per_locus') if _number(total) else None)])),
        ('Loci with at least half the builders', _share(qc.get('loci_with_half_or_more_selected_samples'), total)),
    ])
    lines += _table(['Builders per locus', 'Loci', '% of reference loci'],
                    [[_count(r['samples_with_data']), _count(r['loci']), _percent(r['loci'], total)] for r in data['locus_occupancy']])
    lines += ['  These counts describe builders only, before mapping all samples.', '']
    lines += _section('Similarity Groups', [
        ('All similarity groups', _count(graph.get('components'))),
        ('Groups with multiple sequences from a builder', _count(summary.get('duplicated_components_seen'))),
        ('Groups tested for same-builder merging', _count(summary.get('same_sample_reconciliation_attempted'))),
        ('Groups using same-builder merging', _count(summary.get('components_reconciled'))),
        ('Output loci using same-builder merging', _count(qc.get('reconciled_loci'))),
        ('Output loci from groups with builder duplicates', _count(qc.get('duplicated_component_loci'))),
    ])
    lines += ['  A similarity group (component) connects sequences through similarity links.',
              '  It may contain several loci. Compatible sequences from one builder may be',
              '  combined before splitting (reconciliation). Duplicate sequences do not',
              '  by themselves establish paralogy; these counts are not accuracy estimates.', '']
    audited = audit.get('audited_components')
    lines += _section('Detailed Group Audit', [
        ('Groups recorded in audit', _count(audited)),
        ('Audited groups processed normally', _count(audit.get('processed_components'))),
        ('Audited groups retained as large placeholders', _count(audit.get('oversize_unsplit_components'))),
        ('Largest audited group (input sequences)', _count(audit.get('largest_component_nodes') if _number(audited) else None)),
    ])
    lines += ['  Audit statistics cover only groups with same-builder duplicates or',
              '  oversized groups, not all similarity groups.']
    if _number(audited):
        lines += _table(['Audited-group size', 'Input sequences', 'After simplification'],
                        [[{'p50':'Median','p90':'90th percentile','p99':'99th percentile','max':'Maximum'}[r['quantile']],
                          _float(r['input_nodes'], 3), _float(r['contracted_nodes'], 3)] for r in data['component_node_summary']])
        lines += ['  Simplification combines compatible same-builder sequences; oversized',
                  '  placeholders instead select representatives. Smaller is not proof of resolution.']
    else:
        lines += ['  No groups required a detailed audit; size statistics are not applicable.']
    lines += ['']
    aligned = data['runtime'].get('alignment_mode') != 'none'
    lines += _section('Sequence Construction', [
        ('Single-sequence loci', _count(summary.get('single_sequence_loci') if aligned else None)),
        ('Loci with identical sequences', _count(summary.get('identical_sequence_loci') if aligned else None)),
        ('Loci requiring alignment in this mode', _count(summary.get('mafft_required_loci') if aligned else None)),
        ('Loci with a spacer: joined PE only', _count(summary.get('joined_spacer_loci'))),
        ('Loci with a spacer: mixed PE', _count(summary.get('mixed_reconciled_spacer_loci'))),
        ('Loci without an inserted spacer', _count(summary.get('stripped_output_loci'))),
        ('Sequence memberships per locus: mean / median / max', ' / '.join([
            _float(qc.get('mean_cores_per_locus') if _number(total) else None),
            _float(qc.get('median_cores_per_locus') if _number(total) else None),
            _count(qc.get('max_cores_per_locus') if _number(total) else None)])),
        ('One-builder loci with multiple input sequences', _count(qc.get('multi_core_single_sample_loci'))),
    ])
    lines += ['  Alignment categories are N/A when --no-alignment bypasses classification.',
              '  Spacer counts reflect N padding actually inserted between PE arms.',
              '  Natural ambiguity Ns are not spacers. SE loci do not receive PE padding.',
              '  Joined-only, mixed PE and no-spacer categories partition the output loci.', '']
    params = data['clustering_parameters'];runtime = data['runtime']
    lines += _section('Settings', [
        ('Within-builder similarity', _percent(params.get('within_similarity'))),
        ('Across-builder similarity', _percent(params.get('across_similarity'))),
        ('Minimum query coverage', _percent(params.get('query_cov'))),
        ('Minimum identical-fragment count', _count(params.get('min_derep_size'))),
        ('Minimum sequence length (bp)', _count(params.get('min_length'))),
        ('Minimum PE merge overlap (bp)', _count(params.get('min_merge_overlap')) if inputs.get('paired_mode') == 'paired-end' else 'N/A (single-end)'),
        ('Maximum PE merge differences', _count(params.get('max_merge_diffs')) if inputs.get('paired_mode') == 'paired-end' else 'N/A (single-end)'),
        ('Search reverse complements', {'True':'yes','False':'no'}.get(str(params.get('allow_reverse_complement')), 'N/A')),
        ('Alignment mode', runtime.get('alignment_mode', 'N/A')),
        ('CPU cores requested', _count(runtime.get('cores'))),
        ('VSEARCH jobs / threads per job', f"{_count(runtime.get('vsearch_worker_processes'))} / {_count(runtime.get('vsearch_threads_per_job'))}"),
        ('MAFFT jobs / threads per job', f"{_count(runtime.get('mafft_worker_processes'))} / {_count(runtime.get('mafft_threads_per_job'))}"),
        ('MAFFT timeout (seconds)', _count(runtime.get('mafft_timeout_seconds'))),
        ('Keep intermediate files', {'True':'yes','False':'no'}.get(str(runtime.get('keep_intermediates')), 'N/A')),
    ])
    outputs = data.get('outputs', {})
    labels = {'reference':'Reference FASTA','mapping':'Sequence-to-locus mapping','loci_stats':'Per-locus statistics',
              'sample_graph_summary':'Per-builder graph statistics','run_stats':'Text report','run_stats_json':'JSON report','audit_dir':'Group audit directory'}
    directory = data.get('run', {}).get('output_directory') or str(Path(outputs.get('reference', '.')).parent)
    lines += _section('Outputs', [('Directory', directory), *[(label, Path(outputs[key]).name) for key, label in labels.items() if key in outputs],
                                              ('Intermediate files', outputs.get('intermediates', 'N/A'))])
    return lines


def _execution_sections(data):
    stages = data.get('stages', [])
    if not stages:
        return []
    lines = ['# Time and Memory']
    lines += _table(['Step', 'Status', 'Seconds', 'Peak sampled RSS (MiB)'],
                    [[STAGE_LABELS.get(s['name'], _safe(s['name'])), s['status'], _float(s.get('elapsed_seconds'), 3),
                      _float(s['sampled_peak_tree_rss_bytes']/2**20 if s.get('sampled_peak_tree_rss_bytes') is not None else None, 1)] for s in stages])
    lines += ['  RSS sums the main process and live children; it can double-count shared',
              '  pages and miss short peaks. Stage memory values must not be summed.', '']
    graph = next((s.get('details') for s in stages if s['name'] == 'graph'), None)
    if graph:
        lines += _section('Execution Details: Group Processing', [
            ('Similarity groups / worker batches', f"{_count(graph.get('components'))} / {_count(graph.get('batches_submitted'))}"),
            ('Running + waiting batch limit', _count(graph.get('max_outstanding_batches'))),
            ('Peak writer reorder queue (groups)', _count(graph.get('peak_buffered_components'))),
            ('Prepare component jobs (seconds)', _float(graph.get('job_build_seconds'), 3)),
            ('Dispatch / collect (seconds)', _float(graph.get('dispatch_collect_seconds'), 3)),
            ('Write locus tables (seconds)', _float(graph.get('table_write_seconds'), 3)),
            ('Write group audit (seconds)', _float(graph.get('audit_write_seconds'), 3)),
            ('Flush output streams (seconds)', _float(graph.get('output_flush_seconds'), 3)),
        ])
        lines += ['  The writer queue excludes results held by the process pool.',
                  '  Batch limits do not bound total memory or one indivisible group.', '']
    alignment = next((s.get('details') for s in stages if s['name'] == 'alignment'), None)
    if alignment:
        lines += _section('Execution Details: Alignment', [
            ('MAFFT calls: total / successful / failed', ' / '.join(_count(alignment.get(k)) for k in ['mafft_invocations','mafft_succeeded','mafft_failed'])),
            ('MAFFT execution (summed worker seconds)', _float(alignment.get('mafft_execute_seconds_sum'), 3)),
            ('Peak running + waiting loci / limit', f"{_count(alignment.get('peak_outstanding_loci'))} / {_count(alignment.get('outstanding_locus_limit'))}"),
            ('Peak waiting results', _count(alignment.get('peak_buffered_results'))),
            ('Peak buffered output / budget (bytes)', f"{_count(alignment.get('peak_buffered_payload_bytes'))} / {_count(alignment.get('output_payload_budget_bytes'))}"),
            ('Loci exceeding output budget (run alone)', _count(alignment.get('oversized_loci'))),
            ('Waiting for buffer space (seconds)', _float(alignment.get('backpressure_wait_seconds'), 3)),
        ])
    lines += ['  Dispatch includes waiting, transport and pool startup/shutdown.',
              '  Worker times overlap each other and parent work; do not add them to wall time.',
              '  Output budgets exclude input records, Python overhead and native-tool memory.',
              '  Detailed timings and duration histograms remain in JSON.', '']
    return lines
