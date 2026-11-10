"""Report labels and denominators must reflect the existing scientific counters."""
import copy
import json
import re

import pytest
from ipyrad2.denovo.reporting import render_report
from ipyrad2.denovo.runtime import write_runtime_report


def example():
    return dict(
        inputs=dict(selected_sample_count='4', total_input_sample_count='375', fastq_count='8', paired_mode='paired-end', sample_selection_mode='imap'),
        denovo_summary=dict(loci_written='70',consensus_records='90',duplicated_components_seen='4',same_sample_reconciliation_attempted='3',components_reconciled='2',single_sequence_loci='2',identical_sequence_loci='3',mafft_required_loci='65',joined_spacer_loci='5',mixed_reconciled_spacer_loci='2',stripped_output_loci='63'),
        locus_qc=dict(loci_with_2plus_samples='68',singleton_loci='2',loci_with_all_selected_samples='40',loci_with_half_or_more_selected_samples='68',mean_samples_per_locus='3.0',median_samples_per_locus='4',max_samples_per_locus='4'),
        component_qc=dict(audited_components='4',processed_components='3',oversize_unsplit_components='1',largest_component_nodes='200'),
        component_node_summary=[dict(quantile='max',input_nodes='200',contracted_nodes='4')],
        selected_sample_summary=[dict(sample='s1',consensus_records='120',n_reads_sum='1200',joined_records='40',merged_records='80',single_records='0')],
        locus_occupancy=[dict(samples_with_data='4',loci='40',fraction_of_final_loci='0.571429')],
        clustering_parameters=dict(within_similarity='.95',across_similarity='.85',query_cov='.75',min_derep_size='5'),
        runtime=dict(alignment_mode='mafft'),outputs=dict(reference='/results/denovo_reference.fa'),
        stages=[dict(name='graph',status='COMPLETE',elapsed_seconds=2,sampled_peak_tree_rss_bytes=None,details=dict(components=1000))])


def line(text,label,value):
    return re.search(rf'^\s*{re.escape(label)}\s+{re.escape(value)}$',text,re.M)


def test_distinct_units_audit_subset_and_representation_denominator():
    data=example();before=copy.deepcopy(data);text=render_report(data)
    assert data==before
    assert line(text,'Input representative sequences','120')
    assert line(text,'Sequence-to-locus output memberships','90')
    assert line(text,'Reference sequences (loci)','70')
    assert line(text,'All similarity groups','1,000')
    assert line(text,'Audited groups processed normally','3')
    assert 'not all similarity groups' in text
    assert '40 / 70 (57.14%)' in text  # output loci, not 375 cohort members
    assert 'before mapping all samples' in text
    assert 'not raw-read retention' in text
    assert 'not resolved individual loci' in text
    assert text.index('# Reference Library') < text.index('# Time and Memory')
    assert line(text,'Running + waiting batch limit','N/A')  # old snapshots have no bound
    assert 'denovo_reference.fa' in text and '/results/denovo_reference.fa' not in text


def test_spacer_categories_are_disjoint():
    text=render_report(example())
    assert line(text,'Loci with a spacer: joined PE only','5')
    assert line(text,'Loci with a spacer: mixed PE','2')
    assert line(text,'Loci without an inserted spacer','63')
    assert 'partition the output loci' in text
    assert 'Of these' not in text


def test_no_alignment_does_not_present_uncomputed_categories_as_zero():
    data=example();data['runtime']['alignment_mode']='none'
    text=render_report(data)
    assert line(text,'Single-sequence loci','N/A')
    assert line(text,'Loci requiring alignment in this mode','N/A')
    assert line(text,'Loci with a spacer: joined PE only','5')


def test_zero_loci_and_no_audit_have_undefined_summaries():
    data=example();data['denovo_summary']['loci_written']='0'
    data['component_qc'].update(audited_components='0',processed_components='0',oversize_unsplit_components='0',largest_component_nodes='0')
    text=render_report(data)
    assert ' / 0 (N/A)' in text and ' (0.00%)' not in text
    assert line(text,'Largest audited group (input sequences)','N/A')
    assert '90th percentile' not in text
    assert 'size statistics are not applicable' in text


@pytest.mark.parametrize('status',['FAILED','RECOVERY_REQUIRED'])
def test_failure_and_recovery_do_not_advertise_a_complete_library(status):
    data=example();data['run']=dict(status=status,error=dict(type='OSError',message='disk failed'))
    text=render_report(data,failure_body='Preserve the journal and backups.')
    assert '# Reference Library' not in text
    assert 'disk failed' in text and 'Preserve the journal and backups.' in text
    assert '# Time and Memory' in text


def test_sample_name_controls_cannot_add_fake_report_rows():
    data=example();data['selected_sample_summary'][0]['sample']='a\n# Fake\trow'
    text=render_report(data)
    assert 'a\\n# Fake\\trow' in text
    assert '\n# Fake' not in text


def test_runtime_publisher_preserves_all_existing_json_fields(tmp_path):
    data=example();original=copy.deepcopy(data)
    js=tmp_path/'denovo.stats.json';txt=tmp_path/'denovo.stats.txt'
    js.write_text(json.dumps(data));txt.write_text('old presentation')
    run=dict(status='COMPLETE',started_utc='start',finished_utc='end',elapsed_seconds=2,
             attempt_directory='/attempt',output_directory='/results',invocation={})
    write_runtime_report(dict(run_stats=txt,run_stats_json=js),run,data['stages'])
    observed=json.loads(js.read_text())
    for key,value in original.items():assert observed[key]==value
    assert observed['definitions']['audit_scope']
    assert '# Reference Library' in txt.read_text()
