"""SE reconciliation and spacer accounting must follow completed output."""
import csv
from itertools import combinations

import pytest

from ipyrad2.denovo import align, graph


def member(locus, core, sample, group, kind, seq, boundary=None, form='spaced'):
    boundary = len(seq) if boundary is None else boundary
    return align.LocusMember(locus, f'locus_{locus}', group, core, sample, kind,
                             seq, seq[:boundary], seq[boundary:], 'mixed', form)


def tables(tmp_path, members):
    summary = tmp_path / 'summary.tsv'
    mapping = tmp_path / 'mapping.tsv'
    with summary.open('w') as f:
        f.write('seed\tsample\trecord_type\tcluster_sequence\tarm_boundary\n')
        for m in members:
            f.write(f'{m.core}\t{m.sample}\t{m.record_type}\t{m.cluster_sequence}\t{len(m.left_arm)}\n')
    with mapping.open('w') as f:
        f.write('locus\tlocus_name\tcore\tcontract_group\toutput_form\n')
        for m in members:
            f.write(f'{m.locus_id}\t{m.locus_name}\t{m.core}\t{m.contract_group}\t{m.output_form}\n')
    return mapping, summary


@pytest.mark.parametrize('form', ['spaced', 'stripped', ''])
@pytest.mark.parametrize('duplicates', [1, 2, 4])
def test_se_consensus_gives_each_reconciled_group_one_vote(monkeypatch, form, duplicates):
    monkeypatch.setattr(align, 'mafft_align_one', lambda record, **kw: list(record))
    members = [member(1, f'a{i}', 'A', 'A', 'single', 'AAAA', form=form) for i in range(duplicates)]
    members += [member(1, 'b', 'B', 'B', 'single', 'GGGG', form=form),
                member(1, 'c', 'C', 'C', 'single', 'GGGG', form=form)]
    assert align.worker_build_consensus(1, members, 'test')[2:] == ('GGGG', False)


def test_graph_reconciled_se_output_is_stripped(tmp_path, monkeypatch):
    work = tmp_path / 'work'; work.mkdir()
    fields = ['sample','cluster_id','seed','length','cluster_length','n_unique','n_reads','record_type','cluster_sequence','arm_boundary']
    seq_a = 'ACGT' * 20
    seq_g = seq_a[:20] + 'G' + seq_a[21:]
    with (work/'concat.summary.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t'); w.writeheader()
        for i, (sample, seq) in enumerate([('A',seq_a),('A',seq_a),('B',seq_g),('C',seq_g)]):
            w.writerow(dict(sample=sample, cluster_id=i, seed=f'{sample};S{i}', length=80,
                cluster_length=80, n_unique=1, n_reads=10, record_type='single', cluster_sequence=seq, arm_boundary=80))
    cores = ['A;S0','A;S1','B;S2','C;S3']
    (work/'global_hits.uc.tsv').write_text(''.join(f'{a}\t{b}\t98.75\t+\t100\t80\t80\n' for a,b in combinations(cores,2)))
    g = graph.make_global_tables(work, cores=1)
    with (tmp_path/'denovo.loci.mapping.tsv').open() as f: rows = list(csv.DictReader(f,delimiter='\t'))
    assert g.loci_written == 1 and len(rows) == 4
    assert {r['output_form'] for r in rows} == {'stripped'}
    assert len({r['contract_group'] for r in rows if r['sample']=='A'}) == 1
    monkeypatch.setattr(align, 'mafft_align_one', lambda record, **kw: list(record))
    result = align.write_ordered_consensus_stream_to_file(tmp_path/'denovo.loci.mapping.tsv',
        work/'concat.summary.tsv', tmp_path/'out.fa', 'test')
    assert (tmp_path/'out.fa').read_text().splitlines()[1] == seq_g
    assert (result.joined_spacer_loci,result.mixed_reconciled_spacer_loci,result.stripped_output_loci) == (0,0,1)


@pytest.mark.parametrize('mode', ['mafft','none'])
@pytest.mark.parametrize('cores', [1,2])
def test_counts_partition_written_loci_and_ignore_natural_ns(tmp_path, monkeypatch, mode, cores):
    monkeypatch.setattr(align, 'mafft_align_one', lambda record, **kw: list(record))
    members = [member(1,'a','a','a','single','AAAA'+'N'*50+'TTTT'),
               member(2,'b','b','b','joined','AAAATTTT',4),
               member(3,'c','c','c','joined','AAAATTTT',4),
               member(3,'d','d','d','merged','AAAATTTT'),
               member(4,'e','e','e','merged','AAAATTTT'),
               member(5,'f','f','f','joined','AAAA',4)]
    mapping, summary = tables(tmp_path,members)
    result = align.write_ordered_consensus_stream_to_file(mapping, summary, tmp_path/'out.fa',
        'test', cores=cores, alignment_mode=mode)
    seqs = (tmp_path/'out.fa').read_text().splitlines()[1::2]
    assert seqs == ['AAAA'+'N'*50+'TTTT','AAAA'+'N'*50+'TTTT','AAAA'+'N'*50+'TTTT','AAAATTTT','AAAA']
    assert result.total_loci == 5
    assert (result.joined_spacer_loci,result.mixed_reconciled_spacer_loci,result.stripped_output_loci) == (1,1,3)
    assert result.output_spacer_length == 50


def test_mixed_pe_worker_fallback_counts_no_spacer(tmp_path, monkeypatch):
    monkeypatch.setattr(align, 'mafft_align_one', lambda record, **kw: list(record))
    monkeypatch.setattr(align, '_infer_mixed_locus_boundary', lambda *a: None)
    mapping, summary = tables(tmp_path,[member(1,'a','a','a','joined','AAAATTTT',4),
                                       member(1,'b','b','b','merged','AAAATTTT')])
    result = align.write_ordered_consensus_stream_to_file(mapping,summary,tmp_path/'out.fa','test',cores=2)
    assert (tmp_path/'out.fa').read_text() == '>locus_1\nAAAATTTT\n'
    assert result.mafft_required_loci == 1
    assert (result.joined_spacer_loci,result.mixed_reconciled_spacer_loci,result.stripped_output_loci) == (0,0,1)
    assert result.output_spacer_length == 0


@pytest.mark.parametrize('mode', ['mafft','none'])
def test_zero_length_padding_is_not_an_inserted_spacer(tmp_path, mode):
    mapping, summary = tables(tmp_path,[member(1,'a','a','a','joined','AAAATTTT',4)])
    result = align.write_ordered_consensus_stream_to_file(mapping,summary,tmp_path/'out.fa',
        'unused',alignment_mode=mode,output_spacer_len=0)
    assert (tmp_path/'out.fa').read_text() == '>locus_1\nAAAATTTT\n'
    assert (result.joined_spacer_loci,result.mixed_reconciled_spacer_loci,result.stripped_output_loci) == (0,0,1)
    assert result.output_spacer_length == 0


def test_writer_failure_does_not_return_successful_counts(tmp_path, monkeypatch):
    mapping, summary = tables(tmp_path,[member(1,'a','a','a','single','AAAA')])
    def fail(*a,**kw): raise OSError('write failed')
    monkeypatch.setattr(align,'_write_planned_locus',fail)
    monkeypatch.setattr(align,'_build_alignment_run_summary',lambda *a,**kw: pytest.fail('must not finalize failed output'))
    with pytest.raises(OSError,match='write failed'):
        align.write_ordered_consensus_stream_to_file(mapping,summary,tmp_path/'out.fa','unused')
