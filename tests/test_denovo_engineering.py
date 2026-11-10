"""Workflow preservation, deterministic sparse graphs, and runtime accounting."""
import csv
from itertools import combinations
import json
import random
import subprocess
import sys
from pathlib import Path

import pytest

from ipyrad2.denovo import denovo as runner
from ipyrad2.denovo import graph, graph_split, align
from ipyrad2.denovo.common import reverse_complement
from ipyrad2.denovo.runtime import OutputTransaction, StageRecorder, process_tree_rss
from ipyrad2.denovo import runtime
from ipyrad2.utils.exceptions import IPyradError


def params(tmp):
    fastq = tmp/'sample.fastq'
    fastq.write_text('@r\nACGT\n+\nIIII\n')
    return dict(fastqs=[fastq], outdir=tmp/'OUT', within_similarity=.95, across_similarity=.85,
                query_cov=.75, min_derep_size=2, min_length=35, min_merge_overlap=20,
                max_merge_diffs=4, delim_str=None, delim_idx=1, allow_reverse_complement=False,
                cores=1, threads=1, no_alignment=True, force=True, imap=None,
                use_all_samples=True, keep_intermediates=False, log_level='ERROR')


@pytest.mark.parametrize('failure', ['missing_input', 'missing_binary', 'mixed_layout', 'bad_imap'])
def test_preflight_failure_leaves_previous_results(tmp_path, monkeypatch, failure):
    args = params(tmp_path)
    out = args['outdir'];out.mkdir()
    old = out/'denovo_reference.fa';old.write_bytes(b'>old\nACGT\n')
    if failure == 'missing_input':
        args['fastqs'] = [tmp_path/'absent.fq']
    elif failure == 'missing_binary':
        monkeypatch.setattr(runner, 'BIN_VSEARCH', str(tmp_path/'absent'))
    elif failure == 'mixed_layout':
        monkeypatch.setattr(runner, 'get_name_to_fastq_dict', lambda *a: {'a':(args['fastqs'][0],None), 'b':(args['fastqs'][0],args['fastqs'][0])})
    else:
        args.update(imap=tmp_path/'absent.imap', use_all_samples=False)
    with pytest.raises((IPyradError, FileNotFoundError)):
        runner.run_denovo(**args)
    assert old.read_bytes() == b'>old\nACGT\n'
    assert sorted(p.name for p in out.iterdir()) == ['denovo_reference.fa']


def test_native_failure_keeps_old_run_and_new_failure_evidence(tmp_path, monkeypatch):
    args = params(tmp_path);out=args['outdir'];out.mkdir()
    old=out/'denovo.stats.json';old.write_text('old successful report')
    monkeypatch.setattr(runner, '_validate_required_binaries', lambda: None)
    monkeypatch.setattr(runner, 'tool_version', lambda path: {'path':path,'version':'test'})
    monkeypatch.setattr(runner, 'run_with_pool', lambda *a, **k: (_ for _ in ()).throw(RuntimeError('native test failure')))
    with pytest.raises(RuntimeError, match='native test failure'):
        runner.run_denovo(**args)
    assert old.read_text()=='old successful report'
    attempt=next(out.glob('.denovo-attempt-*'))
    report=json.loads((attempt/'denovo.stats.json').read_text())
    assert report['run']['status']=='FAILED'
    assert report['stages'][-1]['name']=='within_samples'
    assert report['stages'][-1]['status']=='FAILED'
    assert not (out/'.denovo.lock').exists()


def test_transaction_rolls_back_partial_publication(tmp_path, monkeypatch):
    out=tmp_path/'out';out.mkdir();(out/'a').write_text('old a');(out/'b').write_text('old b')
    with OutputTransaction(out, ['a','b'], True) as tx:
        (tx.attempt/'a').write_text('new a');(tx.attempt/'b').write_text('new b')
        rename=Path.rename
        def fail_b(self, target):
            if self==tx.attempt/'b' and target==out/'b':
                raise OSError('publication failed')
            return rename(self,target)
        monkeypatch.setattr(Path, 'rename', fail_b)
        with pytest.raises(OSError, match='publication failed'):
            tx.publish()
        assert (out/'a').read_text()=='old a'
        assert (out/'b').read_text()=='old b'
        assert (tx.attempt/'a').read_text()=='new a'
        assert json.loads((tx.attempt/'PUBLICATION.json').read_text())['status']=='ROLLED_BACK'


def test_transaction_excludes_concurrent_writer_and_preserves_other_files(tmp_path):
    (tmp_path/'unrelated').write_text('keep')
    with OutputTransaction(tmp_path,['a'],False) as tx:
        with pytest.raises(IPyradError, match='locked'):
            with OutputTransaction(tmp_path,['a'],False):
                pass
        (tx.attempt/'a').write_text('new')
        tx.publish()
    assert (tmp_path/'a').read_text()=='new'
    assert (tmp_path/'unrelated').read_text()=='keep'
    assert not list(tmp_path.glob('.denovo*'))


def test_transaction_rolls_back_when_final_journal_write_fails(tmp_path, monkeypatch):
    (tmp_path/'a').write_text('old')
    with OutputTransaction(tmp_path, ['a'], True) as tx:
        (tx.attempt/'a').write_text('new')
        write = runtime.atomic_text
        def fail_commit(path, text):
            if path.name == 'PUBLICATION.json' and json.loads(text)['status'] == 'PUBLISHED':
                raise OSError('journal commit failed')
            return write(path, text)
        monkeypatch.setattr(runtime, 'atomic_text', fail_commit)
        with pytest.raises(OSError, match='journal commit failed'):
            tx.publish()
        assert (tmp_path/'a').read_text() == 'old'
        assert (tx.attempt/'a').read_text() == 'new'
        assert json.loads((tx.attempt/'PUBLICATION.json').read_text())['status'] == 'ROLLED_BACK'


def test_failed_rollback_retains_lock_and_recovery_journal(tmp_path, monkeypatch):
    (tmp_path/'a').write_text('old a')
    (tmp_path/'b').write_text('old b')
    with OutputTransaction(tmp_path, ['a', 'b'], True) as tx:
        (tx.attempt/'a').write_text('new a')
        (tx.attempt/'b').write_text('new b')
        rename = Path.rename
        def fail_publish_and_restore(self, target):
            if self == tx.attempt/'b' or self == tx.attempt/'_previous/a':
                raise OSError('simulated filesystem failure')
            return rename(self, target)
        monkeypatch.setattr(Path, 'rename', fail_publish_and_restore)
        with pytest.raises(IPyradError, match='mixed generations'):
            tx.publish()
    assert tx.lock.exists()
    assert (tx.attempt/'_previous/a').read_text() == 'old a'
    assert json.loads((tx.attempt/'PUBLICATION.json').read_text())['status'] == 'RECOVERY_REQUIRED'
    with pytest.raises(IPyradError, match='locked'):
        with OutputTransaction(tmp_path, ['a', 'b'], True):
            pass


def test_no_force_preserves_dangling_output_symlink(tmp_path):
    output = tmp_path/'denovo_reference.fa'
    output.symlink_to(tmp_path/'absent')
    with pytest.raises(IPyradError, match='already exist'):
        with OutputTransaction(tmp_path, [output.name], False):
            pass
    assert output.is_symlink()
    assert not (tmp_path/'.denovo.lock').exists()


def test_sparse_threshold_matches_legacy_order_for_random_subsets():
    for seed in range(60):
        rng=random.Random(seed);allnodes=[f'n{i}' for i in range(35)]
        nodes=set(rng.sample(allnodes,25));edges={p:(rng.choice([.8,.9,.95]),1.) for p in combinations(sorted(allnodes),2) if rng.random()<.15}
        expected={p:edges[p] for p in combinations(sorted(nodes),2) if p in edges and edges[p][0]>=.9}
        assert list(graph_split.subset_edges_above_threshold(nodes,edges,.9).items())==list(expected.items())


def test_sparse_contraction_matches_all_group_pair_scan():
    for seed in range(40):
        rng=random.Random(seed);nodes={f'n{i}' for i in range(35)};order={n:i for i,n in enumerate(sorted(nodes))}
        # Repeated samples and exact sequences exercise contraction as well as ties.
        records={n:graph.SummaryRecord(sample=f's{i%8}', cluster_id=i,
            n_unique=1,n_reads=5,length=10,cluster_length=10,record_type='single',
            cluster_sequence=['ACGTACGTAA','ACGTACGTCC'][i%2])
            for i,n in enumerate(sorted(nodes))}
        edges={p:(rng.choice([.85,.95,.99]),.9) for p in combinations(sorted(nodes),2) if rng.random()<.15}
        got=graph._build_contracted_component(1,nodes,records,edges,order,.95)
        expected={}
        ordered=sorted(got.node_members,key=got.node_order.get)
        for a,b in combinations(ordered,2):
            scores=[edges[tuple(sorted((x,y)))] for x in got.node_members[a] for y in got.node_members[b]
                    if tuple(sorted((x,y))) in edges]
            if scores:expected[(a,b)]=max(scores)
        assert list(got.edges.items())==list(expected.items())


def test_memory_is_stage_specific_and_missing_is_null():
    values=iter([100,200,None,None])
    recorder=StageRecorder(interval=999,sampler=lambda:next(values))
    with recorder.stage('first'):pass
    with pytest.raises(RuntimeError):
        with recorder.stage('failed'):
            raise RuntimeError('test')
    a,b=recorder.stages
    assert a['sampled_peak_tree_rss_bytes']==200 and a['rss_samples']==2
    assert b['sampled_peak_tree_rss_bytes'] is None and b['unavailable_rss_samples']==2
    assert b['status']=='FAILED' and b['elapsed_seconds']>=0


@pytest.mark.skipif(not Path('/proc/self/statm').exists(),reason='Linux process-tree sampler')
def test_memory_includes_live_native_child():
    before=process_tree_rss()
    child=subprocess.Popen([sys.executable,'-u','-c',"import time; data=bytearray(32*1024*1024); print('ready',flush=True); time.sleep(20)"],stdout=subprocess.PIPE,text=True)
    try:
        assert child.stdout.readline().strip()=='ready'
        during=process_tree_rss()
        assert before is not None and during is not None
        assert during>before+25*1024*1024
    finally:
        child.terminate();child.wait(timeout=5)


def test_orientation_forest_is_deterministic_and_counts_conflicts():
    nodes=['a','b','c'];order={x:i for i,x in enumerate(nodes)}
    edges={('a','b'):(1.,1.),('b','c'):(.99,1.),('a','c'):(.9,1.)}
    strands={('a','b'):'-',('b','c'):'-',('a','c'):'-'}
    signs,conflicts=graph._orient_seeds(nodes,edges,strands,order)
    assert signs=={'a':0,'b':1,'c':0} and conflicts==1
    assert graph._orient_seeds(nodes,dict(reversed(list(edges.items()))),strands,order)==(signs,conflicts)


def test_orientation_path_compression_preserves_parity_in_long_chains():
    nodes=[f'n{i:03}' for i in range(300)]
    order={node:i for i,node in enumerate(nodes)}
    # Process higher-index pairs first, forcing a long chain of earlier roots.
    edges={(a,b):(.8+i/2000,1.) for i,(a,b) in enumerate(zip(nodes,nodes[1:]))}
    strands={edge:'-' for edge in edges}
    signs,conflicts=graph._orient_seeds(nodes,edges,strands,order)
    assert signs=={node:i%2 for i,node in enumerate(nodes)}
    assert conflicts==0


def test_reverse_orientation_swaps_unequal_joined_arms(tmp_path):
    record=align.SummaryRecord('b','s','joined','AAACCCCC','AAA','CCCCC')
    mapping=tmp_path/'map.tsv';mapping.write_text('locus\tcore\torientation\n1\tb\t-\n')
    members=list(align._iter_locus_members(mapping,{'b':record}))[0][2]
    assert members[0].cluster_sequence=='GGGGGTTT'
    assert members[0].left_arm=='GGGGG' and members[0].right_arm=='TTT'


@pytest.mark.skipif(not Path(runner.BIN_VSEARCH).exists(),reason='native VSEARCH unavailable')
def test_native_across_search_recovers_reverse_complements(tmp_path):
    rng=random.Random(17);seq=''.join(rng.choice('ACGT') for _ in range(120))
    table=tmp_path/'summary.tsv';table.write_text('seed\tcluster_sequence\na\t'+seq+'\nb\t'+reverse_complement(seq)+'\n')
    runner.vsearch_cluster_across(tmp_path,table,.85,.75,1)
    assert (tmp_path/'global_hits.uc.tsv').read_text()==''
    runner.vsearch_cluster_across(tmp_path,table,.85,.75,1,allow_reverse_complement=True)
    rows=(tmp_path/'global_hits.uc.tsv').read_text().splitlines()
    assert len(rows)==2 and all(row.split('\t')[3]=='-' for row in rows)


def test_orientation_survives_graph_and_joined_consensus(tmp_path):
    work=tmp_path/'_denovo_work';work.mkdir()
    left,right='ACGTTGCA','GATCCAGTACGTA'
    sequence=left+right
    columns=['sample','cluster_id','seed','length','cluster_length','n_unique','n_reads','record_type','cluster_sequence','arm_boundary']
    with (work/'concat.summary.tsv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,delimiter='\t');writer.writeheader()
        for name,seq,boundary in [('a',sequence,len(left)),('b',reverse_complement(sequence),len(right))]:
            writer.writerow(dict(sample=name,cluster_id=0,seed=name+';J1',length=len(seq)+24,
                                 cluster_length=len(seq),n_unique=1,n_reads=10,record_type='joined',
                                 cluster_sequence=seq,arm_boundary=boundary))
    (work/'global_hits.uc.tsv').write_text('a;J1\tb;J1\t100\t-\t100\t21\t21\nb;J1\ta;J1\t100\t-\t100\t21\t21\n')
    graph.make_global_tables(work,allow_reverse_complement=True)
    mapping=tmp_path/'denovo.loci.mapping.tsv'
    rows=list(csv.DictReader(mapping.open(),delimiter='\t'))
    assert [(row['sample'],row['orientation']) for row in rows]==[('a','+'),('b','-')]
    out=tmp_path/'ref.fa'
    align.write_ordered_consensus_stream_to_file(mapping,work/'concat.summary.tsv',out,
        mafft_binary='unused',alignment_mode='none')
    assert out.read_text().splitlines()[1]==left+'N'*50+right
    audit=json.loads((tmp_path/'denovo.audit/orientation.json').read_text())
    assert audit['conflicting_weaker_edges']==0


def test_report_receives_full_denovo_invocation(monkeypatch):
    from types import SimpleNamespace
    import ipyrad2.cli.cli_main as cli
    argv=['denovo','-d',*[f'data path/{i}.fastq.gz' for i in range(50)],'-o','out']
    args=cli.setup_parsers().parse_args(argv);seen={}
    monkeypatch.setattr(cli.sys,'argv',['/env/bin/ipyrad2',*argv])
    monkeypatch.setattr(cli.importlib,'import_module',lambda *a,**k:SimpleNamespace(run_denovo=lambda **kw:seen.update(kw)))
    cli.run_subcommand(args,_exit=False)
    assert seen['invocation']['argv']==['/env/bin/ipyrad2',*argv]
