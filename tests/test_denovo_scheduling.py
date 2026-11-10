"""Transport batching must preserve graph decisions, order and failure identity."""
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from ipyrad2.denovo import graph, runtime


def transport_job(key, nodes=1, bases=1, edges=0):
    return key, (None, dict(component=tuple(range(nodes)),
        summary_records={0: SimpleNamespace(cluster_sequence='A' * bases)},
        component_edges=dict.fromkeys(range(edges))))


@pytest.mark.parametrize('budget,small,large', [
    ('COMPONENT_BATCH_NODES', dict(nodes=2), dict(nodes=7)),
    ('COMPONENT_BATCH_BASES', dict(bases=2), dict(bases=7)),
    ('COMPONENT_BATCH_EDGES', dict(edges=2), dict(edges=7)),
])
def test_component_batches_bound_payload_and_keep_large_component_intact(monkeypatch, budget, small, large):
    monkeypatch.setattr(graph, budget, 5)
    jobs = [transport_job(i, **(large if i == 3 else small)) for i in range(7)]
    batches = list(graph._iter_component_batches(iter(jobs), 32))
    groups = [[job[0] for job in kwargs['jobs']] for _, (_, kwargs) in batches]
    assert groups == [[0, 1], [2], [3], [4, 5], [6]]
    assert [key for key, _ in batches] == [0, 2, 3, 4, 6]


def test_batches_are_lazy_and_limit_component_count():
    consumed = []
    def jobs():
        for i in range(100):
            consumed.append(i)
            yield transport_job(i)
    batches = graph._iter_component_batches(jobs(), 3)
    _, (_, kwargs) = next(batches)
    assert [j[0] for j in kwargs['jobs']] == [0, 1, 2]
    assert len(consumed) <= 4  # bounded one-job lookahead
    with pytest.raises(ValueError, match='batch_size'):
        list(graph._iter_component_batches([], 0))


def test_batch_error_identifies_component_and_preserves_original_cause():
    def fails():
        raise ValueError('invalid original evidence')
    with pytest.raises(RuntimeError, match='graph component 83 failed') as err:
        graph._process_component_batch(((83, (fails, {})),))
    assert isinstance(err.value.__cause__, ValueError)
    assert str(err.value.__cause__) == 'invalid original evidence'


def write_graph(work):
    work.mkdir()
    header = ['sample', 'cluster_id', 'seed', 'length', 'cluster_length',
              'n_unique', 'n_reads', 'record_type', 'cluster_sequence', 'arm_boundary']
    hits = []
    with (work/'concat.summary.tsv').open('w') as stream:
        w = csv.DictWriter(stream, fieldnames=header, delimiter='\t');w.writeheader()
        for component in range(137):
            count = 3 if component % 7 == 0 else 1
            for member in range(count):
                seed = f'c{component}_{member}'
                w.writerow(dict(sample=f's{member//2}', cluster_id=component,
                    seed=seed, length=8, cluster_length=8, n_unique=1,
                    n_reads=5+member, record_type='single', cluster_sequence='ACGTACGT', arm_boundary=0))
                if member:
                    hits.append(f'c{component}_0\t{seed}\t99\t+\t100\t8\t8\n')
    (work/'global_hits.uc.tsv').write_text(''.join(hits))


def audit_records(path):
    with path.open() as f:
        return sorted(csv.DictReader(f, delimiter='\t'), key=lambda r: int(r['component_id']))


def test_batched_parallel_outputs_match_unbatched_and_serial(tmp_path, monkeypatch):
    roots = [tmp_path/n for n in ['single', 'unbatched', 'batched']]
    summaries = []
    for root, cores, size in zip(roots, [1, 2, 2], [32, 1, 32]):
        root.mkdir();work = root/'work';write_graph(work)
        monkeypatch.setattr(graph, 'COMPONENT_BATCH_SIZE', size)
        summaries.append(graph.make_global_tables(work, cores=cores))
    assert summaries[0] == summaries[1] == summaries[2]
    names = {str(p.relative_to(roots[0])) for p in roots[0].rglob('*') if p.is_file()}
    for root in roots[1:]:
        assert {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()} == names
        for name in names:
            if name == 'denovo.audit/components.summary.tsv':
                assert audit_records(root/name) == audit_records(roots[0]/name)
            else:
                assert (root/name).read_bytes() == (roots[0]/name).read_bytes()
    p = summaries[2].performance
    assert p['components'] == p['components_completed'] == 137
    assert p['audit_components'] == 20 and p['workers'] == 2
    assert p['batches_submitted'] < summaries[1].performance['batches_submitted'] / 4
    assert sum(p['phases_seconds'].values()) <= p['elapsed_seconds']
    assert all(v >= 0 for k, v in p.items()
               if k != 'phases_seconds' and (k.endswith('seconds') or k.endswith('seconds_sum')))
    assert all(v >= 0 for v in p['phases_seconds'].values())
    assert p['peak_buffered_components'] >= 1
    json.dumps(p, allow_nan=False)


def test_open_audit_summary_is_reused_without_changing_fields(tmp_path, monkeypatch):
    import builtins
    root = tmp_path/'run';root.mkdir();work = root/'work';write_graph(work)
    opens = []
    original = builtins.open
    def observed(file, *args, **kwargs):
        if Path(file).name == 'components.summary.tsv':
            opens.append(args[0] if args else kwargs.get('mode', 'r'))
        return original(file, *args, **kwargs)
    monkeypatch.setattr(builtins, 'open', observed)
    summary = graph.make_global_tables(work, cores=1)
    assert summary.performance['audit_components'] == 20
    assert opens == ['wt', 'at']  # header creation, one persistent append handle
    assert len(audit_records(root/'denovo.audit/components.summary.tsv')) == 20


def test_sparse_buckets_transfer_ownership_without_copying_or_reordering():
    components = [(1, {'a'}), (2, {'b', 'c', 'd'}), (3, {'e', 'f'})]
    edges = {('c', 'd'): (.98, .9), ('b', 'c'): (.99, 1.), ('e', 'f'): (.97, .8)}
    buckets = graph._bucket_component_edges(components, edges)
    assert set(buckets) == {2, 3}
    originals = dict(buckets)
    records = {s: SimpleNamespace(cluster_sequence='ACGT') for s in 'abcdef'}
    jobs = graph._iter_component_jobs(components, summary_records=records,
        component_edges=buckets, raw_oversize_component_ids=set(),
        seed_order={s:i for i,s in enumerate('abcdef')}, within_similarity=.95,
        max_component_nodes=40)
    _, (_, singleton) = next(jobs)
    assert singleton['component_edges'] == {} and set(buckets) == {2, 3}
    _, (_, second) = next(jobs)
    assert second['component_edges'] is originals[2]
    assert list(second['component_edges']) == [('c', 'd'), ('b', 'c')]
    assert set(buckets) == {3}
    _, (_, third) = next(jobs)
    assert third['component_edges'] is originals[3] and not buckets
    assert edges == {('c', 'd'): (.98, .9), ('b', 'c'): (.99, 1.), ('e', 'f'): (.97, .8)}


def test_graph_details_in_runtime_report(tmp_path):
    root = tmp_path/'run';root.mkdir();work = root/'work';write_graph(work)
    g = graph.make_global_tables(work, cores=1).performance
    outputs = {'run_stats': root/'stats.txt', 'run_stats_json': root/'stats.json'}
    run = dict(status='COMPLETE', started_utc='start', finished_utc='finish', elapsed_seconds=1.,
               attempt_directory=str(root), output_directory=str(root), invocation={})
    stages = [dict(name='graph', status='COMPLETE', elapsed_seconds=1.,
                   sampled_peak_tree_rss_bytes=None, details=g)]
    runtime.write_runtime_report(outputs, run, stages)
    text = outputs['run_stats'].read_text()
    assert 'Execution Details: Group Processing' in text and '137 /' in text
    assert 'overlap' in text and 'must not be summed' in text
    assert json.loads(outputs['run_stats_json'].read_text())['stages'][0]['details'] == g


def test_stalled_first_graph_batch_bounds_finished_work(tmp_path, monkeypatch):
    from concurrent.futures import Future
    from ipyrad2.utils.parallel import pool
    submitted = []
    delivered = []
    waits = []
    class ControlledGraphPool(pool._ManagedProcessPool):
        def __init__(self, **kwargs):
            pass
        def submit(self, func, kwargs):
            future = Future()
            submitted.append((future, func, kwargs))
            return future
        def collect(self):
            pass
        def close(self, **kwargs):
            pass
        def abort(self, *args, **kwargs):
            raise AssertionError('successful graph replay must not abort')
    def delayed_wait(futures, **kwargs):
        outstanding = len(submitted) - len(delivered)
        waits.append(outstanding)
        assert outstanding <= 4  # two workers, including completed results
        for index, (future, func, job) in enumerate(submitted):
            if future in futures and (index != 0 or len(waits) >= 4):
                future.set_result(func(**job))
        return {f for f in futures if f.done()}, set()
    original = graph.run_with_pool_iter
    def observed(jobs, **kwargs):
        assert kwargs['ordered'] is True and kwargs['max_inflight'] == 4
        for key, result in original(jobs, **kwargs):
            delivered.append(key)
            yield key, result
    monkeypatch.setattr(pool, '_ManagedProcessPool', ControlledGraphPool)
    monkeypatch.setattr(pool, 'wait', delayed_wait)
    monkeypatch.setattr(graph, 'run_with_pool_iter', observed)
    root = tmp_path/'run';root.mkdir();work = root/'work';write_graph(work)
    result = graph.make_global_tables(work, cores=2)
    assert waits[:4] == [4, 4, 4, 4]
    assert result.performance['components_completed'] == 137
    assert result.performance['peak_buffered_components'] == 1
    assert result.performance['max_outstanding_batches'] == 4
