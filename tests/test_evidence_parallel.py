"""Lossless parsing, mask boundaries, ordering, and failed-worker checks."""
import numpy as np
import pytest
import multiprocessing

from ipyrad2.assembler.evidence_hdf5 import _column
from ipyrad2.assembler.evidence_parallel import parse_record_batch, parsed_records
from ipyrad2.assembler.rad_evidence_hdf5 import _mask_window
from ipyrad2.assembler.reporting import EvidenceCounts
from ipyrad2.utils.exceptions import IPyradError
from ipyrad2.utils.parallel.pool import ParallelJobError


@pytest.mark.parametrize('values,width', [
    (['1', '.', '-1', '+23', '2147483647'], 1),
    (['1,2,3', '.', '.,-1,9', '0,0,0'], 3),
    (['.', '.'], 6),
])
def test_fast_integer_parser_matches_scalar_definition(values, width):
    expected = [[-1]*width if v == '.' else [-1 if x == '.' else int(x) for x in v.split(',')]
                for v in values]
    np.testing.assert_array_equal(_column(values, width), expected)
    assert _column(values, width).dtype == np.dtype('i4')


@pytest.mark.parametrize('values,width', [
    (['1,2', '3,4,5,6'], 3),  # equal total values, invalid per-sample widths
    ([''], 1), (['1,,3'], 3),
    (['2147483648'], 1), (['-2'], 1), (['999999999999999999999999'], 1),
])
def test_fast_integer_parser_rejects_malformed_or_lossy_values(values, width):
    with pytest.raises((IPyradError, ValueError)):
        _column(values, width)


def native_records():
    rows = [
        'chr1\t1\t.\tA\t<*>\t.\t.\t.\tDP:AD:PL\t3:3,0:0,9,90\t0:0,0:.',
        'chr1\t2\t.\tA\tC,G\t.\t.\t.\tPL:GQ:AD:DP\t0,1,2,3,4,5:22:3,4,5:12\t.:.:.:.',
        'chr1\t3\t.\tA\t.\t.\t.\t.\tAD:DP\t3:3\t0:0',
        'chr1\t4\t.\tA\t<*>\t.\t.\t.\tDP:AD:PL\t3:3,0:0,9,90\t0:0,0:.',
    ]
    return [((0, i), row.split('\t')) for i, row in enumerate(rows)]


@pytest.mark.parametrize('workers', [1, 2])
def test_mixed_layouts_preserve_order_values_and_missing_fields(workers):
    rows = native_records()
    result = list(parsed_records(iter(rows), np.array([[0, 0, 4]]), workers, 'ERROR', batch_records=2))
    assert [r[0] for r in result] == [(0, i) for i in range(4)]
    for (_, raw), (_, metadata, data) in zip(rows, result):
        assert metadata == raw[:9]
        tags = raw[8].split(':')
        k = 1 if raw[4] == '.' else len(raw[4].split(',')) + 1
        for sample, cell in enumerate(raw[9:]):
            fields = dict(zip(tags, cell.split(':')))
            for key, width in [('DP', 1), ('GQ', 1), ('AD', k), ('PL', k*(k+1)//2)]:
                value = fields.get(key, '.')
                expected = [-1]*width if value == '.' else list(map(int, value.split(',')))
                assert data[key][sample].tolist() == expected


def test_worker_failure_propagates_and_source_errors_are_not_hidden():
    rows = native_records()
    rows[-1][1][-1] = '0:0,0:1,2'  # malformed PL in a later batch
    with pytest.raises(ParallelJobError):
        list(parsed_records(iter(rows), np.array([[0, 0, 4]]), 2, 'ERROR', batch_records=1))

    def source():
        yield native_records()[0]
        raise RuntimeError('late source failure')

    with pytest.raises(RuntimeError, match='late source failure'):
        list(parsed_records(source(), np.array([[0, 0, 1]]), 1, 'ERROR', batch_records=1))


def test_closing_ordered_parser_reaps_its_workers():
    before = {p.pid for p in multiprocessing.active_children()}
    def source():
        for position in range(50):
            _, fields = native_records()[0]
            yield (0, position), fields
    stream = parsed_records(source(), np.array([[0, 0, 50]]), 2, 'ERROR', batch_records=1)
    assert next(stream)[0] == (0, 0)
    stream.close()
    assert {p.pid for p in multiprocessing.active_children()} <= before


def test_excluded_records_are_consumed_but_not_reinterpreted():
    rows = native_records()
    rows[1][1][8] = 'DP:DP'  # excluded coordinate, as in the original writer
    result = list(parsed_records(iter(rows), np.array([[0, 0, 1], [0, 2, 4]]), 1, 'ERROR'))
    assert [r[0][1] for r in result] == [0, 2, 3]
    with pytest.raises(IPyradError, match='Duplicate'):
        parse_record_batch([rows[1]])


def test_window_masks_match_direct_interval_membership():
    entries = [dict(field_idx=9, cursor=0, intervals=[(0, 3), (2, 10), (4, 6), (17, 19)]),
               dict(field_idx=10, cursor=0, intervals=[(3, 5), (11, 12)]),
               dict(field_idx=11, cursor=0, intervals=[(8, 9), (15, 16)])]
    for start, end in [(0, 4), (4, 7), (8, 13), (15, 20)]:
        observed = _mask_window(entries, start, end, 2)
        expected = [[any(e['field_idx'] in (9+i, 11) and any(a <= p < b for a, b in e['intervals'])
                         for e in entries) for p in range(start, end)] for i in range(2)]
        np.testing.assert_array_equal(observed, expected)


def test_ragged_summaries_match_recordwise_summaries():
    records = parse_record_batch(native_records())
    old, new = EvidenceCounts(['a', 'b'], ['PASS', 'Missing']), EvidenceCounts(['a', 'b'], ['PASS', 'Missing'])
    for _, _, data in records:
        data['status'] = np.array([[0], [1]], dtype='u1')
        old.add_record(data)
    block = {key: np.concatenate([r[2][key] for r in records], axis=1) for key in records[0][2]}
    starts = lambda key: np.r_[0, np.cumsum([r[2][key].shape[1] for r in records[:-1]])].astype('intp')
    new.add_ragged(block, starts('AD'), starts('PL'))
    assert new.summary() == old.summary()


@pytest.mark.parametrize('bad', ['1junk', '1+2', '--2', '2147483648', '-2'])
def test_complete_numeric_batch_rejects_invalid_values(bad):
    row = native_records()[0]
    row[1][-1] = '0:0,0:0,0,0'
    row[1][-2] = bad + ':3,0:0,9,90'
    with pytest.raises(IPyradError):
        parse_record_batch([row])
