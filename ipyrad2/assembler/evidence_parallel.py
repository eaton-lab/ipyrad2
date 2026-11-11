"""Bounded, ordered parsing of native evidence; workers never open HDF5 files."""
from collections import defaultdict
from functools import lru_cache
import re

import numpy as np

from .evidence_hdf5 import _column
from ..utils.exceptions import IPyradError
from ..utils.parallel import run_with_pool_iter


@lru_cache(maxsize=128)
def _numeric_sample_pattern(widths):
    integer = r'[+-]?[0-9]+'
    fields = [rf'{integer}(?:,{integer}){{{width-1}}}' for width in widths]
    return re.compile('^' + ':'.join(fields) + '$', re.MULTILINE)


def _numeric_batch(records, indices, tags, widths, samples):
    """Fast path for complete numeric native FORMAT cells (no GT or missing tokens).

    Validate each sample's field boundaries before parsing the contiguous text
    in NumPy. Equal total cardinality cannot hide a malformed individual cell.
    Missing/unsupported layouts retain the general field-by-field parser.
    """
    if any(key not in widths for key in tags):
        return None
    text = '\n'.join(cell for index in indices for cell in records[index][1][9:])
    if '.' in text:
        return None
    shape = tuple(widths[key] for key in tags)
    if _numeric_sample_pattern(shape).sub('', text) != '\n'*(len(indices)*samples-1):
        raise IPyradError('Native evidence has an invalid integer field or allele/genotype cardinality.')
    values = np.fromstring(text.replace(':', ',').replace('\n', ','), sep=',', dtype='i8')
    if values.size != len(indices)*samples*sum(shape) or np.any(values < -1) or np.any(values > np.iinfo(np.int32).max):
        raise IPyradError('Native evidence exceeds the lossless int32 representation.')
    values = values.astype('i4').reshape(len(indices), samples, sum(shape))
    arrays, start = {}, 0
    for key, width in zip(tags, shape):
        arrays[key] = values[:, :, start:start+width]
        start += width
    for key, width in widths.items():
        if key not in arrays:
            arrays[key] = np.full((len(indices), samples, width), -1, dtype='i4')
    return arrays


def parse_record_batch(records, packed=False):
    """Parse equal FORMAT/cardinality layouts together without reordering records."""
    groups = defaultdict(list)
    result = [None] * len(records)
    blocks = []
    for index, (coord, fields) in enumerate(records):
        alleles = 1 if fields[4] == '.' else 2 + fields[4].count(',')
        groups[fields[8], alleles].append(index)
    for (layout, width), indices in groups.items():
        tags = layout.split(':')
        if len(tags) != len(set(tags)):
            raise IPyradError('Duplicate native FORMAT keys.')
        samples = len(records[indices[0]][1]) - 9
        widths = dict(DP=1, GQ=1, AD=width, PL=width*(width+1)//2)
        arrays = _numeric_batch(records, indices, tags, widths, samples)
        if arrays is None:
            fields = [cell.split(':') for index in indices for cell in records[index][1][9:]]
            arrays = {}
            for key, count in widths.items():
                tag = tags.index(key) if key in tags else -1
                if tag < 0:
                    arrays[key] = np.full((len(indices), samples, count), -1, dtype='i4')
                else:
                    arrays[key] = _column([v[tag] if tag < len(v) else '.' for v in fields], count).reshape(len(indices), samples, count)
        if packed:
            blocks.append((indices, arrays))
        else:
            for row, index in enumerate(indices):
                coord, raw = records[index]
                result[index] = (coord, raw[:9], {key: values[row] for key, values in arrays.items()})
    if packed:
        return [(coord, raw[:9]) for coord, raw in records], blocks
    return result


def _unpack_records(batch):
    metadata, blocks = batch
    result = [None] * len(metadata)
    for indices, arrays in blocks:
        for row, index in enumerate(indices):
            coord, raw = metadata[index]
            result[index] = coord, raw, {key: values[row] for key, values in arrays.items()}
    return result


def parsed_records(records, intervals, workers, log_level, batch_records=512, batch_bytes=8*1024**2):
    """Filter to final intervals, then bound both pending and completed batches.

    The byte target includes raw text and decoded int32 arrays. A single larger
    record remains intact. Up to twice the worker count may be outstanding.
    Source iteration still validates every header and coordinate, including
    excluded records and trailing chunks, and checks subprocess exit status.
    """
    def jobs():
        interval = 0
        batch, size, number = [], 0, 0
        for coord, fields in records:
            while interval < len(intervals) and (int(intervals[interval][0]), int(intervals[interval][2])) <= coord:
                interval += 1
            if interval == len(intervals):
                continue
            chrom, start, end = intervals[interval]
            if coord[0] != chrom or not start <= coord[1] < end:
                continue
            k = 1 if fields[4] == '.' else 2 + fields[4].count(',')
            estimate = sum(map(len, fields)) + 4*(len(fields)-9)*(2+k+k*(k+1)//2)
            if batch and (len(batch) >= batch_records or size + estimate > batch_bytes):
                yield number, (parse_record_batch, {'records': batch, 'packed': workers > 1})
                batch, size, number = [], 0, number+1
            batch.append((coord, fields))
            size += estimate
        if batch:
            yield number, (parse_record_batch, {'records': batch, 'packed': workers > 1})

    if workers == 1:
        for _, (_, kwargs) in jobs():
            yield from parse_record_batch(**kwargs)
        return
    results = run_with_pool_iter(jobs(), log_level, max_workers=workers,
                                 max_inflight=2*workers, ordered=True)
    try:
        for _, batch in results:
            yield from _unpack_records(batch)
            # Release the consumed payload before the pool refills its window.
            del batch
    finally:
        results.close()
