"""Mask duplicated overlap evidence in retained anomalous pairs.

HTSlib's native overlap correction only handles PROPER_PAIR records. Preserve
flags/sequence and conservatively count anomalous overlapping mates once, before
both assemble callers. This is not PCR/UMI duplicate removal.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
from itertools import groupby


def aligned_positions(fields: list[str]) -> dict[int, int]:
    """Map reference offsets to query offsets, excluding clips/insertions/gaps."""
    import re
    ref = int(fields[3]); query = 0; result = {}
    for raw, op in re.findall(r'(\d+)([MIDNSHP=X])', fields[5]):
        size = int(raw)
        if op in 'M=X':
            result.update((ref + i, query + i) for i in range(size))
            ref += size; query += size
        elif op in 'IS':
            query += size
        elif op in 'DN':
            ref += size
    return result


def mask_group(records: list[list[str]], counts: Counter) -> None:
    """Modify QUAL only for unambiguous primary mate pairs; leave flags intact."""
    by_rg = {}
    for fields in records:
        flag = int(fields[1])
        if not flag & 1 or flag & 0x90c:
            continue
        rg = next((t[5:] for t in fields[11:] if t.startswith('RG:Z:')), '')
        by_rg.setdefault(rg, []).append(fields)
    for pair in by_rg.values():
        if len(pair) != 2:
            counts['unmatched_or_ambiguous_groups'] += 1
            continue
        a, b = pair
        fa, fb = int(a[1]), int(b[1])
        if bool(fa & 64) == bool(fb & 64) or bool(fa & 128) == bool(fb & 128):
            counts['unmatched_or_ambiguous_groups'] += 1
            continue
        if fa & 2 and fb & 2:
            continue  # leave normal pairs to HTSlib's existing overlap correction
        if a[2] != b[2] or a[9] == '*' or b[9] == '*' or a[10] == '*' or b[10] == '*':
            continue
        # Require reciprocal mate coordinates; duplicate read names must not
        # accidentally combine unrelated templates.
        if int(a[7]) != int(b[3]) or int(b[7]) != int(a[3]):
            counts['nonreciprocal_pairs'] += 1
            continue
        pa, pb = aligned_positions(a), aligned_positions(b)
        overlap = pa.keys() & pb.keys()
        if not overlap:
            continue
        qa, qb = list(a[10]), list(b[10])
        changed_a = changed_b = False
        for pos in overlap:
            i, j = pa[pos], pb[pos]
            if a[9][i].upper() != b[9][j].upper():
                # Discordant mates provide no reliable consensus base.
                qa[i] = qb[j] = '!'; changed_a = changed_b = True
                counts['disagreeing_overlap_bases'] += 1
            elif ord(qa[i]) >= ord(qb[j]):
                qb[j] = '!'; changed_b = True
            else:
                qa[i] = '!'; changed_a = True
        for fields, quality, changed in ((a, qa, changed_a), (b, qb, changed_b)):
            if changed:
                if not any(t.startswith('OQ:Z:') for t in fields[11:]):
                    fields.append('OQ:Z:' + fields[10])
                fields[10] = ''.join(quality)
        counts['anomalous_overlapping_pairs'] += 1
        counts['overlap_positions'] += len(overlap)


def mask_stream(source, target, stats: Path) -> None:
    counts = Counter()
    def records():
        for line in source:
            if line.startswith('@'):
                target.write(line)
            elif line.strip():
                counts['alignment_records'] += 1
                yield line.rstrip('\n').split('\t')
    for _, group in groupby(records(), key=lambda f: f[0]):
        rows = list(group)
        mask_group(rows, counts)
        for row in rows:
            target.write('\t'.join(row) + '\n')
    stats.write_text(json.dumps(dict(counts), indent=2) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stats', type=Path, required=True)
    args = parser.parse_args()
    mask_stream(sys.stdin, sys.stdout, args.stats)
