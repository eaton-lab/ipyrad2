"""Lossless native pileup evidence over the final RAD BED, HDF5 format 4.0.

Records and bases are separate: uncovered bases have no record, and one base
can anchor several native records (e.g. SNP and indel). Ragged allele/PL axes
preserve symbolic alleles and arbitrary allele cardinality without truncation.
"""
from contextlib import ExitStack
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np
from loguru import logger

from .reporting import EvidenceCounts, store_summary
from ..utils.exceptions import IPyradError
from ..utils.parallel import stream_pipeline_lines
from .evidence_parallel import parsed_records
from .genotype_evidence import _read_header
from .variants import BIN_BCF, BIN, _load_mask_entries_from_paths

FORMAT_VERSION = '4.0'
RECORD_STATUSES = ('PASS', 'NoRead', 'NoLikelihood', 'HardMask', 'NonACGT')
SITE_STATUSES = ('NoRecord', 'HasRecord', 'NonACGT')


def read_rad_intervals(bed, contigs):
    """Validate sorted, disjoint, half-open BEDs without silently merging loci."""
    result = []
    for line in Path(bed).read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.split()
        if len(fields) < 3 or fields[0] not in contigs:
            raise IPyradError('RAD evidence BED has an unknown contig or malformed row.')
        chrom, start, end = contigs[fields[0]], int(fields[1]), int(fields[2])
        if start < 0 or end <= start:
            raise IPyradError('RAD evidence BED requires nonempty, nonnegative intervals.')
        if result and (chrom < result[-1][0] or (chrom == result[-1][0] and start < result[-1][2])):
            raise IPyradError('RAD evidence BED must be sorted and nonoverlapping.')
        result.append((chrom, start, end))
    return np.asarray(result, dtype='u8').reshape(-1, 3)


def _sequences(reference, names, intervals, region_file):
    region_file.write_text(''.join(f'{names[c]}:{a+1}-{b}\n' for c,a,b in intervals))
    if not len(intervals):
        return
    with stream_pipeline_lines([[str(BIN/'samtools'), 'faidx', str(reference), '-r', str(region_file)]]) as lines:
        header, chunks = None, []
        for line in lines:
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(chunks).upper()
                header, chunks = line[1:], []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks).upper()


def _records(paths, names, contigs, group):
    previous = None
    headers = group.create_group('native_headers')
    for i, path in enumerate(paths):
        with stream_pipeline_lines([[BIN_BCF, 'view', '--no-version', str(path)]]) as lines:
            meta, header = _read_header(lines)
            if header.split('\t')[9:] != names:
                raise IPyradError('RAD native evidence and database sample orders differ.')
            for key in ['DP', 'AD', 'PL']:
                if sum(m.startswith(f'##FORMAT=<ID={key},') for m in meta) != 1:
                    raise IPyradError(f'RAD native evidence requires unique FORMAT/{key}.')
            raw = ('\n'.join([*meta, header])+'\n').encode()
            ds = headers.create_dataset(str(i), data=np.frombuffer(raw, dtype='u1'), compression='gzip')
            ds.attrs['source'] = str(path)
            for line in lines:
                f = line.split('\t')
                if len(f) != 9+len(names):
                    raise IPyradError('RAD native record has the wrong sample count.')
                if f[0] not in contigs:
                    # Source pileup may contain contigs excluded by final filtering.
                    continue
                coord = (contigs[f[0]], int(f[1])-1)
                if previous is not None and coord < previous:
                    raise IPyradError('RAD native records/chunks are out of coordinate order.')
                previous = coord
                yield coord, f


def _mask_window(entries, start, end, samples):
    """Expand interval masks once per bounded window instead of once per base."""
    result = np.zeros((samples, end-start), dtype='u1')
    for entry in entries:
        spans, cursor = entry['intervals'], entry['cursor']
        while cursor < len(spans) and spans[cursor][1] <= start:
            cursor += 1
        entry['cursor'] = cursor
        row = entry['field_idx'] - 9
        while cursor < len(spans):
            left, right = spans[cursor]
            if left >= end:
                break
            if right > start:
                columns = slice(max(start, left)-start, min(end, right)-start)
                if row == samples:
                    result[:, columns] = 1
                else:
                    result[row, columns] = 1
            cursor += 1
    return result


def write_rad_evidence_hdf5(*, database, bed, reference, native_files, native_regions,
                            hard_masks=(), shared_mask=None, provenance=None,
                            chunk_sites=1024, workers=1, log_level="INFO"):
    """Append optional all-RAD evidence, never editing canonical/SNP arrays.

    Native values, including zero-depth and masked-cell values, are preserved.
    Explicit hard masks and status must be enforced by every reader. No missing
    likelihood is replaced by a reference call or a uniform likelihood vector.
    """
    paths = list(map(Path, native_files))
    regions = list(map(Path, native_regions))
    if len(regions) != len(paths) or any(not p.is_file() for p in regions):
        raise IPyradError("Each native capture requires its original BED domain.")
    if not paths or any(not p.is_file() for p in paths):
        raise IPyradError('All-RAD evidence requires existing pre-call native pileup files.')
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise IPyradError("RAD evidence workers must be a positive integer.")
    if chunk_sites < 1:
        raise IPyradError('RAD evidence chunk size must be positive.')
    with ExitStack() as stack:
        h = stack.enter_context(h5py.File(database, 'a', rdcc_nbytes=8*1024**2, rdcc_nslots=10007))
        if h.attrs.get('format_version') != '3.0' or h.attrs.get('evidence_scope') != 'snps':
            raise IPyradError('All-RAD evidence requires a completed SNP-evidence HDF5 v3 database.')
        if not h['evidence'].attrs.get('complete', False):
            raise IPyradError('SNP evidence is incomplete.')
        if 'rad_evidence' in h or '_rad_v4_staging' in h:
            raise IPyradError('RAD evidence already exists; refusing to replace it.')
        names = [x.decode() if isinstance(x, bytes) else str(x) for x in h['genos'].attrs['names']]
        scaffolds = [x.decode() if isinstance(x, bytes) else str(x) for x in h.attrs['scaffold_names']]
        contigs = {s:i for i,s in enumerate(scaffolds)}
        intervals = read_rad_intervals(bed, contigs)
        domains = {}
        for path in regions:
            for line in path.read_text().splitlines():
                if not line.strip() or line.startswith('#'): continue
                c,a,b,*_ = line.split();a,b=int(a),int(b)
                if a<0 or b<=a: raise IPyradError('Invalid native capture BED domain.')
                if c in contigs: domains.setdefault(contigs[c],[]).append((a,b))
        for c, spans in domains.items():
            merged=[]
            for a,b in sorted(spans):
                if merged and a<=merged[-1][1]:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
                else:merged.append((a,b))
            domains[c]=merged
        domain_cursor={c:0 for c in domains}
        lengths = intervals[:,2]-intervals[:,1]
        nsites, ns = int(lengths.sum()), len(names)
        if ns == 0:
            raise IPyradError('RAD evidence requires observed samples.')
        # Process startup dominates tiny panels. Retain the same batch parser.
        parse_workers = workers if ns*nsites >= 250_000 else 1
        entries = {}
        layers = list(hard_masks)
        sample_index = dict(zip(names, range(ns)))
        if shared_mask is not None:
            sample_index[''] = ns
            layers.append({'': Path(shared_mask)})
        for layer in layers:
            unknown = set(layer)-set(sample_index)
            if unknown or any(not Path(p).is_file() for p in layer.values()):
                raise IPyradError('RAD hard masks must identify known samples and existing BED files.')
            for chrom, vals in _load_mask_entries_from_paths(layer, sample_index).items():
                for entry in vals:
                    entry['intervals'].sort()
                entries.setdefault(chrom, []).extend(vals)
        g = h.create_group('_rad_v4_staging')
        kw = dict(compression='gzip', compression_opts=4, shuffle=True)
        site_chunk = min(chunk_sites, max(1, nsites))
        g.attrs.update(schema_version=FORMAT_VERSION, scope='final_delimited_RAD_intervals',
                       complete=False, coordinate_convention='zero-based half-open BED; site positions zero-based',
                       status_codes=json.dumps(dict(enumerate(RECORD_STATUSES))),
                       site_status_codes=json.dumps(dict(enumerate(SITE_STATUSES))),
                       missing_integer=-1, provenance=json.dumps(provenance or {}, sort_keys=True),
                       masked_values='native values retained for audit; never usable when hard_mask is set',
                       allele_semantics='literal native REF/ALT; symbolic unseen alleles are not specific nucleotides',
                       genotype_semantics='pre-call likelihoods only; no GT or posterior inferred')
        capture=g.create_group('capture_regions')
        for i,path in enumerate(regions):
            ds=capture.create_dataset(str(i),data=np.frombuffer(path.read_bytes(),dtype='u1'),compression='gzip')
            ds.attrs['source']=str(path)
        g.create_dataset('intervals', data=intervals)
        g.create_dataset('interval_offsets', data=np.r_[0,np.cumsum(lengths)].astype('u8'))
        g.create_dataset('samples', data=names, dtype=h5py.string_dtype())
        g.create_dataset('scaffolds', data=scaffolds, dtype=h5py.string_dtype())
        for key, tail, dtype in [('site_map',(3,),'u8'),('reference',(),'u1'),('site_status',(),'u1'),('record_span',(2,),'u8')]:
            g.create_dataset(key, shape=(nsites,*tail), maxshape=(None,*tail), dtype=dtype,
                             chunks=(site_chunk,*tail), **kw)
        g.create_dataset('hard_mask', shape=(ns,nsites), maxshape=(ns,None), dtype='u1',
                         chunks=(ns,site_chunk), **kw)
        for key,dtype in [('record_site','u8'),('allele_ids','u4'),('QUAL','f8')]:
            g.create_dataset(key, shape=(0,), maxshape=(None,), dtype=dtype, chunks=(chunk_sites,), **kw)
        for key in ['allele_offsets','PL_offsets']:
            g.create_dataset(key, data=np.array([0],dtype='u8'), maxshape=(None,), chunks=(chunk_sites,), **kw)
        for key in ['DP','GQ','AD','PL','status']:
            g.create_dataset(key, shape=(ns,0), maxshape=(ns,None), dtype='u1' if key=='status' else 'i4',
                             chunks=(ns,chunk_sites), **kw)
        # Keep handles alive so partial compressed chunks remain cached between appends.
        datasets = {key: value for key, value in g.items() if isinstance(value, h5py.Dataset)}
        counter = EvidenceCounts(names, RECORD_STATUSES)
        symbolic_records = 0
        site_counts = np.zeros(len(SITE_STATUSES), dtype='u8')
        tokens, token_ids = [], {}
        nr = na = ng = 0
        start = 0
        site_buffer, records_buffer = [], []

        def append(key, data, axis=0):
            ds = datasets[key];old = ds.shape[axis]
            ds.resize(old+data.shape[axis], axis=axis)
            if axis: ds[:,old:] = data
            else: ds[old:] = data

        def flush():
            nonlocal start
            if not site_buffer:
                return
            end = start+len(site_buffer)
            for key,col in [('site_map',0),('reference',1),('site_status',2),('record_span',3)]:
                datasets[key][start:end] = np.asarray([s[col] for s in site_buffer])
            datasets['hard_mask'][:,start:end] = np.array([s[4] for s in site_buffer]).T
            if records_buffer:
                for key,col in [('record_site',0),('QUAL',1)]:
                    append(key,np.array([r[col] for r in records_buffer]))
                append('allele_ids',np.concatenate([r[2] for r in records_buffer]))
                append('allele_offsets',np.array([r[3] for r in records_buffer],dtype='u8'))
                append('PL_offsets',np.array([r[4] for r in records_buffer],dtype='u8'))
                block = {key: np.concatenate([r[5][key] for r in records_buffer], axis=1)
                         for key in ['DP','GQ','AD','PL','status']}
                allele_starts = np.r_[0, np.cumsum([r[5]['AD'].shape[1] for r in records_buffer[:-1]])].astype('intp')
                likelihood_starts = np.r_[0, np.cumsum([r[5]['PL'].shape[1] for r in records_buffer[:-1]])].astype('intp')
                counter.add_ragged(block, allele_starts, likelihood_starts)
                for key, values in block.items():
                    append(key, values, 1)
            start = end
            site_buffer.clear();records_buffer.clear()

        try:
            source_records = _records(paths, names, contigs, g)
            stack.callback(source_records.close)
            records = parsed_records(source_records, intervals, parse_workers, log_level)
            stack.callback(records.close)
            current = next(records, None)
            temp = Path(stack.enter_context(TemporaryDirectory(prefix='.rad-evidence-', dir=Path(database).parent)))
            sequences = _sequences(reference, scaffolds, intervals, temp/'regions.txt')
            stack.callback(sequences.close)
            for li,(chrom,a,b) in enumerate(intervals):
                header, sequence = next(sequences)
                if header != f'{scaffolds[chrom]}:{a+1}-{b}' or len(sequence) != b-a:
                    raise IPyradError('Reference FASTA does not match RAD evidence coordinates.')
                mask_end = 0
                for off,base in enumerate(sequence):
                    pos = int(a)+off;site = start+len(site_buffer);coord = (int(chrom),pos)
                    if base in 'ACGT':
                        spans=domains.get(int(chrom),());cursor=domain_cursor.get(int(chrom),0)
                        while cursor<len(spans) and spans[cursor][1]<=pos:cursor+=1
                        domain_cursor[int(chrom)]=cursor
                        if cursor==len(spans) or not spans[cursor][0]<=pos<spans[cursor][1]:
                            raise IPyradError('A RAD reference base lies outside the captured BED domains.')
                    if off == mask_end:
                        mask_end = min(len(sequence), off + chunk_sites)
                        mask_start = off
                        masks = _mask_window(entries.get(scaffolds[chrom], ()), int(a)+off, int(a)+mask_end, ns)
                    hard = masks[:, off-mask_start]
                    while current is not None and current[0] < coord:
                        current = next(records,None)
                    first = nr;seen = set()
                    while current is not None and current[0] == coord:
                        f = current[1];alleles = [f[3]] + ([] if f[4]=='.' else f[4].split(','))
                        if tuple(alleles) in seen or len(set(alleles))!=len(alleles):
                            raise IPyradError('Duplicate native RAD allele records are ambiguous.')
                        seen.add(tuple(alleles))
                        if not alleles[0] or alleles[0][0].upper()!=base:
                            raise IPyradError('Native RAD REF does not match the indexed reference base.')
                        ids = []
                        for token in alleles:
                            if token not in token_ids:
                                token_ids[token]=len(tokens);tokens.append(token)
                            ids.append(token_ids[token])
                        k=len(alleles);width=k*(k+1)//2;data=current[2]
                        status=np.zeros(ns,dtype='u1')
                        status[(data['PL']<0).any(axis=1)]=2
                        status[data['DP'][:,0]<=0]=1
                        status[hard.astype(bool)]=3
                        if base not in 'ACGT':status[:]=4
                        data['status']=status[:,None]
                        symbolic_records += any(a == '*' or a.startswith('<') or '[' in a or ']' in a for a in alleles)
                        na+=k;ng+=width;nr+=1
                        records_buffer.append((site,float(f[5]) if f[5]!='.' else np.nan,np.array(ids,'u4'),na,ng,data))
                        current=next(records,None)
                    status=2 if base not in 'ACGT' else int(nr>first)
                    site_counts[status]+=1
                    site_buffer.append(((chrom,pos,li),ord(base),status,(first,nr-first),hard))
                    if len(site_buffer)>=chunk_sites:flush()
                if li%1000==0:logger.info('RAD evidence: {}/{} intervals, {} native records',li+1,len(intervals),nr)
            if next(sequences,None) is not None:
                raise IPyradError('Unexpected extra reference intervals.')
            # Consume remaining source records so subprocess errors cannot be hidden.
            for _ in records: pass
            flush()
            if start!=nsites:
                raise IPyradError('RAD evidence site count mismatch.')
            counts = counter.counts.sum(axis=0)
            g.create_dataset('allele_tokens',data=tokens,dtype=h5py.string_dtype())
            g.attrs.update(capture_domains_validated=True,sites=nsites,records=nr,status_counts=json.dumps(dict(zip(RECORD_STATUSES,map(int,counts)))),
                           site_status_counts=json.dumps(dict(zip(SITE_STATUSES,map(int,site_counts)))),complete=True)
            summary = dict(counter.summary(), sites=nsites, records=nr,
                           site_status_counts=dict(zip(SITE_STATUSES,map(int,site_counts))),
                           symbolic_allele_records=symbolic_records, unit='sample_record_cell',
                           scope='all_bases_in_final_retained_locus_BED')
            store_summary(g, summary)
            h.move('_rad_v4_staging','rad_evidence')
            h.attrs.update(version=4.0,format_version=FORMAT_VERSION,evidence_scope='rad_sites')
        except BaseException:
            if '_rad_v4_staging' in h:del h['_rad_v4_staging']
            raise
    return summary
