"""Stream native SNP evidence into the self-contained HDF5 v3 database."""
from contextlib import ExitStack
import json
from pathlib import Path

import h5py
import numpy as np
from loguru import logger

from .. import __version__
from ..utils.exceptions import IPyradError
from ..utils.parallel import stream_pipeline_lines
from .genotype_evidence import _read_header, _coordinate_groups, _masked_fields
from .variants import BIN_BCF, _load_mask_entries_from_paths
from .reporting import EvidenceCounts, store_summary


FORMAT_VERSION = "3.0"
STATUSES = ("PASS", "NoHardCall", "NoRead", "HardMask", "NoMatch", "NoLikelihood")


def _ints(value, width):
    if value == ".":
        return [-1] * width
    fields = value.split(",")
    if len(fields) != width:
        raise IPyradError("Native SNP evidence has an invalid allele/genotype cardinality.")
    values = [-1 if v == "." else int(v) for v in fields]
    if any(v < -1 or v > np.iinfo(np.int32).max for v in values):
        raise IPyradError("Native SNP evidence exceeds the lossless int32 representation.")
    return values


def _column(values, width):
    """Parse a sample column in NumPy, keeping cardinality/missingness checks."""
    missing = ','.join(['-1'] * width)
    if any(v != '.' and v.count(',') != width-1 for v in values):
        raise IPyradError('Native SNP evidence has an invalid allele/genotype cardinality.')
    text = ','.join(missing if v == '.' else v for v in values).replace('.', '-1')
    data = np.fromstring(text, sep=',', dtype=np.int64)
    if data.size != len(values)*width or np.any(data < -1) or np.any(data > np.iinfo(np.int32).max):
        raise IPyradError('Native SNP evidence exceeds the lossless int32 representation.')
    return data.reshape(len(values), width).astype('i4')


def write_snp_evidence_hdf5(*, database, final_vcf, raw_vcf=None,
                            hard_masks=(), shared_mask=None, provenance=None,
                            chunk_sites=1024):
    """Add ordered alleles and optional native evidence without changing calls.

    Only SNP rows already present in genos/snpsmap are included. Exact native
    coordinate AND ordered-allele matches are mandatory for evidence. Unmatched
    rows remain in the database with NoMatch rather than silently disappearing.
    Hard masks and zero-read cells never acquire usable likelihoods.
    """
    enabled = raw_vcf is not None
    with ExitStack() as stack:
        io = stack.enter_context(h5py.File(database, "a", rdcc_nbytes=256*1024**2,
                                         rdcc_nslots=1000003))
        if any(k in io for k in ("alleles", "evidence", "_snp_v3_staging")):
            raise IPyradError("HDF5 SNP metadata/evidence already exists; refusing to replace it.")
        names = [str(x.decode() if isinstance(x, bytes) else x) for x in io['genos'].attrs['names']]
        scaffolds = [str(x.decode() if isinstance(x, bytes) else x) for x in io.attrs['scaffold_names']]
        ns, nv = io['genos'].shape[:2]
        # Keep handles open so the large existing call chunks remain cached
        # while the smaller evidence chunks are streamed.
        genos_ds, map_ds, ref_ds = io['genos'], io['snpsmap'], io['reference']
        final = stack.enter_context(stream_pipeline_lines([[BIN_BCF, 'view', '--no-version', str(final_vcf)]]))
        final_meta, final_header = _read_header(final)
        if final_header.split('\t')[9:] != names:
            raise IPyradError('Final VCF and HDF5 sample orders differ.')
        raw_meta = []
        raw_groups = iter(())
        contigs = {line.split('=', 2)[2].split(',', 1)[0].rstrip('>'): i
                   for i, line in enumerate(final_meta) if line.startswith('##contig=<ID=')}
        if enabled:
            raw = stack.enter_context(stream_pipeline_lines([[BIN_BCF, 'view', '--no-version', str(raw_vcf)]]))
            raw_meta, raw_header = _read_header(raw)
            if raw_header.split('\t')[9:] != names:
                raise IPyradError('Native VCF and HDF5 sample orders differ.')
            contigs = {line.split('=', 2)[2].split(',', 1)[0].rstrip('>'): i
                       for i, line in enumerate(raw_meta) if line.startswith('##contig=<ID=')}
            for tag in ('DP', 'AD', 'GQ', 'PL'):
                if sum(line.startswith(f'##FORMAT=<ID={tag},') for line in raw_meta) != 1:
                    raise IPyradError(f'Native VCF lacks a unique FORMAT/{tag} declaration.')
            raw_groups = _coordinate_groups(raw, contigs)
        current = next(raw_groups, None)
        sample_index = {name: i for i, name in enumerate(names)}
        masks = list(hard_masks)
        if shared_mask is not None:
            sample_index[''] = ns
            masks.append({'': Path(shared_mask)})
        entries = {}
        for layer in masks:
            existing = {name: Path(path) for name, path in layer.items() if Path(path).is_file()}
            for chrom, values in _load_mask_entries_from_paths(existing, sample_index).items():
                for entry in values:
                    entry['intervals'].sort()
                entries.setdefault(chrom, []).extend(values)
        stage = io.create_group('_snp_v3_staging')
        ncols = max(1, min(chunk_sites, max(1, nv)))
        kw = dict(compression='gzip', compression_opts=4, shuffle=True)
        alleles_ds = stage.create_dataset('alleles', shape=(nv, 4), maxshape=(None, 4),
                                         dtype='u1', chunks=(ncols, 4), fillvalue=0, **kw)
        alleles_ds.attrs['description'] = 'Ordered REF then ALT ASCII bases; zero pads unused slots.'
        evidence = None
        if enabled:
            evidence = stage.create_group('evidence')
            evidence.attrs['schema_version'] = FORMAT_VERSION
            evidence.attrs['scope'] = 'retained_snps'
            evidence.attrs['status_codes'] = json.dumps(dict(enumerate(STATUSES)))
            evidence.attrs['missing_integer'] = -1
            evidence.attrs['provenance'] = json.dumps(provenance or {}, sort_keys=True)
            evidence.create_dataset('native_header', data='\n'.join(raw_meta), dtype=h5py.string_dtype())
            evidence.create_dataset('final_header', data='\n'.join(final_meta), dtype=h5py.string_dtype())
            for key, width in (('DP', 0), ('AD', 4), ('GQ', 0), ('PL', 10), ('status', 0)):
                extra = (width,) if width else ()
                dtype, fill = ('u1', 4) if key == 'status' else ('i4', -1)
                evidence.create_dataset(key, shape=(ns, nv, *extra), maxshape=(ns, None, *extra),
                                        chunks=(ns, ncols, *extra), dtype=dtype, fillvalue=fill, **kw)
        index = 0
        start = end = 0
        counts = np.zeros(len(STATUSES), dtype=np.int64)
        matched = 0
        counter = EvidenceCounts(names, STATUSES)

        def flush():
            if end <= start:
                return
            alleles_ds[start:end] = allele_buffer
            if evidence is not None:
                counter.add(buffers, np.count_nonzero(allele_buffer, axis=1))
                for key, value in buffers.items():
                    evidence[key][:, start:end] = value

        try:
            for coord, group in _coordinate_groups(final, contigs):
                if index == nv:
                    break
                while current is not None and current[0] < coord:
                    current = next(raw_groups, None)
                originals = current[1] if current is not None and current[0] == coord else []
                excluded = set(_masked_fields(entries.get(group[0][0], ()), coord[1] - 1))
                for called in group:
                    alleles = [called[3], *called[4].split(',')]
                    if called[6] != 'PASS' or not (2 <= len(alleles) <= 4) or not all(a in 'ACGT' and len(a) == 1 for a in alleles):
                        continue
                    if index == nv:
                        continue
                    if index == end:
                        flush()
                        start, end = index, min(nv, index + chunk_sites)
                        maps = map_ds[start:end]
                        calls = genos_ds[:, start:end, :2]
                        refs = ref_ds[start:end]
                        allele_buffer = np.zeros((end-start, 4), dtype='u1')
                        buffers = {key: np.full((ns, end-start, *ds.shape[2:]),
                                                4 if key == 'status' else -1, dtype=ds.dtype)
                                   for key, ds in (evidence.items() if evidence is not None else ())
                                   if key in ('DP', 'AD', 'GQ', 'PL', 'status')}
                    j = index - start
                    expected = (contigs[scaffolds[int(maps[j, 3])]], int(maps[j, 4]) + 1)
                    if expected != coord:
                        if expected < coord:
                            raise IPyradError(f'Final VCF lacks HDF5 SNP row {index}.')
                        continue
                    if len(set(alleles)) != len(alleles) or ord(alleles[0]) != refs[j]:
                        raise IPyradError('HDF5 SNP REF/allele mismatch.')
                    allele_buffer[j, :len(alleles)] = list(map(ord, alleles))
                    tags = called[8].split(':')
                    gtidx = tags.index('GT')
                    gt_text = ','.join(x.split(':')[gtidx] for x in called[9:])
                    gt_text = gt_text.replace('.', '255').replace('/', ',').replace('|', ',')
                    gt = np.fromstring(gt_text, sep=',', dtype=np.int64)
                    if gt.size != ns*2 or not np.array_equal(gt.reshape(ns,2), calls[:,j]):
                        raise IPyradError('HDF5 SNP GT does not match final VCF.')
                    if enabled:
                        statuses = np.full(ns, 4, dtype='u1')
                        matches = [x for x in originals if x[3:5] == called[3:5]]
                        unique = len(matches) == 1 and sum(x[3:5] == called[3:5] for x in group) == 1
                        if unique:
                            matched += 1
                            native = matches[0]
                            native_tags = native[8].split(':')
                            fields = [x.split(':') for x in native[9:]]
                            na = len(alleles)
                            ng = na*(na+1)//2
                            for key, width in (('DP',1),('AD',na),('GQ',1),('PL',ng)):
                                tagidx = native_tags.index(key) if key in native_tags else -1
                                values = [x[tagidx] if 0 <= tagidx < len(x) else '.' for x in fields]
                                data = _column(values, width)
                                if width == 1:
                                    buffers[key][:,j] = data[:,0]
                                else:
                                    buffers[key][:,j,:width] = data
                            statuses[:] = np.where(np.any(calls[:,j] == 255, axis=1), 1, 0)
                            statuses[np.any(buffers['PL'][:,j,:ng] < 0, axis=1)] = 5
                            no_read = buffers['DP'][:,j] <= 0
                            statuses[no_read] = 2
                            for key in ('AD','GQ','PL'):
                                buffers[key][no_read,j] = -1
                        hard = np.array([i+9 in excluded or ns+9 in excluded for i in range(ns)])
                        statuses[hard] = 3
                        for key in ('DP','AD','GQ','PL'):
                            buffers[key][hard,j] = -1
                        buffers['status'][:,j] = statuses
                        counts += np.bincount(statuses, minlength=len(STATUSES))
                    index += 1
                    if index % 10000 == 0:
                        logger.info('HDF5 SNP evidence progress: {}/{} records', index, nv)
            if index != nv:
                raise IPyradError(f'Only {index}/{nv} HDF5 SNPs matched the final VCF.')
            flush()
            if evidence is not None:
                evidence.attrs['matched_native_records'] = matched
                evidence.attrs['status_counts'] = json.dumps(dict(zip(STATUSES, map(int, counts))))
                summary = dict(counter.summary(), sites=nv, matched_native_records=matched,
                               unit='sample_snp_cell', scope='retained_snps')
                store_summary(evidence, summary)
                evidence.attrs['complete'] = True
            io.move('_snp_v3_staging/alleles', 'alleles')
            if enabled:
                io.move('_snp_v3_staging/evidence', 'evidence')
            del io['_snp_v3_staging']
            io.attrs['version'] = 3.0
            io.attrs['format_version'] = FORMAT_VERSION
            io.attrs['evidence_scope'] = 'snps' if enabled else 'none'
            io.attrs['ipyrad2_version'] = __version__
        except BaseException:
            if '_snp_v3_staging' in io:
                del io['_snp_v3_staging']
            raise
    logger.info('HDF5 v3: {} SNP rows; evidence {}; {} exact native matches', nv,
                'enabled' if enabled else 'disabled', matched)

    return summary if enabled else None
