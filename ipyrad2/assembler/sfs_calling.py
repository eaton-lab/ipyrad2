"""Bounded, cross-fitted SFS calling on the unchanged native pileup stream.

Native variants remain a separate input for unsupported records. Acceptance
BEDs cover reference decisions as well as emitted variants, so a failed SFS
decision cannot silently become a reference base during consensus generation.
"""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack, nullcontext
import json
from pathlib import Path
import time

import numpy as np
from loguru import logger

from ..analysis.methods.popgen.sfs_prior import folded_saf, fit_folded_sfs, sfs_posterior_moments
from ..utils.exceptions import IPyradError
from ..utils.parallel import run_pipeline, stream_pipeline_lines
from .evidence_parallel import parse_record_batch
from .genotype_evidence import _read_header, _coordinate_groups
from .rad_evidence_hdf5 import read_rad_intervals
from .variants import BIN_BCF, BIN_BED

TRAINING_BYTES = 256 * 1024**2


def phred(probability):
    return -10 * np.log10(np.maximum(1 - np.asarray(probability), 1e-12))


def allele_axis(records):
    """Return kind/triad indices only for the validated single-record model."""
    if len(records) != 1:
        return 0, None, None
    row = records[0]
    alleles = [row[3], *row[4].split(',')]
    if row[3] not in ('A', 'C', 'G', 'T') or len(set(alleles)) != len(alleles):
        return 0, None, None
    concrete = [i for i, a in enumerate(alleles) if a in ('A', 'C', 'G', 'T')]
    symbolic = [i for i, a in enumerate(alleles) if a == '<*>']
    if len(concrete) + len(symbolic) != len(alleles) or len(symbolic) > 1:
        return 0, None, None
    if len(concrete) == 1 and len(symbolic) == 1:
        a, b = concrete[0], symbolic[0]
        kind = 1
    elif len(concrete) == 2:
        a, b = concrete
        kind = 2
    else:
        return 0, None, None
    return kind, (a, b), (a*(a+1)//2+a, b*(b+1)//2+a, b*(b+1)//2+b)


def reference_fallback(records, allowed, min_depth, geno_q):
    """Qualify REF only on a single unsupported multiallelic nucleotide record.

    The caller must additionally require absence of a native variant record.
    Confidence normalizes *all* literal native diploid PL entries with a uniform
    genotype prior; it is not an SFS posterior or an empirical accuracy score.
    Missing evidence, indels and multiple records never qualify here.
    """
    keep = np.zeros(len(allowed), dtype=bool)
    quality = np.full(len(allowed), np.nan)
    if len(records) != 1:
        return keep, quality
    row = records[0]
    alleles = [row[3], *row[4].split(',')]
    concrete = sum(a in ('A', 'C', 'G', 'T') for a in alleles)
    if (row[3] not in ('A', 'C', 'G', 'T') or concrete not in (3, 4)
            or len(set(alleles)) != len(alleles)
            or any(a not in ('A', 'C', 'G', 'T', '<*>') for a in alleles)):
        return keep, quality
    if len(row) - 9 != len(allowed):
        raise IPyradError('Reference fallback sample count differs from assembled samples.')
    _, _, values = parse_record_batch([((0, 0), row)])[0]
    pl, ad, dp = values['PL'], values['AD'], values['DP'][:, 0]
    valid = (np.asarray(allowed, dtype=bool) & (dp >= max(5, min_depth))
             & (pl >= 0).all(axis=1) & (ad >= 0).all(axis=1) & (ad[:, 0] > 0))
    ids = np.flatnonzero(valid)
    if len(ids):
        likelihoods = 10. ** (-(pl[ids].astype(float) - pl[ids].min(axis=1)[:, None]) / 10)
        probability = likelihoods[:, 0] / likelihoods.sum(axis=1)
        quality[ids] = phred(probability)
        keep[ids] = (quality[ids] >= geno_q) & (pl[ids, 0, None] < pl[ids, 1:]).all(axis=1)
    return keep, quality


class NativeSource:
    """Repeatable bounded stream; no permanent additional observation database."""
    def __init__(self, paths, names, bed, ref_info, retained_beds=None, block_sites=2048):
        self.paths, self.names = list(map(Path, paths)), list(names)
        self.contigs = {line.split()[0]: i for i, line in enumerate(Path(ref_info).read_text().splitlines())}
        self.intervals = read_rad_intervals(bed, self.contigs).astype(np.int64)
        self.block_sites = block_sites
        self.retained = {}
        for name, path in (retained_beds or {}).items():
            spans = {}
            for line in Path(path).read_text().splitlines():
                if line.strip() and not line.startswith('#'):
                    c, a, b, *_ = line.split()
                    spans.setdefault(c, []).append((int(a), int(b)))
            self.retained[name] = {c: np.array(sorted(v), dtype=np.int64) for c, v in spans.items()}

    def lines(self):
        for path in self.paths:
            with stream_pipeline_lines([[BIN_BCF, 'view', '--no-version', str(path)]]) as lines:
                meta, header = _read_header(lines)
                if header.split('\t')[9:] != self.names:
                    raise IPyradError('SFS capture sample order differs from assembled samples.')
                for key in ('DP', 'AD', 'PL'):
                    if sum(x.startswith(f'##FORMAT=<ID={key},') for x in meta) != 1:
                        raise IPyradError(f'SFS capture requires unique FORMAT/{key}.')
                for line in lines:
                    if len(line.split('\t')) != 9 + len(self.names):
                        raise IPyradError('SFS capture has an inconsistent sample count.')
                    yield line

    def allows(self, chrom, positions):
        result = np.ones((len(self.names), len(positions)), dtype=bool)
        for i, name in enumerate(self.names):
            if name not in self.retained:
                continue
            spans = self.retained[name].get(chrom)
            if spans is None or not len(spans):
                result[i] = False
            else:
                ix = np.searchsorted(spans[:, 0], positions, side='right') - 1
                result[i] = (ix >= 0) & (positions < spans[np.maximum(ix, 0), 1])
        return result

    def decode(self, records):
        ns, size = len(self.names), len(records)
        pl = np.full((ns, size, 3), -1, dtype=np.int32)
        dp = np.full((ns, size), -1, dtype=np.int32)
        ad = np.full((ns, size, 2), -1, dtype=np.int32)
        valid = np.zeros((ns, size), dtype=bool)
        kinds = np.zeros(size, dtype=np.uint8)
        singles, columns = [], []
        for j, (key, locus, rows) in enumerate(records):
            kind, axis, triad = allele_axis(rows)
            kinds[j] = kind
            if kind:
                singles.append((key, rows[0])); columns.append((j, axis, triad))
        for (_, _, values), (j, axis, triad) in zip(parse_record_batch(singles), columns):
            dp[:, j] = values['DP'][:, 0]
            ad[:, j] = values['AD'][:, list(axis)]
            pl[:, j] = values['PL'][:, list(triad)]
            valid[:, j] = (dp[:, j] > 0) & (values['PL'] >= 0).all(axis=1)
        for c in {r[2][0][0] for r in records}:
            js = [j for j, r in enumerate(records) if r[2][0][0] == c]
            valid[:, js] &= self.allows(c, np.array([records[j][0][1]-1 for j in js]))
        safe = np.where(valid[:, :, None], pl, 0).astype(float)
        safe -= safe.min(axis=2, keepdims=True)
        likelihoods = 10. ** (-safe / 10)
        return dict(records=records, loci=np.array([r[1] for r in records]),
                    kind=kinds, x=likelihoods, valid=valid, dp=dp, ad=ad)

    def blocks(self):
        locus, batch = 0, []
        for key, rows in _coordinate_groups(self.lines(), self.contigs):
            coord = (key[0], key[1]-1)
            while locus < len(self.intervals) and tuple(self.intervals[locus, [0, 2]]) <= coord:
                locus += 1
            if locus == len(self.intervals):
                continue
            c, a, b = self.intervals[locus]
            if coord[0] != c or not a <= coord[1] < b:
                continue
            batch.append((key, locus, rows))
            if len(batch) == self.block_sites:
                yield self.decode(batch); batch = []
        if batch:
            yield self.decode(batch)


def fit_models(source, groups, *, folds, fit_loci, max_iter, cores):
    """Match the frozen per-group/per-fold eligible-locus selection procedure."""
    if len(source.intervals) < folds:
        raise IPyradError(f'SFS needs at least {folds} retained loci for {folds} folds.')
    available = {g: set() for g in groups}
    for block in source.blocks():
        for group, ix in groups.items():
            available[group].update(block['loci'][block['valid'][ix].any(axis=0)].tolist())
    selected = {}
    for group in groups:
        loci = np.array(sorted(available[group]), dtype=int)
        for fold in range(folds):
            pool = loci[loci % folds != fold]
            if not len(pool):
                raise IPyradError(f'SFS group {group}, fold {fold}: no usable training loci.')
            selected[group, fold] = pool[np.unique(np.linspace(0, len(pool)-1, min(fit_loci, len(pool)), dtype=int))]
    unions = {g: np.unique(np.concatenate([selected[g, f] for f in range(folds)])) for g in groups}
    parts = {g: [] for g in groups}
    labels = {g: [] for g in groups}
    nbytes = 0
    for block in source.blocks():
        for g, ix in groups.items():
            take = np.flatnonzero(np.isin(block['loci'], unions[g]) & block['valid'][ix].any(axis=0))
            if not len(take):
                continue
            matrix = folded_saf(block['x'][ix][:, take], block['valid'][ix][:, take])
            nbytes += matrix.nbytes + take.size*8
            if nbytes > TRAINING_BYTES:
                raise IPyradError('SFS training exceeds 256 MiB; reduce --sfs-fit-loci or use smaller explicit cohorts.')
            parts[g].append(matrix); labels[g].append(block['loci'][take])
    matrices = {g: np.concatenate(v) for g, v in parts.items()}
    labels = {g: np.concatenate(v) for g, v in labels.items()}
    parts.clear()
    def fit(key):
        g, f = key
        try:
            weights, diagnostic = fit_folded_sfs(matrices[g][np.isin(labels[g], selected[key])],
                                                pseudocount=.5, max_iter=max_iter)
        except ValueError as exc:
            raise IPyradError(f'SFS group {g}, held-out fold {f}: {exc}') from exc
        return key, weights, dict(group=g, fold=f, individuals=len(groups[g]),
                                 training_loci=selected[key].tolist(), weights=weights.tolist(), **diagnostic)
    try:
        from threadpoolctl import threadpool_limits
        limits = threadpool_limits(limits=1)
        workers = min(max(1, cores), 4, len(selected))
    except ImportError:
        limits, workers = nullcontext(), 1
    models, diagnostics = {}, []
    with limits, ThreadPoolExecutor(max_workers=workers) as pool:
        for key, weights, diag in pool.map(fit, selected):
            models[key] = weights; diagnostics.append(diag)
            logger.info('SFS {} fold {}: {} informative sites, {} iterations', key[0], key[1],
                        diag['informative_sites'], diag['iterations'])
    return models, diagnostics, nbytes


def call_block(block, groups, models, folds, min_depth, geno_q, site_q, wgs_indices):
    """Pure posterior/acceptance calculation; masks never change fitting input."""
    valid = block['valid']
    gp = np.full((*valid.shape, 3), np.nan)
    for group, ix in groups.items():
        for fold in np.unique(block['loci'] % folds):
            js = np.flatnonzero(block['loci'] % folds == fold)
            values = sfs_posterior_moments(block['x'][ix][:, js], valid[ix][:, js], models[group, int(fold)])['GP']
            gp[ix[:, None], js[None, :]] = values
    gt = np.argmax(np.where(np.isfinite(gp), gp, -1), axis=2)
    score = phred(np.max(np.where(valid[:, :, None], gp, 0), axis=2))
    # A conservative lower bound on P(any nonreference), avoiding an
    # independence assumption among the joint posterior's sample marginals.
    support = phred(np.max(np.where(valid, 1-gp[:, :, 0], 0), axis=0))
    keep = valid & (block['dp'] >= min_depth) & (score >= geno_q)
    keep &= (gt == 0) | (support >= site_q)[None, :]
    keep &= (block['kind'] == 2)[None, :] | (gt == 0)
    for i in wgs_indices:
        ad = block['ad'][i]
        total = ad.sum(axis=1)
        balance = np.divide(ad[:, 1], total, out=np.full(len(total), .5),
                            where=(ad >= 0).all(axis=1) & (total > 0))
        keep[i] &= (gt[i] != 1) | ((balance >= .2) & (balance <= .8))
    return np.where(keep, gt, -1), score, support, gp


class AcceptedSpans:
    """Coalesce accepted spans online; at most one pending span per sample."""
    def __init__(self, handles):
        self.handles, self.pending = handles, [None]*len(handles)

    def add(self, sample, chrom, start, end):
        old = self.pending[sample]
        if old is not None and old[0] == chrom and start <= old[2]:
            self.pending[sample] = (chrom, old[1], max(end, old[2]))
        else:
            self.flush(sample)
            self.pending[sample] = (chrom, start, end)

    def flush(self, sample):
        span = self.pending[sample]
        if span is not None:
            self.handles[sample].write('\t'.join(map(str, span))+'\n')
            self.pending[sample] = None


def run_sfs_calling(*, tmpdir, names, group_file, retained_beds, min_depth,
                    geno_q, site_q, fit_loci=1024, folds=5, max_iter=10000,
                    cores=1, threads=1, wgs_samples=()):
    """Replace filtered variant decisions, retaining native inputs and masks."""
    started = time.monotonic()
    tmpdir = Path(tmpdir); vcfs = tmpdir/'vcfs'; beds = tmpdir/'beds'
    groups = {}
    if group_file is None:
        groups = {s: np.array([i]) for i, s in enumerate(names)}
    else:
        assignments = dict(line.split() for line in Path(group_file).read_text().splitlines() if line.strip())
        if set(assignments) != set(names):
            raise IPyradError('SFS population assignments must cover every assembled sample exactly.')
        for i, sample in enumerate(names):
            groups.setdefault(assignments[sample], []).append(i)
        groups = {g: np.array(ix) for g, ix in groups.items()}
    if any(len(ix) > 128 for ix in groups.values()):
        raise IPyradError('SFS currently supports at most 128 individuals per explicit group.')
    max_group = max(map(len, groups.values()))
    # Forward/backward polynomial arrays scale quadratically in group size.
    block_sites = max(1, min(2048, (32*1024**2)//(64*(max_group+1)**2),
                             (8*1024**2)//(128*len(names))))
    paths = [vcfs/(p.stem+'.pileup.bcf') for p in sorted(beds.glob('chunk-*.bed'), key=lambda p: int(p.stem.split('-')[1]))]
    source = NativeSource(paths, names, beds/'loci.bed', tmpdir/'REF_info.txt', retained_beds, block_sites)
    models, fits, training_bytes = fit_models(source, groups, folds=folds, fit_loci=fit_loci,
                                              max_iter=max_iter, cores=cores)
    fitting_seconds = time.monotonic()-started
    native = vcfs/'loci.filtered.vcf.gz'
    plain = vcfs/'loci.sfs.vcf'
    counts = Counter(sfs_bases=0, unseen_bases=0, native_fallback_records=0,
                     accepted_sfs_cells=0, rejected_sfs_cells=0, accepted_native_cells=0)
    acceptance = {s: beds/f'{s}.sfs.accepted.bed' for s in names}
    wgs_indices = [names.index(s) for s in wgs_samples]
    post_start = time.monotonic()
    with ExitStack() as stack:
        lines = stack.enter_context(stream_pipeline_lines([[BIN_BCF, 'view', '--no-version', str(native)]]))
        meta, header = _read_header(lines)
        if header.split('\t')[9:] != names:
            raise IPyradError('Native fallback sample order differs from SFS sample order.')
        output = stack.enter_context(plain.open('x'))
        meta = [x for x in meta if not x.startswith('##FORMAT=<ID=GQ,')]
        meta += ['##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality: posterior Phred confidence for SFS; native GQ for native fallback">',
                 '##INFO=<ID=IPYRAD_CALLER,Number=1,Type=String,Description="Genotype decision model: sfs or native fallback">',
                 '##ipyrad2_sfs_site_quality=Conservative marginal lower bound on any nonreference posterior; native QUAL for fallback']
        if not any(x.startswith('##FORMAT=<ID=AD,') for x in meta):
            meta.append('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Native allele depths on the named allele axis">')
        output.write('\n'.join([*meta, header])+'\n')
        accepted = AcceptedSpans([stack.enter_context(p.open('x')) for p in acceptance.values()])
        fallback = _coordinate_groups(lines, source.contigs)
        current = next(fallback, None)
        for block in source.blocks():
            gt, score, support, _ = call_block(block, groups, models, folds, min_depth, geno_q, site_q, wgs_indices)
            for j, (key, locus, rows) in enumerate(block['records']):
                while current is not None and current[0] < key:
                    # A called indel may be shifted by bcftools call. Refuse
                    # silent loss rather than guessing an unsupported mapping.
                    raise IPyradError('Native fallback contains a coordinate absent from its captured evidence.')
                original = current[1] if current is not None and current[0] == key else []
                if original:
                    current = next(fallback, None)
                chrom, pos = rows[0][0], key[1]-1
                if block['kind'][j]:
                    counts['sfs_bases'] += 1
                    counts['unseen_bases'] += int(block['kind'][j] == 1)
                    keep = gt[:, j] >= 0
                    counts['accepted_sfs_cells'] += int(keep.sum())
                    counts['rejected_sfs_cells'] += int((~keep).sum())
                    for i in np.flatnonzero(keep):
                        accepted.add(i, chrom, pos, pos+1)
                    if block['kind'][j] != 2:
                        continue
                    _, axis, _ = allele_axis(rows)
                    alleles = [rows[0][3], *rows[0][4].split(',')]
                    cells = []
                    for i in range(len(names)):
                        g = int(gt[i, j]); call = ('0/0', '0/1', '1/1')[g] if g >= 0 else './.'
                        depth = str(int(block['dp'][i, j])) if block['dp'][i, j] >= 0 else '.'
                        ad = ','.join(str(int(v)) if v >= 0 else '.' for v in block['ad'][i, j])
                        q = str(int(np.floor(score[i, j]))) if block['valid'][i, j] else '.'
                        cells.append(':'.join([call, depth, ad, q]))
                    output.write('\t'.join([chrom, str(pos+1), '.', alleles[axis[0]], alleles[axis[1]],
                                             str(float(support[j])), 'PASS', 'IPYRAD_CALLER=sfs',
                                             'GT:DP:AD:GQ', *cells])+'\n')
                else:
                    allowed = source.allows(chrom, np.array([pos]))[:, 0]
                    if not original:
                        keep, _ = reference_fallback(rows, allowed, min_depth, geno_q)
                        counts['accepted_native_reference_cells'] += int(keep.sum())
                        for i in np.flatnonzero(keep):
                            accepted.add(i, chrom, pos, pos+1)
                    for row in original:
                        row = list(row)
                        tags = row[8].split(':')
                        for i, cell in enumerate(row[9:]):
                            fields = cell.split(':'); values = dict(zip(tags, fields))
                            g = values.get('GT', './.')
                            try:
                                native_ok = (int(values.get('DP', '.')) >= max(5, min_depth)
                                             and float(values.get('GQ', '.')) >= geno_q
                                             and float(row[5]) >= site_q)
                            except ValueError:
                                native_ok = False
                            if (not allowed[i] or not native_ok) and 'GT' in tags:
                                fields[tags.index('GT')] = './.'; g = './.'
                                row[9+i] = ':'.join(fields)
                            if row[6] == 'PASS' and '.' not in g:
                                accepted.add(i, chrom, pos, pos+len(row[3]))
                                counts['accepted_native_cells'] += 1
                        row[7] = (row[7]+';' if row[7] != '.' else '')+'IPYRAD_CALLER=native'
                        output.write('\t'.join(row)+'\n')
                        counts['native_fallback_records'] += 1
        if current is not None or next(fallback, None) is not None:
            raise IPyradError('Unconsumed native fallback records after SFS capture evaluation.')
        for i in range(len(names)):
            accepted.flush(i)
    staged = vcfs/'loci.sfs.vcf.gz'
    run_pipeline([[BIN_BCF, '+fill-tags', str(plain), '-Oz', '-o', str(staged), '--', '-t', 'AC,AN,AF,MAF,F_MISSING']])
    run_pipeline([[BIN_BCF, 'index', '-c', str(staged)]])
    native.rename(vcfs/'loci.native.filtered.vcf.gz')
    native.with_suffix('.gz.csi').rename(vcfs/'loci.native.filtered.vcf.gz.csi')
    staged.replace(native); staged.with_suffix('.gz.csi').replace(native.with_suffix('.gz.csi'))
    # Complement accepted spans across the complete shared locus domain,
    # including positions absent from the capture (zero reads/non-ACGT).
    for sample, path in acceptance.items():
        run_pipeline([[BIN_BED, 'subtract', '-a', str(beds/'loci.bed'), '-b', str(path),
                       '-sorted', '-g', str(tmpdir/'REF_info.txt')]], beds/f'{sample}.sfs.mask.bed')
    result = dict(model='folded cross-fitted SFS confidence', counts=dict(counts),
                  groups={g: [names[i] for i in ix] for g, ix in groups.items()},
                  folds=folds, max_training_loci_per_fold=fit_loci, max_iterations=max_iter,
                  pseudocount=.5, native_min_depth=max(5, min_depth), fits=fits,
                  training_bytes=training_bytes, posterior_block_sites=block_sites,
                  fitting_domain='Post-paralog shared RAD loci; usable DP>=1, independent of output depth/confidence thresholds.',
                  quality='GQ=-10 log10(1-max GP); QUAL conservative lower bound using max_i P(G_i != REF). Native variant fallback retains native GQ/QUAL. Reference-only multiallelic fallback uses all native PL with a uniform genotype prior, requires a unique best REF genotype and DP>=max(5,min_depth), and never overrides an existing native variant record.',
                  uncertainty='Conditional on fitted priors; scores are not empirical accuracy guarantees.',
                  timing_seconds=dict(fitting_and_training_reads=fitting_seconds,
                                      posterior_and_acceptance=time.monotonic()-post_start))
    (tmpdir/'sfs.calling.json').write_text(json.dumps(result, indent=2)+'\n')
    return result
