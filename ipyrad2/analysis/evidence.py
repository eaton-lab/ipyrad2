"""Shared, bounded access to native SNP evidence in HDF5 formats 3.0/4.0."""
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import h5py
import numpy as np

from ..utils.exceptions import IPyradError
from ..utils.parallel import run_pipeline


def _names(values):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in values]


class GenotypeEvidenceReader:
    """Read selected SNPs without loading the entire evidence matrix.

    AD/PL use the exact allele order in `alleles`. Negative entries are missing
    or padding, never zero likelihood. `usable_PL` enforces coverage and masks.
    Returned arrays retain the requested sample and SNP orders.
    """
    def __init__(self, path):
        self.io = h5py.File(path, 'r', rdcc_nbytes=256*1024**2, rdcc_nslots=1000003)
        try:
            if self.io.attrs.get('format_version') not in ('3.0', '4.0'):
                raise IPyradError('Native SNP evidence requires HDF5 format 3.0 or 4.0.')
            if self.io.attrs.get('evidence_scope') not in ('snps', 'rad_sites') or 'evidence' not in self.io:
                raise IPyradError('This HDF5 has no SNP evidence; assemble with --evidence snps or --evidence loci.')
            self.group = self.io['evidence']
            if not self.group.attrs.get('complete', False):
                raise IPyradError('HDF5 SNP evidence is incomplete.')
            self.names = _names(self.io['genos'].attrs['names'])
            self.scaffolds = _names(self.io.attrs['scaffold_names'])
            self.status_codes = {int(k): v for k, v in json.loads(self.group.attrs['status_codes']).items()}
            ns, nv = self.io['genos'].shape[:2]
            for key, shape in [('DP', (ns,nv)), ('AD',(ns,nv,4)), ('GQ',(ns,nv)),
                               ('PL',(ns,nv,10)), ('status',(ns,nv))]:
                if self.group[key].shape != shape:
                    raise IPyradError(f'HDF5 evidence/{key} shape does not match the SNP database.')
            if self.io['alleles'].shape != (nv,4):
                raise IPyradError('HDF5 allele catalogue shape does not match the SNP database.')
            self.datasets = {k: self.group[k] for k in ('DP','AD','GQ','PL','status')}
            self.genos = self.io['genos']
            self.alleles = self.io['alleles']
            self.snpsmap = self.io['snpsmap']
        except BaseException:
            self.io.close()
            raise

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.io.close()

    def read(self, sites, samples=None):
        raw_sites = np.asarray(sites)
        sites = np.asarray(sites, dtype=np.int64)
        if not np.array_equal(raw_sites, sites):
            raise IPyradError('SNP evidence indices must be integers.')
        if sites.ndim != 1 or np.any(sites < 0) or np.any(sites >= self.io['genos'].shape[1]):
            raise IPyradError('Invalid SNP evidence indices.')
        samples = self.names if samples is None else list(samples)
        if len(set(samples)) != len(samples) or any(s not in self.names for s in samples):
            raise IPyradError('Evidence requires unique, observed sample names; synthetic reference rows have no read evidence.')
        sidx = [self.names.index(s) for s in samples]
        unique, inverse = np.unique(sites, return_inverse=True)
        data = {key: self.datasets[key][:, unique][sidx][:, inverse]
                for key in ('DP', 'AD', 'GQ', 'PL', 'status')}
        data['alleles'] = self.alleles[unique][inverse]
        data['GT'] = self.genos[:, unique, :2][sidx][:, inverse]
        data['allele_count'] = (data['alleles'] != 0).sum(axis=1)
        lengths = data['allele_count'] * (data['allele_count'] + 1) // 2
        present = np.arange(10)[None,:] < lengths[:,None]
        data['usable_PL'] = ((data['DP'] > 0) & np.isin(data['status'], (0,1)) &
                             np.all((data['PL'] >= 0) | ~present[None,:,:], axis=2))
        return data

    def indices_for_map(self, selected):
        """Resolve an ordered snpex map, checking all five coordinate columns."""
        selected = np.asarray(selected)
        result = []
        target = 0
        for start in range(0, len(self.snpsmap), 16384):
            for offset, row in enumerate(self.snpsmap[start:start+16384]):
                if target == len(selected):
                    return np.asarray(result, dtype=np.int64)
                if np.array_equal(row, selected[target]):
                    result.append(start+offset)
                    target += 1
        if target != len(selected):
            raise IPyradError('Selected SNP map is not an ordered subset of the evidence database.')
        return np.asarray(result, dtype=np.int64)


def export_evidence_bcf(*, database, output, samples, snpsmap, genotypes, force=False):
    """Write indexed BCF with snpex's calls and unchanged native evidence.

    This exports the same filtered/subsampled SNPs as other snpex outputs.
    Native GQ is retained as caller provenance, not recalculated export quality.
    """
    output = Path(output)
    index = Path(str(output) + '.csi')
    if not force and (output.exists() or index.exists()):
        raise IPyradError(f'BCF output already exists: {output}')
    with GenotypeEvidenceReader(database) as reader:
        if any(s not in reader.names for s in samples):
            raise IPyradError('BCF evidence export cannot include a synthetic reference sample.')
        indices = reader.indices_for_map(snpsmap)
        if genotypes.shape != (len(samples), len(indices)):
            raise IPyradError('BCF calls and selected SNP evidence have different shapes.')
        output.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(prefix='.snpex-bcf-', dir=output.parent) as tmp:
            plain = Path(tmp)/'export.vcf'
            with plain.open('w') as handle:
                handle.write(f"##fileformat=VCFv4.3\n##source=ipyrad2-HDF5-{reader.io.attrs['format_version']}\n")
                for name, length in zip(reader.scaffolds, reader.io.attrs['scaffold_lengths']):
                    handle.write(f'##contig=<ID={name},length={int(length)}>\n')
                for key, number, kind, desc in (
                    ('GT','1','String','Exported hard genotype'),
                    ('DP','1','Integer','Native caller read depth'),
                    ('AD','R','Integer','Native caller allele depths in REF ALT order'),
                    ('GQ','1','Integer','Native caller genotype quality; not recalculated after export filtering'),
                    ('PL','G','Integer','Native phred-scaled genotype likelihoods'),
                    ('FT','1','String','Evidence status; ExportMask means snpex masked the hard call')):
                    handle.write(f'##FORMAT=<ID={key},Number={number},Type={kind},Description="{desc}">\n')
                handle.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)+'\n')
                for start in range(0, len(indices), 512):
                    stop = min(len(indices), start+512)
                    data = reader.read(indices[start:stop], samples)
                    quals = reader.io['site_qual'][indices[start:stop]]
                    for j in range(stop-start):
                        alleles = ''.join(map(chr, data['alleles'][j,:data['allele_count'][j]]))
                        ng = len(alleles)*(len(alleles)+1)//2
                        fields = []
                        for i in range(len(samples)):
                            dosage = int(genotypes[i,start+j])
                            gt = {0:'0/0',1:'0/1',2:'1/1',255:'./.'}.get(dosage)
                            if gt is None:
                                raise IPyradError('BCF export encountered an unsupported genotype dosage.')
                            status = reader.status_codes[int(data['status'][i,j])]
                            if status == 'HardMask':
                                gt = './.'
                            elif gt == './.' and np.all(data['GT'][i,j] != 255):
                                status = 'ExportMask'
                            def value(key, width=None):
                                vals = np.atleast_1d(data[key][i,j])
                                if width is not None:
                                    vals = vals[:width]
                                if np.all(vals < 0):
                                    return '.'
                                return ','.join('.' if x < 0 else str(int(x)) for x in vals)
                            fields.append(':'.join((gt,value('DP'),value('AD',len(alleles)),
                                                    value('GQ'),value('PL',ng),status)))
                        row = snpsmap[start+j]
                        qual = '.' if not np.isfinite(quals[j]) else str(float(quals[j]))
                        handle.write('\t'.join((reader.scaffolds[int(row[3])],str(int(row[4])+1),'.',
                                                alleles[0],','.join(alleles[1:]),qual,'PASS','.',
                                                'GT:DP:AD:GQ:PL:FT',*fields))+'\n')
            staged = Path(tmp)/'export.bcf'
            bcf = str(Path(sys.prefix)/'bin/bcftools')
            run_pipeline([[bcf,'view','--no-version','-Ob','-o',str(staged),str(plain)]])
            run_pipeline([[bcf,'index','-c',str(staged)]])
            staged.replace(output)
            Path(str(staged)+'.csi').replace(index)
    return output
