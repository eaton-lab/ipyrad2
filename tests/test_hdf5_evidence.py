from pathlib import Path
import subprocess

import h5py
import numpy as np
import pandas as pd
import pytest

from ipyrad2.assembler.evidence_hdf5 import write_snp_evidence_hdf5
from ipyrad2.assembler.genotype_evidence import BIN_BCF
from ipyrad2.assembler.write_snps import _chunk_query_to_gt_arrays
from ipyrad2.analysis.evidence import GenotypeEvidenceReader, export_evidence_bcf
from ipyrad2.analysis.extracters.snps_extracter import run_snps_extracter
from ipyrad2.utils.exceptions import IPyradError


HEADER = '''##fileformat=VCFv4.3
##contig=<ID=chr1,length=200>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Genotype likelihoods">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tRAD\tWGS
'''


def fixture(tmp_path):
    raw = tmp_path/'raw.vcf'
    final = tmp_path/'final.vcf'
    lines = [
        'chr1\t10\t.\tA\tC\t40\tPASS\t.\tGT:DP:AD:GQ:PL\t0/1:8:4,4:30:30000,0,30000\t0/0:1:1,0:3:0,3,30',
        'chr1\t20\t.\tA\tC,G\t40\tPASS\t.\tGT:DP:AD:GQ:PL\t1/2:6:0,2,4:20:60,20,30,30,0,20\t./.:0:0,0,0:.:.',
        'chr1\t30\t.\tA\tG,C\t40\tPASS\t.\tGT:DP:AD:GQ:PL\t1/2:6:0,2,4:20:60,20,30,30,0,20\t./.:0:0,0,0:.:.',
        'chr1\t40\t.\tA\tC\t40\tPASS\t.\tGT:DP:AD:GQ:PL\t0/1:8:4,4:30:30,0,30\t0/0:1:1,0:.:.',
    ]
    raw.write_text(HEADER+'\n'.join(lines)+'\n')
    rows = [x.split('\t') for x in lines]
    rows[0][10] = './.:1'
    rows[2][4] = 'C,G'
    final.write_text(HEADER+'\n'.join('\t'.join(x) for x in rows)+'\n')
    query = pd.DataFrame([[x[0],x[1],x[5],x[3],x[4],*[':'.join(v.split(':')[:2]) for v in x[9:]]] for x in rows])
    genos, ref, _, dp, qual = _chunk_query_to_gt_arrays(query)
    db = tmp_path/'test.hdf5'
    with h5py.File(db,'w') as io:
        io.attrs['version'] = 2.0
        io.attrs['names'] = ['RAD','WGS']
        io.attrs['scaffold_names'] = ['chr1']
        io.attrs['scaffold_lengths'] = [200]
        io.attrs['nsnps'] = 4
        io.create_dataset('genos',data=genos).attrs['names'] = ['RAD','WGS']
        io.create_dataset('reference',data=ref)
        io.create_dataset('sample_dp',data=dp)
        io.create_dataset('site_qual',data=qual)
        io.create_dataset('snpsmap',data=np.array([[i,0,0,0,int(x[1])-1] for i,x in enumerate(rows)],dtype='u4'))
    return db,raw,final


def test_lossless_evidence_and_standard_arrays(tmp_path):
    db,raw,final = fixture(tmp_path)
    before = final.read_bytes(),raw.read_bytes()
    with h5py.File(db) as io:
        arrays = {k:io[k][:] for k in io}
    write_snp_evidence_hdf5(database=db,raw_vcf=raw,final_vcf=final,chunk_sites=2)
    with h5py.File(db) as io:
        for k,v in arrays.items():
            np.testing.assert_array_equal(v,io[k][:])
        assert io.attrs['format_version']=='3.0'
    with GenotypeEvidenceReader(db) as reader:
        data = reader.read([1,0,3,2],['WGS','RAD'])
        assert data['PL'][1,1,:3].tolist()==[30000,0,30000]
        assert data['AD'][1,0,:3].tolist()==[0,2,4]
        assert data['PL'][1,0,:6].tolist()==[60,20,30,30,0,20]
        assert data['GT'][0,1].tolist()==[255,255]
        assert data['usable_PL'][0,1]  # Covered evidence survives missing GT.
        assert data['status'][0,0]==2  # Zero reads do not produce a posterior.
        assert not data['usable_PL'][0,0]
        assert data['status'][0,2]==5  # Positive DP but no usable PL.
        assert (data['status'][:,3]==4).all()  # Reordered ALT is not silently permuted.
        assert (data['PL'][:,3]==-1).all()
    assert before==(final.read_bytes(),raw.read_bytes())


def test_masks_optout_and_invalid_evidence(tmp_path):
    db,raw,final=fixture(tmp_path)
    mask=tmp_path/'hard.bed';mask.write_text('chr1\t9\t10\n')
    write_snp_evidence_hdf5(database=db,raw_vcf=raw,final_vcf=final,hard_masks=({'RAD':mask},))
    with GenotypeEvidenceReader(db) as reader:
        d=reader.read([0]);assert d['status'][0,0]==3
        assert np.all(d['PL'][0,0]==-1) and d['DP'][0,0]==-1
        assert not d['usable_PL'][0,0]
    other=tmp_path/'optout';other.mkdir();db,raw,final=fixture(other)
    write_snp_evidence_hdf5(database=db,final_vcf=final)
    with h5py.File(db) as io:
        assert io.attrs['evidence_scope']=='none' and 'evidence' not in io
        assert io['alleles'][2,:3].tolist()==list(map(ord,'ACG'))
    with pytest.raises(IPyradError,match='--evidence snps or --evidence loci'):
        GenotypeEvidenceReader(db)
    broken=tmp_path/'broken';broken.mkdir();db,raw,final=fixture(broken)
    raw.write_text(raw.read_text().replace('30000,0,30000','30000,0'))
    with pytest.raises(IPyradError,match='cardinality'):
        write_snp_evidence_hdf5(database=db,raw_vcf=raw,final_vcf=final)
    with h5py.File(db) as io:
        assert 'evidence' not in io and '_snp_v3_staging' not in io


def test_bcf_selection_masks_and_evidence_roundtrip(tmp_path):
    db,raw,final=fixture(tmp_path)
    write_snp_evidence_hdf5(database=db,raw_vcf=raw,final_vcf=final)
    with h5py.File(db) as io:
        selected=io['snpsmap'][[0,3]]
    bcf=tmp_path/'out.bcf'
    export_evidence_bcf(database=db,output=bcf,samples=['WGS','RAD'],snpsmap=selected,
                        genotypes=np.array([[255,255],[1,1]],dtype='u1'))
    assert Path(str(bcf)+'.csi').exists()
    result=subprocess.check_output([BIN_BCF,'query','-f',r'%POS[\t%SAMPLE=%GT:%AD:%PL:%FT]\n',str(bcf)],text=True)
    assert result.splitlines()==[
        '10\tWGS=./.:1,0:0,3,30:NoHardCall\tRAD=0/1:4,4:30000,0,30000:PASS',
        '40\tWGS=./.:1,0:.:ExportMask\tRAD=0/1:4,4:30,0,30:PASS',
    ]
    with pytest.raises(IPyradError,match='already exists'):
        export_evidence_bcf(database=db,output=bcf,samples=['WGS','RAD'],snpsmap=selected,
                            genotypes=np.array([[255,255],[1,1]],dtype='u1'))
    kwargs=dict(data=db,name='snps',outdir=tmp_path/'snpex',min_sample_coverage=1,
                max_sample_missing=1.,min_minor_allele_frequency=0.,imap=None,minmap=None)
    with pytest.raises(IPyradError,match='does not support imputation'):
        run_snps_extracter(**kwargs,write_bcf=True,impute_method='sample')


def test_empty_database_has_valid_empty_evidence(tmp_path):
    db,raw,final=fixture(tmp_path)
    with h5py.File(db,'a') as io:
        for key in list(io.keys()):
            arr=io[key][:]
            del io[key]
            if key in ('genos','sample_dp'):
                arr=arr[:,:0]
            else:
                arr=arr[:0]
            ds=io.create_dataset(key,data=arr)
            if key=='genos':ds.attrs['names']=['RAD','WGS']
    write_snp_evidence_hdf5(database=db,raw_vcf=raw,final_vcf=final)
    with GenotypeEvidenceReader(db) as reader:
        assert reader.read([])['PL'].shape==(2,0,10)


@pytest.mark.parametrize('subsample', [False, True])
def test_snpex_bcf_matches_selected_matrices_and_native_likelihoods(tmp_path, subsample):
    db, raw, final = fixture(tmp_path)
    write_snp_evidence_hdf5(database=db, raw_vcf=raw, final_vcf=final)
    with h5py.File(db, 'a') as io:
        # Put both biallelic sites in one linkage block to exercise subsampling.
        io['snpsmap'][3, 0] = 0
    args = dict(data=db, name='snps', min_sample_coverage=1,
                max_sample_missing=1., min_minor_allele_frequency=0.,
                imap={'pop': ['R*', 'W*']}, minmap=None,
                subsample=subsample, random_seed=73, min_genotype_depth=2)
    run_snps_extracter(**args, outdir=tmp_path/'plain')
    run_snps_extracter(**args, outdir=tmp_path/'bcf', write_bcf=True)
    for suffix in ('genos.npy', 'snps.npy', 'snpsmap.tsv', 'samples.txt',
                   'sample_data_summary.tsv'):
        assert (tmp_path/'plain'/f'snps.{suffix}').read_bytes() == (
            tmp_path/'bcf'/f'snps.{suffix}').read_bytes()
    prefix = tmp_path/'bcf/snps'
    names = prefix.with_suffix('.samples.txt').read_text().splitlines()
    coords = pd.read_csv(prefix.with_suffix('.snpsmap.tsv'), sep='\t')
    genos = np.load(prefix.with_suffix('.genos.npy'))
    bcf = prefix.with_suffix('.bcf')
    assert Path(str(bcf)+'.csi').exists()
    assert subprocess.check_output([BIN_BCF, 'query', '-l', str(bcf)], text=True).splitlines() == names
    rows = subprocess.check_output(
        [BIN_BCF, 'query', '-f', r'%POS[\t%GT:%PL:%FT]\n', str(bcf)], text=True).splitlines()
    assert len(rows) == len(coords) == (1 if subsample else 2)
    for j, row in enumerate(rows):
        pos, *cells = row.split('\t')
        assert int(pos) == coords.iloc[j]['pos'] + 1
        for i, cell in enumerate(cells):
            gt, pl, status = cell.split(':')
            assert gt == {0: '0/0', 1: '0/1', 2: '1/1', 255: './.'}[genos[i, j]]
            if names[i] == 'RAD':
                assert pl == ('30000,0,30000' if pos == '10' else '30,0,30')
            elif pos == '10':
                assert pl == '0,3,30' and status == 'NoHardCall'
            else:
                assert pl == '.' and status == 'ExportMask'
    with pytest.raises(IPyradError, match='already exists'):
        run_snps_extracter(**args, outdir=tmp_path/'bcf', write_bcf=True)


def test_four_alleles_and_duplicate_native_records(tmp_path):
    db,raw,final=fixture(tmp_path)
    old='A\tC,G\t40\tPASS\t.\tGT:DP:AD:GQ:PL\t1/2:6:0,2,4:20:60,20,30,30,0,20\t./.:0:0,0,0:.:.'
    new='A\tC,G,T\t40\tPASS\t.\tGT:DP:AD:GQ:PL\t1/2:6:0,2,4,0:20:60000,20000,30000,30000,0,20000,40000,50000,60000,70000\t./.:0:0,0,0,0:.:.'
    raw.write_text(raw.read_text().replace(old,new))
    final.write_text(final.read_text().replace(old,new))
    lines=raw.read_text().splitlines()
    i=next(i for i,x in enumerate(lines) if x.startswith('chr1\t40\t'))
    lines.insert(i,lines[i]);raw.write_text('\n'.join(lines)+'\n')
    write_snp_evidence_hdf5(database=db,raw_vcf=raw,final_vcf=final)
    with GenotypeEvidenceReader(db) as reader:
        data=reader.read([1,3])
        assert data['allele_count'][0]==4
        assert data['AD'][0,0].tolist()==[0,2,4,0]
        assert data['PL'][0,0].tolist()==[60000,20000,30000,30000,0,20000,40000,50000,60000,70000]
        assert (data['status'][:,1]==4).all()
        with pytest.raises(IPyradError,match='integers'):
            reader.read([0.5])
