from pathlib import Path
import hashlib
import subprocess

import h5py
import numpy as np
import pytest

from ipyrad2.assembler.rad_evidence_hdf5 import write_rad_evidence_hdf5
from ipyrad2.assembler.evidence_hdf5 import write_snp_evidence_hdf5
from ipyrad2.analysis.rad_evidence import RadSiteEvidenceReader
from ipyrad2.analysis.evidence import GenotypeEvidenceReader
from ipyrad2.assembler.variants import BIN
from ipyrad2.utils.exceptions import IPyradError
from test_hdf5_evidence import fixture, HEADER


def setup(tmp_path):
    db,raw,final=fixture(tmp_path)
    write_snp_evidence_hdf5(database=db,raw_vcf=raw,final_vcf=final)
    ref=tmp_path/'ref.fa';ref.write_text('>chr1\nAAAAANAA'+'A'*192+'\n')
    subprocess.run([str(BIN/'samtools'),'faidx',str(ref)],check=True)
    bed=tmp_path/'rad.bed';bed.write_text('chr1\t0\t8\n')
    p=tmp_path/'pileup.vcf'
    rows=[
        'chr1\t1\t.\tA\t<*>\t0\t.\t.\tDP:AD:PL\t1:1,0:0,3,30\t0:0,0:0,0,0',
        'chr1\t2\t.\tA\tC,<*>\t0\t.\t.\tDP:AD:PL\t8:4,4,0:60,0,60,90,90,120\t4:4,0,0:0,12,120,12,120,120',
        'chr1\t3\t.\tA\tC,G,T,<*>\t0\t.\t.\tDP:AD:PL\t10:2,3,4,1,0:'+','.join(str(i*5000) for i in range(15))+'\t0:0,0,0,0,0:.',
        'chr1\t4\t.\tA\t<*>\t0\t.\t.\tDP:AD:PL\t10:10,0:0,30,255\t4:4,0:0,12,120',
        'chr1\t4\t.\tA\tAA\t0\t.\t.\tDP:AD:PL\t10:8,2:20,0,80\t4:4,0:0,12,120',
        'chr1\t5\t.\tA\t<*>\t0\t.\t.\tDP:AD:PL\t0:0,0:0,0,0\t4:.:.',
        'chr1\t8\t.\tA\t.\t0\t.\t.\tDP:AD:PL\t4:4:.\t4:4:.',
    ]
    p.write_text(HEADER+'\n'.join(rows)+'\n')
    mask=tmp_path/'hard.bed';mask.write_text('chr1\t1\t2\n')
    return db,ref,bed,p,mask


def snapshot(db):
    result={}
    with h5py.File(db) as h:
        def visit(name,x):
            if isinstance(x,h5py.Dataset):result[name]=x[()]
        h.visititems(visit)
    return result


def test_all_bases_native_vectors_masks_and_standard_datasets(tmp_path):
    db,ref,bed,raw,mask=setup(tmp_path);before=snapshot(db)
    inputs={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ref,bed,raw,mask]}
    metrics=write_rad_evidence_hdf5(database=db,reference=ref,bed=bed,native_files=[raw],native_regions=[bed],hard_masks=({'RAD':mask},),chunk_sites=2)
    after=snapshot(db)
    for k,v in before.items():np.testing.assert_array_equal(v,after[k])
    assert all(hashlib.sha256(p.read_bytes()).hexdigest()==h for p,h in inputs.items())
    assert metrics['sites']==8 and metrics['records']==7
    with RadSiteEvidenceReader(db) as r:
        d=r.read([3,0,1,2,4,5,6,7],['WGS','RAD'])
        assert d['record_site'].tolist()==[3,3,0,1,2,4,7]
        assert d['alleles'][0]==('A','<*>') and d['alleles'][1]==('A','AA')
        assert d['allele_count'][4]==5 and d['PL_count'][4]==15
        assert d['PL'][1,4,14]==70000
        assert d['PL'][1,2,:3].tolist()==[0,3,30]
        assert d['GQ'].min()==d['GQ'].max()==-1
        assert d['status'][1,3]==3 and not d['usable_PL'][1,3]
        assert d['PL'][1,3,:6].tolist()==[60,0,60,90,90,120] # preserved behind explicit hard mask
        assert d['DP'][0,2]==0 and d['PL'][0,2,:3].tolist()==[0,0,0]
        assert not d['usable_PL'][0,2] # zero PL never makes zero-read evidence usable
        assert d['status'][0,5]==2 and not d['usable_PL'][0,5]
        assert d['site_status'][5:7].tolist()==[2,0] # N and absent record remain distinct
        assert not d['is_base_alleles_only'][:6].any()
        assert r.read([])['PL'].shape==(2,0,0)
    with GenotypeEvidenceReader(db) as r:
        assert r.read([0])['PL'][0,0,:3].tolist()==[30000,0,30000]
    with pytest.raises(IPyradError,match='completed SNP-evidence|already exists'):
        write_rad_evidence_hdf5(database=db,reference=ref,bed=bed,native_files=[raw],native_regions=[bed])


@pytest.mark.parametrize('error',['PL','duplicate','order','reference','sample'])
def test_bad_native_evidence_fails_without_publishing(tmp_path,error):
    db,ref,bed,raw,mask=setup(tmp_path);before=snapshot(db);s=raw.read_text()
    if error=='PL':s=s.replace('0,3,30','0,3')
    if error=='duplicate':s+=s.splitlines()[-1]+'\n'
    if error=='order':s+='chr1\t1\t.\tA\t<*>\t0\t.\t.\tDP:AD:PL\t1:1,0:0,3,30\t0:0,0:0,0,0\n'
    if error=='reference':s=s.replace('chr1\t1\t.\tA','chr1\t1\t.\tG')
    if error=='sample':s=s.replace('\tRAD\tWGS','\tWGS\tRAD')
    raw.write_text(s)
    with pytest.raises(IPyradError):write_rad_evidence_hdf5(database=db,reference=ref,bed=bed,native_files=[raw],native_regions=[bed],chunk_sites=2)
    with h5py.File(db) as h:
        assert 'rad_evidence' not in h and '_rad_v4_staging' not in h
        assert h.attrs['format_version']=='3.0'
    for k,v in before.items():np.testing.assert_array_equal(v,snapshot(db)[k])


def test_empty_and_disjoint_intervals_and_shared_mask(tmp_path):
    db,ref,bed,raw,mask=setup(tmp_path);bed.write_text('chr1\t0\t2\nchr1\t6\t8\n');mask.write_text('chr1\t0\t8\n')
    write_rad_evidence_hdf5(database=db,reference=ref,bed=bed,native_files=[raw],native_regions=[bed],shared_mask=mask)
    with RadSiteEvidenceReader(db) as r:
        d=r.read([0,1,2,3]);assert d['site_map'][:,1].tolist()==[0,1,6,7]
        assert d['hard_mask'].all() and not d['usable_PL'].any()
        for bad in [[1,1],[.1],[-1],[4]]:
            with pytest.raises(IPyradError):r.read(bad)
        with pytest.raises(IPyradError):r.read([0],['RAD','RAD'])
    other=tmp_path/'empty';other.mkdir();db,ref,bed,raw,mask=setup(other);bed.write_text('')
    write_rad_evidence_hdf5(database=db,reference=ref,bed=bed,native_files=[raw],native_regions=[bed])
    with RadSiteEvidenceReader(db) as r:assert r.nsites==r.nrecords==0


def test_invalid_bed_is_not_silently_merged(tmp_path):
    db,ref,bed,raw,mask=setup(tmp_path);bed.write_text('chr1\t0\t8\nchr1\t4\t6\n')
    with pytest.raises(IPyradError,match='nonoverlapping'):
        write_rad_evidence_hdf5(database=db,reference=ref,bed=bed,native_files=[raw],native_regions=[bed])


def test_uncaptured_bases_are_not_labeled_as_uncovered(tmp_path):
    db,ref,bed,raw,mask=setup(tmp_path)
    partial=tmp_path/'partial.bed';partial.write_text('chr1\t0\t2\n')
    with pytest.raises(IPyradError,match='outside the captured'):
        write_rad_evidence_hdf5(database=db,reference=ref,bed=bed,native_files=[raw],native_regions=[partial])
    with h5py.File(db) as h:assert 'rad_evidence' not in h
