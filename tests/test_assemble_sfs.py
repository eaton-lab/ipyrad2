"""Scientific and native-format integration checks for opt-in assembly SFS."""
import shutil
import subprocess

import numpy as np
import pytest

from ipyrad2.assembler.calling_options import resolve_calling_options
from ipyrad2.assembler.sfs_calling import NativeSource, fit_models, call_block, run_sfs_calling
from ipyrad2.assembler.variants import BIN_BCF, get_filtered_vcf
from ipyrad2.analysis.methods.popgen.sfs_prior import folded_saf, fit_folded_sfs
from ipyrad2.utils.exceptions import IPyradError
from test_popgen_sfs_prior import enumerate_states


def command(args):
    return subprocess.run(list(map(str, args)), check=True, capture_output=True, text=True).stdout


def fixture(root, block_sites=31):
    root.mkdir(exist_ok=True); beds=root/'beds'; beds.mkdir(); vcfs=root/'vcfs'; vcfs.mkdir()
    (root/'REF_info.txt').write_text('chr\t1000\n')
    bed=beds/'loci.bed';bed.write_text(''.join(f'chr\t{l*100}\t{l*100+40}\n' for l in range(10)))
    shutil.copyfile(bed,beds/'chunk-0.bed')
    header='''##fileformat=VCFv4.2
##contig=<ID=chr,length=1000>
##FILTER=<ID=PASS,Description="Passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tRAD\tWGS\tRAD2
'''
    rows=[]
    for l in range(10):
        for j in range(40):
            pos=l*100+j+1
            alt='<*>';cells=['1:1,0:0,3,40']*3
            if j%10 == 1:
                alt='C,<*>';cells=['12:6,6,0:60,0,60,90,90,120','12:6,6,0:60,0,60,90,90,120','8:8,0,0:0,24,240,24,240,240']
            if j == 20:
                cells=['1:1,0:0,0,0']*3  # uninformative evidence must not get fabricated certainty
            if j == 21:
                alt='<*>'
                cells=['0:0,0:0,0,0']*3
            if j == 22:
                continue  # no captured record; consensus must remain missing
            if j == 23:
                alt='C,G';cells=['10:4,3,3:70,50,40,30,0,80']*3
            if j == 24:
                alt='AA';cells=['10:5,5:80,0,80']*3
            rows.append(f'chr\t{pos}\t.\tA\t{alt}\t60\tPASS\t.\tDP:AD:PL\t'+'\t'.join(cells))
    raw=vcfs/'pileup.vcf';raw.write_text(header+'\n'.join(rows)+'\n')
    capture=vcfs/'chunk-0.pileup.bcf'
    command([BIN_BCF,'view','-Ob','-o',capture,raw])
    # Explicit native fallback calls: strong indel/multiallelic, one low-depth
    # sample to check that fallback does not inherit SFS's permissive floor.
    native=vcfs/'native.vcf'
    native.write_text(header+''.join(f'chr\t{l*100+j+1}\t.\tA\t{alt}\t60\tPASS\t.\tGT:DP:AD:GQ\t1/2:10:4,3,3:40\t1/2:2:0,1,1:40\t1/2:10:4,3,3:40\n' if j==23 else
        f'chr\t{l*100+j+1}\t.\tA\t{alt}\t60\tPASS\t.\tGT:DP:AD:GQ\t0/1:10:5,5:40\t0/1:2:1,1:40\t0/1:10:5,5:40\n'
        for l in range(10) for j,alt in [(23,'C,G'),(24,'AA')]))
    command([BIN_BCF,'view','-Oz','-o',vcfs/'loci.raw.vcf.gz',native])
    get_filtered_vcf(root,5,13,13,1)
    return NativeSource([capture],['RAD','WGS','RAD2'],bed,root/'REF_info.txt',block_sites=block_sites)


@pytest.mark.parametrize('model,depth,locus,expected',[
    ('native',None,None,(5,5)),('native',2,None,(2,2)),('native',2,8,(2,8)),
    ('sfs',None,None,(1,5)),('sfs',7,None,(7,5)),('sfs',2,3,(2,3))])
def test_depth_resolution(model,depth,locus,expected):
    assert resolve_calling_options(model,depth,locus,13,13)==expected


@pytest.mark.parametrize('depth,locus',[(0,None),(-1,None),(None,0),(True,None)])
def test_invalid_depth_rejected(depth,locus):
    with pytest.raises(IPyradError):resolve_calling_options('sfs',depth,locus,13,13)


@pytest.mark.parametrize('n',[1,2,3])
def test_crossfit_matches_frozen_recipe_and_independent_enumeration(tmp_path,n):
    source=fixture(tmp_path)
    groups={'cohort':np.arange(n),**{f's{i}':np.array([i]) for i in range(n,3)}}
    models,fits,_=fit_models(source,groups,folds=5,fit_loci=4,max_iter=10000,cores=2)
    blocks=list(source.blocks());x=np.concatenate([b['x'] for b in blocks],axis=1)
    valid=np.concatenate([b['valid'] for b in blocks],axis=1);loci=np.concatenate([b['loci'] for b in blocks])
    # The experiment selected eligible loci separately for every cohort/fold.
    for g,ix in groups.items():
        available=np.unique(loci[valid[ix].any(axis=0)])
        for f in range(5):
            pool=available[available%5!=f]
            selected=pool[np.unique(np.linspace(0,len(pool)-1,min(4,len(pool)),dtype=int))]
            ids=np.flatnonzero(np.isin(loci,selected)&valid[ix].any(axis=0))
            weights,_=fit_folded_sfs(folded_saf(x[ix][:,ids],valid[ix][:,ids]),pseudocount=.5,max_iter=10000)
            np.testing.assert_allclose(models[g,f],weights,atol=1e-13)
    for block in blocks:
        gt,score,support,gp=call_block(block,groups,models,5,1,13,13,[1])
        for j in range(len(block['loci'])):
            for g,ix in groups.items():
                _,expected,_,_=enumerate_states(block['x'][ix,j],block['valid'][ix,j],models[g,int(block['loci'][j]%5)])
                np.testing.assert_allclose(gp[ix,j],expected,atol=1e-12)
        assert not np.any((gt>=0)&~block['valid'])
        assert not np.any((gt>0)&(block['kind'][None,:]!=2))
        # Raising acceptance depth must not alter posterior inference.
        high=call_block(block,groups,models,5,5,13,13,[1])
        np.testing.assert_equal(gp,high[3]);assert not np.any((high[0]>=0)&(block['dp']<5))
    assert all(d['converged'] and not any(l%5==d['fold'] for l in d['training_loci']) for d in fits)


def test_native_formats_reference_masks_and_fallback(tmp_path):
    source=fixture(tmp_path)
    group=tmp_path/'populations.tsv';group.write_text('RAD\tP\nWGS\tP\nRAD2\tP\n')
    result=run_sfs_calling(tmpdir=tmp_path,names=source.names,group_file=group,
                          retained_beds=None,min_depth=1,geno_q=13,site_q=13,
                          fit_loci=4,cores=2,wgs_samples=['WGS'])
    assert result['counts']['native_fallback_records']==20
    assert result['counts']['accepted_native_cells']==40
    lines=command([BIN_BCF,'view','-H',tmp_path/'vcfs/loci.filtered.vcf.gz']).splitlines()
    native=[x.split('\t') for x in lines if 'IPYRAD_CALLER=native' in x]
    assert all(x[10].startswith('./.') for x in native)
    assert all(x[9].split(':')[0] in ('0/1','1/2') for x in native)
    masks=[line.split() for line in (tmp_path/'beds/RAD.sfs.mask.bed').read_text().splitlines()]
    def masked(p):return any(int(a)<=p<int(b) for c,a,b in masks)
    assert not masked(0)  # one-read reference accepted on learned SFS
    assert masked(21) and masked(22)  # no observations cannot inherit REF
    assert not masked(23) and not masked(24)  # conservative native fallback
    assert (tmp_path/'vcfs/loci.native.filtered.vcf.gz').is_file()


def test_stream_chunking_and_training_failures(tmp_path,monkeypatch):
    source=fixture(tmp_path,block_sites=7);groups={'all':np.arange(3)}
    models,_,_=fit_models(source,groups,folds=5,fit_loci=4,max_iter=10000,cores=1)
    source.block_sites=113
    other,_,_=fit_models(source,groups,folds=5,fit_loci=4,max_iter=10000,cores=2)
    for k in models:np.testing.assert_array_equal(models[k],other[k])
    with pytest.raises(IPyradError,match='group all.*fold.*converge'):
        fit_models(source,groups,folds=5,fit_loci=4,max_iter=1,cores=1)
    monkeypatch.setattr('ipyrad2.assembler.sfs_calling.TRAINING_BYTES',1)
    with pytest.raises(IPyradError,match='256 MiB'):
        fit_models(source,groups,folds=5,fit_loci=4,max_iter=10000,cores=1)


def test_retained_domain_controls_training_and_calls(tmp_path):
    source=fixture(tmp_path)
    keep=tmp_path/'rad.keep.bed';keep.write_text('chr\t0\t1\n')
    source=NativeSource(source.paths,source.names,tmp_path/'beds/loci.bed',tmp_path/'REF_info.txt',{'RAD':keep})
    blocks=list(source.blocks())
    assert sum(int(b['valid'][0].sum()) for b in blocks)==1
    with pytest.raises(IPyradError,match='no usable training loci'):
        fit_models(source,{s:np.array([i]) for i,s in enumerate(source.names)},folds=5,fit_loci=4,max_iter=10000,cores=1)


def test_discovery_depth_preserves_low_coverage_evidence(tmp_path):
    """A depth-one flank is retained in coverage but cannot define a locus."""
    from ipyrad2.assembler.beds import get_coverage_bed_graphs, BIN_SAM
    sam = tmp_path/'reads.sam'
    sam.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr\tLN:100\n'+
        ''.join(f'r{i}\t0\tchr\t{1 if i == 0 else 6}\t60\t20M\t*\t0\t0\t'+
                'A'*20+'\t'+'I'*20+'\n' for i in range(5)))
    bam = tmp_path/'reads.bam'
    command([BIN_SAM,'view','-b','-o',bam,sam])
    (tmp_path/'REF_info.txt').write_text('chr\t100\n')
    merged = get_coverage_bed_graphs('RAD',bam,False,tmp_path/'unused.fa',tmp_path,
        min_map_q=40,min_sample_depth=1,min_merge_distance=0,threads=1,min_locus_depth=5)
    assert merged.read_text() == 'chr\t5\t20\n'
    assert (tmp_path/'beds/RAD.fragments.bedgraph').read_text() == 'chr\t0\t5\t1\nchr\t5\t20\t5\nchr\t20\t25\t4\n'


def test_no_prior_or_depth_only_reference_fallback(tmp_path):
    source = fixture(tmp_path)
    groups = {'all':np.arange(3)}
    models,_,_ = fit_models(source,groups,folds=5,fit_loci=4,max_iter=10000,cores=1)
    block = next(source.blocks())
    gt,_,_,_ = call_block(block,groups,models,5,1,100,13,[1])
    # Even covered reference bases become missing when confidence fails.
    assert gt[0,0] == -1 and block['dp'][0,0] == 1
    block['valid'][:] = False
    gt,_,_,gp = call_block(block,groups,models,5,1,0,0,[1])
    assert (gt == -1).all() and np.isnan(gp).all()
