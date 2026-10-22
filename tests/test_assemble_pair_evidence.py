"""Real-tool regression for retained non-PROPER_PAIR evidence in both callers."""
from pathlib import Path
import random
import subprocess
import sys
import pytest
from ipyrad2.assembler.read_filters import prepare_filtered_analysis_bam
from ipyrad2.assembler.paralogs import call_vcf_from_bam
from ipyrad2.assembler.variants import get_group_called_variants_in_vcf_chunks

BIN = Path(sys.prefix) / 'bin'
pytestmark = pytest.mark.skipif(
    not all((BIN/x).exists() for x in ('samtools','bcftools')),
    reason='real samtools/bcftools required',
)

def run(*args):
    return subprocess.check_output(list(map(str,args)),text=True)

def make_bam(tmp, name, refseq, records):
    sam=tmp/f'{name}.sam'
    sam.write_text('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:500\n'
                   f'@RG\tID:{name}\tSM:{name}\n'+'\n'.join(records)+'\n')
    bam=tmp/f'{name}.bam'
    run(BIN/'samtools','sort','-o',bam,sam)
    run(BIN/'samtools','index',bam)
    return bam

def record(qname,flag,start,seq,mate,cigar='120M',mq=60,nm=0,sample='s'):
    return '\t'.join(map(str,[qname,flag,'chr1',start,mq,cigar,'=',mate,160 if start==41 else -160,seq,'I'*len(seq),f'NM:i:{nm}',f'RG:Z:{sample}']))

def query(vcf):
    out=run(BIN/'bcftools','query','-f','%POS[\t%GT:%DP:%AD]\n',vcf)
    return {int(line.split('\t')[0]):line.split('\t')[1:] for line in out.splitlines()}

@pytest.mark.parametrize('paired', [True,False,'proper'])
def test_both_callers_recover_retained_evidence_and_keep_overlap_correction(tmp_path,paired):
    rng=random.Random(81);seq=''.join(rng.choices('ACGT',k=500))
    ref=tmp_path/'ref.fa';ref.write_text('>chr1\n'+seq+'\n')
    run(BIN/'samtools','faidx',ref)
    bed=tmp_path/'loci.bed';bed.write_text('chr1\t0\t500\n')
    pos=101;alt=next(a for a in 'ACGT' if a!=seq[pos-1])
    records=[]
    for i in range(20):
        s1=list(seq[40:160]);s1[pos-41]=alt;s1=''.join(s1)
        records.append(record(f'good{i}',(99 if paired=='proper' else 97) if paired else 0,41,s1,81 if paired else 0,nm=1))
        if paired:
            s2=list(seq[80:200]);s2[pos-81]=alt
            records.append(record(f'good{i}',147 if paired=='proper' else 145,81,''.join(s2),41,nm=1))
    # Excluded records support the opposite allele; neither caller may count them.
    for flag in (0x200,0x400,0x800):
        for i in range(20):
            records.append(record(f'excluded{flag}_{i}',(97 if paired else 0)|flag,41,seq[40:160],81 if paired else 0))
    bam=make_bam(tmp_path,'s',seq,records)
    filtered=prepare_filtered_analysis_bam(sname='s',bam_file=bam,is_paired=paired,
        tmpdir=tmp_path/'filter',min_map_q=30,max_tlen=1000,max_softclip=5,
        max_nm=5,min_aligned_len=100,threads=1).bam_path
    stage=tmp_path/'final';(stage/'vcfs').mkdir(parents=True)
    final=get_group_called_variants_in_vcf_chunks(tmpdir=stage,reference=ref,
        bam_files=[filtered],locus_chunk=bed,min_map_q=30,min_base_q=13,threads=1)
    paralog=call_vcf_from_bam(filtered,bed,ref,tmp_path/'paralog.vcf.gz',min_map_q=30,min_base_q=13)
    for vcf in (final,paralog):
        call=query(vcf)[pos][0].split(':')
        assert call[0]=='1/1'
        assert int(call[1])==20  # overlapping mates count once, not twice
        assert call[2].split(',')[:2]==['0','20']
    if paired is True:
        # The pre-fix mpileup default really excludes this otherwise eligible pair.
        baseline=run(BIN/'bcftools','mpileup','-f',ref,'-q','30','-Q','13','-a','FORMAT/DP','-Ov',filtered)
        assert not any(not l.startswith('#') for l in baseline.splitlines())


def test_existing_read_filters_reject_short_poor_or_unsupported_pairs(tmp_path):
    seq='ACGT'*125;records=[]
    cases=[('good',97,60,'120M',0,81,160),
           ('proper',99,60,'120M',0,81,160),
           ('short',97,60,'60M',0,81,160),
           ('clipped',97,60,'60M60S',0,81,160),
           ('lowmap',97,10,'120M',0,81,160),
           ('edits',97,60,'120M',20,81,160),
           ('mate_unmapped',105,60,'120M',0,0,0),
           ('too_long',97,60,'120M',0,81,2000)]
    for name,flag,mq,cigar,nm,mate,tlen in cases:
        read=seq[40:100] if name=='short' else seq[40:160]
        fields=record(name,flag,41,read,mate,cigar,mq,nm).split('\t');fields[8]=str(tlen)
        records.append('\t'.join(fields))
    bam=make_bam(tmp_path,'s',seq,records)
    result=prepare_filtered_analysis_bam(sname='s',bam_file=bam,is_paired=True,
        tmpdir=tmp_path/'filter',min_map_q=30,max_tlen=1000,max_softclip=5,
        max_nm=5,min_aligned_len=100,threads=1)
    names={l.split('\t')[0] for l in run(BIN/'samtools','view',result.bam_path).splitlines()}
    assert names=={'good','proper'}


def test_overlap_mask_preserves_flags_sequence_and_handles_disagreeing_mates():
    from collections import Counter
    from ipyrad2.assembler.pair_overlaps import mask_group
    a=record('pair',97,41,'A'*120,81).split('\t')
    b=record('pair',145,81,'A'*120,41).split('\t')
    b[9]='C'+b[9][1:]
    original=[(x[1],x[9],x[5]) for x in (a,b)]
    stats=Counter();mask_group([a,b],stats)
    assert [(x[1],x[9],x[5]) for x in (a,b)]==original
    assert a[10][40]=='!' and b[10][0]=='!'  # conflicting overlapping bases masked
    assert all((a[10][40+i]=='!') != (b[10][i]=='!') for i in range(1,80))
    assert a[10][:40]=='I'*40 and b[10][80:]=='I'*40
    assert any(t=='OQ:Z:'+'I'*120 for t in a[11:])
    assert stats['anomalous_overlapping_pairs']==1


def test_overlap_mask_leaves_proper_pairs_and_different_read_groups_unchanged():
    from collections import Counter
    from ipyrad2.assembler.pair_overlaps import mask_group
    for flags,groups in [((99,147),('s','s')),((97,145),('s','other'))]:
        a=record('pair',flags[0],41,'A'*120,81,sample=groups[0]).split('\t')
        b=record('pair',flags[1],81,'A'*120,41,sample=groups[1]).split('\t')
        mask_group([a,b],Counter())
        assert a[10]==b[10]=='I'*120


def test_overlap_position_mapping_uses_cigar_gaps_and_clips():
    from ipyrad2.assembler.pair_overlaps import aligned_positions
    a=record('pair',97,41,'A'*12,81,cigar='2S3M1I2M2D4M').split('\t')
    assert aligned_positions(a)=={41:2,42:3,43:4,44:6,45:7,48:8,49:9,50:10,51:11}


@pytest.mark.parametrize('excluded_flag', [0x200,0x400,0x800])
def test_excluded_mate_cannot_mask_eligible_primary(excluded_flag):
    from collections import Counter
    from ipyrad2.assembler.pair_overlaps import mask_group
    a=record('pair',97,41,'A'*120,81).split('\t')
    b=record('pair',145|excluded_flag,81,'A'*120,41).split('\t')
    b[10]='J'*120
    mask_group([a,b],Counter())
    assert a[10]=='I'*120
    assert b[10]=='J'*120
