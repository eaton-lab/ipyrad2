"""Streaming concatenation agrees with the original typed pandas serialization."""
from pathlib import Path
import io
import pandas as pd
import pytest
from ipyrad2.denovo.cluster import concat_summaries

DTYPES=dict(sample='string',cluster_id='string',seed='string',length='Int64',cluster_length='Int64',n_unique='Int64',n_reads='Int64',record_type='string',cluster_sequence='string',arm_boundary='Int64')

@pytest.mark.parametrize('n',[0,1,4095,4096,4097,8200])
def test_chunk_boundaries_preserve_nullable_and_string_values(tmp_path,n):
    paths=[]
    for name,number in [('b',n),('a',2)]:
        rows=[]
        for i in range(number):
            rows.append(dict(sample='001' if i%2 else 'sample\tquoted',cluster_id=f'{i:05}',seed=f'{name};{i}',length=12,cluster_length=12,n_unique=None if i%3 else 1,n_reads=5,record_type='single',cluster_sequence='ACGTNNACGTAA',arm_boundary=12))
        path=tmp_path/f'{name}.summary.tsv';pd.DataFrame(rows,columns=DTYPES).to_csv(path,sep='\t',index=False);paths.append(path)
    expected=pd.concat([pd.read_csv(p,sep='\t',dtype=DTYPES) for p in sorted(paths)],ignore_index=True).to_csv(sep='\t',index=False)
    assert concat_summaries(tmp_path)==n+2
    assert (tmp_path/'concat.summary.tsv').read_text()==expected
    # A previous combined output cannot be included as another input sample.
    assert concat_summaries(tmp_path)==n+2
    assert (tmp_path/'concat.summary.tsv').read_text()==expected


def test_only_empty_summaries_keep_one_header(tmp_path):
    for n in ['a','b']:pd.DataFrame(columns=DTYPES).to_csv(tmp_path/f'{n}.summary.tsv',sep='\t',index=False)
    assert concat_summaries(tmp_path)==0
    assert (tmp_path/'concat.summary.tsv').read_text()=='\t'.join(DTYPES)+'\n'


def test_missing_and_inconsistent_inputs_fail(tmp_path):
    with pytest.raises(ValueError,match='No sample summaries'):concat_summaries(tmp_path)
    pd.DataFrame(columns=DTYPES).to_csv(tmp_path/'a.summary.tsv',sep='\t',index=False)
    pd.DataFrame(columns=['sample','seed']).to_csv(tmp_path/'b.summary.tsv',sep='\t',index=False)
    with pytest.raises(ValueError,match='Inconsistent'):concat_summaries(tmp_path)
