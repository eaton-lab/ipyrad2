"""All-base numerator/denominator and native symbolic-likelihood controls."""
import hashlib
import itertools
import math

import h5py
import numpy as np
import pytest

from ipyrad2.analysis.methods.popgen.locus_evidence_backend import run_locus_evidence_popgen, biallelic_locus_chunk
from ipyrad2.analysis.methods.popgen.runner import run_popgen_method
from ipyrad2.analysis.rad_evidence import RadSiteEvidenceReader
from ipyrad2.assembler.rad_evidence_hdf5 import write_rad_evidence_hdf5
from ipyrad2.utils.exceptions import IPyradError
from test_rad_site_evidence import setup
from test_popgen_evidence import kwargs, brute_force


def make(tmp_path):
    db,ref,bed,raw,mask=setup(tmp_path)
    write_rad_evidence_hdf5(database=db,reference=ref,bed=bed,native_files=[raw],native_regions=[bed],hard_masks=({'RAD':mask},))
    return db


def test_explicit_base_classification_native_projection_and_masks(tmp_path):
    db=make(tmp_path)
    with RadSiteEvidenceReader(db) as reader:
        pl,good,kind,reasons=biallelic_locus_chunk(reader,np.arange(8),reader.names)
        assert kind.tolist()==[1,2,0,0,1,0,0,0]
        assert dict(reasons)==dict(multiallelic_nucleotides=1,multiple_records=1,non_acgt_reference=1,no_record=1,no_alternative_likelihood=1)
        np.testing.assert_array_equal(pl[0,0],[0,3,30])
        np.testing.assert_array_equal(pl[1,1],[0,12,120])
        assert good.sum()==2 and good[0,0] and good[1,1]
        assert not good[0,1]  # native PL survives, but hard mask excludes it
        assert not good[:,4].any()  # DP zero and missing PL do not become invariant zero
        # Sparse and contiguous reader paths must preserve exact fields.
        x=reader.read(np.arange(8)); y=reader.read(np.arange(7,-1,-1))
        np.testing.assert_array_equal(x['hard_mask'][:,::-1],y['hard_mask'])
        for site in range(8):
            one=reader.read([site]); ix=np.where(x['record_site']==site)[0]
            np.testing.assert_array_equal(one['DP'],x['DP'][:,ix])


@pytest.mark.parametrize('prior',['uniform','population'])
def test_per_base_means_match_exhaustive_genotypes_and_leave_input_unchanged(tmp_path,prior):
    db=make(tmp_path);before=hashlib.sha256(db.read_bytes()).hexdigest()
    args=kwargs(db);args['requested_stats']=['pi','heterozygosity']
    result=run_locus_evidence_popgen(**args,evidence_prior=prior)
    p=result.population_stats.iloc[0]
    assert p.callable_bases==2 and p.callable_sample_cells==2
    assert p.callable_unseen_alt_bases==p.callable_observed_alt_bases==1
    means=[]
    for pl in [[0,3,30],[0,12,120]]:
        x=10.**(-np.array(pl)/10)
        # With one individual both priors integrate to uniform genotype weights.
        means.append(x[1]/x.sum())
    assert p.posterior_pi==pytest.approx(np.mean(means))
    assert p.posterior_pairwise_differences_sum==pytest.approx(sum(means))
    assert p.posterior_pi != pytest.approx(sum(means)/8) # uncovered/excluded bases aren't zeros
    assert result.sample_stats.bases_with_likelihood.tolist()==[1,1]
    assert hashlib.sha256(db.read_bytes()).hexdigest()==before


def test_population_covariance_certain_genotypes_and_depth_prior_sensitivity(tmp_path):
    db=make(tmp_path)
    with h5py.File(db,'a') as h:
        g=h['rad_evidence'];g['DP'][:,0]=10;g['status'][:,0]=0
        # Analytically certain REF/REF and ALT/ALT at a reference/unseen site.
        g['PL'][:,0:3]=[[0,70000,70000],[70000,70000,0]]
    args=kwargs(db);args['min_sample_coverage']=2
    for prior in ['uniform','population']:
        result=run_locus_evidence_popgen(**args,evidence_prior=prior)
        p=result.population_stats.iloc[0]
        assert p.callable_bases==1 and p.posterior_pi==pytest.approx(2/3)
    with h5py.File(db,'a') as h:h['rad_evidence/PL'][:,0:3]=0
    u=run_locus_evidence_popgen(**args)
    p=run_locus_evidence_popgen(**args,evidence_prior='population')
    assert u.population_stats.iloc[0].posterior_pi==pytest.approx(4/9)
    assert p.population_stats.iloc[0].posterior_pi==pytest.approx(1/3)
    # Flat positive-depth likelihoods retain uncertainty; they aren't certain REF.
    assert u.population_stats.iloc[0].posterior_pi>0
    with h5py.File(db,'a') as h:h['rad_evidence/PL'][:,0:3]=[[0,70000,70000]]*2
    assert run_locus_evidence_popgen(**args).population_stats.iloc[0].posterior_pi==0


def test_filters_empty_denominators_individuals_and_runner(tmp_path):
    db=make(tmp_path);args=kwargs(db)
    for update in [dict(min_minor_allele_frequency=.01),dict(min_site_qual=1),dict(subsample_unlinked=True),dict(random_seed=42)]:
        with pytest.raises(IPyradError,match='all callable bases'):run_locus_evidence_popgen(**dict(args,**update))
    with pytest.raises(IPyradError,match='one row per individual'):
        run_locus_evidence_popgen(**args,individual_map={'RAD':'same','WGS':'same'})
    result=run_locus_evidence_popgen(**dict(args,min_sample_coverage=2))
    assert result.population_stats.callable_bases.sum()==0
    assert result.population_stats.posterior_pi.isna().all()
    result=run_locus_evidence_popgen(**dict(args,min_genotype_depth=2))
    assert result.population_stats.iloc[0].callable_bases==1
    with pytest.raises(IPyradError,match='No samples remain'):
        run_locus_evidence_popgen(**dict(args,max_sample_missing=.5))
    runner=dict(args);runner.pop('requested_stats');runner.update(name='out',outdir=tmp_path/'result',stats='all',force=False,genotype_mode='likelihoods',evidence_scope='loci')
    run_popgen_method(**runner)
    assert 'posterior_pi' in (tmp_path/'result/out.population_stats.tsv').read_text()
    assert 'biallelic' in (tmp_path/'result/out.manifest.txt').read_text()


def test_hard_calls_do_not_determine_likelihood_callability(tmp_path):
    db=make(tmp_path)
    with h5py.File(db,'a') as h:
        h['genos'][:]=255
        g=h['rad_evidence'];g['hard_mask'][:,1]=0;g['status'][:,1]=0
    args=kwargs(db)
    result=run_locus_evidence_popgen(**args,evidence_prior='population')
    x=np.array([[60,0,60],[0,12,120]])
    expected=brute_force(10.**(-x/10))[1]
    first=(10**(-.3))/(1+10**(-.3)+10**(-3))
    assert result.population_stats.iloc[0].posterior_pi==pytest.approx((first+expected)/2)


def test_observed_pair_projection_uses_native_allele_indices(tmp_path):
    db=make(tmp_path)
    with RadSiteEvidenceReader(db) as r:
        before=biallelic_locus_chunk(r,np.arange(8),r.names)
    with h5py.File(db,'a') as h:
        g=h['rad_evidence'];a,b=g['allele_offsets'][1:3];x,y=g['PL_offsets'][1:3]
        ids=g['allele_ids'][a:b];g['allele_ids'][a:b]=ids[[0,2,1]]
        old=g['PL'][:,x:y];g['PL'][:,x:y]=old[:,[0,3,5,1,4,2]]
    with RadSiteEvidenceReader(db) as r:
        after=biallelic_locus_chunk(r,np.arange(8),r.names)
    for a,b in zip(before[:3],after[:3]):np.testing.assert_array_equal(a,b)
