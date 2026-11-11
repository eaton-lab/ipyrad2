import hashlib
import itertools
import math

import h5py
import numpy as np
import pandas as pd
import pytest

from ipyrad2.analysis.methods.popgen.likelihoods import likelihoods_from_pl, posterior_moments
from ipyrad2.analysis.methods.popgen.evidence_backend import run_evidence_popgen
from ipyrad2.analysis.methods.popgen.runner import run_popgen_method
from ipyrad2.assembler.evidence_hdf5 import write_snp_evidence_hdf5
from ipyrad2.utils.exceptions import IPyradError
from test_hdf5_evidence import fixture


def brute_force(x):
    n = len(x)
    states = np.array(list(itertools.product(range(3), repeat=n)))
    k = states.sum(axis=1)
    # Integral p^K(1-p)^(2n-K) dp, including each heterozygote's factor 2.
    weights = np.array([math.exp(math.lgamma(int(a)+1)+math.lgamma(2*n-int(a)+1)-math.lgamma(2*n+2)) for a in k])
    weights *= 2.**(states == 1).sum(axis=1)
    weights *= np.prod(np.array([x[i, states[:, i]] for i in range(n)]), axis=0)
    weights /= weights.sum()
    gp = np.array([[weights[states[:, i] == g].sum() for g in range(3)] for i in range(n)])
    diversity = weights @ (2*k*(2*n-k)/(2*n*(2*n-1)))
    return gp, diversity, weights @ k/(2*n), (weights @ k**2-(weights @ k)**2)/(2*n)**2


def test_joint_moments_match_exhaustive_genotype_integration():
    rng = np.random.default_rng(27)
    pl = rng.integers(0, 80, (4, 5, 3))
    good = np.ones((4, 5), bool);good[3, 2] = False
    before = pl.copy()
    x, good = likelihoods_from_pl(pl, good)
    result = posterior_moments(x, good, 'population')
    for site in range(5):
        expected = brute_force(x[good[:, site], site])
        np.testing.assert_allclose(result['GP'][good[:, site], site], expected[0], atol=1e-12)
        np.testing.assert_allclose([result[k][site] for k in ['gene_diversity','allele_frequency','allele_frequency_variance']], expected[1:], atol=1e-12)
    np.testing.assert_array_equal(pl, before)
    assert np.isnan(result['GP'][3, 2]).all()


def test_symmetry_missingness_and_certain_calls():
    pl = np.array([[[0,70000,70000],[0,0,0],[-1,-1,-1]],
                   [[70000,0,70000],[0,0,0],[0,0,0]],
                   [[70000,70000,0],[0,0,0],[0,0,0]]], dtype=np.int32)
    valid = np.ones((3,3), bool);valid[:,2] = False
    for prior in ('uniform','population'):
        x, v = likelihoods_from_pl(pl, valid)
        a = posterior_moments(x, v, prior)
        b = posterior_moments(x[..., ::-1], v, prior)
        np.testing.assert_allclose(a['GP'], b['GP'][..., ::-1], equal_nan=True, atol=1e-12)
        np.testing.assert_allclose(a['allele_frequency'], 1-b['allele_frequency'], equal_nan=True)
        assert a['gene_diversity'][0] == pytest.approx(.6)
        assert np.isnan(a['gene_diversity'][2]) and np.isnan(a['GP'][:,2]).all()
        assert np.allclose(a['GP'][:,0], np.eye(3))
    # Shared-frequency covariance matters, even when marginals are identical.
    x = np.ones((2,1,3));v = np.ones((2,1),bool)
    assert posterior_moments(x,v,'population')['gene_diversity'][0] == pytest.approx(1/3)
    assert posterior_moments(x,v,'uniform')['gene_diversity'][0] == pytest.approx(4/9)


def kwargs(db):
    return dict(data=db, requested_stats=['heterozygosity'], min_sample_coverage=1,
                max_sample_missing=1., min_minor_allele_frequency=0., imap={'p':['RAD','WGS']},
                minmap=None, exclude=None, include_reference=False, subsample_unlinked=False,
                random_seed=None, cores=1, log_level='WARNING')


def test_database_backend_uses_missing_gt_evidence_without_modifying_database(tmp_path):
    db,raw,final = fixture(tmp_path)
    write_snp_evidence_hdf5(database=db, raw_vcf=raw, final_vcf=final)
    digest = hashlib.sha256(db.read_bytes()).hexdigest()
    a = run_evidence_popgen(**kwargs(db))
    rows = a.sample_stats.set_index('sample')
    assert rows.loc['WGS','sites_with_likelihood'] == 1  # missing GT but positive DP/PL
    assert rows.loc['WGS','canonical_called_sites_on_same_cells'] == 0
    assert rows.loc['RAD','sites_with_likelihood'] == 2
    assert a.summary['exported_snps'] == 2  # multi-allelic and NoMatch never enter
    assert 0 < rows.loc['WGS','posterior_heterozygosity'] < 1
    b = run_evidence_popgen(**kwargs(db), evidence_prior='population')
    assert np.isfinite(b.population_stats['posterior_snp_gene_diversity']).all()
    assert hashlib.sha256(db.read_bytes()).hexdigest() == digest
    with pytest.raises(IPyradError, match='one row per individual'):
        run_evidence_popgen(**kwargs(db), individual_map={'RAD':'same','WGS':'same'})
    bad = kwargs(db);bad['imap'] = None
    with pytest.raises(IPyradError, match='explicit --imap'):
        run_evidence_popgen(**bad, evidence_prior='population')
    bad = kwargs(db);bad['requested_stats']=['pi']
    with pytest.raises(IPyradError, match='not pi'):
        run_evidence_popgen(**bad)
    bad = kwargs(db);bad['min_minor_allele_frequency']=.01
    with pytest.raises(IPyradError, match='MAF filtering'):
        run_evidence_popgen(**bad)


def test_hard_masks_sample_coverage_and_empty_catalogue(tmp_path):
    db,raw,final = fixture(tmp_path)
    mask=tmp_path/'mask.bed';mask.write_text('chr1\t9\t40\n')
    write_snp_evidence_hdf5(database=db, raw_vcf=raw, final_vcf=final, hard_masks=({'RAD':mask},))
    a=run_evidence_popgen(**kwargs(db))
    assert a.sample_stats.set_index('sample').loc['RAD','sites_with_likelihood']==0
    k=kwargs(db);k['max_sample_missing']=.9
    b=run_evidence_popgen(**k)
    assert b.summary['samples_final']==['WGS'] and b.summary['exported_snps']==1
    k=kwargs(db);k['min_sample_coverage']=2
    c=run_evidence_popgen(**k)
    assert c.summary['exported_snps']==0
    assert c.sample_stats['posterior_heterozygosity'].isna().all()


def test_runner_keeps_default_call_outputs_identical(tmp_path):
    db,raw,final=fixture(tmp_path)
    write_snp_evidence_hdf5(database=db, raw_vcf=raw, final_vcf=final)
    args=kwargs(db);args.pop('requested_stats');args.update(stats='heterozygosity',name='out',force=False)
    run_popgen_method(**args, outdir=tmp_path/'default')
    run_popgen_method(**args, outdir=tmp_path/'calls', genotype_mode='calls')
    for suffix in ['sample_stats.tsv','population_stats.tsv']:
        assert (tmp_path/'default'/('out.'+suffix)).read_bytes()==(tmp_path/'calls'/('out.'+suffix)).read_bytes()
    run_popgen_method(**args, outdir=tmp_path/'evidence', genotype_mode='likelihoods', evidence_prior='population')
    table=pd.read_csv(tmp_path/'evidence/out.population_stats.tsv',sep='\t')
    assert 'posterior_snp_gene_diversity' in table
    assert 'not per-base diversity' in (tmp_path/'evidence/out.manifest.txt').read_text()


def test_bad_likelihoods_fail_without_silent_repair():
    with pytest.raises(ValueError, match='positive mass'):
        posterior_moments(np.zeros((2,1,3)), np.ones((2,1),bool))
    with pytest.raises(ValueError, match='shape'):
        likelihoods_from_pl(np.ones((3,3)))


def test_population_groups_are_isolated_and_patterns_expand(tmp_path):
    db, raw, final = fixture(tmp_path)
    write_snp_evidence_hdf5(database=db, raw_vcf=raw, final_vcf=final)
    args = kwargs(db)
    args['imap'] = {'rad': ['R*'], 'wgs': ['W*']}
    grouped = run_evidence_popgen(**args, evidence_prior='population')
    for group, sample in [('rad', 'RAD'), ('wgs', 'WGS')]:
        alone = dict(args, imap={group: [sample]})
        expected = run_evidence_popgen(**alone, evidence_prior='population')
        observed = grouped.sample_stats.set_index('sample').loc[sample]
        # Global site selection may include a missing cell, but cannot change
        # that individual's usable-site posterior or denominator.
        for field in ('sites_with_likelihood', 'posterior_heterozygosity',
                      'mean_posterior_alt_dosage'):
            assert observed[field] == pytest.approx(expected.sample_stats.iloc[0][field])
    args['imap'] = {'one': ['RAD', 'WGS'], 'two': ['RAD']}
    with pytest.raises(IPyradError, match='assigns sample.*multiple times'):
        run_evidence_popgen(**args, evidence_prior='population')


def test_native_depth_and_quality_filters_use_likelihood_denominators(tmp_path):
    db, raw, final = fixture(tmp_path)
    write_snp_evidence_hdf5(database=db, raw_vcf=raw, final_vcf=final)
    result = run_evidence_popgen(**kwargs(db), min_genotype_depth=2)
    rows = result.sample_stats.set_index('sample')
    assert rows.loc['RAD', 'sites_with_likelihood'] == 2
    assert rows.loc['WGS', 'sites_with_likelihood'] == 0
    assert np.isnan(rows.loc['WGS', 'posterior_heterozygosity'])
    result = run_evidence_popgen(**kwargs(db), min_site_qual=41)
    assert result.summary['exported_snps'] == 0
    assert result.population_stats['posterior_snp_gene_diversity'].isna().all()
