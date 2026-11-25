"""Independent enumeration of joint count posteriors and aggregate statistics."""
import itertools
import math

import h5py
import numpy as np
import pytest

from ipyrad2.analysis.methods.popgen.likelihoods import posterior_moments
from ipyrad2.analysis.methods.popgen.sfs_prior import sfs_posterior_moments
from ipyrad2.analysis.methods.popgen.sfs_fitting import CrossFittedSFS
from ipyrad2.analysis.methods.popgen.posterior_stats import PosteriorStats
from ipyrad2.analysis.methods.popgen.models import resolve_likelihood_stats
from ipyrad2.analysis.methods.popgen.runner import run_popgen_method
from ipyrad2.utils.exceptions import IPyradError
from test_popgen_locus_evidence import make
from test_popgen_evidence import kwargs


def exhaustive(x, good, prior, folded):
    n = len(x)
    values = np.zeros(2*n+1)
    for state in itertools.product(range(3), repeat=n):
        k = sum(state)
        multiplicity = 2**state.count(1)
        if prior == 'uniform':
            p = 3.**(-n)
        elif prior == 'population':
            p = multiplicity*math.factorial(k)*math.factorial(2*n-k)/math.factorial(2*n+1)
        else:
            p = folded[min(k, 2*n-k)]/(1 if k == n else 2)
            p *= multiplicity/math.comb(2*n, k)
        weight = p*np.prod([x[i, g] if good[i] else 1. for i, g in enumerate(state)])
        values[sum(g for g, use in zip(state, good) if use)] += weight
    return values/values.sum()


@pytest.mark.parametrize('prior', ['uniform', 'population', 'sfs'])
@pytest.mark.parametrize('n', [1, 2, 3, 4])
def test_joint_counts_against_exhaustive_states(prior, n):
    rng = np.random.default_rng(920+n)
    x = rng.uniform(.01, 1, (n, 11, 3)); good = rng.random((n, 11)) > .3
    good[:, -1] = False
    w = rng.dirichlet(np.ones(n+1))
    fn = (lambda counts: sfs_posterior_moments(x, good, w, include_counts=counts)) if prior == 'sfs' else (lambda counts: posterior_moments(x, good, prior, include_counts=counts))
    old = fn(False); result = fn(True)
    for key in old:
        np.testing.assert_array_equal(result[key], old[key])
    for j in range(11):
        if not good[:, j].any():
            assert np.isnan(result['count_posterior'][j]).all()
            continue
        expected = exhaustive(x[:, j], good[:, j], prior, w)
        np.testing.assert_allclose(result['count_posterior'][j], expected, atol=1e-13)
        m = 2*good[:, j].sum(); k = np.arange(2*n+1)
        assert result['gene_diversity'][j] == pytest.approx(np.sum(expected*2*k*(m-k)/(m*(m-1))), abs=1e-13)
    agg = PosteriorStats(); agg.add(result, good.any(axis=0))
    strata, spectrum = agg.tables('p')
    for m, row in strata.set_index('chromosome_count').iterrows():
        assert spectrum.loc[spectrum.chromosome_count == m, 'posterior_expected_sites'].sum() == pytest.approx(row.callable_bases)
    expected_theta = []
    for j in np.flatnonzero(good.any(axis=0)):
        m = 2*good[:, j].sum(); p = exhaustive(x[:, j], good[:, j], prior, w)
        expected_theta.append(p[1:m].sum()/sum(1/i for i in range(1, m)))
    assert agg.summary()['posterior_theta_w'] == pytest.approx(np.mean(expected_theta))


def test_counts_keep_population_dependence():
    x = np.ones((2, 1, 3)); good = np.ones((2, 1), bool)
    p = posterior_moments(x, good, 'population', include_counts=True)
    np.testing.assert_allclose(p['count_posterior'][0], .2)
    independent = posterior_moments(x, good, 'uniform', include_counts=True)
    assert not np.allclose(p['count_posterior'], independent['count_posterior'])


@pytest.mark.parametrize('prior', ['uniform', 'population', 'sfs'])
def test_deterministic_genotypes_missingness_and_negative_fis(prior):
    # Both heterozygous, opposite homozygotes, invariant, singleton, absent.
    gt = np.array([[1, 0, 0, 1, 0], [1, 2, 0, 0, 0]])
    good = np.ones_like(gt, bool); good[1, 3] = False; good[:, 4] = False
    x = np.eye(3)[gt]
    p = sfs_posterior_moments(x, good, [.8,.1,.1], include_counts=True) if prior == 'sfs' else posterior_moments(x, good, prior, include_counts=True)
    agg = PosteriorStats(); agg.add(p, good.any(axis=0)); s = agg.summary()
    assert s['posterior_theta_w'] == pytest.approx((2/(1+1/2+1/3)+1)/4)
    assert s['posterior_segregating_sites'] == pytest.approx(3)
    assert s['fis_callable_bases'] == 3
    assert s['posterior_fis'] == pytest.approx(1-1/(4/3))
    negative = PosteriorStats(); negative.add(p, np.array([True, False, False, False, False]))
    assert negative.summary()['posterior_fis'] == pytest.approx(-.5)
    invariant = PosteriorStats(); invariant.add(p, np.array([False,False,True,False,False]))
    assert invariant.summary()['posterior_theta_w'] == 0
    assert np.isnan(invariant.summary()['posterior_fis'])
    empty = PosteriorStats(); empty.add(p, np.zeros(5, bool))
    assert np.isnan(empty.summary()['posterior_theta_w'])


def test_crossfold_count_array_scatter():
    obj = CrossFittedSFS(np.arange(8), {'p': np.array([0, 1])}, fit_loci=8, folds=2)
    obj.models = {('p', 0): np.array([.8,.1,.1]), ('p', 1): np.array([.3,.3,.4])}
    x = np.random.default_rng(5).uniform(.01, 1, (2,8,3)); good = np.ones((2,8), bool)
    good[1,1] = False
    result = obj.moments('p', x, good, np.arange(8), include_counts=True)
    for j in range(8):
        expected = sfs_posterior_moments(x[:,j:j+1], good[:,j:j+1], obj.models['p',j%2], include_counts=True)
        np.testing.assert_allclose(result['count_posterior'][j], expected['count_posterior'][0])


def test_capabilities_all_and_explicit_rejection_before_fitting(tmp_path, monkeypatch):
    stats, caps = resolve_likelihood_stats('all', 'loci')
    assert set(stats) == {'pi','heterozygosity','theta_w','fis','sfs'}
    assert set(caps['omitted']) == {'dxy','fst','fit','tajima_d'}
    assert resolve_likelihood_stats('all','snps')[0] == ['heterozygosity']
    db = make(tmp_path)
    args = kwargs(db); args.pop('requested_stats')
    args.update(name='out',outdir=tmp_path/'result',force=False,genotype_mode='likelihoods',evidence_scope='loci')
    def fail(*a, **k):
        raise AssertionError('unsupported request reached fitting')
    monkeypatch.setattr(CrossFittedSFS, 'fit', fail)
    for stat in ['dxy','fst','fit','tajima_d']:
        with pytest.raises(IPyradError, match='No hard-call fallback'):
            run_popgen_method(**args, stats='pi,'+stat, evidence_prior='sfs')
        assert not args['outdir'].exists()
    run_popgen_method(**args, stats='all')
    manifest = (args['outdir']/'out.manifest.txt').read_text()
    assert 'statistic_capabilities' in manifest and 'tajima_d' in manifest
    assert (args['outdir']/'out.posterior_sfs.tsv').is_file()
    assert (args['outdir']/'out.chromosome_counts.tsv').is_file()
    assert not (args['outdir']/'out.sfs.tsv').exists()  # No fitted prior under uniform.
