"""Independent exhaustive checks of SFS likelihoods and posterior moments."""
import itertools
from math import comb

import numpy as np
import pytest

from ipyrad2.analysis.methods.popgen.sfs_prior import (
    folded_saf, fit_folded_sfs, sfs_posterior_moments, unfold_sfs)
from ipyrad2.analysis.methods.popgen.sfs_fitting import CrossFittedSFS


def enumerate_states(x, good, weights):
    n = len(x)
    p = unfold_sfs(weights)
    gp = np.zeros_like(x)
    likelihood = np.zeros(2*n+1)
    terms = []
    for state in itertools.product(range(3), repeat=n):
        k = sum(state)
        multiplicity = 2**state.count(1)/comb(2*n, k)
        like = np.prod([x[i, g] if good[i] else 1. for i, g in enumerate(state)])
        likelihood[k] += multiplicity*like
        terms.append((state, p[k]*multiplicity*like))
    normalizer = sum(t[1] for t in terms)
    e = e2 = 0.
    for state, weight in terms:
        weight /= normalizer
        count = sum(g for g, use in zip(state, good) if use)
        e += weight*count
        e2 += weight*count*count
        for i, g in enumerate(state):
            gp[i, g] += weight
    gp[~good] = np.nan
    m = 2*good.sum()
    return likelihood, gp, e/m if m else np.nan, 2*(m*e-e2)/(m*(m-1)) if m else np.nan


@pytest.mark.parametrize('n', [1, 2, 3, 4])
def test_exhaustive_saf_and_moments(n):
    rng = np.random.default_rng(n)
    x = rng.uniform(.001, 1, (n, 13, 3))
    good = rng.random((n, 13)) > .25
    w = rng.dirichlet(np.ones(n+1))
    saf = folded_saf(x, good)
    result = sfs_posterior_moments(x, good, w)
    for j in range(13):
        lk, gp, freq, diversity = enumerate_states(x[:, j], good[:, j], w)
        folded = (lk[:n+1]+lk[n:][::-1])/2
        np.testing.assert_allclose(saf[j], folded/folded.max(), atol=1e-13)
        np.testing.assert_allclose(result['GP'][:, j], gp, atol=1e-13)
        np.testing.assert_allclose(result['allele_frequency'][j], freq, atol=1e-13)
        np.testing.assert_allclose(result['gene_diversity'][j], diversity, atol=1e-13)
    swapped = sfs_posterior_moments(x[..., ::-1], good, w)
    np.testing.assert_allclose(swapped['GP'], result['GP'][..., ::-1], atol=1e-13)
    np.testing.assert_allclose(swapped['gene_diversity'], result['gene_diversity'], atol=1e-13)


@pytest.mark.parametrize('n', [2, 3])
def test_small_group_known_counts(n):
    # Certain genotype configurations identify the folded count exactly.
    counts = np.arange(1, n+2)*20
    f = np.repeat(np.eye(n+1), counts, axis=0)
    w, d = fit_folded_sfs(f)
    np.testing.assert_allclose(w, (counts+.5)/(counts.sum()+.5*(n+1)), atol=1e-10)
    assert d['converged'] and d['informative_sites'] == counts.sum()
    with pytest.raises(ValueError, match='No informative'):
        fit_folded_sfs(np.ones((20, n+1)))
    noisy = np.random.default_rng(n).uniform(.01, 1, (100, n+1))
    with pytest.raises(ValueError, match='did not converge'):
        fit_folded_sfs(noisy, max_iter=1)


def test_crossfit_loci_excluded_and_budget():
    obj = CrossFittedSFS(np.arange(20), {'a': np.array([0, 1])}, fit_loci=10, folds=2)
    rng = np.random.default_rng(4)
    x = rng.uniform(.01, 1, (2, 20, 3))
    good = np.ones((2, 20), bool)
    obj.add(x, good, np.arange(20))
    obj.fit(2)
    for fold in range(2):
        chosen = obj.fit_loci[obj.fold_ids(obj.fit_loci) != fold]
        expected, _ = fit_folded_sfs(folded_saf(x[:, chosen], good[:, chosen]))
        np.testing.assert_allclose(obj.models['a', fold], expected)
    result = obj.moments('a', x, good, np.arange(20))
    assert result['GP'].shape == x.shape
    np.testing.assert_allclose(result['GP'].sum(axis=-1), 1.)


def test_missing_flat_and_invalid():
    x = np.ones((3, 4, 3))
    good = np.zeros((3, 4), bool)
    result = sfs_posterior_moments(x, good, [.8, .1, .05, .05])
    assert np.isnan(result['GP']).all() and np.isnan(result['gene_diversity']).all()
    np.testing.assert_allclose(folded_saf(x, good), 1.)
    with pytest.raises(ValueError):
        sfs_posterior_moments(x, good, [1., 0.])


@pytest.mark.parametrize('scope', ['snps', 'loci'])
def test_native_backend_crossfit_and_read_only(tmp_path, scope):
    import hashlib
    import h5py
    from test_popgen_evidence import kwargs
    from test_hdf5_evidence import fixture
    from test_rad_site_evidence import setup
    from ipyrad2.assembler.evidence_hdf5 import write_snp_evidence_hdf5
    from ipyrad2.assembler.rad_evidence_hdf5 import write_rad_evidence_hdf5
    from ipyrad2.analysis.methods.popgen.runner import run_popgen_method
    if scope == 'snps':
        db, raw, final = fixture(tmp_path)
        write_snp_evidence_hdf5(database=db, raw_vcf=raw, final_vcf=final)
    else:
        db, ref, bed, raw, mask = setup(tmp_path)
        bed.write_text('chr1\t0\t1\nchr1\t1\t8\n')
        write_rad_evidence_hdf5(database=db, reference=ref, bed=bed,
            native_files=[raw], native_regions=[bed], hard_masks=({'RAD': mask},))
    before = hashlib.sha256(db.read_bytes()).hexdigest()
    args = kwargs(db)
    args.pop('requested_stats')
    args.update(name='out', outdir=tmp_path/'out', stats='all', force=False,
        genotype_mode='likelihoods', evidence_scope=scope, evidence_prior='sfs',
        sfs_fit_loci=2, sfs_folds=2)
    run_popgen_method(**args)
    assert (tmp_path/'out/out.sfs.tsv').is_file()
    assert 'whole-locus cross-fitting' in (tmp_path/'out/out.manifest.txt').read_text()
    assert hashlib.sha256(db.read_bytes()).hexdigest() == before
    # Canonical hard calls do not fit or select the SFS likelihood domain.
    with h5py.File(db, 'a') as h:
        h['genos'][:] = 255
    args['outdir'] = tmp_path/'missing_calls'
    run_popgen_method(**args)
    assert (tmp_path/'out/out.sfs.tsv').read_bytes() == (tmp_path/'missing_calls/out.sfs.tsv').read_bytes()
    assert (tmp_path/'out/out.population_stats.tsv').read_bytes() == (tmp_path/'missing_calls/out.population_stats.tsv').read_bytes()
