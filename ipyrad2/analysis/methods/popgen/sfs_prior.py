"""Folded finite-sample SFS estimation and posterior chromosome moments.

SAF(k) is the likelihood averaged over exchangeable assignments of k alternate
chromosomes, not a likelihood for an estimated allele frequency. Missing
individuals have constant likelihood. Reference is never assumed ancestral.
"""
from functools import lru_cache
from math import comb

import numpy as np


def _inputs(likelihoods, valid):
    x = np.asarray(likelihoods, float)
    good = np.asarray(valid, bool)
    if x.ndim != 3 or x.shape[-1] != 3 or good.shape != x.shape[:2] or not len(x):
        raise ValueError('Expected nonempty individual x site x 3 likelihoods and validity mask.')
    if len(x) > 128:
        raise ValueError('SFS mode currently supports at most 128 individuals per group.')
    if np.any(~np.isfinite(x[good])) or np.any(x[good] < 0) or np.any(x[good].max(axis=-1) <= 0):
        raise ValueError('Usable likelihoods must be finite, nonnegative and nonzero.')
    x = np.where(good[..., None], x, 1.)
    return x / x.max(axis=-1, keepdims=True), good


def _forward(x):
    """Scaled coefficients of product(L0 + 2 L1 z + L2 z**2)."""
    a = x * np.array([1., 2., 1.])
    forward = [np.ones((x.shape[1], 1))]
    for row in a:
        old = forward[-1]
        new = np.zeros((len(row), old.shape[1] + 2))
        for g in range(3):
            new[:, g:g+old.shape[1]] += old * row[:, g, None]
        new /= new.max(axis=1, keepdims=True)
        forward.append(new)
    return a, forward


def folded_saf(likelihoods, valid):
    """Return site x (N+1) folded SAF, scaled independently at each site."""
    x, _ = _inputs(likelihoods, valid)
    _, forward = _forward(x)
    n = len(x)
    saf = forward[-1] / np.array([comb(2*n, k) for k in range(2*n+1)], float)
    folded = (saf[:, :n+1] + saf[:, n:][:, ::-1]) / 2
    folded /= folded.max(axis=1, keepdims=True)
    return folded


def fit_folded_sfs(saf, *, max_iter=10000, tolerance=1e-7, pseudocount=0.5):
    """Regularized EM; Dirichlet alpha = 1 + pseudocount per folded bin.

    Flat SAF rows contain no information and are omitted, not counted as
    monomorphic observations. No convergence or uninformative-fit fallback.
    """
    f = np.asarray(saf, float)
    if f.ndim != 2 or f.shape[1] < 2 or not len(f):
        raise ValueError('SFS fitting requires a nonempty site x folded-count matrix.')
    if np.any(~np.isfinite(f)) or np.any(f < 0) or np.any(f.max(axis=1) <= 0):
        raise ValueError('Invalid SAF values.')
    if max_iter < 1 or tolerance <= 0 or pseudocount <= 0:
        raise ValueError('SFS iterations, tolerance and pseudocount must be positive.')
    f = f / f.max(axis=1, keepdims=True)
    informative = np.ptp(f, axis=1) > 1e-12
    total = len(f)
    f = f[informative]
    if not len(f):
        raise ValueError('No informative sites for SFS fitting.')
    # Exact deduplication is useful for recurrent low-depth likelihood patterns.
    f, counts = np.unique(f, axis=0, return_counts=True)
    w = np.full(f.shape[1], 1/f.shape[1])
    objective = -np.inf
    for iteration in range(1, max_iter+1):
        marginal = f @ w
        updated = w * (f.T @ (counts / marginal)) + pseudocount
        updated /= updated.sum()
        new_objective = float(counts @ np.log(f @ updated) + pseudocount*np.log(updated).sum())
        delta = float(np.max(np.abs(updated-w)))
        if new_objective < objective - 1e-8*max(1., abs(objective)):
            raise ValueError('SFS EM objective decreased unexpectedly.')
        w = updated
        if delta <= tolerance:
            return w, dict(converged=True, iterations=iteration, weight_delta=delta,
                penalized_log_likelihood=new_objective, sites=total,
                informative_sites=int(informative.sum()), unique_saf_rows=len(f),
                tolerance=tolerance, pseudocount=pseudocount)
        objective = new_objective
    raise ValueError(f'SFS fit did not converge in {max_iter} iterations (weight delta {delta:.3g}).')


def unfold_sfs(weights):
    w = np.asarray(weights, float)
    if w.ndim != 1 or len(w) < 2 or np.any(~np.isfinite(w)) or np.any(w < 0) or w.sum() <= 0:
        raise ValueError('Invalid folded SFS weights.')
    n = len(w)-1
    p = np.zeros(2*n+1)
    p[:n+1] += w/2
    p[n:] += w[::-1]/2
    return p/p.sum()


@lru_cache(maxsize=128)
def _projection(n, m):
    return np.array([[comb(k, j)*comb(2*n-k, 2*m-j)/comb(2*n, 2*m)
        if j <= k and 0 <= 2*m-j <= 2*n-k else 0.
        for j in range(2*m+1)] for k in range(2*n+1)])


def sfs_posterior_moments(likelihoods, valid, weights):
    """Joint posterior with hypergeometric projection onto usable individuals.

    Missing cells receive NaN, not imputed posteriors. Count variance retains
    dependence between individuals induced by the finite-count prior.
    """
    x, good = _inputs(likelihoods, valid)
    ns, nv = good.shape
    if len(weights) != ns+1:
        raise ValueError('SFS sample count differs from likelihood group size.')
    prior = unfold_sfs(weights)
    gp = np.full_like(x, np.nan)
    n = good.sum(axis=0)
    ek = np.zeros(nv)
    ek2 = np.zeros(nv)
    for m in np.unique(n[n > 0]):
        sites = np.flatnonzero(n == m)
        order = np.argsort(~good[:, sites], axis=0, kind='stable')[:m]
        xx = x[order, sites[None, :]]
        a, forward = _forward(xx)
        p = prior @ _projection(ns, int(m))
        end = p / np.array([comb(2*int(m), k) for k in range(2*int(m)+1)], float)
        post = forward[-1] * end
        mass = post.sum(axis=1, keepdims=True)
        if np.any(mass <= 0):
            raise ValueError('SFS prior and likelihoods have no common support.')
        post /= mass
        k = np.arange(2*m+1)
        ek[sites] = post @ k
        ek2[sites] = post @ (k*k)
        backward = np.broadcast_to(end, post.shape).copy()
        for i in range(int(m)-1, -1, -1):
            width = forward[i].shape[1]
            marginal = np.stack([(forward[i]*backward[:, g:g+width]).sum(axis=1)*a[i, :, g]
                                 for g in range(3)], axis=-1)
            marginal /= marginal.sum(axis=1, keepdims=True)
            gp[order[i], sites] = marginal
            backward = sum(a[i, :, g, None]*backward[:, g:g+width] for g in range(3))
            backward /= backward.max(axis=1, keepdims=True)
    chromosomes = 2*n
    with np.errstate(divide='ignore', invalid='ignore'):
        frequency = ek/chromosomes
        variance = np.maximum(0., ek2-ek*ek)/chromosomes**2
        diversity = 2*(chromosomes*ek-ek2)/(chromosomes*(chromosomes-1))
    return dict(GP=gp, dosage=gp[..., 1]+2*gp[..., 2], n_individuals=n,
        allele_frequency=frequency, allele_frequency_variance=variance,
        gene_diversity=np.where(n > 0, np.maximum(0., diversity), np.nan))
