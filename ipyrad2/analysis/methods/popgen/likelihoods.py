"""Diploid biallelic posterior moments, conditional on a fixed SNP catalogue.

The population model integrates a Beta(1, 1) allele-frequency prior. Its
joint count variance includes covariance induced by the shared frequency;
plugging expected dosages into a hard-call estimator would omit this term.
No dropout correction or missing-cell imputation is performed here.
"""
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss


@lru_cache(maxsize=32)
def _quadrature(nindividuals):
    # The marginal likelihood is a polynomial of degree 2*n. Gauss-Legendre
    # with n+1 nodes integrates its moments exactly up to floating error.
    p, w = leggauss(max(2, nindividuals + 1))
    p = (p + 1) / 2
    return np.stack(((1-p)**2, 2*p*(1-p), p*p), axis=1), w / 2


def likelihoods_from_pl(pl, valid=None):
    """Return scaled native likelihoods and a validity mask, without mutation."""
    pl = np.asarray(pl, dtype=float)
    if pl.ndim != 3 or pl.shape[-1] != 3:
        raise ValueError('Expected PL with shape (individual, site, 3).')
    good = np.isfinite(pl).all(axis=-1) & (pl >= 0).all(axis=-1)
    if valid is not None:
        if np.shape(valid) != pl.shape[:2]:
            raise ValueError('PL and validity mask shapes differ.')
        good &= np.asarray(valid, bool)
    safe = np.where(good[..., None], pl, 0.)
    # Subtract first in float space to avoid signed-integer overflow.
    loglike = -(safe - safe.min(axis=-1, keepdims=True)) * (np.log(10.)/10.)
    return np.exp(loglike), good


def posterior_moments(likelihoods, valid, prior='uniform'):
    """Return marginal GP and joint sampled-chromosome diversity moments.

    Each row must represent a distinct individual. Under `population`,
    integrating all individual likelihoods once gives the same marginal GP
    as multiplying focal evidence by a leave-focal-out partner prior.
    Missing individuals contribute no evidence and receive no posterior.
    Gene diversity is E[2*K*(M-K)/(M*(M-1)) | reads], M=2*n covered samples.
    It is conditional on these sites/samples, not per-base nucleotide pi.
    """
    x = np.asarray(likelihoods, dtype=float)
    good = np.asarray(valid, bool)
    if x.ndim != 3 or x.shape[-1] != 3 or good.shape != x.shape[:2]:
        raise ValueError('Likelihoods must be individual x site x 3 with matching validity.')
    if prior not in ('uniform', 'population'):
        raise ValueError('Prior must be uniform or population.')
    if np.any(~np.isfinite(x[good])) or np.any(x[good] < 0) or np.any(x[good].sum(axis=-1) <= 0):
        raise ValueError('Usable likelihoods must be finite, nonnegative, and have positive mass.')
    x = np.where(good[..., None], x, 1.)
    x = x / x.max(axis=-1, keepdims=True)
    n = good.sum(axis=0)
    ns, nv = good.shape
    gp = np.zeros_like(x)
    ek = np.zeros(nv)
    ek2 = np.zeros(nv)
    if prior == 'uniform':
        gp = x / x.sum(axis=-1, keepdims=True)
        ds = (gp[..., 1] + 2*gp[..., 2]) * good
        second = (gp[..., 1] + 4*gp[..., 2]) * good
        ek = ds.sum(axis=0)
        ek2 = ek**2 + (second-ds**2).sum(axis=0)
    elif ns:
        h, w = _quadrature(ns)
        logweights = np.empty((len(w), nv))
        for j, hh in enumerate(h):
            marginal = (x * hh).sum(axis=-1)
            logweights[j] = np.log(w[j]) + np.where(good, np.log(marginal), 0).sum(axis=0)
        logweights -= logweights.max(axis=0, keepdims=True)
        weights = np.exp(logweights)
        weights /= weights.sum(axis=0, keepdims=True)
        for hh, weight in zip(h, weights):
            conditional = x * hh
            conditional /= conditional.sum(axis=-1, keepdims=True)
            gp += conditional * weight[None, :, None]
            ds = (conditional[..., 1] + 2*conditional[..., 2]) * good
            second = (conditional[..., 1] + 4*conditional[..., 2]) * good
            mean = ds.sum(axis=0)
            ek += weight * mean
            ek2 += weight * (mean**2 + (second-ds**2).sum(axis=0))
    gp[~good] = np.nan
    m = 2*n
    with np.errstate(divide='ignore', invalid='ignore'):
        frequency = ek/m
        frequency_variance = np.maximum(0., ek2-ek**2)/m**2
        diversity = 2*(m*ek-ek2)/(m*(m-1))
    return dict(GP=gp, dosage=gp[..., 1]+2*gp[..., 2], n_individuals=n,
                allele_frequency=frequency, allele_frequency_variance=frequency_variance,
                gene_diversity=np.where(n > 0, np.maximum(0., diversity), np.nan))
