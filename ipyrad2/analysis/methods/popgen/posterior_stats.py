"""Streaming all-base summaries from joint posterior allele counts."""
import numpy as np
import pandas as pd


class PosteriorStats:
    """Keep chromosome-count strata separate; never impute missing individuals."""

    def __init__(self):
        self.strata = {}
        self.spectra = {}
        self.fis_sites = 0
        self.ho = self.he = 0.

    def add(self, moments, eligible):
        n = moments['n_individuals']
        post = moments['count_posterior']
        for size in np.unique(n[eligible & (n > 0)]):
            use = eligible & (n == size)
            m = 2*int(size)
            p = post[use, :m+1]
            # Sum interior bins rather than subtracting near-one endpoints.
            seg = float(p[:, 1:m].sum())
            harmonic = float(np.sum(1./np.arange(1, m)))
            row = self.strata.setdefault(m, dict(callable_bases=0,
                posterior_segregating_sites=0., theta_w_sum=0.))
            row['callable_bases'] += int(use.sum())
            row['posterior_segregating_sites'] += seg
            row['theta_w_sum'] += seg/harmonic
            folded = p[:, :size+1].sum(axis=0)
            folded[:size] += p[:, size+1:][:, ::-1].sum(axis=0)
            self.spectra.setdefault(m, np.zeros(size+1))[:] += folded
        use = eligible & (n >= 2)
        self.fis_sites += int(use.sum())
        self.ho += float(np.sum(np.nansum(moments['GP'][:, use, 1], axis=0)/n[use]))
        self.he += float(np.sum(moments['gene_diversity'][use]))

    def summary(self):
        bases = sum(r['callable_bases'] for r in self.strata.values())
        theta = sum(r['theta_w_sum'] for r in self.strata.values())
        return dict(theta_w_callable_bases=bases,
            posterior_segregating_sites=sum(r['posterior_segregating_sites'] for r in self.strata.values()),
            posterior_theta_w=theta/bases if bases else np.nan,
            fis_callable_bases=self.fis_sites,
            fis_observed_heterozygosity=self.ho/self.fis_sites if self.fis_sites else np.nan,
            fis_expected_heterozygosity=self.he/self.fis_sites if self.fis_sites else np.nan,
            posterior_fis=1-self.ho/self.he if self.fis_sites and self.he > 0 else np.nan)

    def tables(self, population):
        strata = [dict(population=population, chromosome_count=m, **r,
                       posterior_theta_w=r['theta_w_sum']/r['callable_bases'])
                  for m, r in sorted(self.strata.items())]
        spectra = [dict(population=population, chromosome_count=m,
                        minor_allele_count=k, posterior_expected_sites=float(value),
                        callable_bases=self.strata[m]['callable_bases'])
                   for m, values in sorted(self.spectra.items()) for k, value in enumerate(values)]
        return pd.DataFrame(strata), pd.DataFrame(spectra)
