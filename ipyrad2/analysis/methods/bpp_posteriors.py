#!/usr/bin/env python

"""Deprecated posterior-transformation helpers for legacy BPP workflows."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from ...utils.exceptions import IPyradError


warnings.warn(
    "ipyrad2.analysis.methods.bpp_posteriors is deprecated and retained only for "
    "legacy notebook workflows.",
    DeprecationWarning,
    stacklevel=2,
)


try:
    import scipy.stats as ss
except ImportError as exc:  # pragma: no cover - optional dependency
    raise IPyradError(
        "You are missing required packages to use legacy BPP posterior helpers.\n"
        "First run: conda install scipy -c conda-forge"
    ) from exc


def _coerce_positive_number(value, label: str) -> float:
    """Parse one positive float."""
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise IPyradError(f"{label} must be numeric.") from exc
    if parsed <= 0:
        raise IPyradError(f"{label} must be > 0.")
    return parsed


def draw_gamma_from_range(min_value: float, max_value: float) -> tuple[float, float]:
    """Approximate a bounded uncertainty range with a gamma(a, b-rate) prior."""
    min_value = _coerce_positive_number(min_value, "minimum value")
    max_value = _coerce_positive_number(max_value, "maximum value")
    if max_value <= min_value:
        raise IPyradError("maximum value must be greater than minimum value.")
    mean = (max_value + min_value) / 2.0
    var = ((max_value - min_value) ** 2) / 16.0
    return mean ** 2 / var, mean / var


class Transformer:
    """Transform posterior theta/tau samples into Ne and divergence-time units."""

    def __init__(self, df, gentime_min, gentime_max, mutrate_min, mutrate_max, seed=123):
        self.df = df
        if self.df is None or self.df.empty:
            raise IPyradError("Cannot transform an empty BPP posterior table.")

        self.seed = int(seed)
        self._rng = np.random.default_rng(self.seed)
        self.gentime_a, self.gentime_b = draw_gamma_from_range(gentime_min, gentime_max)
        self.mutrate_a, self.mutrate_b = draw_gamma_from_range(mutrate_min, mutrate_max)
        self._sample_gentime_rvs()
        self._sample_mutrate_rvs()

    def _sample_gentime_rvs(self):
        """Sample generation times from the user-provided uncertainty range."""
        self.gentime_rvs = ss.gamma.rvs(
            self.gentime_a,
            scale=1 / self.gentime_b,
            random_state=self._rng,
            size=self.df.shape[0],
        )

    def _sample_mutrate_rvs(self):
        """Sample mutation rates from the user-provided uncertainty range."""
        self.mutrate_rvs = ss.gamma.rvs(
            self.mutrate_a,
            scale=1 / self.mutrate_b,
            random_state=self._rng,
            size=self.df.shape[0],
        )

    def _get_parameter_values(self, colname):
        """Return one posterior parameter column plus the aligned uncertainty draws."""
        if colname not in self.df.columns:
            raise IPyradError(f"posterior column not found: {colname}")
        series = pd.to_numeric(self.df[colname], errors="coerce")
        valid = series.notna().to_numpy()
        if not np.any(valid):
            raise IPyradError(f"posterior column has no numeric values: {colname}")
        values = series.to_numpy(dtype=float, na_value=np.nan)[valid]
        return values, self.gentime_rvs[valid], self.mutrate_rvs[valid]

    def transform(self, colname):
        """Transform one posterior column to divergence-time or Ne units."""
        values, gentime, mutrate = self._get_parameter_values(colname)
        if "tau" in colname:
            return (values * gentime) / mutrate
        if "theta" in colname:
            return values / (mutrate * 4)
        raise IPyradError(f"Unsupported BPP posterior parameter: {colname}")
