#!/usr/bin/env python

"""Impute missing SNP genotypes from population allele frequencies."""

from __future__ import annotations

import numpy as np

from ...utils.exceptions import IPyradError


_MISSING_GENO = 255


class SNPImputer(object):
    """
    Impute missing diploid SNP genotypes from within-group allele frequencies.

    This class is used internally by analysis tools such as PCA and SNP export.
    The active implementation supports two behaviors only:

    - `"sample"`: sample missing genotypes from each group's inferred derived
      allele frequency.
    - `"zero-fill"` and null-like compatibility aliases (`None`, `False`,
      `"none"`): replace missing genotypes with `0`.
    """

    def __init__(
        self,
        data,
        names,
        imap=None,
        impute_method="sample",
        quiet=False,
    ):

        self.quiet = quiet
        self.snps = self._validate_data(data)
        self.names = self._validate_names(names, self.snps.shape[0])
        self.name_to_index = {name: idx for idx, name in enumerate(self.names)}
        self.imap = self._normalize_imap(imap)
        self.impute_method = self._normalize_impute_method(impute_method)
        self._mvals = int(np.sum(self.snps == _MISSING_GENO))

    def _print(self, msg):
        if not self.quiet:
            print(msg)

    @staticmethod
    def _validate_data(data) -> np.ndarray:
        """Return a writable 2D genotype array copy."""
        snps = np.asarray(data)
        if snps.ndim != 2:
            raise IPyradError("SNPImputer data must be a 2D genotype matrix.")
        if not np.issubdtype(snps.dtype, np.integer):
            raise IPyradError("SNPImputer data must have an integer genotype dtype.")
        return snps.copy()

    @staticmethod
    def _validate_names(names, nsamples: int) -> list[str]:
        """Validate and normalize the ordered sample names."""
        normalized = [str(name) for name in names]
        if len(normalized) != nsamples:
            raise IPyradError(
                "SNPImputer names length must match the number of genotype rows."
            )
        return normalized

    def _normalize_imap(self, imap) -> dict[str, list[str]]:
        """Normalize IMAP input and validate group/sample assignments."""
        if imap is None:
            return {"1": list(self.names)}
        if not isinstance(imap, dict):
            raise IPyradError("SNPImputer imap must be a dict or None.")

        normalized = {}
        seen = set()
        for group, samples in imap.items():
            group_name = str(group)
            sample_names = [str(name) for name in samples]
            if not sample_names:
                raise IPyradError(f"SNPImputer imap group {group_name!r} is empty.")
            if len(set(sample_names)) != len(sample_names):
                raise IPyradError(
                    f"SNPImputer imap group {group_name!r} contains duplicate sample names."
                )
            missing = sorted(set(sample_names).difference(self.name_to_index))
            if missing:
                raise IPyradError(
                    "SNPImputer imap contains sample names not present in the genotype matrix: "
                    + ", ".join(missing)
                )
            duplicates = sorted(name for name in sample_names if name in seen)
            if duplicates:
                dupes = ", ".join(sorted(set(duplicates)))
                raise IPyradError(
                    "SNPImputer imap assigns a sample to multiple groups: " + dupes
                )
            seen.update(sample_names)
            normalized[group_name] = sample_names
        return normalized

    @staticmethod
    def _normalize_impute_method(impute_method):
        """Normalize supported imputation method names."""
        if impute_method == "sample":
            return "sample"
        if impute_method is None or impute_method is False:
            return "zero-fill"
        if isinstance(impute_method, str) and impute_method.lower() in {"none", "zero", "zero-fill"}:
            return "zero-fill"
        raise IPyradError(
            "Unsupported SNPImputer impute_method. Use 'sample' or 'zero-fill'."
        )

    def run(self):
        """Return an imputed copy of the genotype matrix."""
        if not self._mvals:
            self._print("Imputation: no missing genotypes; returning input unchanged.")
            return self.snps

        if self.impute_method == "sample":
            self.snps = self._impute_sample()
        else:
            self.snps = self._fill_missing_with_zero()
        return self.snps

    def _fill_missing_with_zero(self) -> np.ndarray:
        """Replace missing genotypes with homozygous reference calls."""
        newdata = self.snps.copy()
        missing = newdata == _MISSING_GENO
        newdata[missing] = 0
        imputed = newdata[missing]
        self._print(self._format_imputation_stats("zero-fill", imputed))
        return newdata

    def _impute_sample(self, imap=None):
        """
        Sample genotypes from each population's derived-allele frequency.

        Sites that are entirely missing within a population get derived-allele
        frequency zero, which deterministically imputes genotype `0`.
        """
        groups = self.imap if imap is None else self._normalize_imap(imap)
        newdata = self.snps.copy()

        for samps in groups.values():
            sidxs = sorted(self.name_to_index[name] for name in samps)
            data = newdata[sidxs, :].copy()

            # Convert genotype counts to allele counts while ignoring missing.
            nalleles = np.sum(data != _MISSING_GENO, axis=0) * 2
            tmp = data.copy()
            tmp[tmp == _MISSING_GENO] = 0
            with np.errstate(divide="ignore", invalid="ignore"):
                fderived = tmp.sum(axis=0) / nalleles
            fderived[np.isnan(fderived)] = 0

            sampled = np.random.binomial(n=2, p=fderived, size=data.shape)
            missing = data == _MISSING_GENO
            data[missing] = sampled[missing]
            newdata[sidxs, :] = data

        imputed = newdata[self.snps == _MISSING_GENO]
        self._print(self._format_imputation_stats("sample", imputed))
        return newdata

    @staticmethod
    def _format_imputation_stats(label: str, imputed: np.ndarray) -> str:
        """Render a compact summary of imputed genotype frequencies."""
        if imputed.size == 0:
            return f"Imputation: '{label}'; no genotypes were imputed."
        return (
            "Imputation: '{}'; (0, 1, 2) = {:.1f}%, {:.1f}%, {:.1f}%".format(
                label,
                100 * np.sum(imputed == 0) / imputed.size,
                100 * np.sum(imputed == 1) / imputed.size,
                100 * np.sum(imputed == 2) / imputed.size,
            )
        )

__all__ = ["SNPImputer"]
