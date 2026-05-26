#!/usr/bin/env python

"""Load, filter, subsample, and export SNP data from HDF5 databases."""

from __future__ import annotations

import gzip
import itertools
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple

import h5py
import numpy as np
import pandas as pd
from loguru import logger

from ...utils.exceptions import IPyradError
from ...utils.parallel import run_with_pool
from ...utils.pops import expand_imap_patterns, parse_imap, parse_minmap, parse_pops_file

# Value of missing data in the snps matrix
_MISSING_GENO = 255
_MISSING_SNP = 78

# how many cols of SNPs to load in at once from snps, genos, snpsmap
CHUNKSIZE = 10_000
REFERENCE_SAMPLE_NAME = "assembly_reference_sequence"
SNPSMAP_COLUMNS = ["loc", "loc_idx", "loc_pos", "scaff", "pos"]
STATS_HEADER = [
    "samples",
    "pre_filter_snps",
    "pre_filter_percent_missing",
    "masked_genotypes_by_min_depth",
    "filter_by_indels_present",
    "filter_by_non_biallelic",
    "filter_by_mincov",
    "filter_by_minmap",
    "filter_by_min_site_qual",
    "filter_by_invariant_after_subsampling",
    "filter_by_minor_allele_frequency",
    "post_filter_snps",
    "post_filter_snp_containing_linkage_blocks",
    "post_filter_percent_missing",
]

_IUPAC_TO_BASES = {
    "R": ("A", "G"),
    "K": ("G", "T"),
    "S": ("G", "C"),
    "Y": ("C", "T"),
    "W": ("A", "T"),
    "M": ("A", "C"),
}
_BASES_TO_IUPAC = {
    tuple(sorted(("A", "C"))): "M",
    tuple(sorted(("A", "G"))): "R",
    tuple(sorted(("A", "T"))): "W",
    tuple(sorted(("C", "G"))): "S",
    tuple(sorted(("C", "T"))): "Y",
    tuple(sorted(("G", "T"))): "K",
}
NEXHEADER = """#nexus
begin data;
  dimensions ntax={} nchar={};
  format datatype=dna missing=N gap=- interleave=yes;
  matrix
"""


def _decode_h5_names(values) -> list[str]:
    """Decode an HDF5 string array into a plain Python string list."""
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


class SNPExportView(NamedTuple):
    """Aligned linked or unlinked SNP data ready to write."""

    genos: np.ndarray
    snps: np.ndarray
    snpsmap: np.ndarray
    reference: np.ndarray | None


class PreparedSNPExport(NamedTuple):
    """Prepared SNP export payload shared across all output formats."""

    view: SNPExportView
    genos: np.ndarray
    snps: np.ndarray
    imputation_summary: object | None
    sample_data_summary: pd.DataFrame


class SNPsExtracter:
    """Extract and subsample SNP data from HDF5 after filtering."""

    def __init__(
        self,
        data: Path,
        min_sample_coverage: float,
        max_sample_missing: float,
        min_minor_allele_frequency: float,
        imap: Path | str | Dict | None,
        minmap: Path | Dict | None,
        min_genotype_depth: int = 0,
        min_site_qual: float = 0.0,
        exclude: Path | str | List | None = None,
        include_reference: bool = False,
        cores: int = 1,
    ):
        imap, minmap = self._parse_imap_minmap_inputs(imap, minmap)

        # store params
        self.data = Path(data)
        self.user_imap = bool(imap)
        self.imap = imap if imap else {}
        self.minmap = minmap if minmap else {}
        self.mincov = min_sample_coverage
        self.max_sample_missing = min(1.0, max(0.0, max_sample_missing))
        self.maf = min_minor_allele_frequency
        self.min_genotype_depth = max(0, int(min_genotype_depth))
        self.min_site_qual = max(0.0, float(min_site_qual))
        self.exclude = exclude if exclude else []
        self.include_reference = include_reference
        self.cores = cores

        # attributes to be filled
        self.nsnps: int | None = None
        self.dbnames: List[str] | None = None
        self.dbname_to_index: Dict[str, int] = {}
        self.dbname_to_genos_index: Dict[str, int | None] = {}
        self.snames: List[str] | None = None
        self.initial_snames: List[str] = []
        self.dropped_samples_by_missing: List[str] = []
        self.sidxs: List[int] | None = None
        self.mask: np.ndarray | None = None
        self.snps: np.ndarray | None = None
        self.genos: np.ndarray | None = None
        self.snpsmap: np.ndarray | None = None
        self.reference: np.ndarray | None = None
        self.has_reference_dataset = False
        self.has_sample_dp_dataset = False
        self.has_site_qual_dataset = False
        self.synthetic_reference_row = False
        self.stats: pd.Series | None = None
        self.sample_missing: pd.Series | None = None

        self._set_nsnps_and_check_h5_file_format()
        self._check_requested_h5_capabilities()
        self._get_exclude(self.exclude)
        self._get_snames_and_sidxs_subset(self.imap)
        self._get_imap_minmap(self.imap, self.minmap)
        self.initial_snames = list(self.snames)
        self.names = self.snames

    def _get_exclude(self, exclude):
        if isinstance(exclude, str):
            exclude = Path(exclude)
        if isinstance(exclude, Path):
            if not exclude.exists():
                raise IPyradError(f"sample exclude file does not exist: {exclude}")
            with open(exclude, encoding="utf-8") as infile:
                self.exclude = [x.strip() for x in infile.readlines()]
        elif isinstance(exclude, list):
            self.exclude = exclude
        self.exclude = set(self.exclude)
        logger.debug(f"Excluding samples: {self.exclude}")

    def _parse_imap_minmap_inputs(self, imap, minmap):
        """Normalize imap/minmap inputs from dicts or files."""
        if isinstance(imap, str):
            imap = Path(imap)
        if isinstance(minmap, str):
            minmap = Path(minmap)

        if imap is not None and not isinstance(imap, (Path, Dict)):
            raise IPyradError("imap must be one of Path, str, or Dict")
        if minmap is not None and not isinstance(minmap, (Path, Dict)):
            raise IPyradError("minmap must be one of Path, str, or Dict")

        if isinstance(minmap, Path):
            minmap = parse_minmap(minmap)

        if isinstance(imap, Path):
            parsed_minmap = None
            try:
                imap, parsed_minmap = parse_pops_file(imap)
            except IPyradError:
                logger.info(
                    "imap file doesn't include minmap info, parsing standard imap file format."
                )
                imap = parse_imap(imap)
            if minmap is None:
                minmap = parsed_minmap

        return imap, minmap

    def _get_snames_and_sidxs_subset(self, imap) -> None:
        if imap:
            imap, _unmatched = expand_imap_patterns(
                imap,
                self.snames,
                mapping_name="IMAP",
                available_name=f"data file {self.data}",
            )
            self.imap = imap
            imapset = set(itertools.chain(*imap.values()))
            badnames = imapset.difference(self.snames)
            if badnames:
                badlist = ", ".join(sorted(badnames))
                raise IPyradError(f"Samples {badlist} are not in data file: {self.data}")
            if (
                self.include_reference
                and REFERENCE_SAMPLE_NAME not in self.exclude
                and REFERENCE_SAMPLE_NAME not in imapset
            ):
                raise IPyradError(
                    "assembly_reference_sequence was requested with -R, "
                    "but it must also be assigned to an IMAP group."
                )
        else:
            imapset = set()

        dbnames = set(self.snames)
        if (
            REFERENCE_SAMPLE_NAME in dbnames
            and REFERENCE_SAMPLE_NAME not in self.exclude
            and not self.include_reference
            and REFERENCE_SAMPLE_NAME not in imapset
        ):
            self.exclude.add(REFERENCE_SAMPLE_NAME)

        if imap:
            self.exclude.update(set(self.snames).difference(imapset))
            logger.debug(
                "dropping samples that are either not in the imap dict, "
                f"or are in the exclude list: {self.exclude}"
            )

        self.sidxs = [i for (i, name) in enumerate(self.snames) if name not in self.exclude]
        self.snames = [name for (i, name) in enumerate(self.snames) if i in self.sidxs]
        if (
            self.synthetic_reference_row
            and REFERENCE_SAMPLE_NAME in self.snames
            and not self.has_reference_dataset
        ):
            raise IPyradError(
                "assembly_reference_sequence was selected for SNP analysis, but this HDF5 "
                "cannot synthesize its genotype row because the `reference` dataset is missing."
            )

    def _get_imap_minmap(self, imap, minmap):
        """Set _imap and _minmap for seqarr filtering."""
        if not imap:
            self.imap = {"all": list(self.snames)}
            self.minmap = {"all": self.mincov}
        else:
            if not minmap:
                logger.info(
                    "No minmap specified. The global `-m` filter still applies; "
                    "defaulting per-population minimums to 0 so `-g` has no effect "
                    "unless it is provided."
                )
                minmap = {group: 0 for group in imap}

            if set(minmap) != set(imap):
                raise IPyradError("imap and minmap keys must match.")
            included_names = set(self.snames)
            self.imap = {
                key: [name for name in names if name in included_names]
                for key, names in imap.items()
            }
            self.minmap = minmap.copy()
            logger.debug(f"loaded imap = {self.imap}")
            logger.debug(f"sample coverage minmap = {self.minmap}")

    def _set_nsnps_and_check_h5_file_format(self) -> None:
        """Check input data is proper format, get nsnps, and set snames."""
        if self.data.suffix in [".vcf", ".vcf.gz"]:
            raise TypeError("input should be hdf5, see the vcf_to_hdf5 tool.")
        with h5py.File(self.data, "r") as io5:
            self.nsnps = int(io5.attrs["nsnps"])
            if "genos" not in io5:
                raise IPyradError(f"SNP HDF5 is missing the required `genos` dataset: {self.data}")
            self.has_reference_dataset = "reference" in io5
            self.has_sample_dp_dataset = "sample_dp" in io5
            self.has_site_qual_dataset = "site_qual" in io5
            self.dbnames, self.dbname_to_genos_index = self._resolve_snp_sample_names(io5)
            self.snames = list(self.dbnames)
            self.dbname_to_index = {name: idx for idx, name in enumerate(self.dbnames)}

    def _check_requested_h5_capabilities(self) -> None:
        """Reject legacy HDF5 files when new optional SNP filters are requested."""
        if self.min_genotype_depth > 0 and not self.has_sample_dp_dataset:
            raise IPyradError(
                "The `--min-genotype-depth` filter requires the HDF5 `sample_dp` dataset. "
                "Rebuild the SNP HDF5 with a current assemble or `ipyrad2 vcf2hdf5` run."
            )
        if self.min_site_qual > 0 and not self.has_site_qual_dataset:
            raise IPyradError(
                "The `--min-site-qual` filter requires the HDF5 `site_qual` dataset. "
                "Rebuild the SNP HDF5 with a current assemble or `ipyrad2 vcf2hdf5` run."
            )

    def _resolve_snp_sample_names(
        self, io5: h5py.File
    ) -> tuple[list[str], dict[str, int | None]]:
        """Return sample names aligned to rows of the SNP genotype matrix."""
        genos = io5["genos"]
        genos_rows = int(genos.shape[0])
        root_names = _decode_h5_names(io5.attrs["names"])
        dataset_names_raw = genos.attrs.get("names")

        if dataset_names_raw is not None:
            dataset_names = _decode_h5_names(dataset_names_raw)
            if len(dataset_names) != genos_rows:
                raise IPyradError(
                    "SNP HDF5 sample metadata mismatch: "
                    f"`genos`.attrs['names'] has {len(dataset_names)} entries but "
                    f"`genos` has {genos_rows} sample rows in {self.data}."
                )
            return dataset_names, {name: idx for idx, name in enumerate(dataset_names)}

        if len(root_names) == genos_rows:
            return root_names, {name: idx for idx, name in enumerate(root_names)}

        filtered_root_names = [
            name for name in root_names if name != REFERENCE_SAMPLE_NAME
        ]
        if (
            len(root_names) == genos_rows + 1
            and len(filtered_root_names) == genos_rows
            and REFERENCE_SAMPLE_NAME in root_names
        ):
            self.synthetic_reference_row = True
            logger.debug(
                "top-level sample metadata includes assembly_reference_sequence, "
                "but the SNP genotypes do not; synthesizing a homozygous-reference row "
                "when the reference sample is selected"
            )
            mapping: dict[str, int | None] = {}
            actual_index = 0
            for name in root_names:
                if name == REFERENCE_SAMPLE_NAME:
                    mapping[name] = None
                else:
                    mapping[name] = actual_index
                    actual_index += 1
            return root_names, mapping

        raise IPyradError(
            "SNP HDF5 sample metadata mismatch: "
            f"attrs['names'] has {len(root_names)} entries but `genos` has "
            f"{genos_rows} sample rows in {self.data}."
        )

    def _calculate_missing_frequencies(self, genos: np.ndarray) -> pd.Series:
        """Return per-sample missing-data frequencies on the current matrix."""
        if genos.size == 0 or genos.shape[1] == 0:
            values = np.zeros(genos.shape[0], dtype=float)
        else:
            values = np.mean(genos == _MISSING_GENO, axis=1)
        return pd.Series(values, index=self.snames, dtype=float)

    def _sync_imap_after_sample_drop(self) -> None:
        """Update IMAP after dropping samples above max missingness once."""
        if not self.user_imap:
            self.imap = {"all": list(self.snames)}
            self.minmap = {"all": self.mincov}
            return

        selected = set(self.snames)
        new_imap = {
            group: [name for name in names if name in selected]
            for group, names in self.imap.items()
        }
        empty = sorted(group for group, names in new_imap.items() if not names)
        if empty:
            raise IPyradError(
                "IMAP group(s) became empty after max_sample_missing filtering: "
                + ", ".join(empty)
            )
        self.imap = new_imap

    def _drop_samples_by_missingness_once(self, log_level: str = "INFO") -> bool:
        """Drop high-missing samples once, then rerun SNP filtering."""
        self.sample_missing = self._calculate_missing_frequencies(self.genos)
        if self.max_sample_missing >= 1.0:
            return False
        if self.genos is None or self.genos.shape[1] == 0:
            return False

        keep_mask = self.sample_missing.to_numpy() <= self.max_sample_missing
        if np.all(keep_mask):
            return False

        keep_names = [name for name, keep in zip(self.snames, keep_mask, strict=False) if keep]
        dropped = [
            name for name, keep in zip(self.snames, keep_mask, strict=False) if not keep
        ]
        if not keep_names:
            raise IPyradError("No samples passed max_sample_missing filter.")

        logger.log(
            log_level,
            "dropping {} sample(s) above max_sample_missing {}: {}",
            len(dropped),
            self.max_sample_missing,
            ", ".join(dropped),
        )

        self.dropped_samples_by_missing = dropped
        self.snames = keep_names
        self.sidxs = [self.dbname_to_index[name] for name in self.snames]
        self._sync_imap_after_sample_drop()
        return True

    def _run_filter_pass(
        self, log_level: str = "INFO"
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None, pd.Series]:
        """Parse genotype calls from HDF5 and apply filters once."""
        stats = pd.Series(index=STATS_HEADER, dtype=float)

        mask_arrs = []
        genos_arrs = []
        snpsmap_arrs = []
        snps_arrs = []
        reference_arrs = []
        nmissing = 0
        ntotal = 0
        depth_masked_genotypes = 0

        if self.cores == 1:
            results = {}
            for start in range(0, self.nsnps, CHUNKSIZE):
                results[start] = self._get_masks_chunk(start)
        else:
            jobs = {}
            for start in range(0, self.nsnps, CHUNKSIZE):
                jobs[start] = (self._get_masks_chunk, {"start": start})

            results = run_with_pool(jobs, log_level, self.cores, msg="Filtering SNPs")

        for job in results:
            try:
                snpsmap, snps, genos, reference, masks, nmiss, ntot, ndp = results[job].get()
            except (NameError, AttributeError):
                snpsmap, snps, genos, reference, masks, nmiss, ntot, ndp = results[job]
            snpsmap_arrs.append(snpsmap)
            snps_arrs.append(snps)
            genos_arrs.append(genos)
            mask_arrs.append(masks)
            if reference is not None:
                reference_arrs.append(reference)
            nmissing += nmiss
            ntotal += ntot
            depth_masked_genotypes += ndp

        if not mask_arrs:
            raise ValueError("No SNPs found.")

        mask = np.concatenate(mask_arrs)
        genos = np.concatenate(genos_arrs, axis=1).astype(np.uint8)
        snps = np.concatenate(snps_arrs, axis=1)
        snpsmap = np.concatenate(snpsmap_arrs)
        reference = None
        if reference_arrs:
            reference = np.concatenate(reference_arrs).astype(np.uint8, copy=False)

        if genos.size:
            missing_cells = np.sum(genos == _MISSING_GENO)
            missing_percent = missing_cells / genos.size
        else:
            missing_percent = 1.0

        stats.samples = len(self.snames)
        stats.pre_filter_snps = self.nsnps
        stats.pre_filter_percent_missing = 100 * (nmissing / ntotal) if ntotal else 100.0
        stats.masked_genotypes_by_min_depth = depth_masked_genotypes
        stats.filter_by_indels_present = mask[:, 0].sum()
        stats.filter_by_non_biallelic = mask[:, 1].sum()
        stats.filter_by_mincov = mask[:, 2].sum()
        stats.filter_by_minmap = mask[:, 3].sum()
        stats.filter_by_min_site_qual = mask[:, 4].sum()
        stats.filter_by_invariant_after_subsampling = mask[:, 5].sum()
        stats.filter_by_minor_allele_frequency = mask[:, 6].sum()
        stats.post_filter_snps = snpsmap.shape[0]
        stats.post_filter_snp_containing_linkage_blocks = np.unique(snpsmap[:, 0]).size
        stats.post_filter_percent_missing = 100 * missing_percent
        return snpsmap, snps, genos, reference, stats

    def run(self, log_level: str = "INFO"):
        """Run SNP filtering and apply max-sample-missing pruning once."""
        self.dropped_samples_by_missing = []
        self.snpsmap, self.snps, self.genos, self.reference, self.stats = self._run_filter_pass(
            log_level=log_level
        )
        if self._drop_samples_by_missingness_once(log_level=log_level):
            self.snpsmap, self.snps, self.genos, self.reference, self.stats = self._run_filter_pass(
                log_level=log_level
            )

        self.sample_missing = self._calculate_missing_frequencies(self.genos)
        self.mask = np.zeros(self.snpsmap.shape[0], dtype=bool)
        self.names = self.snames

        logger.info("SNP extraction summary")
        pretty = self.stats.map(lambda value: f"{value:.3f}".rstrip("0").rstrip("."))
        for key, value in pretty.items():
            if key == "pre_filter_percent_missing":
                logger.info(
                    "filter statistic {}: {} (linked genotype cells missing before site filtering)",
                    key,
                    value,
                )
                continue
            if key == "post_filter_percent_missing":
                logger.info(
                    "filter statistic {}: {} (linked post-filter genotype cells missing before optional subsampling)",
                    key,
                    value,
                )
                continue
            logger.info("filter statistic {}: {}", key, value)

    def _get_masks_chunk(
        self, start: int
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None, np.ndarray, int, int, int]:
        """Load one chunk from HDF5 and calculate filters."""
        chunkslice = slice(start, start + CHUNKSIZE)
        with h5py.File(self.data, "r") as io5:
            snps_ds = io5["genos"][:, :, 2]
            genos_ds = io5["genos"][:, :, :2]
            snpsmap = io5["snpsmap"]
            reference = io5["reference"] if "reference" in io5 else None
            sample_dp_ds = io5["sample_dp"] if "sample_dp" in io5 else None
            site_qual_ds = io5["site_qual"] if "site_qual" in io5 else None

            start = chunkslice.start
            end = min(chunkslice.stop, genos_ds.shape[1])
            nsnps = end - start

            snpsmap = snpsmap[start:end, :]
            ref = reference[start:end].astype(np.uint8) if reference is not None else None
            site_qual = (
                site_qual_ds[start:end].astype(np.float32, copy=False)
                if site_qual_ds is not None
                else None
            )
            snps_rows = []
            genos_rows = []
            dp_rows = []
            for name in self.snames:
                row_index = self.dbname_to_genos_index[name]
                if row_index is None:
                    if ref is None:
                        raise IPyradError(
                            "assembly_reference_sequence was selected for SNP analysis, "
                            "but the HDF5 `reference` dataset is missing."
                        )
                    snps_rows.append(ref.copy())
                    genos_rows.append(np.zeros((nsnps, 2), dtype=np.uint8))
                    dp_rows.append(np.full(nsnps, np.iinfo(np.uint32).max, dtype=np.uint32))
                else:
                    snps_rows.append(snps_ds[row_index, start:end])
                    genos_rows.append(genos_ds[row_index, start:end, :].astype(np.uint8))
                    if sample_dp_ds is not None:
                        dp_rows.append(sample_dp_ds[row_index, start:end].astype(np.uint32))
            snps = np.stack(snps_rows, axis=0).astype(np.uint8, copy=False)
            genos = np.stack(genos_rows, axis=0).astype(np.uint8, copy=False)
            sample_dp = (
                np.stack(dp_rows, axis=0).astype(np.uint32, copy=False)
                if sample_dp_ds is not None
                else None
            )

            nmissing = np.sum(genos == _MISSING_GENO)
            ntotal = genos.size
            depth_masked_genotypes = 0

            if self.min_genotype_depth > 0 and sample_dp is not None:
                called_mask = np.any(genos != _MISSING_GENO, axis=2)
                depth_mask = called_mask & (sample_dp < self.min_genotype_depth)
                if np.any(depth_mask):
                    genos = genos.copy()
                    snps = snps.copy()
                    genos[depth_mask] = _MISSING_GENO
                    snps[depth_mask] = _MISSING_SNP
                    depth_masked_genotypes = int(np.count_nonzero(depth_mask))

            masks, diplos = self._masks_filter(nsnps, snps, genos, site_qual=site_qual)
            flat_mask = np.invert(masks.sum(axis=1).astype(bool))

            snpsmap = snpsmap[flat_mask]
            snps = snps[:, flat_mask]
            diplos = diplos[:, flat_mask]
            if ref is not None:
                ref = ref[flat_mask]

        return snpsmap, snps, diplos, ref, masks, nmissing, ntotal, depth_masked_genotypes

    def _masks_filter(
        self,
        nsnps: int,
        snps: np.ndarray,
        genos: np.ndarray,
        *,
        site_qual: np.ndarray | None = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return arrays with filter masks and diploid genotypes."""
        masks = np.zeros((nsnps, 7), dtype=bool)

        masks[:, 0] = np.any(snps == 45, axis=0)
        masks[:, 1] = np.sum(genos == 2, axis=2).sum(axis=0).astype(bool)
        masks[:, 1] += np.sum(genos == 3, axis=2).sum(axis=0).astype(bool)

        genomask = np.ma.array(data=genos, mask=(genos == _MISSING_GENO))
        nhaplos = (~genomask.mask).sum(axis=2).sum(axis=0)

        if isinstance(self.mincov, int):
            masks[:, 2] = nhaplos < (2 * self.mincov)
        elif isinstance(self.mincov, float):
            masks[:, 2] = nhaplos < (2 * self.mincov * len(self.sidxs))
        else:
            raise ValueError("mincov should be an int or float.")

        for pop, samps in self.imap.items():
            mincov = self.minmap[pop]
            imap_sidxs = [self.snames.index(sample) for sample in samps]
            subarr = genomask[imap_sidxs, :, :]
            nhaplos = (~subarr.mask).sum(axis=2).sum(axis=0)
            if isinstance(mincov, int):
                masks[:, 3] += nhaplos < (2 * mincov)
            elif isinstance(mincov, float):
                masks[:, 3] += nhaplos < (2 * mincov * len(imap_sidxs))
            else:
                raise ValueError("minmap dictionary malformed.")

        if self.min_site_qual > 0:
            if site_qual is None:
                raise IPyradError(
                    "The `--min-site-qual` filter requires the HDF5 `site_qual` dataset."
                )
            masks[:, 4] = np.isnan(site_qual) | (site_qual < self.min_site_qual)

        diplo_common = (
            genomask.sum(axis=2).mean(axis=0).round().astype(int).data
        )
        diplos = genomask.sum(axis=2).data
        masks[:, 5] = np.all(diplo_common == diplos, axis=0)

        called_0 = (genomask == 0).sum(axis=2).sum(axis=0).data
        called_1 = (genomask == 1).sum(axis=2).sum(axis=0).data

        with np.errstate(divide="ignore", invalid="ignore"):
            if isinstance(self.maf, int):
                freqs = called_1
            else:
                freqs = called_1 / (called_0 + called_1)
                freqs[freqs > 0.5] = 1 - freqs[freqs > 0.5]
            masks[:, 6] = freqs < self.maf

        diplos[snps == _MISSING_SNP] = _MISSING_GENO
        return masks, diplos

    def _require_run(self) -> None:
        if self.genos is None or self.snps is None or self.snpsmap is None:
            raise IPyradError("SNPsExtracter.run() must be called before accessing data.")

    def subsample_column_indices(
        self, random_seed: Optional[int] = None, log_level: str = "INFO"
    ) -> np.ndarray:
        """Return filtered column indices for one SNP sampled per linkage block."""
        self._require_run()
        if self.snpsmap.shape[0] == 0:
            return np.zeros(0, dtype=np.int64)

        rng = np.random.default_rng(random_seed)
        locs = self.snpsmap[:, 0].astype(np.int64, copy=False)
        _, starts, counts = np.unique(locs, return_index=True, return_counts=True)
        indices = np.empty(starts.size, dtype=np.int64)
        for idx, (start, count) in enumerate(zip(starts, counts, strict=False)):
            indices[idx] = start + int(rng.integers(count))
        logger.log(
            log_level,
            "subsampled {} unlinked SNPs from {} filtered SNPs.",
            indices.size,
            self.snpsmap.shape[0],
        )
        return indices

    def get_view(
        self,
        *,
        subsample: bool = False,
        random_seed: Optional[int] = None,
        log_level: str = "INFO",
    ) -> SNPExportView:
        """Return linked or unlinked aligned outputs from the filtered matrix."""
        self._require_run()
        if not subsample:
            return SNPExportView(self.genos, self.snps, self.snpsmap, self.reference)

        indices = self.subsample_column_indices(random_seed=random_seed, log_level=log_level)
        reference = None if self.reference is None else self.reference[indices]
        return SNPExportView(
            self.genos[:, indices],
            self.snps[:, indices],
            self.snpsmap[indices],
            reference,
        )

    def subsample_snps(
        self, random_seed: Optional[int] = None, log_level: str = "INFO"
    ) -> np.ndarray:
        """Return one SNP character column sampled per linkage block."""
        return self.get_view(
            subsample=True, random_seed=random_seed, log_level=log_level
        ).snps

    def subsample_genos(
        self, random_seed: Optional[int] = None, log_level: str = "INFO"
    ) -> np.ndarray:
        """Return one diploid genotype column sampled per linkage block."""
        return self.get_view(
            subsample=True, random_seed=random_seed, log_level=log_level
        ).genos

    def subsample_loci(
        self,
        random_seed: Optional[int] = None,
        return_sites: bool = False,
        log_level: str = "INFO",
    ) -> np.ndarray:
        """Return loci resampled with replacement from the filtered matrix."""
        self._require_run()
        if self.snpsmap.shape[0] == 0:
            return self.snps[:, :0] if return_sites else self.genos[:, :0]

        rng = np.random.default_rng(random_seed)
        locs = self.snpsmap[:, 0].astype(np.int64, copy=False)
        unique_locs, starts, counts = np.unique(
            locs, return_index=True, return_counts=True
        )
        spans = {
            int(loc): np.arange(start, start + count, dtype=np.int64)
            for loc, start, count in zip(unique_locs, starts, counts, strict=False)
        }
        sampled_locs = rng.choice(unique_locs, size=unique_locs.size, replace=True)
        indices = np.concatenate([spans[int(loc)] for loc in sampled_locs]).astype(np.int64)
        subarr = self.snps[:, indices] if return_sites else self.genos[:, indices]
        logger.log(
            log_level,
            "subsampled {} SNPs from {} variable loci w/ replacement.",
            subarr.shape[1],
            unique_locs.size,
        )
        return subarr

    def get_population_geno_counts(
        self,
        subsample: bool = False,
        random_seed: Optional[int] = None,
        genos: np.ndarray | None = None,
        imap: Dict[str, List[str]] | None = None,
        log_level="INFO",
    ):
        """Return dataframe with genotypes in treemix count format."""
        if imap is None:
            if self.user_imap:
                imap = self.imap.copy()
            else:
                imap = {name: [name] for name in self.names}

        data = pd.DataFrame(columns=list(imap))
        if genos is not None:
            pass
        elif subsample:
            genos = self.subsample_genos(random_seed=random_seed, log_level=log_level)
        else:
            genos = self.genos

        for pop in data.columns:
            samp = [self.names.index(name) for name in imap[pop]]
            ances = np.sum(genos[samp, :] == 0, axis=0) * 2
            deriv = np.sum(genos[samp, :] == 2, axis=0) * 2
            heter = np.sum(genos[samp, :] == 1, axis=0)
            ances += heter
            deriv += heter
            data.loc[:, pop] = [f"{i},{j}" for i, j in zip(ances, deriv, strict=False)]
        return data

    def get_population_geno_frequency(
        self,
        subsample: bool = False,
        random_seed=None,
        genos: np.ndarray | None = None,
        imap: Dict[str, List[str]] | None = None,
        log_level="INFO",
    ):
        """Return a dataframe with genotype frequencies as in construct format."""
        if imap is None:
            if self.user_imap:
                imap = self.imap.copy()
            else:
                imap = {name: [name] for name in self.names}

        data = pd.DataFrame(columns=list(imap))
        if genos is not None:
            pass
        elif subsample:
            genos = self.subsample_genos(random_seed=random_seed, log_level=log_level)
        else:
            genos = self.genos

        for pop in data.columns:
            samp = [self.names.index(name) for name in imap[pop]]
            ances = np.sum(genos[samp, :] == 0, axis=0) * 2
            deriv = np.sum(genos[samp, :] == 2, axis=0) * 2
            heter = np.sum(genos[samp, :] == 1, axis=0)
            ances += heter
            deriv += heter
            with np.errstate(divide="ignore", invalid="ignore"):
                data.loc[:, pop] = deriv / (deriv + ances)
        return data.T

    def write_plink(
        self,
        prefix: Path,
        view: SNPExportView | None = None,
        *,
        genos: np.ndarray | None = None,
        snps: np.ndarray | None = None,
        snpsmap: np.ndarray | None = None,
        reference: np.ndarray | None = None,
        impute_method: str | None = None,
    ) -> dict[str, Path]:
        """Write PLINK BED/BIM/FAM files for the selected SNP view."""
        if view is not None:
            if any(value is not None for value in (genos, snps, snpsmap, reference)):
                raise IPyradError(
                    "write_plink accepts either `view` or explicit genos/snps/snpsmap/reference data, not both."
                )
            genos = view.genos
            snps = view.snps
            snpsmap = view.snpsmap
            reference = view.reference
            if impute_method is not None:
                from ..methods.common import normalize_impute_method
                from ..methods.snps_imputer import SNPsImputer

                method = normalize_impute_method(impute_method)
                genos = SNPsImputer(
                    genos,
                    self.snames,
                    imap=self.imap,
                    impute_method=method,
                    quiet=True,
                ).run()

        if any(value is None for value in (genos, snps, snpsmap)):
            raise IPyradError(
                "write_plink requires either `view` or explicit genos/snps/snpsmap data."
            )
        if reference is None:
            raise IPyradError(
                "PLINK export requires the HDF5 `reference` dataset. "
                "Rebuild the SNP HDF5 with a current assemble or `ipyrad2 vcf2hdf5` run."
            )

        paths = {
            "bed": prefix.with_suffix(".bed"),
            "bim": prefix.with_suffix(".bim"),
            "fam": prefix.with_suffix(".fam"),
        }

        with open(paths["fam"], "w", encoding="utf-8") as out:
            for sample in self.snames:
                out.write(f"{sample}\t{sample}\t0\t0\t0\t-9\n")

        with open(paths["bim"], "w", encoding="utf-8") as out:
            for idx in range(snpsmap.shape[0]):
                loc, _loc_idx, _loc_pos, scaff, pos = [
                    int(value) for value in snpsmap[idx]
                ]
                ref_base = _decode_base(reference[idx])
                alt_base = _infer_alt_allele(ref_base, snps[:, idx])
                chrom = str(scaff + 1)
                snp_id = f"loc{loc}_pos{pos + 1}"
                out.write(
                    f"{chrom}\t{snp_id}\t0\t{pos + 1}\t{ref_base}\t{alt_base}\n"
                )

        with open(paths["bed"], "wb") as out:
            out.write(bytes([0x6C, 0x1B, 0x01]))
            for idx in range(genos.shape[1]):
                out.write(_encode_plink_site(genos[:, idx]))

        return paths


def _decode_base(value: np.uint8) -> str:
    """Convert one uint8 ASCII code to a base string."""
    base = chr(int(value))
    if base not in {"A", "C", "G", "T"}:
        raise IPyradError(f"PLINK export requires A/C/G/T reference alleles, found {base!r}.")
    return base


def _infer_alt_allele(ref_base: str, calls: np.ndarray) -> str:
    """Infer the non-reference allele from filtered site calls."""
    for value in calls:
        if int(value) == _MISSING_SNP:
            continue
        base = chr(int(value))
        if base in {"A", "C", "G", "T"} and base != ref_base:
            return base
        if base in _IUPAC_TO_BASES:
            alleles = set(_IUPAC_TO_BASES[base])
            if ref_base in alleles:
                alt = sorted(alleles.difference({ref_base}))
                if alt:
                    return alt[0]
    raise IPyradError(
        f"Could not infer ALT allele for PLINK export at reference allele {ref_base!r}."
    )


def _encode_plink_site(genotypes: np.ndarray) -> bytes:
    """Encode one SNP's diploid genotype column into PLINK BED SNP-major bytes."""
    mapping = {
        0: 0b00,  # homozygous allele1 (reference)
        1: 0b10,  # heterozygous
        2: 0b11,  # homozygous allele2 (alternate)
        _MISSING_GENO: 0b01,  # missing
    }
    nbytes = (genotypes.shape[0] + 3) // 4
    encoded = bytearray(nbytes)
    for idx, genotype in enumerate(genotypes):
        code = mapping.get(int(genotype), 0b01)
        encoded[idx // 4] |= code << ((idx % 4) * 2)
    return bytes(encoded)


def _encode_snp_char(genotype: int, ref_base: str, alt_base: str) -> str:
    """Return one SNP character from a diploid genotype and reference/alt alleles."""
    if genotype == _MISSING_GENO:
        return "N"
    if genotype == 0:
        return ref_base
    if genotype == 2:
        return alt_base
    if genotype == 1:
        code = _BASES_TO_IUPAC.get(tuple(sorted((ref_base, alt_base))))
        if code is None:
            raise IPyradError(
                f"Could not encode heterozygous SNP for alleles {ref_base!r}/{alt_base!r}."
            )
        return code
    raise IPyradError(f"Unsupported diploid genotype for SNP export: {genotype}.")


def _reconstruct_snp_chars(
    *,
    genos: np.ndarray,
    snps: np.ndarray,
    reference: np.ndarray | None,
) -> np.ndarray:
    """Return an SNP character matrix aligned to one genotype matrix."""
    if reference is None:
        raise IPyradError(
            "SNP character reconstruction requires the HDF5 `reference` dataset. "
            "Rebuild the SNP HDF5 with a current assemble or `ipyrad2 vcf2hdf5` run."
        )
    out = np.empty(genos.shape, dtype=np.uint8)
    for idx in range(genos.shape[1]):
        ref_base = _decode_base(reference[idx])
        alt_base = _infer_alt_allele(ref_base, snps[:, idx])
        for sidx in range(genos.shape[0]):
            out[sidx, idx] = ord(
                _encode_snp_char(int(genos[sidx, idx]), ref_base, alt_base)
            )
    return out


def _format_phy_alignment(samples: list[str], snps: np.ndarray) -> str:
    """Return a PHYLIP alignment string from one SNP character matrix."""
    ntaxa = len(samples)
    nsites = snps.shape[1]
    longname = max(len(name) for name in samples)
    padded = [name.ljust(longname + 5) for name in samples]
    rows = []
    for idx, name in enumerate(padded):
        seq = snps[idx].tobytes().decode("utf-8")
        rows.append(f"{name} {seq}")
    return f"{ntaxa} {nsites}\n" + "\n".join(rows) + "\n"


def _format_nexus_alignment(samples: list[str], snps: np.ndarray) -> str:
    """Return a NEXUS alignment string from one SNP character matrix."""
    ntaxa = len(samples)
    nsites = snps.shape[1]
    longname = max(len(name) for name in samples)
    padded = [name.ljust(longname + 5) for name in samples]
    lines = [NEXHEADER.format(ntaxa, nsites)]
    for block in range(0, snps.shape[1], 100):
        stop = min(block + 100, snps.shape[1])
        for idx, name in enumerate(padded):
            seq = snps[idx, block:stop].tobytes().decode("utf-8")
            lines.append(f"  {name}{seq}\n")
        lines.append("\n")
    lines.append("  ;\nend;")
    return "".join(lines)


def _format_fasta_alignment(samples: list[str], snps: np.ndarray) -> str:
    """Return a FASTA alignment string from one SNP character matrix."""
    records = []
    for idx, sample in enumerate(samples):
        seq = snps[idx].tobytes().decode("utf-8")
        records.append(f">{sample}\n{seq}")
    return "\n".join(records) + "\n"


def _write_text(path: Path, contents: str) -> None:
    with open(path, "w", encoding="utf-8") as out:
        out.write(contents)


def _eems_dissimilarity_matrix(genos: np.ndarray, *, use_mean_fill: bool) -> np.ndarray:
    """Return the EEMS pairwise-difference matrix from one genotype matrix."""
    matrix = genos.astype(float, copy=True)
    if use_mean_fill:
        missing = matrix == _MISSING_GENO
        with np.errstate(invalid="ignore"):
            col_means = np.nanmean(np.where(missing, np.nan, matrix), axis=0)
        col_means = np.where(np.isnan(col_means), 0.0, col_means)
        matrix[missing] = np.broadcast_to(col_means, matrix.shape)[missing]
    similarities = matrix @ matrix.T / matrix.shape[1]
    self_similarities = np.diag(similarities)
    diffs = (
        self_similarities[:, np.newaxis]
        + self_similarities[np.newaxis, :]
        - (2.0 * similarities)
    )
    np.fill_diagonal(diffs, 0.0)
    return diffs


def _prepare_snp_export(
    *,
    extracter: SNPsExtracter,
    view: SNPExportView,
    impute_method: str | None,
    random_seed: int | None,
) -> PreparedSNPExport:
    """Return one prepared export bundle shared across all SNP writers."""
    from ..methods.common import (
        calculate_sample_missing_fraction,
        build_imputed_sample_data_summary,
        build_sample_data_summary,
        impute_genotype_matrix,
        normalize_impute_method,
        summarize_imputation,
    )

    method = normalize_impute_method(impute_method)
    if method is None:
        missing = calculate_sample_missing_fraction(view.genos, extracter.snames)
        sample_summary = build_sample_data_summary(
            samples=extracter.snames,
            missing_fraction=missing,
            post_imputation_missing_fraction=missing,
            imputation_algorithm="not-imputed",
        )
        return PreparedSNPExport(
            view=view,
            genos=view.genos,
            snps=view.snps,
            imputation_summary=None,
            sample_data_summary=sample_summary,
        )

    if view.reference is None:
        raise IPyradError(
            "Global SNP imputation for snpex requires the HDF5 `reference` dataset so all "
            "written outputs, including the SNP character matrix, can be reconstructed. "
            "Rebuild the SNP HDF5 with a current assemble or `ipyrad2 vcf2hdf5` run."
        )

    imputation = summarize_imputation(view.genos, method)
    genos = impute_genotype_matrix(
        view.genos,
        extracter,
        impute_method=method,
        random_seed=random_seed,
    )
    snps = _reconstruct_snp_chars(genos=genos, snps=view.snps, reference=view.reference)
    sample_summary = build_imputed_sample_data_summary(
        samples=extracter.snames,
        matrix=view.genos,
        impute_method=method,
    )
    return PreparedSNPExport(
        view=view,
        genos=genos,
        snps=snps,
        imputation_summary=imputation,
        sample_data_summary=sample_summary,
    )


def _write_snps_stats(
    path: Path,
    extracter: SNPsExtracter,
    prepared: PreparedSNPExport,
    *,
    subsample: bool,
    random_seed: int | None,
    impute_method: str | None,
    written_formats: list[str],
) -> None:
    with open(path, "w", encoding="utf-8") as out:
        out.write("Summary\n")
        out.write("-------\n")
        out.write(f"infile: {extracter.data}\n")
        out.write(f"samples_selected_initial: {extracter.initial_snames}\n")
        out.write(f"samples_dropped_by_max_missing: {extracter.dropped_samples_by_missing}\n")
        out.write(f"samples_final: {extracter.snames}\n")
        out.write(f"imap: {extracter.imap}\n")
        out.write(f"minmap: {extracter.minmap}\n")
        out.write(f"include_reference: {extracter.include_reference}\n")
        out.write(f"max_sample_missing: {extracter.max_sample_missing}\n")
        out.write(f"min_genotype_depth: {extracter.min_genotype_depth}\n")
        out.write(f"min_site_qual: {extracter.min_site_qual}\n")
        out.write(f"subsample: {subsample}\n")
        out.write(f"random_seed: {random_seed}\n")
        out.write(f"impute_method: {impute_method if impute_method is not None else 'none'}\n")
        out.write(
            "imputation_algorithm: "
            f"{prepared.sample_data_summary['imputation_algorithm'].iloc[0]}\n"
        )
        if prepared.imputation_summary is None:
            out.write("imputed_snp_count: 0\n")
            out.write("imputed_snp_fraction: 0.0\n")
            out.write("imputed_genotype_count: 0\n")
            out.write("imputed_genotype_fraction: 0.0\n")
        else:
            out.write(
                f"imputed_snp_count: {prepared.imputation_summary.imputed_snp_count}\n"
            )
            out.write(
                f"imputed_snp_fraction: {prepared.imputation_summary.imputed_snp_fraction}\n"
            )
            out.write(
                "imputed_genotype_count: "
                f"{prepared.imputation_summary.imputed_genotype_count}\n"
            )
            out.write(
                "imputed_genotype_fraction: "
                f"{prepared.imputation_summary.imputed_genotype_fraction}\n"
            )
        out.write(f"written_formats: {', '.join(written_formats)}\n")
        out.write(
            "linked_post_filter_snps: "
            f"{int(extracter.stats['post_filter_snps'])}\n"
        )
        out.write(
            "linked_post_filter_snp_containing_linkage_blocks: "
            f"{int(extracter.stats['post_filter_snp_containing_linkage_blocks'])}\n"
        )
        out.write(f"exported_snps: {prepared.view.snpsmap.shape[0]}\n")
        out.write(
            "exported_snp_containing_linkage_blocks: "
            f"{np.unique(prepared.view.snpsmap[:, 0]).size if prepared.view.snpsmap.size else 0}\n"
        )
        out.write("\n")
        out.write("Filter statistics\n")
        out.write("-----------------\n")
        for key in extracter.stats.index:
            out.write(f"{key}: {extracter.stats[key]}\n")


def run_snps_extracter(
    *,
    data: Path | str,
    name: str,
    outdir: Path | str,
    min_sample_coverage: float,
    max_sample_missing: float,
    min_minor_allele_frequency: float,
    imap: Path | str | Dict | None,
    minmap: Path | Dict | None,
    min_genotype_depth: int = 0,
    min_site_qual: float = 0.0,
    exclude: Path | str | List | None = None,
    include_reference: bool = False,
    cores: int = 1,
    force: bool = False,
    log_level: str = "INFO",
    subsample: bool = True,
    random_seed: int | None = None,
    write_plink: bool = False,
    write_phylip: bool = False,
    write_nexus: bool = False,
    write_fasta: bool = False,
    write_treemix: bool = False,
    write_eems: bool = False,
    impute_method: str | None = None,
) -> None:
    from ..methods.common import (
        log_snp_imputation_summary,
        log_snp_view_summary,
        normalize_impute_method,
        summarize_prepared_snp_view,
        write_sample_data_summary,
    )

    outdir = Path(outdir).expanduser().absolute()
    prefix = outdir / name
    paths = {
        "genos": outdir / f"{name}.genos.npy",
        "snps": outdir / f"{name}.snps.npy",
        "snpsmap": outdir / f"{name}.snpsmap.tsv",
        "samples": outdir / f"{name}.samples.txt",
        "sample_data_summary": outdir / f"{name}.sample_data_summary.tsv",
        "stats": outdir / f"{name}.stats.txt",
    }
    if write_plink:
        paths.update(
            {
                "plink_bed": prefix.with_suffix(".bed"),
                "plink_bim": prefix.with_suffix(".bim"),
                "plink_fam": prefix.with_suffix(".fam"),
            }
        )
    if write_phylip:
        paths["phylip"] = prefix.with_suffix(".phy")
    if write_nexus:
        paths["nexus"] = prefix.with_suffix(".nex")
    if write_fasta:
        paths["fasta"] = prefix.with_suffix(".fa")
    if write_treemix:
        paths["treemix"] = outdir / f"{name}.treemix.gz"
    if write_eems:
        paths["eems"] = outdir / f"{name}.eems"

    normalized_impute = normalize_impute_method(impute_method)

    existing = next((path for path in paths.values() if path.exists()), None)
    if existing is not None and not force:
        raise IPyradError(
            f"Output file already exists: {existing}. Use --force to overwrite."
        )

    outdir.mkdir(parents=True, exist_ok=True)
    tool = SNPsExtracter(
        data=Path(data),
        min_sample_coverage=min_sample_coverage,
        max_sample_missing=max_sample_missing,
        min_minor_allele_frequency=min_minor_allele_frequency,
        min_genotype_depth=min_genotype_depth,
        min_site_qual=min_site_qual,
        imap=imap,
        minmap=minmap,
        exclude=exclude,
        include_reference=include_reference,
        cores=cores,
    )
    tool.run(log_level=log_level)
    view = tool.get_view(
        subsample=subsample,
        random_seed=random_seed,
        log_level="DEBUG",
    )
    prepared = _prepare_snp_export(
        extracter=tool,
        view=view,
        impute_method=normalized_impute,
        random_seed=random_seed,
    )
    log_snp_imputation_summary("snpex", prepared.imputation_summary, subsample=subsample)
    log_snp_view_summary(
        "snpex",
        summarize_prepared_snp_view(tool, prepared.view, subsample=subsample),
        view_label="exported",
    )

    np.save(paths["genos"], prepared.genos, allow_pickle=False)
    np.save(paths["snps"], prepared.snps, allow_pickle=False)
    pd.DataFrame(prepared.view.snpsmap, columns=SNPSMAP_COLUMNS).to_csv(
        paths["snpsmap"],
        sep="\t",
        index=False,
    )
    with open(paths["samples"], "w", encoding="utf-8") as out:
        out.write("\n".join(tool.snames) + "\n")
    write_sample_data_summary(paths["sample_data_summary"], prepared.sample_data_summary)

    written_formats = ["genos", "snps", "snpsmap", "samples", "sample_data_summary"]

    if write_plink:
        plink_paths = tool.write_plink(
            prefix,
            genos=prepared.genos,
            snps=prepared.snps,
            snpsmap=prepared.view.snpsmap,
            reference=prepared.view.reference,
        )
        logger.info("wrote PLINK BED to {}", plink_paths["bed"])
        logger.info("wrote PLINK BIM to {}", plink_paths["bim"])
        logger.info("wrote PLINK FAM to {}", plink_paths["fam"])
        written_formats.append("plink")

    if write_phylip:
        _write_text(paths["phylip"], _format_phy_alignment(tool.snames, prepared.snps).rstrip("\n"))
        logger.info(
            "wrote PHYLIP alignment ({}, {}) to {}",
            len(tool.snames),
            prepared.snps.shape[1],
            paths["phylip"],
        )
        written_formats.append("phylip")

    if write_nexus:
        _write_text(paths["nexus"], _format_nexus_alignment(tool.snames, prepared.snps))
        logger.info(
            "wrote NEXUS alignment ({}, {}) to {}",
            len(tool.snames),
            prepared.snps.shape[1],
            paths["nexus"],
        )
        written_formats.append("nexus")

    if write_fasta:
        _write_text(paths["fasta"], _format_fasta_alignment(tool.snames, prepared.snps))
        logger.info(
            "wrote FASTA alignment ({}, {}) to {}",
            len(tool.snames),
            prepared.snps.shape[1],
            paths["fasta"],
        )
        written_formats.append("fasta")

    if write_treemix:
        counts = tool.get_population_geno_counts(
            genos=prepared.genos,
            imap=None if not tool.user_imap else tool.imap,
        )
        with gzip.open(paths["treemix"], "wt", encoding="utf-8") as out:
            out.write(" ".join(counts.columns) + "\n")
            for row in counts.itertuples(index=False, name=None):
                out.write(" ".join(str(val) for val in row) + "\n")
        logger.info("wrote TreeMix counts to {}", paths["treemix"])
        written_formats.append("treemix")

    if write_eems:
        diffs = _eems_dissimilarity_matrix(
            prepared.genos,
            use_mean_fill=(normalized_impute is None),
        )
        pd.DataFrame(diffs, index=tool.snames, columns=tool.snames).to_csv(
            paths["eems"],
            sep="\t",
            header=False,
            index=False,
            float_format="%.12g",
        )
        logger.info("wrote EEMS genetic matrix to {}", paths["eems"])
        written_formats.append("eems")

    _write_snps_stats(
        paths["stats"],
        tool,
        prepared,
        subsample=subsample,
        random_seed=random_seed,
        impute_method=normalized_impute,
        written_formats=written_formats,
    )

    logger.info("wrote exported diploid genotypes to {}", paths["genos"])
    logger.info("wrote exported SNP character matrix to {}", paths["snps"])
    logger.info("wrote exported SNP map to {}", paths["snpsmap"])
    logger.info("wrote filtered sample order to {}", paths["samples"])
    logger.info("wrote sample data summary to {}", paths["sample_data_summary"])
    logger.info("wrote stats/log to: {}", paths["stats"])


__all__ = [
    "CHUNKSIZE",
    "REFERENCE_SAMPLE_NAME",
    "SNPSMAP_COLUMNS",
    "SNPExportView",
    "SNPsExtracter",
    "_MISSING_GENO",
    "_MISSING_SNP",
    "run_snps_extracter",
]
