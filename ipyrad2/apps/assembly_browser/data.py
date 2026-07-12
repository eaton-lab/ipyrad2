"""Data access helpers for the Streamlit assembly browser."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterator

import h5py
import numpy as np
import pandas as pd


REFERENCE_SAMPLE_NAME = "assembly_reference_sequence"
SNPSMAP_COLUMNS = ["loc", "loc_idx", "loc_pos", "scaff", "pos"]


def _hdf5_access_error(path: Path, exc: Exception) -> RuntimeError:
    """Return an actionable HDF5 read error."""
    return RuntimeError(
        f"Could not read HDF5 datasets from {path}. If the error mentions a bad "
        "layout message version, open the file with an environment using an HDF5 "
        "library new enough for the file that created it."
    )


def _decode_h5_strings(values) -> list[str]:
    """Decode HDF5 string-like arrays into plain Python strings."""
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


@dataclass(frozen=True)
class AssemblyOutputs:
    """Known files belonging to one assembly output prefix."""

    outdir: Path
    prefix: str
    hdf5: Path
    stats_json: Path | None = None
    stats_txt: Path | None = None
    vcf: Path | None = None
    loci: Path | None = None
    bed: Path | None = None


def discover_outputs(path: str | Path) -> AssemblyOutputs:
    """Find the primary ipyrad2 output files from an output directory or HDF5 path."""
    path = Path(path).expanduser().resolve()
    if path.is_file():
        hdf5 = path
        outdir = path.parent
        prefix = path.stem
    else:
        outdir = path
        hdf5_files = sorted(outdir.glob("*.hdf5"))
        if not hdf5_files:
            raise FileNotFoundError(f"No .hdf5 assembly file found in {outdir}")
        hdf5 = hdf5_files[0]
        prefix = hdf5.stem

    def existing(name: str) -> Path | None:
        candidate = outdir / name
        return candidate if candidate.exists() else None

    return AssemblyOutputs(
        outdir=outdir,
        prefix=prefix,
        hdf5=hdf5,
        stats_json=existing(f"{prefix}.stats.json"),
        stats_txt=existing(f"{prefix}.stats.txt"),
        vcf=existing(f"{prefix}.vcf.gz"),
        loci=existing(f"{prefix}.loci.gz"),
        bed=existing(f"{prefix}.bed"),
    )


def load_stats_json(path: Path | None) -> dict:
    """Load stats JSON if present."""
    if path is None:
        return {}
    with path.open(encoding="utf-8") as infile:
        return json.load(infile)


class AssemblyStore:
    """Lazy HDF5-backed access to assembly data and static stats."""

    def __init__(self, outputs: AssemblyOutputs):
        self.outputs = outputs
        self.stats = load_stats_json(outputs.stats_json)

    @classmethod
    def from_path(cls, path: str | Path) -> "AssemblyStore":
        """Create a store from an output directory or an HDF5 path."""
        return cls(discover_outputs(path))

    def _open(self) -> h5py.File:
        return h5py.File(self.outputs.hdf5, "r")

    @property
    def summary(self) -> dict:
        """Return top-level assembly summary stats."""
        return dict(self.stats.get("summary", {}))

    @property
    def baseline_sample_summary(self) -> pd.DataFrame:
        """Return the assembly-time sample summary table."""
        rows = self.stats.get("sample_summary", [])
        return pd.DataFrame(rows)

    @property
    def baseline_locus_occupancy(self) -> pd.DataFrame:
        """Return the assembly-time locus occupancy table."""
        rows = self.stats.get("locus_occupancy", [])
        return pd.DataFrame(rows)

    def metadata(self) -> dict:
        """Return lightweight HDF5 metadata for display and UI defaults."""
        with self._open() as h5:
            try:
                attrs = h5.attrs
                genos = h5["genos"]
                snpsmap = h5["snpsmap"]
            except Exception as exc:
                raise _hdf5_access_error(self.outputs.hdf5, exc)
            dataset_names = genos.attrs.get("names")
            sample_names = (
                _decode_h5_strings(dataset_names)
                if dataset_names is not None
                else self.root_sample_names()
            )
            scaffold_names = _decode_h5_strings(attrs.get("scaffold_names", []))
            scaffold_lengths = [int(value) for value in attrs.get("scaffold_lengths", [])]
            return {
                "hdf5": str(self.outputs.hdf5),
                "prefix": self.outputs.prefix,
                "nsnps": int(attrs.get("nsnps", snpsmap.shape[0])),
                "samples": sample_names,
                "sample_count": int(genos.shape[0]),
                "scaffold_count": len(scaffold_names),
                "scaffold_names": scaffold_names,
                "scaffold_lengths": scaffold_lengths,
                "has_sample_dp": "sample_dp" in h5,
                "has_site_qual": "site_qual" in h5,
                "has_reference": "reference" in h5,
            }

    def root_sample_names(self) -> list[str]:
        """Return root-level sample names, including reference if recorded there."""
        with self._open() as h5:
            try:
                return _decode_h5_strings(h5.attrs.get("names", []))
            except Exception as exc:
                raise _hdf5_access_error(self.outputs.hdf5, exc)

    def genotype_sample_names(self) -> list[str]:
        """Return sample names aligned to rows of the `genos` dataset."""
        with self._open() as h5:
            try:
                genos = h5["genos"]
            except Exception as exc:
                raise _hdf5_access_error(self.outputs.hdf5, exc)
            dataset_names = genos.attrs.get("names")
            if dataset_names is not None:
                return _decode_h5_strings(dataset_names)
            root_names = _decode_h5_strings(h5.attrs.get("names", []))
            if len(root_names) == genos.shape[0]:
                return root_names
            return [name for name in root_names if name != REFERENCE_SAMPLE_NAME]

    def sample_name_to_index(self) -> dict[str, int]:
        """Return genotype-row indices keyed by sample name."""
        return {name: idx for idx, name in enumerate(self.genotype_sample_names())}

    def site_chunks(
        self,
        sample_indices: list[int],
        *,
        chunk_size: int = 50_000,
    ) -> Iterator[dict[str, np.ndarray]]:
        """Yield SNP data chunks for selected samples."""
        with self._open() as h5:
            try:
                genos_ds = h5["genos"]
                snpsmap_ds = h5["snpsmap"]
                sample_dp_ds = h5["sample_dp"] if "sample_dp" in h5 else None
                site_qual_ds = h5["site_qual"] if "site_qual" in h5 else None
            except Exception as exc:
                raise _hdf5_access_error(self.outputs.hdf5, exc)
            nsnps = int(genos_ds.shape[1])
            row_index = np.asarray(sample_indices, dtype=np.int64)
            read_order = np.argsort(row_index)
            sorted_index = row_index[read_order]
            restore_order = np.argsort(read_order)
            for start in range(0, nsnps, chunk_size):
                end = min(start + chunk_size, nsnps)
                genos = genos_ds[sorted_index, start:end, :2].astype(np.uint8)
                snps = genos_ds[sorted_index, start:end, 2].astype(np.uint8)
                if restore_order.size:
                    genos = genos[restore_order]
                    snps = snps[restore_order]
                sample_dp = None
                if sample_dp_ds is not None:
                    sample_dp = sample_dp_ds[sorted_index, start:end].astype(np.uint32)
                    if restore_order.size:
                        sample_dp = sample_dp[restore_order]
                yield {
                    "start": np.asarray(start),
                    "genos": genos,
                    "snps": snps,
                    "snpsmap": snpsmap_ds[start:end, :].astype(np.int64),
                    "sample_dp": sample_dp,
                    "site_qual": (
                        site_qual_ds[start:end].astype(np.float32)
                        if site_qual_ds is not None
                        else None
                    ),
                }
