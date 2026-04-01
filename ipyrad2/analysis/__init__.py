"""Canonical public surface for phase-1 and phase-2 analysis tools and helpers."""

from importlib import import_module

from .converters.vcf_to_hdf5 import VCFToHDF5, run_vcf_to_hdf5
from .extractors.locus_extractor import LocusExtractor, run_locus_extractor
from .extractors.snp_extractor import SNPExtractor, run_snp_extractor
from .extractors.window_extractor import WindowExtractor, run_window_extractor

__all__ = [
    "Bpp",
    "LocusExtractor",
    "PCAFamilyResult",
    "SNPExtractor",
    "SNPImputer",
    "VCFToHDF5",
    "WindowExtractor",
    "run_admixture_method",
    "run_dapc_method",
    "run_locus_extractor",
    "run_pca_analysis",
    "run_pca_method",
    "run_popgen_method",
    "run_snp_extractor",
    "run_snmf_method",
    "run_tsne_analysis",
    "run_umap_analysis",
    "run_vcf_to_hdf5",
    "run_window_extractor",
]


_LAZY_EXPORTS = {
    "SNPImputer": (".methods.snp_imputer", "SNPImputer"),
    "Bpp": (".methods.bpp", "Bpp"),
    "run_admixture_method": (".methods.admixture", "run_admixture_method"),
    "run_dapc_method": (".methods.dapc", "run_dapc_method"),
    "run_popgen_method": (".methods.popgen", "run_popgen_method"),
    "run_snmf_method": (".methods.snmf", "run_snmf_method"),
    "PCAFamilyResult": (".methods.pca", "PCAFamilyResult"),
    "run_pca_analysis": (".methods.pca", "run_pca_analysis"),
    "run_pca_method": (".methods.pca", "run_pca_method"),
    "run_tsne_analysis": (".methods.pca", "run_tsne_analysis"),
    "run_umap_analysis": (".methods.pca", "run_umap_analysis"),
}


def __getattr__(name: str):
    """Load heavier or optional analysis helpers lazily."""
    if name in _LAZY_EXPORTS:
        module_name, attr_name = _LAZY_EXPORTS[name]
        module = import_module(module_name, __name__)
        value = getattr(module, attr_name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
