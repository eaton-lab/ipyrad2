"""Canonical public surface for phase-1 and phase-2 analysis tools and helpers."""

from importlib import import_module

__all__ = [
    "Bpp",
    "LocusExtracter",
    "PCAFamilyResult",
    "SNPsExtracter",
    "SNPsImputer",
    "VCFToHDF5",
    "WindowExtracter",
    "run_admixture_method",
    "run_bpp_method",
    "run_dapc_method",
    "run_locus_extracter",
    "run_pca_analysis",
    "run_pca_method",
    "run_popgen_method",
    "run_snps_extracter",
    "run_snmf_method",
    "run_tsne_analysis",
    "run_umap_analysis",
    "run_vcf_to_hdf5",
    "run_window_extracter",
]


_LAZY_EXPORTS = {
    "VCFToHDF5": (".converters.vcf_to_hdf5", "VCFToHDF5"),
    "run_vcf_to_hdf5": (".converters.vcf_to_hdf5", "run_vcf_to_hdf5"),
    "LocusExtracter": (".extracters.locus_extracter", "LocusExtracter"),
    "run_locus_extracter": (".extracters.locus_extracter", "run_locus_extracter"),
    "SNPsExtracter": (".extracters.snps_extracter", "SNPsExtracter"),
    "run_snps_extracter": (".extracters.snps_extracter", "run_snps_extracter"),
    "WindowExtracter": (".extracters.window_extracter", "WindowExtracter"),
    "run_window_extracter": (".extracters.window_extracter", "run_window_extracter"),
    "SNPsImputer": (".methods.snps_imputer", "SNPsImputer"),
    "Bpp": (".methods.bpp", "Bpp"),
    "run_admixture_method": (".methods.admixture", "run_admixture_method"),
    "run_bpp_method": (".methods.bpp", "run_bpp_method"),
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
