"""Downstream analysis methods built on extraction tools."""

from importlib import import_module

__all__ = [
    "Bpp",
    "PCA",
    "PCAFamilyResult",
    "SNPsImputer",
    "run_admixture_method",
    "run_baba_method",
    "run_bpp_method",
    "run_dapc_method",
    "run_pca_analysis",
    "run_pca_method",
    "run_popgen_method",
    "run_snmf_method",
    "run_tsne_analysis",
    "run_umap_analysis",
]


_LAZY_EXPORTS = {
    "Bpp": (".bpp", "Bpp"),
    "SNPsImputer": (".snps_imputer", "SNPsImputer"),
    "run_admixture_method": (".admixture", "run_admixture_method"),
    "run_baba_method": (".baba", "run_baba_method"),
    "run_bpp_method": (".bpp", "run_bpp_method"),
    "run_dapc_method": (".dapc", "run_dapc_method"),
    "run_popgen_method": (".popgen", "run_popgen_method"),
    "run_snmf_method": (".snmf", "run_snmf_method"),
    "PCA": (".pca", "PCA"),
    "PCAFamilyResult": (".pca", "PCAFamilyResult"),
    "run_pca_analysis": (".pca", "run_pca_analysis"),
    "run_pca_method": (".pca", "run_pca_method"),
    "run_tsne_analysis": (".pca", "run_tsne_analysis"),
    "run_umap_analysis": (".pca", "run_umap_analysis"),
}


def __getattr__(name: str):
    """Load method runners lazily so package import stays lightweight."""
    if name in _LAZY_EXPORTS:
        module_name, attr_name = _LAZY_EXPORTS[name]
        module = import_module(module_name, __name__)
        value = getattr(module, attr_name)
        globals()[name] = value
        return value
    raise AttributeError(name)
