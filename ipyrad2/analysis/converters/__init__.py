"""Data-conversion tools for the analysis workflow."""

from importlib import import_module

__all__ = ["VCFToHDF5", "run_vcf_to_hdf5"]


_LAZY_EXPORTS = {
    "VCFToHDF5": (".vcf_to_hdf5", "VCFToHDF5"),
    "run_vcf_to_hdf5": (".vcf_to_hdf5", "run_vcf_to_hdf5"),
}


def __getattr__(name: str):
    """Load converter helpers lazily so parser imports stay lightweight."""
    if name in _LAZY_EXPORTS:
        module_name, attr_name = _LAZY_EXPORTS[name]
        module = import_module(module_name, __name__)
        value = getattr(module, attr_name)
        globals()[name] = value
        return value
    raise AttributeError(name)
