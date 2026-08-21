"""Sequence and SNP extracter tools."""

from importlib import import_module

__all__ = [
    "SeqexEngine",
    "SNPsExtracter",
    "run_seqex",
    "run_snps_extracter",
]


_LAZY_EXPORTS = {
    "SeqexEngine": (".seqex", "SeqexEngine"),
    "run_seqex": (".seqex", "run_seqex"),
    "SNPsExtracter": (".snps_extracter", "SNPsExtracter"),
    "run_snps_extracter": (".snps_extracter", "run_snps_extracter"),
}


def __getattr__(name: str):
    """Load extracter helpers lazily so parser imports stay lightweight."""
    if name in _LAZY_EXPORTS:
        module_name, attr_name = _LAZY_EXPORTS[name]
        module = import_module(module_name, __name__)
        value = getattr(module, attr_name)
        globals()[name] = value
        return value
    raise AttributeError(name)
