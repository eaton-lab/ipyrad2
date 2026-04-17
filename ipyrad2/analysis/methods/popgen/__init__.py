"""Population-genetic analysis helpers."""

from importlib import import_module

__all__ = ["run_popgen_method"]


_LAZY_EXPORTS = {
    "run_popgen_method": (".runner", "run_popgen_method"),
}


def __getattr__(name: str):
    """Load popgen helpers lazily so parser imports stay lightweight."""
    if name in _LAZY_EXPORTS:
        module_name, attr_name = _LAZY_EXPORTS[name]
        module = import_module(module_name, __name__)
        value = getattr(module, attr_name)
        globals()[name] = value
        return value
    raise AttributeError(name)
