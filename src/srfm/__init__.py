"""Top-level SRFM package."""

from importlib import import_module
import warnings

version = "0.0.1"

__all__ = [
    "units",
    "readers",
    "plotting",
    "utilities",
    "disort_functions",
    "forward_model",
    "rfm_functions",
    "DISORT",
    "DISORT_dbl",
    "ARIA_module",
    "optical_properties",
    "size_distribution",
    "quadrature",
    "mie_module",
    "layer",
    "orography",
    "iasi_main",
    "rfm_helper",
    "RFM",
    "main",
    "inputs",
    "oxharp_iasi_main",
]


def _safe_import(name: str, optional: bool = False):
    try:
        module = import_module(f".{name}", __name__)
        globals()[name] = module
    except ModuleNotFoundError:
        if optional:
            warnings.warn(
                f"Optional extension '{name}' is unavailable. "
                "Run 'python build_extensions.py' to compile the native modules.",
                stacklevel=2,
            )
            globals()[name] = None
        else:
            raise


_OPTIONAL_MODULES = {"mie_module"}

for module_name in __all__:
    _safe_import(module_name, optional=module_name in _OPTIONAL_MODULES)
