"""Top-level SRFM package."""

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version as _distribution_version
import warnings

try:
    __version__ = _distribution_version("SRFM")
except PackageNotFoundError:
    # Source trees are not necessarily installed. Distribution metadata remains
    # authoritative; this sentinel cannot silently masquerade as a release version.
    __version__ = "0+unknown"

# Backwards-compatible public alias.
version = __version__

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
    "input_schema",
    "inputs",
    "oxharp_main",
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
