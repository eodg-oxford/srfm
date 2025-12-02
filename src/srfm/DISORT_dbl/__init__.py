"""DISORT double-precision extension."""

from importlib import import_module
import warnings

__all__ = ["disort_module_d"]

try:
    disort_module_d = import_module(".disort_module_d", __name__)
except ModuleNotFoundError:
    disort_module_d = None
    warnings.warn(
        "The DISORT double-precision extension is not available. "
        "Run 'python build_extensions.py disort_double' to build it.",
        stacklevel=2,
    )
