"""DISORT single-precision extension."""

from importlib import import_module
import warnings

__all__ = ["disort_module_s"]

try:
    disort_module_s = import_module(".disort_module_s", __name__)
except ModuleNotFoundError:
    disort_module_s = None
    warnings.warn(
        "The DISORT single-precision extension is not available. "
        "Run 'python build_extensions.py disort_single' to build it.",
        stacklevel=2,
    )
