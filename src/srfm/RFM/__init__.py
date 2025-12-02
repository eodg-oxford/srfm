"""Python bindings for the RFM code."""

from importlib import import_module
import warnings

__all__ = ["rfm_py"]

try:
    rfm_py = import_module(".rfm_py", __name__)
except ModuleNotFoundError:
    rfm_py = None
    warnings.warn(
        "The RFM extension is not available. "
        "Run 'python build_extensions.py rfm' to build it.",
        stacklevel=2,
    )
