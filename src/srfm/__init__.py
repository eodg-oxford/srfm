version = "0.0.1"
print(f"Welcome to srfm v.{version}")

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
]

"""Import srfm modules."""
from . import units
from . import readers
from . import plotting
from . import utilities
from . import disort_functions
from . import forward_model
from . import rfm_functions
from . import DISORT
from . import DISORT_dbl
from . import ARIA_module
from . import optical_properties
from . import quadrature
from . import size_distribution
from . import mie_module
from . import layer
from . import orography
from . import iasi_main
