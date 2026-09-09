"""
Run SRFM with a driver table
============================

This example calculates brightness-temperature spectra for an atmosphere with
sulphuric-acid aerosol, ash, and a water cloud. It demonstrates the traditional
SRFM workflow in which a small runner loads a separate, editable driver table.

External spectroscopy data
---------------------------

SRFM does not distribute HITRAN line data or RFM cross sections. Download the
`complete basic example bundle <../../_downloads/basic_example.zip>`_, extract
it, and set ``SRFM_HITRAN_FILE`` and ``SRFM_XSC_DIR`` as described in its
``README.rst``. No machine-specific paths are stored in this example.

The driver table
----------------

The complete configuration is kept in a normal Python file so it is easy to
adapt. In particular, review the spectral grids, atmosphere, gases, scattering
layers, and output geometry before running a scientific calculation.

.. literalinclude:: driver_table.py
   :language: python
   :caption: driver_table.py

Example output
--------------

The figure below is a precomputed result from this configuration. The hosted
documentation does not rerun the model because its external spectroscopy data
is not available on Read the Docs.

.. image:: example_output.png
   :alt: Example SRFM brightness-temperature spectrum
   :width: 100%
"""

import os
from pathlib import Path

from srfm import inputs, main

# Use the script location rather than the caller's working directory so that
# the downloaded example can be run from anywhere.
EXAMPLE_DIR = Path(__file__).resolve().parent

required_variables = ("SRFM_HITRAN_FILE", "SRFM_XSC_DIR")
missing_variables = [name for name in required_variables if not os.environ.get(name)]
if missing_variables:
    names = ", ".join(missing_variables)
    raise RuntimeError(f"Set {names} before running this example; see README.rst.")

srfm_inputs = inputs.Inputs()
srfm_inputs.read_srfm_drv(EXAMPLE_DIR / "driver_table.py")

result = main.run_srfm(srfm_inputs)

# Use ``result`` for further analysis in the same Python process. Files and
# plots requested by the driver table are written beneath ``results/``.

# sphinx_gallery_thumbnail_path = '../../examples/basic_example/example_output.png'
