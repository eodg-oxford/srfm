Reference optical depth
=======================

``compare_layers.py`` generates a reproducible random particle layer from number
concentration, then reconstructs its optical properties using the optical depth
at one wavelength. It compares concentration, mass loading, extinction, optical
depth, single-scattering albedo, phase functions, and Legendre coefficients.

Run from an installed SRFM environment with the compiled Mie extension::

   python examples/optical_depth_input/compare_layers.py

For a source checkout, use ``PYTHONPATH=src`` with that command. The script prints
its physical inputs and numerical comparison errors and writes no files. No
HITRAN data, RFM calculation, or DISORT calculation is required.
