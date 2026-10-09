Input files
===========

SRFM's driver table is executable Python and may be replaced by an equivalent
in-memory dictionary. Any relative path in that dictionary is resolved from the
process's current working directory, not from the directory containing a driver
table.

Driver and atmospheric files
----------------------------

``Inputs.read_srfm_drv`` reads a Python file defining one top-level ``inputs``
dictionary. OXHARP driver files instead define ``STATE`` and ``ANCILLARY``
dictionaries, which are recursively merged. These files are imported as Python,
so values may be scalars, containers, NumPy arrays, paths, or constructed RFM
helper objects. An equivalent dictionary can be supplied directly without a
file.

RFM atmosphere, spectroscopy, cross-section, and optional driver files are named
in ``driver_inputs`` and ``rfm_config``. Atmosphere files use the RFM ATM format:
comments begin with ``!``, the first data item is the number of profile levels,
each ``*NAME [units]`` record is followed by that many values, and ``*END``
terminates the file. HITRAN input is the binary line database expected by RFM;
cross-sections use RFM's XSC format and may be selected with path wildcards.
Generated or existing RFM driver files use its ``*HDR``, ``*FLG``, and related
section records. The complete native definitions are in the `RFM documentation
<https://eodg.atm.ox.ac.uk/RFM/index.html>`_. The basic example contains usable
atmosphere samples.

.. _aria-refractive-indices:

ARIA refractive indices
-----------------------

SRFM bundles ARIA beneath ``srfm/data/ARIA``. Set each Mie layer's ``comp`` to
an exact ``.ri`` basename; subdirectory names are not part of the identifier.
The October 2026 replacement contains 677 datasets, including all 450 datasets
from the previous bundle. All 371 renamed basenames remain accepted so existing
driver tables continue to select the same samples. Current names include more
explicit temperatures, concentrations, relative humidities, and sample labels.

The generic compositions retain their original datasets:

.. list-table:: Generic composition selection
   :header-rows: 1
   :widths: 15 40 45

   * - ``comp``
     - Previous filename (still accepted)
     - Current filename
   * - ``ash``
     - ``eyjafjallajokull-ash_Reed.ri``
     - ``eyjafjallajokull_ash_58.5%SiO2_Reed_2018.ri``
   * - ``ice``
     - ``ICE_Warren_2008.ri``
     - ``ice_266K_Warren_2008.ri``
   * - ``sulphuric acid``
     - ``H2SO4_75_Palmer_1975.ri``
     - ``H2SO4_75%_300K_Palmer_1975.ri``

For water, specify a dataset such as ``H2O_263K_Rowe_2020.ri``, as the basic
example does. There is no generic ``water`` selection. Filename matching is
case-sensitive; punctuation such as ``%`` is part of the filename.

To list current filenames or resolve a legacy name::

   from importlib.resources import files
   from pathlib import Path
   from srfm.ARIA_module import find_ri_files, get_ri_filepathname

   filenames = sorted(
       Path(path).name for path in find_ri_files(files("srfm.data") / "ARIA")
   )
   acid_path = get_ri_filepathname("H2SO4_75_Palmer_1975.ri")

``ARIA_module.RI`` reads the file's ``FORMAT`` header and interpolates linearly
on the requested wavelength or wavenumber grid. The scattering calculation
uses wavelengths in micrometres and converts ARIA's extinction coefficient to
the Mie ``n - ik`` convention. Both ``n`` and ``k`` must be present, and the
dataset must cover the requested scattering grid. As before, archive entries
containing only one index component cannot be used directly for Mie scattering.
To read an external file directly, use ``ARIA_module.read_ri_file(path)``.

Numerical compatibility
~~~~~~~~~~~~~~~~~~~~~~~

All stored numerical columns are identical for 448 of the 450 previous
datasets. The supplied replacement changes the following values in two
Peterson quartz files; SRFM retains these new values. Calculations using these
datasets can therefore change where interpolation includes the listed samples.

.. list-table:: Numerical changes in the supplied quartz datasets
   :header-rows: 1

   * - Filename
     - Wavelength (micrometres)
     - Component
     - Old value
     - New value
   * - ``quartz_E_Peterson_1969.ri``
     - 17.8
     - k
     - 0.07604
     - 0.07904
   * - ``quartz_O_Peterson_1969.ri``
     - 9.0
     - n
     - 2.59701
     - 0.17463
   * - ``quartz_O_Peterson_1969.ri``
     - 9.4
     - n
     - 1.37131
     - 6.38644
   * - ``quartz_O_Peterson_1969.ri``
     - 9.4
     - k
     - 4.87845
     - 1.37131
   * - ``quartz_O_Peterson_1969.ri``
     - 22.3
     - k
     - 3.03182
     - 3.03132
   * - ``quartz_O_Peterson_1969.ri``
     - 25.4
     - n
     - 5.99961
     - 5.99862

Regression tests compare all original numerical data with the replacements,
check interpolation in both spectral units for previously readable unchanged
datasets, and compare six small example-derived RFM/Mie/DISORT runs against
original ARIA files using exact equality for radiance, brightness temperature,
and fluxes.

Instrument and IASI files
-------------------------

The optional ``iasi_ils`` field and the IASI runner's ``ils`` field name an RFM
ILS text file. After three comment/header lines, its fourth line contains the
number of points, first spectral offset, and offset spacing. The following
whitespace-separated values are the line-shape samples. See the example
``iasi.ils`` file and the `RFM ILS specification
<https://eodg.atm.ox.ac.uk/RFM/sum/ilsfil.html>`_.

The IASI runner's ``nedt`` file is whitespace-delimited text with three header
rows followed by ``wavenumber_cm-1 NEDT_K`` pairs. Its ``iasi_spc_fldr`` and
``iasi_fl`` fields select a Python pickle containing per-pixel arrays. Required
keys are ``spec_bbt``, ``zen``, ``azi``, ``sza``, ``saa``, ``ecmwf_z``,
``ecmwf_T``, ``ecmwf_p``, ``ecmwf_O3``, ``ecmwf_CO``, ``ecmwf_N2O``,
``ecmwf_CO2``, ``ecmwf_CH4``, and ``ecmwf_H2O``. The leading dimension indexes
pixels; ``spec_bbt`` uses the fixed 645--2760 cm\ :sup:`-1` IASI grid at
0.25 cm\ :sup:`-1`, and every ECMWF field contains one vertical profile per
pixel. Because pickle can execute code while loading, use only trusted files.

Common spectral text format
---------------------------

Albedo, custom solar irradiance, prescribed optical depth, prescribed
single-scattering albedo, HG asymmetry, and normalized Legendre moments can be
read from whitespace-delimited text files. Blank lines and lines beginning with
``#`` are accepted by NumPy's text reader. ``skiprows`` counts physical header
rows. Column numbers are zero-based.

A scalar spectral field uses::

   {
       "file": "field.txt",
       "grid_column": 0,
       "value_column": 1,
       "skiprows": 0,
       "grid_units": "cm-1",
   }

``grid_column`` defaults to 0, ``value_column`` defaults to 1, and ``skiprows``
defaults to 0. The grid and value columns must differ. The corresponding
in-memory representation replaces the file and column keys with ``grid`` and
``values`` arrays.

Legendre moments have several value columns::

   {
       "file": "moments.txt",
       "grid_column": 0,
       "value_columns": [1, 2, 3, 4],
       "skiprows": 1,
       "grid_units": "cm-1",
   }

Each row is ``coordinate beta_0 beta_1 ... beta_n``. ``value_columns`` is
required and its order defines moment order zero through ``n``. In memory,
``values`` has shape ``(spectral_points, moments)``.

All grids must contain at least two finite, positive, strictly monotonic points.
Accepted ``grid_units`` spellings are exactly ``"cm-1"``, ``"um"``, and
``"nm"``. SRFM converts the source coordinate to cm :sup:`-1`, orders it by
increasing wavenumber, and performs linear interpolation in wavenumber. The
source must bracket the complete computational grid; endpoint clamping and
extrapolation are not allowed.

Surface albedo files
--------------------

The top-level ``albedo`` field accepts the common scalar-field file mapping.
Every value must be finite and in the inclusive interval [0, 1]. At each
computational wavenumber the interpolated value is passed to DISORT's Lambertian
``ALBEDO`` input. The historical scalar form remains supported and uses the
one-time setter fast path.

Custom solar-spectrum files
---------------------------

The top-level ``solar_spectrum`` mapping uses the common scalar-field format and
adds a required ``value_units`` key. Values are beam-normal DISORT ``FBEAM``
spectral irradiance and must be finite and non-negative. Accepted spellings and
conversion to ``W m-2 (cm-1)-1`` are:

* ``"W m-2 (cm-1)-1"``: unchanged;
* ``"W m-2 um-1"``: multiply by :math:`\lambda_{um}^2/10^4`;
* ``"W m-2 nm-1"``: multiply by :math:`\lambda_{um}^2/10`;
* ``"mW cm-2 um-1"``: multiply by :math:`\lambda_{um}^2/10^3`.

Conversion is performed once at source points before interpolation. A custom
spectrum is never normalized or scaled by date, target distance, solar zenith
angle, or its cosine. Solar zenith controls only DISORT ``UMU0``. When ``sun``
is false, ``FBEAM`` is zero even if a custom file was provided. Omitting
``solar_spectrum`` preserves the Gueymard-2018 interpolation and date correction.

Prescribed optical-property files
---------------------------------

For tabulated column optical depth, put the common scalar-field keys beside its
type::

   "optical_depth": {
       "type": "tabulated",
       "file": "optical_depth.txt",
       "grid_column": 0,
       "value_column": 1,
       "skiprows": 0,
       "grid_units": "cm-1",
   }

Optical depth must be finite and non-negative. ``ssalb`` uses the common scalar
file mapping without a ``type`` key and must lie in [0, 1]. A finite scalar SSA
is also accepted and broadcast. A bare array is rejected because it has no
interpolation coordinate.

Spectral HG asymmetry is nested under ``asymmetry``::

   "phase_function": {
       "type": "henyey_greenstein",
       "asymmetry": {
           "file": "asymmetry.txt",
           "grid_column": 0,
           "value_column": 1,
           "skiprows": 0,
           "grid_units": "cm-1",
       },
   }

The asymmetry factor must lie in [-1, 1]. A scalar asymmetry factor is also
supported and broadcast. SRFM evaluates normalized HG moments analytically as
:math:`\beta_l=g^l` for orders zero through ``nmom``.

Tabulated phase moments use the matrix file form directly inside
``phase_function`` and require ``"type": "legendre_moments"`` plus
``"convention": "normalised"``. Moment zero must equal one within
``1e-5`` wherever ``optical_depth * ssalb`` is positive. An absorption-only
layer may omit ``phase_function``; SRFM supplies harmless moments internally.
