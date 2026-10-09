Basic SRFM example
==================

This directory contains a complete SRFM driver-table example. The atmospheric
profiles and optional IASI files are included, but the spectroscopy inputs must
be supplied separately because they are not distributed with SRFM.

Requirements
------------

Install SRFM and obtain:

* a HITRAN line database converted to the binary RFM format with ``hitbin``;
* a directory containing any required RFM cross-section (``.xsc``) files.

Run
---

From this directory, run::

   python run_srfm.py

The output is written beneath ``results/``. Edit ``driver_table.py`` to change
the grids, output geometry, atmospheric inputs, scattering layers, or viewing
angles. The table also demonstrates a coarse spectral Lambertian albedo, a
custom beam-normal solar spectrum that is passed to DISORT without amplitude or
date scaling, and a prescribed Angstrom/Henyey--Greenstein layer separated from
the Mie and grey-body layers.

Refractive indices
------------------

The ``comp`` field in each scattering layer selects bundled ARIA data.
``"sulphuric acid"`` selects ``H2SO4_75%_300K_Palmer_1975.ri`` and ``"ash"``
selects ``eyjafjallajokull_ash_58.5%SiO2_Reed_2018.ri``. The water cloud uses
``H2O_263K_Rowe_2020.ri``. Existing generic names and legacy filenames remain
valid with the updated database, so these example selections are unchanged.
For an ice calculation, use ``"ice"`` or ``"ice_266K_Warren_2008.ri"`` and
adjust the particle and layer parameters for the intended scene.
