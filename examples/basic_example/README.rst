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

Configure external data
-----------------------

Set two environment variables before running the example.

Linux and macOS::

   export SRFM_HITRAN_FILE=/path/to/hitran.bin
   export SRFM_XSC_DIR=/path/to/xsc

Windows PowerShell::

   $env:SRFM_HITRAN_FILE = "C:\path\to\hitran.bin"
   $env:SRFM_XSC_DIR = "C:\path\to\xsc"

Run
---

From this directory, run::

   python run_srfm.py

The output is written beneath ``results/``. Edit ``driver_table.py`` to change
the grids, output geometry, atmospheric inputs, scattering layers, or viewing
angles.
