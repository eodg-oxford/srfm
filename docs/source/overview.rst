Overview
========

Welcome to SRFM, a package designed to manage satellite data and perform retrievals. It combines Python tooling with compiled Fortran extensions for radiative transfer.

Components
----------
This is a Python package that ships three external Fortran codes compiled via ``numpy.f2py``:

- **DISORT** — Radiative transfer solver (single and double precision) with minor tweaks to printing routines. See the `DISORT site <http://www.rtatmocn.com/disort/>`_.
- **RFM** — Reference Forward Model for gas absorption with small adjustments (6-digit filenames, Python wrappers). `RFM home <https://eodg.atm.ox.ac.uk/RFM/index.html>`_.
- **mie_ewp** — Mie scattering on particles; outputs extinction, single-scatter albedo, phase function, and Legendre expansion.

Installation
------------
``pip install srfm`` builds the Fortran extensions from source. Only sdists are published (no prebuilt wheels).

If installation fails, download the sdist and rebuild locally. Fortran compilation depends on ``numpy`` plus ``meson`` and ``ninja`` (installed automatically when missing).

Supported platforms and prerequisites
-------------------------------------
- CPython 3.10–3.13 on Linux and Windows.
- A Fortran compiler (``gfortran`` works cross-platform), Meson, and Ninja. Helper scripts bootstrap Meson/Ninja via ``pip`` if absent.

Building and packaging
----------------------
Linux / macOS quick start::

   make native          # compile the Fortran extensions
   make dist            # rebuild wheel + sdist under dist/
   make install         # install from the freshly built artifacts

Or run::

   python build_extensions.py      # builds Fortran extensions

Windows (PowerShell / cmd)::

   py tools/build_package.py native
   py tools/build_package.py dist
   py tools/build_package.py install

Optional arguments let you rebuild only certain components, e.g. ``COMPONENTS="mie disort_double"`` for ``make`` or ``--components mie disort_double`` for the Python helper. Compilation logs are stored in ``src/srfm/build.log`` and ``python clean_srfm.py`` removes compiled artifacts.

The ``.f2py_f2cmap`` file in ``src/srfm/RFM`` must remain in place to enforce correct types during f2py; deleting it will cause wrong precision and allocation errors.

Installing from dist/
---------------------
After building locally, install from ``dist/`` with::

   python -m pip install --no-index --find-links dist SRFM

The ``tools/install_from_dist.py`` wrapper runs the same command and bootstraps Meson/Ninja first.

Verification
------------
If installation succeeded, ``import srfm`` works. Failures often mean the Fortran extensions were not compiled correctly.

Usage
-----
To run the packaged example as a traditional radiative transfer model:

1. Install SRFM.
2. Obtain external dependencies (line parameters, cross sections) as needed.
3. Adjust paths in the driver table (output locations, external data, etc.).
4. Run ``run_srfm.py``.
5. Results are written to the paths you configured.

You can also import modules directly::

   from srfm import *
   import srfm

External dependencies
---------------------
Python runtime deps (installed automatically): ``numpy``, ``pandas``, ``scipy``, ``netCDF4``, ``matplotlib``, ``cartopy``, ``numba``, ``psutil``, ``mergedeep``.

Fortran build deps: ``gfortran``, ``meson``, ``ninja`` (Meson/Ninja auto-installed by helpers).

Domain data requirements:

- **Spectroscopy** — HITRAN line parameters (download separately and convert to binary with `hitbin <https://eodg.atm.ox.ac.uk/RFM/hitbin.html>`_). Cross sections may be used directly; specify paths in the driver table.
- **Refractive indices** — ARIA database is bundled (`ARIA <https://eodg.atm.ox.ac.uk/ARIA/>`_). Alternatives must keep the same format; specify in the driver table.

Documentation
-------------
Sphinx-generated docs are built from function docstrings. Use the navigation sidebar for API details.

References
----------
The package publication link will be added when available.
