# TODO
This README is outdated and need updating (22 Oct 2025, AK).

# SRFM
Welcome to SRFM, a package designed to manage satelitte data and perform retrievals.

## Components
The package contains three external codes - RFM, DISORT and mie_ewp.
The general idea is to run RFM first, mie_ewp after that and DISORT at the end.
Other models can be used, of course. 
This code is designed to be as general-purpose as possible.

### DISORT
The package contains DISORT code as and external addition.
DISORT is a code that calculates radiative transfer (with scattering).
This code is written if Fortran.
Two versions - single and double precision - are present.

### RFM
The code contains the Reference Forward Model (RFM) which calculates gas absorption in the atmosphere.

### mie_ewp
The mie_ewp module calculates Mie scattering on particles.
The main outputs are the extinction coefficient, single scatter albedo, the phas function and its Legendre polynomial expansion.
The outputs are used by DISORT.

## Supported platforms and prerequisites
SRFM now targets every currently supported CPython release (3.10–3.13) on
Linux and Windows. The Fortran components require a Fortran compiler
(`gfortran` works well cross-platform), plus [Meson](https://mesonbuild.com)
and [Ninja](https://ninja-build.org). The helper scripts install Meson and
Ninja automatically through `pip` if they are missing, and the build metadata
in `pyproject.toml` lists them as build requirements so `pip install` also
bootstraps them in isolated build environments.

The runtime Python dependencies are declared in the `pyproject.toml`. Building
the Fortran extensions requires `numpy` (for `numpy.f2py`), `meson`, and
`ninja`. You only need to install those manually when working outside of the
supplied helpers.

## Installation
The package was built in Python 3.13 originally, but is now actively tested on
Python 3.12 and 3.13 and is expected to work on 3.10+.
The package incorporates RFM, DISORT and Mie scattering codes, which are in Fortran.
RFM and Mie are in Fortran 90, DISORT is a mixture of both Fortran 77 and 90.

## Building and packaging (multi-platform)
The repository root now contains both a GNU Makefile (for Linux/macOS) and a
Python driver (`tools/build_package.py`) that offer the same targets so that
Windows users do not have to rely on `make`. Every workflow below automatically
installs Meson/Ninja if they are missing.

### Linux / macOS quick start
```
make native          # compile the Fortran extensions
make dist            # rebuild wheel + sdist under dist/
make install         # install from the freshly built artifacts
```

### Windows (PowerShell or Command Prompt)
```
py tools/build_package.py native
py tools/build_package.py dist
py tools/build_package.py install
```

Both entry points accept optional arguments such as
`COMPONENTS="mie disort_double"` for `make` or
`--components mie disort_double` for the Python helper when you want to rebuild
only a subset of extensions. Compilation logs continue to be stored in
`src/srfm/build.log`. `python clean_srfm.py` removes every compiled artifact.

### Installing from dist/
When you run `pip` with the local `dist/` folder as an extra index, it
automatically inspects the platform/Python tags exposed by each wheel and
falls back to the source distribution when there is no suitable wheel.
```
python -m pip install --no-index --find-links dist SRFM
```
The `tools/install_from_dist.py` wrapper runs the same command and performs the
Meson/Ninja bootstrap step first. This works unchanged on Linux and Windows.

The package comes with an example test script, which you should be able to just run,
and also use as a template to create your own codes.
At this stage, users are recommended to develop their own script similar to the test script.

If required, the package can be easily imported as a whole by typing
    `from srfm import *`
    or
    `import srfm`
Modules can also be imported individually.

RFM also needs to be compiled.
This can be done manually or from python after the required inputs are set (recommended).
There is a class method that does that from python and is shown in the example code (srfm.forward\_model.RFM.compile\_rfm()).

# Developer instructions
The full documentation and instructions can be found on the links below:

Full documentation (working version) can be found here: https://www.overleaf.com/6287434921cxhkjptrnpvm#4c7cfa

And some underlying science here: https://www.overleaf.com/8669653368vnvvvhdgsyvs#04edb8

## Developer notes
The RFM bindings are compiled as a python module via `build_extensions.py`
(or the Makefile, which delegates to the same helper). 
To succeed, the compiler needs the ".f2py_f2cmap" file in src/srfm/RFM. Do not delete that!!!
That file ensures that correct types are enforced during the compilation.
For example, it forces things such as real(kind=r8)/double precision to emit C doubles during the f2py run.
If the file is deleted, the C emits single precision values and underallocates arrays, leading to malloc() errors.
