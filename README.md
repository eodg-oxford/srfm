# SRFM
Welcome to SRFM, a package designed to manage satelitte data and perform retrievals.

## Components
This is a Python package.
Apart from pure python routines, the package contains three external 
codes - RFM, DISORT and mie_ewp.
These codes are in Fortran and are compiled as python modules with numpy-f2py.

### DISORT
The package contains the [DISORT](http://www.rtatmocn.com/disort/) code as 
an external module.
DISORT is a code that solves the radiative transfer equation (with scattering).
The original code, which is in single precision and mostly fixed-form Fortran,
has been updated to double precision.
Both versions - single and double precision - are present.
Slight modifications to DISORT were performed in its printing routines.

### RFM
The code contains the Reference Forward Model 
([RFM](https://eodg.atm.ox.ac.uk/RFM/index.html)) which calculates gas 
absorption in the atmosphere.
Small changes with respect to the original RFM were made,
such as 6-digit filenames for better precision and a few scripts for 
compilation into Python modules and capturing the output.

### mie_ewp
The mie_ewp module calculates Mie scattering on particles.
The main outputs are the extinction coefficient, single scatter albedo, 
the scattering phase function and its Legendre polynomial expansion for a given
particle size distribution and composition in a given wavenumber/wavelength
range.

## Installation
Installation through pip should compile the Fortran extenstions with numpy on
the user's system.
To install, `pip install srfm` should work.
Only source (sdist) is provided. Platform specific wheels are currently not 
shipped.
If installation fails, the second option is to download the sdist and manually
recompile the Fortran routines and debug - see below.

### Supported platforms and prerequisites
SRFM targets every currently supported CPython release (3.10–3.13) on
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

### Building and packaging (multi-platform)
The repository root contains both a GNU Makefile (for Linux/macOS) and a
Python driver (`tools/build_package.py`) that offer the same targets so that
Windows users do not have to rely on `make`. Every workflow below automatically
installs Meson/Ninja if they are missing.

### Linux / macOS quick start
```
make native          # compile the Fortran extensions
make dist            # rebuild wheel + sdist under dist/
make install         # install from the freshly built artifacts
```

or

```
python build_extensions.py      # builds Fortran extensions 
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

To succeed, the compiler needs the ".f2py_f2cmap" file in src/srfm/RFM. 
Do not delete that!!!
That file ensures that correct types are enforced during the compilation.
For example, it forces things such as real(kind=r8)/double precision to emit 
C doubles during the f2py run.
If the file is deleted, the C emits single precision values and underallocates
arrays, leading to malloc() errors.

### Installing from dist/
If installing manually, first generate the sdist or wheel on your machine 
(as described above) and then install from the local `dist/` directory by
`pip install sdist-or-wheel-filename`
When you run `pip` with the local `dist/` folder as an extra index, it
automatically inspects the platform/Python tags exposed by each wheel and
falls back to the source distribution when there is no suitable wheel.
Alternatively, run 
```
python -m pip install --no-index --find-links dist SRFM
```
The `tools/install_from_dist.py` wrapper runs the same command and performs the
Meson/Ninja bootstrap step first. This works unchanged on Linux and Windows.

## Verifcation
If the srfm was successfully installed, you should be able to import it in 
Python (`import srfm`).
If this import fails, then so has your installation.
A common error is not having correctly compiled the Fortran extensions, as a 
result of which they can't be imported. 
A hopefully helpful error message will be printed.

## Usage
To use this code as a traditional radiative transfer model program, the package
has with an example.
The example can be found at:
**INSERT LINK**

To run this example, you need to:
1. Successfully install srfm (see above).
2. If required (probably yes), obtain the external dependencies as described 
below.
3. Modify any paths in the driver table (i.e. where to save results, where your
external dependencies are on your computer, etc.
Note that you need line parameters or cross-sections for all gases listed in 
the gases section of the RFM inputs.
4. Run the program run_srfm.py
5. Your results will be stored in the folder you specified in the driver table.

Users are recommended to use this example as a template to develop their own.

Alternatively, the package can be easily imported by
    `from srfm import *`
    or
    `import srfm`
Modules can also be imported individually.
In this case, the package works and should therefore be treated as a standalone
python package with its modules and functions, which the user may use as 
required.

## External dependencies
All current package dependencies are standard packages. 
Installation requires `gfortran, meson and ninja`, the latter two being 
installed automatically at installation.
Usage requires packages listed in pyproject.toml:
```
numpy, pandas, scipy, netCDF4, matplotlib, cartopy, numba, psutil and 
mergedeep
```
These should be installed automatically by pip at installation.

Besides those, there are several external dependencies.
The prime among those is line transition parameters for the line-by-line 
radiative transfer calculation.
The default option is to use line by line parameters from 
[HITRAN](https://hitran.org/).
The transition parameters need to be downloaded separately and converted to 
a binary format that the code can efficiently use.
A code to convert hitran database data to the binary file can be downloaded 
[here](https://eodg.atm.ox.ac.uk/RFM/hitbin.html).
For molecules not in the HITRAN line parameters database, absorption
cross sections may be used.
These can also be downloaded from the HITRAN website (and do not need to be 
converted).
The location of both the line parameters file and the cross-sections folder
are specified in the driver table when running the program as paths (see 
example).
Other sources of line parameters are possible, but in that case, care should be
taken to convert them to the appropriate format.

Next, the scattering calculations require complex particle refractive indices.
The [ARIA](https://eodg.atm.ox.ac.uk/ARIA/) database of refractive indices is 
shipped with the code and is the default option (no need to do anything by the 
user).
Alternative sources are permitted if they keep the same format.
Refractive indices to use are again specified in the driver table 
(see exmaple).

## Documentation
The full documentation can be found at:
**INSERT LINK**
This is a Sphinx-generated documentation from function docstrings, i.e. 
equivalent information may be found in the docstrings within the code.

## References
The package has been published at:
**INSERT LINK**
