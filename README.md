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

## Prerequisites
The required python packages are:
> numpy, matplotlib.pyplot, pandas, pickle, os, sys, pathlib, psutil, numba, meson, datetime, time, warnings, multiprocessing, scipy, abc, bisect.

Another required package is 
> f2py

However, this should be part of numpy if installed correctly, in full and up-to-date.
Next, a fortran compiler if required.
This code was tested with the *gfortran* compiler on Linux.
The code was also observed to _NOT_ work with ifort.

## Installation
The package was built in Python3.13 and was not tested on earlier versions.
The package incorporates RFM, DISORT and Mie scattering codes, which are in Fortran.
RFM and Mie are in Fortran 90, DISORT is a mixture of both Fortran 77 and 90.

## Usage
First navigate to the folder srfm.
> cd srfm

Once there, run the provided Makefile by typing
> make

To clean the build, type
> make clean

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
the RFM is now compiled as a python module. This is done through the Makefile. 
To succeed, the compiler needs the ".f2py_f2cmap" file in src/srfm/RFM. Do not delete that!!!
That file ensures that correct types are enforced during the compilation.
For example, it forces things such as real(kind=r8)/double precision to emit C doubles during the f2py run.
If the file is deleted, the C emits single precision values and underallocates arrays, leading to malloc() errors.

