# SRFM
Welcome to SRFM, a package designed to manage satelitte data and perform retrievals.

## Components
The package contains three external codes - RFM, DISORT and mie_ewp.
The general idea is to run RFM first, mie_ewp after that and DISORT at the end.

### DISORT
The package contains DISORT code as and external addition.
DISORT is a code that calculates radiative transfer (with scattering).

### RFM
The code contains RFM as the reference forward model which calculates the gas absorption in the atmosphere. 
This code should be run before DISORT is run, which can be done either manually or from python (see examples).

### mie_ewp
The mie_ewp module calculates mie scattering on particles.
The main outputs are the extinction coefficient, single scatter albedo, the phas function and its Legendre polynomial expansion.
The outputs are used by DISORT.

## Prerequisites
The required python packages are:
> numpy, matplotlib.pyplot, pandas, pickle, os, sys, pathlib, psutil and meson. 
Another required package is 
> f2py
However, this should be part of numpy if installed correctly, in full and up-to-date.
Next, a fortran compiler if required. This code was tested with the *gfortran* compiler on Linux.

## Installation
The package was built in Python3.13 and was not tested on earlier versions.
The package incorporates RFM, DISORT and mie scattering codes, which are in Fortran.
RFM and mie are in Fortran 90, DISORT is a mixture of both Fortran 77 and 90.

## Usage
First navigate to the folder srfm.
Once there, run the provided bash script to compile required fortran modules by typing
> bash prepare_all.sh

After that, the package can be easily imported as a whole by typing
    `from srfm import *`
    or
    `import srfm`
Modules can also be imported individually.

RFM also needs to be compiled.
This can be done manually or from python after the required inputs are set (recommended).
There is a class method that does that from python and is shown in the example code.

The package comes with an example test script, which you should be able to just run,
and also use as a template to create your own codes.

## Particle scattering optical properties
1. Scattering code uses

- optical_properties.py
- ARIA_module.py
- size_distribution.py
- quadrature.py
- mie_ewp.f90

2. The routines you want to use primarily in size_distribution.py to create a size distribution and optical_properties.py to get the extinction, 
single scatter albedo and Legendre coefficients.  The other modules support these calculations.

The calls are
sd = szd.create_distribution("log_normal", n=x, r=r, s=s)
or
sd = szd.create_distribution("log_normal", surface_area_density=x, r=r, s=s)
or
sd = szd.create_distribution("log_normal", volume_density=x, r=r, s=s)

where
n                     the total concentration in 1/cm^3
r                     median radius in um
s                     spread is typically 1.5 and should be \ge 1 
surface_area_density  surface area of the particles in um^2/cm^3
volume_density        volume of the particles in um^3/cm^3

NOTE: The concentration is set directly using n or indirectly using surface_area_density or volume_density.

3. Once the size distribution object is formed the optical properties can be found through a call to 

e, w, p, l = op.ewp_hs(wavelength, composition, sd, legendre_coefficients_flag=True)

where
wavelength           in um (must be an array)
composition          is one of these strings 'ice' 'ash' 'sulphuric acid'
sd                   size distribution
p                    phase function 
legendre_coefficients_flag True to get Legendre coefficeints

The outputs are
e is the extinction in 1/m
w is the single scatter albedo
p is the phase functions
l are the Legendre coefficients

# Developer instructions
The full documentation and instructions can be found on the links below:

Full documentation (working version) can be found here: https://www.overleaf.com/6287434921cxhkjptrnpvm#4c7cfa

And some underlying science here: https://www.overleaf.com/8669653368vnvvvhdgsyvs#04edb8

