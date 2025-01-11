1) Scattering code uses

optical_properties.py
ARIA_module.py
size_distribution.py
quadrature.py
mie_ewp.f90

Compile the fortran  module to python using
f2py -c mie_ewp.f90 -m mie_module

This creates
mie_module.mod
and other files which you don't need to worry about

2) The routines you want to use primarily in size_distribution.py to create a size distribution and optical_properties.py to get the extinction, 
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

3) Once the size distribution object is formed the optical properties can be found through a call to 

e, w, p, l = op.ewp_hs(wavelength, composition, sd, legendre_coefficients_flag=True)

where
wavelength           in um (can be an array)
composition          is one of these strings 'ice' 'ash' 'sulphuric acid'
sd                   size distribution
p                    phase function 
legendre_coefficients_flag True to get Legendre coefficeints

The outputs are
e is the extinction in 1/m
w is the single scatter albedo
p is the phase functions
l are the Legendre coefficients





the particle composition is set using one of these strings

python3 mie_test.py

Notes on Using Fortran with Python

Example

Start with Fortran test programme

mie_test.f90

which calls

mie_ewp.f90

then build mie_test using

gfortran -O2 -ffixed-line-length-132 -fcheck=all mie_ewp.f90 mie_test.f90 -o mie_test

note parameters
O2                     Enables level 2 optimization, 
fixed-line-length-132  Extends the default fixed line length in fixed-format Fortran (like .f files) from 72 to 132 characters.
check=all              Enables runtime checks for various conditions, such as array bounds, uninitialized variables, and type mismatches.




Code for calculating scattering
size_distribution.py : functions descibing the particle size distribution
Use as
import size_distribution as sd 

Test Code
test_size_distibution.py

optical_properties.py



