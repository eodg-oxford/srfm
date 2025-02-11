#!/bin/sh
#This code compiles mie_module and DISORT, which can then be imported as modules within
# srfm. Always run this before working with srfm

# compile mie_module
f2py -c mie_ewp.f90 -m mie_module --f90flags='-O3 -g -fcheck=all -fdump-core -fbounds-check -Wall'

# compile disort_module
cd ./DISORT
cat DISOTESTAUX.f DISORT.f BDREF.f DISOBRDF.f ERRPACK.f LINPAK.f LAPACK.f RDI1MACH.f > code.f
f2py -c code.f -m disort_module --backend meson --f90flags='-O3 -g -fcheck=all -fdump-core -fbounds-check -Wall'
