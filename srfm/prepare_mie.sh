#!/bin/sh
#This bash code compiles mie_ewp.f90 for f2py
f2py -c mie_ewp.f90 -m mie_module --f90flags='-O3 -g -fcheck=all -fdump-core -fbounds-check -Wall'


cat DISOTESTAUX.f DISORT.f BDREF.f DISOBRDF.f ERRPACK.f LINPAK.f LAPACK.f RDI1MACH.f > code.f
f2py -c code.f -m disort_module --backend meson --f90flags='-O3 -g -fcheck=all -fdump-core -fbounds-check -Wall'
