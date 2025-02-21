#!/bin/sh
#This shell code compiles DISORT for f2py
cat DISOTESTAUX.f DISORT.f BDREF.f DISOBRDF.f ERRPACK.f LINPACK_D.f LAPACK.f RDI1MACH.f > code.f
f2py -c code.f -m disort_module --backend meson --f90flags='-O3 -g -fcheck=all -fdump-core -fbounds-check -Wall'

