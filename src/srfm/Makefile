#This Makefile compiles mie_module and DISORT and rfm_py, which can then be imported as modules within
# srfm. Always run this before working with srfm

# for performance compile
FFLAGS = -O3 -g -Wall

#for debug compile
#FFLAGS = -O3 -g -fcheck=all -fdump-core -fbounds-check -Wall

LOGFILE = build.log

all: mie disort_single disort_double rfm_py

.PHONY: all mie disort_single disort_double clean rfm_py

mie:
	@echo "Building mie_module..." | tee -a $(LOGFILE)
	@set -e; \
	if [ ! -f mie_ewp.f90 ]; then \
		echo "Error: mie_ewp.f90 not found." | tee -a $(LOGFILE); exit 1; \
	fi; \
	f2py -c mie_ewp.f90 -m mie_module --f90flags='$(FFLAGS)' >> $(LOGFILE) 2>&1

disort_single:
	@echo "Building disort_module_s (single precision)..." | tee -a $(LOGFILE)
	@set -e; cd DISORT; \
	for f in DISOTESTAUX.f DISORT.f BDREF.f DISOBRDF.f ERRPACK.f LINPAK.f LAPACK.f RDI1MACH.f; do \
		if [ ! -f $$f ]; then echo "Error: $$f not found in DISORT." | tee -a ../$(LOGFILE); exit 1; fi; \
	done; \
	cat DISOTESTAUX.f DISORT.f BDREF.f DISOBRDF.f ERRPACK.f LINPAK.f LAPACK.f RDI1MACH.f > code.f; \
	f2py -c code.f -m disort_module_s --backend meson --f90flags='$(FFLAGS)' >> ../$(LOGFILE) 2>&1

disort_double:
	@echo "Building disort_module_d (double precision)..." | tee -a $(LOGFILE)
	@set -e; cd DISORT_dbl; \
	for f in DISOTESTAUX.f DISORT.f BDREF.f DISOBRDF.f ERRPACK.f LINPACK_D.f LAPACK.f RDI1MACH.f; do \
		if [ ! -f $$f ]; then echo "Error: $$f not found in DISORT_dbl." | tee -a ../$(LOGFILE); exit 1; fi; \
	done; \
	cat DISOTESTAUX.f DISORT.f BDREF.f DISOBRDF.f ERRPACK.f LINPACK_D.f LAPACK.f RDI1MACH.f > code.f; \
	f2py -c code.f -m disort_module_d --backend meson --f90flags='$(FFLAGS)' >> ../$(LOGFILE) 2>&1

rfm_py:
	@echo "Building rfm_py..." | tee -a $(LOGFILE)
	@set -e; cd RFM; \
	python3 -m numpy.f2py -c -m rfm_py rfm_wrapper.pyf $$(find source -maxdepth 1 -name '*.f90' ! -name 'rfm.f90' ! -name 'combined.f90') >> ../$(LOGFILE) 2>&1

clean:
	@echo "Cleaning build artifacts..." | tee -a $(LOGFILE)
	@rm -f DISORT/code.f DISORT_dbl/code.f RFM/source/combined.f90 *.so
	@rm -rf build __pycache__
	@rm -f $(LOGFILE)
	@echo "Clean complete."
