PYTHON ?= python3
BUILD_SCRIPT := $(abspath $(CURDIR)/../../build_extensions.py)

all: mie disort_single disort_double rfm_py

.PHONY: all mie disort_single disort_double clean rfm_py

mie:
	$(PYTHON) $(BUILD_SCRIPT) mie

disort_single:
	$(PYTHON) $(BUILD_SCRIPT) disort_single

disort_double:
	$(PYTHON) $(BUILD_SCRIPT) disort_double

rfm_py:
	$(PYTHON) $(BUILD_SCRIPT) rfm

clean:
	$(PYTHON) $(abspath $(CURDIR)/../../clean_srfm.py)
