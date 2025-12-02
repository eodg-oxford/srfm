PYTHON ?= python3
DIST_DIR ?= dist
COMPONENTS ?=
PIP_ARGS ?=

ifneq ($(strip $(COMPONENTS)),)
COMPONENT_FLAGS := --components $(COMPONENTS)
endif

ifneq ($(strip $(PIP_ARGS)),)
PIP_FLAGS := --pip-args "$(PIP_ARGS)"
endif

.PHONY: native wheel sdist dist install clean clean-dist bootstrap

native:
	$(PYTHON) tools/build_package.py native $(COMPONENT_FLAGS)

wheel:
	$(PYTHON) tools/build_package.py wheel

sdist:
	$(PYTHON) tools/build_package.py sdist

dist:
	$(PYTHON) tools/build_package.py dist

install:
	$(PYTHON) tools/build_package.py install $(PIP_FLAGS)

clean:
	$(PYTHON) clean_srfm.py

clean-dist:
	$(PYTHON) tools/build_package.py clean-dist

bootstrap:
	$(PYTHON) tools/build_package.py native $(COMPONENT_FLAGS)
