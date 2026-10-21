# SRFM test suite

This directory contains the manual automated test suite for SRFM. It protects
scientific calculations, input validation, file formats, and the boundaries to
the compiled RFM, Mie, and DISORT components. The suite does not use the
network, machine-specific paths, interactive plotting, or release automation.

## Installing the test requirements

From a local checkout, install SRFM with its test dependency:

```bash
python -m pip install -e '.[test]'
```

The native tests also require SRFM's compiled extensions. A source checkout can
build them with:

```bash
python build_extensions.py
```

The package currently supports Python 3.10 through 3.13. A Fortran compiler and
the build requirements declared in `pyproject.toml` are needed when compiling
the native extensions. No separately installed RFM executable is required: the
suite exercises the packaged f2py RFM wrapper.

## Running the tests

Run the complete suite from the repository root:

```bash
pytest
```

Select one test level with its registered marker:

```bash
pytest -m unit
pytest -m integration
pytest -m e2e
```

Normal pytest output controls work for the whole suite or any marker. For
example:

```bash
pytest -q
pytest -v
pytest -m integration -v
pytest tests/unit/test_quadrature.py -q
pytest tests/unit/test_quadrature.py::test_shift_quadrature_preserves_integrals_and_reference -v
```

Plain `pytest` includes all three levels. If a compiled extension is absent,
the test that needs it is skipped with a message identifying
`python build_extensions.py` as the remedy. Such a skip is an environmental
result, not evidence that the native calculation passed.

### Reusing compiled code

The two Numba Legendre-expansion functions cache compiled code between Python
processes. Tests default `NUMBA_CACHE_DIR` to `/tmp/srfm-numba`; set that variable
before running pytest to use a different cache directory. The numerical
calculations and comparisons still run on every invocation.

The DISORT validation driver executable is stored under pytest's cache directory
in `d/srfm-disort`. Its cache key includes all nine source files' contents and
paths, the compiler path and version, compiler flags, and platform. A changed
build input triggers recompilation. Every test run executes the driver and
checks all 57 cases and its 47 reference comparisons. Use `pytest --cache-clear`
to force this executable to be rebuilt; this does not clear the Numba cache.

## Test levels

### Unit (`tests/unit/`)

The unit tests are fast, deterministic, and isolate Python-side behaviour.
They cover:

- unit conversions and conversion round trips;
- Gaussian, Radau, Lobatto, and trapezoidal quadrature, including the numerical
  reference values formerly embedded in `quadrature.py`;
- Gaussian and log-normal size distributions, parameter equivalence,
  normalisation, and the reference values formerly in the module's example
  block;
- strict parsing, interpolation, ascending/descending spectra, and range
  handling for tiny synthetic ARIA `.ri` files;
- every current ARIA filename, all legacy filename mappings, full numeric
  fingerprints of the 450 original datasets, and unchanged wavelength and
  wavenumber interpolation (with explicit assertions for the supplied Peterson
  quartz changes);
- RFM output/profile/atmosphere parsing and driver, grid, level, atmosphere,
  and HITRAN cross-section writers;
- the generic, OXHARP, and processed-IASI input schemas, including the unchanged
  public example, runner-specific fields, cross-field constraints, aggregated
  dotted-path errors, and validation before execution side effects;
- `Inputs`, structured RFM helper dataclasses and driver generation;
- layer construction, extent/grid/concentration/mass-loading preparation, and
  size-distribution configuration without invoking Mie;
- optional reference optical-depth inputs in cm-1, um, and nm, density
  requirements, conflict detection, reference index interpolation, and repeated
  layer calculations;
- scalar and spectral albedo, custom beam-normal solar spectra, common
  in-memory/file-backed spectral fields, prescribed Angstrom or tabulated
  optics, analytic HG moments, normalized tabulated moments, and cross-type
  vertical-layer separation;
- Legendre expansion and other deterministic optical-property mathematics;
- DISORT setter, dimension, validation, and model-integrity behaviour without
  invoking the native solver;
- deterministic utility functions and validation; and
- agreement between `srfm.__version__` and installed distribution metadata.

Temporary files are produced with pytest's `tmp_path`; no user's filesystem is
assumed.

Altitude reconciliation regressions cover RFM's single precision relative
tolerance, sequential insertions, grid subdivision, effective particle geometry
and concentration conversion. A native RFM regression checks the 3.000/3.001 km
collision, and all three complete runners compare a near-grid layer against its
explicitly snapped equivalent while checking input preservation and metadata.
Run these focused cases with:

```bash
pytest tests/unit/test_rfm_levels.py tests/integration/test_native_rfm.py -q
pytest tests/e2e/test_srfm_end_to_end.py -k reconciles_collision -q
```

### Integration (`tests/integration/`)

Integration coverage includes:

- size distribution -> quadrature -> optical properties -> compiled Mie;
- synthetic ARIA refractive indices -> compiled Mie optical calculation;
- reproducible random number-concentration layers -> reference optical depth ->
  recovered concentration, mass loading, extinction, albedo, phase functions,
  and Legendre coefficients, including gamma distributions, refractive-index
  arrays, ARIA composition lookup, optional density, and multiprocessing;
- `Inputs` -> `MieLayer` Python-side preparation;
- structured RFM driver -> compiled RFM -> parsed optical-depth output, using
  a three-level atmosphere and three-point F11 cross section;
- every RFM option flag, incompatible flag pair, dependency rule, and malformed
  flag path in the native ``*FLG`` parser, including repeated calls after
  rejected configurations;
- a genuine malformed native RFM run, including stderr log forwarding and
  transient-log removal;
- Python DISORT configuration -> compiled single- and double-precision DISORT;
- analytic Henyey--Greenstein moments -> native ``getmom(iphas=3)`` convention
  oracle, with all eight configurations checked in one isolated subprocess;
- all 17 families and 57 cases from the authoritative DISORT Fortran driver,
  including its 47 embedded reference comparisons;
  and
- each specialized top-level runner -> its own schema before filesystem side
  effects.

The RFM reference optical depths use explicit numerical tolerances because they
are scientific floating-point outputs, not serialized text contracts. Mie and
DISORT calls run in dedicated Python subprocesses so a Fortran `STOP` or native
crash becomes an ordinary pytest failure instead of terminating the suite.

### End to end (`tests/e2e/`)

Self-contained cases follow the public user paths:

```text
Inputs -> main/oxharp_main/iasi_main.run_srfm -> RFM -> Mie -> DISORT
       -> radiance and brightness temperature
```

All use a tiny synthetic F11 atmosphere/cross section and aerosol layer. The
matrix comprises:

- a direct `Inputs` run with no file or plotting output;
- matching radiance and brightness temperature from number-concentration and
  reference-optical-depth Mie inputs without density, including retained
  reference coordinates and recovered concentration in NetCDF metadata;
- a complete basic-example-style `driver_table.py` loaded by
  `Inputs.read_srfm_drv()`;
- an OXHARP-tailored run using precomputed cosine/secant geometry;
- a processed-IASI run using a synthetic one-pixel spectrum, ECMWF profiles,
  observation geometry, NEDT table, and instrument line shape;
- single-precision DISORT with text brightness-temperature/radiance output and
  a non-interactive radiance plot; and
- double-precision DISORT with solar illumination, IASI convolution, NetCDF
  brightness-temperature/radiance output, and a non-interactive brightness-
  temperature plot.
- a prescribed HG particle layer with spectral albedo and an unscaled custom
  FBEAM spectrum, including NetCDF optical provenance.

Together these cover `None`, text, and NetCDF output modes; default and custom
filenames; both plot quantities; solar on/off; IASI convolution on/off; and
single/double DISORT. They verify output contents as well as file existence.
Every complete native pathway runs in a dedicated subprocess. The driver-table
path also verifies that native `rfm.log` content is forwarded to standard output
and its transient file is cleaned up. Unit tests verify that failed-run logs go
to standard error and are removed after emission.

Six additional cases load `examples/basic_example/driver_table.py`, retain its
particle parameters, and reduce the grid and atmosphere to synthetic inputs.
They exercise ash, ice, sulphuric acid, water, quartz, and malic acid. Each case
runs with an archived original ARIA file and its replacement, including legacy
and generic names where applicable. Radiance, brightness temperature, and all
three retained flux fields must agree exactly. This exact equality requirement
checks database substitution within the same environment, without assuming
native outputs are bitwise portable across compilers or platforms.

Each material's comparisons share one outer subprocess to reuse imports and
Mie compilation. Each complete model calculation retains its own RFM worker
subprocess. Named batch progress messages identify the active comparison if a
native call crashes or stops before returning a complete result.

Run the ARIA migration checks with:

```bash
pytest tests/unit/test_aria_module.py tests/unit/test_aria_database.py -q
pytest tests/e2e/test_srfm_end_to_end.py -k aria_replacement -v
```

## Fixtures and scientific reference data

Reusable fixture factories live in `tests/conftest.py` and create their files
inside pytest temporary directories. `tests/fixtures/` documents the unit,
integration, and E2E fixture boundaries. Static data is kept only when it adds
scientific provenance that is not better expressed by a small generator.

The previous top-level `tests/` files duplicated `examples/basic_example/` and
the old driver/output contained private machine paths. Those duplicates and
generated artifacts were removed. The original runnable example and its input
files remain under `examples/basic_example/`. A small, non-executable record of
selected legacy brightness-temperature values is retained in
`tests/fixtures/e2e/legacy_bbt_reference.tsv`; it is not asserted by the new
E2E test because reproducing it requires the external HITRAN/IASI inputs named
by that historical run.

## Adding tests

Put deterministic Python tests in `unit/`, component-boundary tests in
`integration/`, and reserve `e2e/` for the complete application pathway. Mark
every module or test explicitly with `unit`, `integration`, or `e2e`. Prefer
small generated files and physical/mathematical invariants. For floating-point
results use `numpy.testing.assert_allclose` or `pytest.approx` with a tolerance
chosen for the calculation; use exact equality only for discrete contracts.

Do not make routine tests download databases or depend on local absolute paths.
Large native test matrices belong in specialised validation studies rather than
this developer feedback suite.

## Testing driver-table changes

The three runner schemas are centralized in `src/srfm/input_schema.py`, while
the user formats remain their existing dictionaries. Run the focused tests
while changing driver fields or validation:

```bash
pytest tests/unit/test_input_schema.py -q
pytest tests/unit/test_input_schema.py -v
```

These tests load `examples/basic_example/driver_table.py` and assert that
validation retains its complete top-level, RFM, and scattering-layer field
structure. They also check that the nested RFM schema stays synchronized with
the arguments accepted by `run_rfm_with_parameters()`. The E2E driver-table
tests then confirm that a separately generated, self-contained complete table
loads and executes through the real RFM, Mie, and DISORT components.
