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
- RFM output/profile/atmosphere parsing and driver, grid, level, atmosphere,
  and HITRAN cross-section writers;
- the generic, OXHARP, and processed-IASI input schemas, including the unchanged
  public example, runner-specific fields, cross-field constraints, aggregated
  dotted-path errors, and validation before execution side effects;
- `Inputs`, structured RFM helper dataclasses and driver generation;
- layer construction, extent/grid/concentration/mass-loading preparation, and
  size-distribution configuration without invoking Mie;
- Legendre expansion and other deterministic optical-property mathematics;
- DISORT setter, dimension, validation, and model-integrity behaviour without
  invoking the native solver;
- deterministic utility functions and validation; and
- agreement between `srfm.__version__` and installed distribution metadata.

Temporary files are produced with pytest's `tmp_path`; no user's filesystem is
assumed.

### Integration (`tests/integration/`)

Nine deliberately small component-boundary cases cover:

- size distribution -> quadrature -> optical properties -> compiled Mie;
- synthetic ARIA refractive indices -> compiled Mie optical calculation;
- `Inputs` -> `MieLayer` Python-side preparation;
- structured RFM driver -> compiled RFM -> parsed optical-depth output, using
  a three-level atmosphere and three-point F11 cross section;
- a genuine malformed native RFM run, including stderr log forwarding and
  transient-log removal;
- Python DISORT configuration -> compiled single- and double-precision DISORT;
  and
- each specialized top-level runner -> its own schema before filesystem side
  effects.

The RFM reference optical depths use explicit numerical tolerances because they
are scientific floating-point outputs, not serialized text contracts. Mie and
DISORT calls run in dedicated Python subprocesses so a Fortran `STOP` or native
crash becomes an ordinary pytest failure instead of terminating the suite.

### End to end (`tests/e2e/`)

Six self-contained cases follow the public user paths:

```text
Inputs -> main/oxharp_main/iasi_main.run_srfm -> RFM -> Mie -> DISORT
       -> radiance and brightness temperature
```

All use a tiny synthetic F11 atmosphere/cross section and aerosol layer. The
matrix comprises:

- a direct `Inputs` run with no file or plotting output;
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

Together these cover `None`, text, and NetCDF output modes; default and custom
filenames; both plot quantities; solar on/off; IASI convolution on/off; and
single/double DISORT. They verify output contents as well as file existence.
Every complete native pathway runs in a dedicated subprocess. The driver-table
path also verifies that native `rfm.log` content is forwarded to standard output
and its transient file is cleaned up. Unit tests verify that failed-run logs go
to standard error and are removed after emission.

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
