# SRFM memory benchmark

Measured on 2026-09-04 with Python 3.13 using the driver and atmospheric data
from `srfm_main_backup_pre_test_suite/tests`. The matched calculation covered
848–1002 cm-1 at 0.1 cm-1 computational spacing (1,541 points), three particle
layers, and one DISORT output level. The optimized run retained brightness
temperature only and used a 10,000-point scattering block.

| Revision | Peak process-tree RSS | Peak runner RSS | Runtime |
|---|---:|---:|---:|
| Backup implementation | 673.36 MiB | 401.92 MiB | 17.62 s |
| Memory-efficient implementation | 611.70 MiB | 338.85 MiB | 19.89 s |
| Difference | -61.66 MiB (-9.16%) | -63.07 MiB (-15.69%) | +2.27 s (+12.90%) |

The final 601-point brightness-temperature spectra agree to a maximum absolute
difference of 1.49e-9 K and a maximum relative difference of 5.62e-12. Particle
optical depths are identical.

## Production-grid projection

The requested 800–1500 cm-1 range at 0.0005 cm-1 spacing contains 1,400,001
points. Array-size extrapolation using the example's three 200-angle particle
layers, 158 total retained Legendre coefficients, 52 atmospheric layers, and
one output level gives at least 8.59 GiB less live array payload:

| Removed or bounded allocation | Estimated saving |
|---|---:|
| Full-grid phase functions | 6.26 GiB |
| Full-grid Legendre coefficients (10,000-point blocks) | 1.64 GiB |
| Integrated RFM optical-depth matrix | 0.54 GiB |
| Retained DISORT numeric result payload | 0.08 GiB |
| Unrequested SRFM numeric output arrays | 0.07 GiB |

The DISORT history saving excludes Python dictionary and per-array object
overhead, so its realized benefit will be larger than the numeric-payload value.
The compact RFM estimate likewise excludes the legacy DataFrame's one-string-
per-wavenumber column metadata. Consequently 8.59 GiB is a conservative
projection, not a peak-RSS prediction. RFM's own native peak during its worker
process is unchanged.

Reproduce a measurement with:

```bash
python tools/benchmark_memory.py --output memory.json -- \
    python examples/basic_example/run_srfm.py
```

## Production-sized 650–2760 cm-1 run

A completed run over 650–2760 cm-1 at 0.001 cm-1 spacing contained 2,110,001
spectral points. It used the example's three scattering layers, retained BBT
only, and calculated four output levels (10, 15, and 20 km plus TOA).

| Measurement | Result |
|---|---:|
| Peak process-tree RSS | 9.52 GiB |
| Peak runner RSS | 1.46 GiB |
| Peak RFM child RSS | 9.13 GiB |
| End-to-end wall time | 2,347.30 s (39 min 7.30 s) |
| `run_srfm` time | 2,345.77 s (39 min 5.77 s) |

The process-tree peak occurred 957.44 seconds after launch during native RFM.
At that instant the runner used 399.36 MiB and its child processes used
9.13 GiB. The resulting NetCDF grid spans exactly 650–2760 cm-1 and has shape
``(2110001, 1, 4, 1)``.
