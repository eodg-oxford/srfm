# End-to-end fixtures

The active E2E atmospheres, cross section, driver tables, small IASI-like
instrument line shape, processed one-pixel IASI observation, and NEDT table are
generated in pytest temporary directories. The TSV in this directory records a
few values from the removed legacy
`tests/results/bbt.nc` artifact solely as scientific provenance. It is not an
active regression fixture because the old run depended on external HITRAN and
IASI files.

`aria_legacy_samples.zip` contains byte-for-byte copies of six original ARIA
files from commit `019c2b2404892befaa2fe95b37813e74fd08c6f0`:

- `eyjafjallajokull-ash_Reed.ri`
- `ICE_Warren_2008.ri`
- `H2SO4_75_Palmer_1975.ri`
- `H2O_263K_Rowe_2020.ri`
- `quartz100_Henning_1997.ri`
- `malic-acid_Laskina_2014.ri`

Their source headers and numerical data are preserved, and their byte hashes
are checked against `../unit/aria_legacy_manifest.json`. Pytest extracts the
archive to a temporary directory. The ARIA E2E tests run the old files through
the same reader, interpolation, RFM, Mie, and DISORT pipeline as their
replacements, using small synthetic variants of the public example. They
require exact equality of the returned spectra and fluxes. The original
database is not installed alongside the replacement.
