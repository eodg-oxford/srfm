# End-to-end fixtures

The active E2E atmospheres, cross section, driver tables, small IASI-like
instrument line shape, processed one-pixel IASI observation, and NEDT table are
generated in pytest temporary directories. The TSV in this directory records a
few values from the removed legacy
`tests/results/bbt.nc` artifact solely as scientific provenance. It is not an
active regression fixture because the old run depended on external HITRAN and
IASI files.
