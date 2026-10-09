# Unit fixtures

Unit-test `.ri`, atmosphere, profile, grid, and output files are intentionally
generated with `tmp_path` by the tests or `tests/conftest.py`. Keeping their
contents beside the assertion makes malformed-format cases easier to review and
prevents accidental dependence on a developer's filesystem.

`aria_legacy_manifest.json` records the 450 ARIA files bundled at commit
`019c2b2404892befaa2fe95b37813e74fd08c6f0`, before the October 2026 replacement.
Each record contains the old path, verified replacement basename, SHA-256 of
the original bytes, FORMAT columns, numeric array shape, and SHA-256 of every
stored numeric value as a C-order little-endian float64 array. Headers and
whitespace do not affect the numeric fingerprint.

For the 386 files readable by the previous parser, the manifest also records
SHA-256 fingerprints of `(n, k)` returned by `RI.select` on seven evenly spaced
points from the minimum to maximum wavelength and wavenumber, serialized in
the same array format. The regression excludes the two changed Peterson quartz
datasets from the interpolation comparison and instead asserts each approved
changed cell and the unchanged remainder of both original numeric arrays.

These fingerprints were captured from the original files, not regenerated from
the replacements. Of the remaining 64 originals, 59 lack n or k; five had
non-UTF-8 headers, which the replacement fixes. Numeric fingerprints cover all
450 files, including those 64.
