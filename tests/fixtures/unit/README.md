# Unit fixtures

Unit-test `.ri`, atmosphere, profile, grid, and output files are intentionally
generated with `tmp_path` by the tests or `tests/conftest.py`. Keeping their
contents beside the assertion makes malformed-format cases easier to review and
prevents accidental dependence on a developer's filesystem.
