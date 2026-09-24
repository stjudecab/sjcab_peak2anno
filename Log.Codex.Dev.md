# Development Log

## 2026-09-24

- Standardized the package/import name as `sjcab_peak2anno` and retained the
  shorter name only for CLI compatibility.
- Added Python and optional explicit bedtools backends, removed pybedtools and
  bedtools runtime dependencies, and made Python the automatic default.
- Added reproducible bedtools review scripts with one-stage and two-stage
  modes, plus matching output and boundary behavior across backends.
- Improved peak-to-gene interval indexing, sorting, nearest-gene tie handling,
  feature order-list support, worker options, and related tests.
- Updated README, API, install, configuration, command, and changelog
  documentation.
