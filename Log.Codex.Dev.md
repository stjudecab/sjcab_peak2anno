# Development Log

## 2026-09-25

- Added comma-separated and `.lst`/`.list` multi-input support with parallel
  processing through `--workers`.
- Changed multi-input outputs to use per-input `.anno` or `-o` suffix files;
  `stdout` and `/dev/stdout` merge results with one text header.
- Changed headerless BED text output to assume BED3 and name extra columns
  `field4`, `field5`, and so on.
- Corrected `list-db` to hide compatibility `def` folders and select concrete
  installed versions.
- Added RC support for `SJCAB_PEAK2ANNO_WRITE_README`.

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
- Made missing database installation automatic by default, limited installs to
  the requested species, and updated the database CLI subcommands.
- Corrected `list-db` species/default detection and aligned its output columns;
  switched gene-reference discovery to `DB_PATH/genebed`.
- Normalized text output regions to `chr:start-end`, corrected strand-aware
  TSS and BED half-open distances, and removed gene-type options from
  feature-only commands.
- Fixed feature/state BED loading so feature-only commands do not apply the
  gene-type filter.
- Removed the database-package version constraint from pip, conda, and
  automatic database installation.
