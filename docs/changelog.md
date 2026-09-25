# Changelog

## Unreleased

- Changed multi-input file output to write one result per input using an
  `.anno` default or the `-o` suffix; `stdout` and `/dev/stdout` retain merged
  output. Headerless BED text output now assumes BED3 and names extra columns
  `field4`, `field5`, and so on.

- Added comma-separated and `.lst`/`.list` multi-input support; independent
  input files can be processed in parallel with `--workers`.

- Updated `list-db` to hide compatibility `def` directories when a concrete
  installed version is available and mark selected gene/feature entries as
  default.

- Added `--write-readme`/`--readme` to write command-specific methods and
  output-column notes, and simplified `.run.log` to expanded command-only
  entries.
- Added `SJCAB_PEAK2ANNO_WRITE_README` to generated rc files; it can be set in
  the RC file or environment, with the environment taking precedence. An
  existing README is not changed.
- Simplified `list-db` to list gene BED files and feature folders with columns
  `species`, `version`, `annotation`, `default`, and `path`.
- Replaced setuptools-scm with versioningit for automatic Git-based package
  versioning in pip and conda builds.
- Removed the `sjcab_peak2anno_db` version constraint from pip, conda, and
  automatic database installation.
- Added `SJCAB_PEAK2ANNO_AUTO_INSTALL_DB` rc/environment configuration for
  the automatic database-installation default.

## 0.2.0 / 2026-09-24

- Updated automatic database installation to use the current
  `install-genebed` and `install-feature` subcommands.
- Restricted automatic installation to the requested species and improved
  `list-db` discovery for installed current-layout BED files.
- Updated gene-reference discovery to use `DB_PATH/genebed`, corrected
  species/default reporting, and aligned `list-db` output into columns.
- Corrected closest-gene distances to use standard BED half-open gap lengths;
  for example, a gap from peak end 36,936,544 to TSS start 36,974,543 is
  reported as 37,999 bp.
- Corrected minus-strand TSS conversion to use the BED end coordinate, so
  strand-aware closest-gene distances match the GENCODE/voom convention.
- Removed the gene-type option from feature-only commands.
- Fixed feature/state BED loading to bypass gene-type filtering entirely.
- Enabled automatic installation of missing implicit database references by
  default; use `--no-auto-install-db` to opt out.
- Changed the default `SJCAB_PEAK2ANNO_GENE_TYPE` and all related API/CLI
  defaults from `all` to `nomicro`.
- Reworked the interval backend and documentation for reproducible Python and
  optional bedtools execution, including reviewable one-stage and two-stage
  bedtools scripts and matching peak-to-gene results.
- Improved indexed peak-to-gene queries, deterministic sorting and nearest-gene
  tie handling, feature order-list support, worker options, and regression
  coverage.
- Renamed the Python package and source directory to `sjcab_peak2anno`; the
  shorter `peak2anno` name is retained only as a CLI alias.
- Added `docs/api.md` with Python API examples using canonical
  `sjcab_peak2anno` imports.
- Added selectable interval backends through the
  `SJCAB_PEAK2ANNO_BACKEND` configuration setting. `auto` prefers `bedtools`
  when available and records the resolved backend in `.run.log`.
- Implemented batched single-window `bedtools` queries for fixed-cutoff
  peak-to-gene annotation while retaining exact Python-compatible results.
- Removed the `pybedtools` backend and CLI backend options. The backend remains
  configurable through `SJCAB_PEAK2ANNO_BACKEND`; selected `bedtools` commands
  are recorded in `bedtools-peak2anno.sh` for review.
- Changed the default `auto` backend to Python and removed bedtools from pip
  and conda runtime requirements. Bedtools is now an explicit optional
  environment-selected backend.
- Aligned promoter-first and single-wide-window cutoff boundaries so Python
  and bedtools approaches produce identical peak-to-gene assignments.
- Generated a standalone raw-BED bedtools review script that emits a
  package-compatible table with selectable one-window and two-stage modes.
- Added indexed interval-window and nearest-gene queries to reduce repeated
  full-reference scans.
- Made `peak2gene` output voom-compatible: rows are chromosome/start sorted,
  grouped names and IDs are lexical, promoter and distal assignments are
  mutually exclusive, and closest-gene ties follow reference order.
- Added `--order-lst`/`-O` for feature priority ordering, including database
  defaults and the `utr` shorthand for `order.utr.lst`.
- Made default feature discovery independent of numeric cutoff prefixes and
  added two-column order-list feature-name mappings.
- Added 5'UTR and 3'UTR to the default feature set.
- Added default-valued `-n/--workers`; loop anchors can now annotate in
  parallel and combined steps use the same worker setting.
- Added short command-line aliases and documented the legacy `voom2anno.sh`
  equivalence for peak-to-gene runs.
- Added regression coverage for backend resolution, distance compatibility,
  and command-line parsing.

## 0.1.8 / 2026-09-11

- Updated the package and conda dependency to `sjcab_peak2anno_db` 0.1.8.

## 2026-09-10

- Added separate strand-aware upstream/downstream promoter cutoffs for gene annotation.
- Added rc-file and environment configuration for database, species versions, and gene cutoffs.
- Added rc-file and environment defaults for gene type and isoform set.
- Added configurable max/percent output modes for feature and state annotations.
- Renamed `loop2anno` to `loop2gene` and replaced legacy gene cutoff options with `--prom-enha-cutoffs`.

## 2026-09-04

- Renamed the distribution to `sjcab_peak2anno`.
- Added the `sjcab-peak2anno` CLI alias and pinned `sjcab_peak2anno_db` to 0.1.5.
- Previously added bedtools and pybedtools to the conda runtime dependencies;
  these are no longer required.
- Lowered the declared Python minimum from 3.9 to 3.7.

## 2026-07-14

- Added the `peak2anno` Python package scaffold with `peak2anno` and `sjcab_peak2anno` console entry points.
- Added `peak2gene`, `narrow2feature`, `broad2feature`, `peak2state`, and `list-db` subcommands.
- Added a noarch conda recipe depending on `sjcab_peak2anno_db`.
- Added package documentation and offline unit tests for toy BED/database inputs.
