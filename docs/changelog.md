# Changelog

## 2026-09-04

- Renamed the distribution to `sjcab_peak2anno`.
- Added the `sjcab-peak2anno` CLI alias and pinned `sjcab_peak2anno_db` to 0.1.5.
- Added `bedtools` and `pybedtools` to the conda runtime dependencies.
- Lowered the declared Python minimum from 3.9 to 3.7.

## 2026-07-14

- Added the `peak2anno` Python package scaffold with `peak2anno` and `sjcab_peak2anno` console entry points.
- Added `peak2gene`, `narrow2context`, `broad2context`, `peak2state`, and `list-db` subcommands.
- Added a noarch conda recipe depending on `sjcab_peak2anno_db`.
- Added package documentation and offline unit tests for toy BED/database inputs.
