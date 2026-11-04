# Package Development Guide

This package is a reference layout for creating and maintaining small Python
pip/conda bioinformatics command-line packages.

## Repository layout

- `src/<package_name>/`: importable implementation; keep the Python package
  name identical to the distribution name.
- `src/<package_name>/cli.py`: one parser and dispatcher for all subcommands;
  expose a short CLI alias only through entry points, not as the import name.
- `src/<package_name>/intervals.py`: shared BED, BEDPE, region-text parsing,
  format detection, interval indexing, and table output.
- `src/<package_name>/config.py`: built-in defaults, rc-file loading, and
  environment-variable overrides.
- `src/<package_name>/db.py`: database-root discovery and reference resolution;
  keep database layout logic out of command implementations.
- `src/<package_name>/loopops.py`: indexed paired-loop operations such as
  anchor-level intersection, paired-loop intersection, and merge.
- `src/<package_name>/<command>.py`: command-specific configuration dataclasses
  and annotation functions that can be called from both the CLI and Python API.
- `tests/test_package.py`: standard pytest tests for public APIs, CLI parsing,
  output formats, database discovery, and representative end-to-end commands.
- `tests_data/`: minimal BED, BEDPE, TSS, feature, and state fixtures.
- `docs/`: organized MkDocs pages such as installation, commands,
  input/output, configuration, API, and changelog.
- `mkdocs.yml`: Read the Docs-compatible MkDocs configuration.
- `conda-recipe/meta.yaml`: noarch conda recipe with package entry points and
  metadata links.

## Design conventions

- Use dataclasses for command configuration and keep annotation functions
  deterministic and directly testable.
- Support BED3 as the minimum input. Do not assume BED6 metadata; generic
  headerless extra columns should be named `field4`, `field5`, and so on.
- Auto-detect BED, headered text, and headerless region text. Support a
  configurable region delimiter through an environment variable and rc file.
- Accept one input, comma-separated inputs, or `.lst`/`.list` files. Ignore
  blank/comment lines in list files.
- For multiple inputs, process independent files with `--workers`. Merge only
  for stdout; file output should use one output per input with `.anno` or the
  requested suffix.
- Make output format default to the input format. Support BED/BEDPE, TXT with
  one header, and TXT without a header.
- Keep shared options short and long (`-i/--input`, `-o/--output`,
  `-n/--workers`, etc.) and show resolved defaults in `-h`.
- Resolve references from explicit CLI paths first, then configured database
  paths. Automatic database installation should be opt-out, species-specific,
  and protected by a species/version lock.
- Use a Python interval backend by default. If an optional native tool backend
  is selected, detect the executable, warn when it is unavailable, and write a
  reviewable shell script containing the equivalent commands.
- Keep methods and output-column explanations available through an optional
  command-generated README, without overwriting an existing README.
- Record the command and expanded input/reference paths in `.run.log`.

## Configuration conventions

- Resolve settings in this order: CLI, environment, first existing rc file,
  built-in defaults.
- Support an explicit rc path, then
  `$XDG_CONFIG_HOME/<package>/.<package>.rc`, then `~/.<package>.rc`.
- When the selected rc file does not exist, create it with every supported
  setting as a commented example. Environment-only settings should be added
  to the template if they are supported configuration variables.
- Use package-prefixed variables, for example
  `<PACKAGE>_DB_PATH`, `<PACKAGE>_SPECIES_VERSIONS`, and
  `<PACKAGE>_AUTO_INSTALL_DB`.

## Packaging and release workflow

- Use a `src/` layout and `versioningit` for Git-tag-derived versions. The
  generated version file must be inside the actual package.
- Keep the pip dependency on the database package unpinned when the database
  is independently versioned; conda runtime requirements should match.
- Include Homepage, Repository, Issue, and documentation URLs in both pip and
  conda metadata.
- Keep pip and conda entry points aligned and test all aliases with `--help`.
- Run `ruff`, `git diff --check`, and the complete pytest suite before release.
- Update `docs/changelog.md` and `Log.Codex.Dev.md` before pushing. Move the
  release notes from `Unreleased` into a dated version section.
- Before pushing, squash empty or automatic checkpoint commits, create one
  descriptive release commit, create an annotated version tag, and push both
  the branch and tag.

## Learned User Preferences

- 2026-09-25: Prefer reusable package templates with a `src/` layout, shared
  CLI/config/database modules, standard pytest coverage, MkDocs documentation,
  and synchronized pip/conda release metadata.
- 2026-09-25: Prefer concise CLI options with explicit defaults in help,
  automatic reference discovery/install, reproducible run logs, and changelog
  updates before Git pushes.
