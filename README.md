# sjcab_peak2anno

`sjcab_peak2anno` annotates genomic peak BED files to nearby genes, genomic
contexts, and chromatin states.

Python 3.7 or newer is supported. Python 3.6 is not supported because the
implementation uses standard-library features introduced in Python 3.7.

The package provides both `peak2anno` and `sjcab-peak2anno` command-line
entry points. It requires `sjcab_peak2anno_db==0.1.5` and an external
`bedtools` executable available in `PATH` when installed with pip.

For conda installations, use the St. Jude CAB channel:

```bash
conda install -c stjudecab sjcab_peak2anno sjcab_peak2anno_db=0.1.5 bedtools pybedtools
```

See [docs/peak2anno.md](docs/peak2anno.md) for command examples.

## Install with pip

Install `bedtools` separately and make sure it is on `PATH`:

```bash
bedtools --version
python -m pip install sjcab_peak2anno
```

To install from a source checkout with test support:

```bash
python -m pip install -e ".[test]"
python -m pytest -q
```

## Install with conda

```bash
conda install -c stjudecab sjcab_peak2anno sjcab_peak2anno_db=0.1.5 bedtools pybedtools
```

## Test from a source checkout

The minimal test command is:

```bash
python -m pytest -q
```

If the package is not installed, use the source tree directly:

```bash
PYTHONPATH=src python -m pytest -q
```

## Smoke test after installation

The repository includes the minimal BED fixtures `tests_data/peaks.bed` and
`tests_data/tss.bed`. After installing the package, run:

```bash
bedtools --version
peak2anno peak2gene tests_data/peaks.bed \
  --tss-bed tests_data/tss.bed \
  --promoter-cutoff 100bp \
  --enhancer-cutoff 3kb \
  --output /tmp/peak2anno-smoke.tsv
cat /tmp/peak2anno-smoke.tsv
```

The same command is available through the alternate CLI name:

```bash
sjcab-peak2anno peak2gene tests_data/peaks.bed \
  --tss-bed tests_data/tss.bed \
  --output /tmp/sjcab-peak2anno-smoke.tsv
```

## Minimal files

For a source installation, keep these files and directories:

```text
pyproject.toml
README.md
src/peak2anno/
tests/
tests_data/peaks.bed
tests_data/tss.bed
```

Runtime annotation data comes from the separately installed
`sjcab_peak2anno_db==0.1.5` package or from `--db-path`/`SJCAB_PEAK2ANNO_DB_PATH`.

## Output and run log

`-o/--output` is optional. If it is omitted, the main annotation table is
written to stdout, so it can be piped to another command. Use
`--summary summary.tsv` if a context/state summary is also needed without a
main output file.

Every CLI run appends its command line and resolved input/reference files to
`.run.log` in the current directory. For example:

```bash
peak2anno peak2gene tests_data/peaks.bed \
  --tss-bed tests_data/tss.bed > /tmp/peak2anno.tsv
cat .run.log
```

## Input formats and columns

All single-region commands accept BED-like input or region-text input.

BED input defaults to columns 0, 1, and 2 for chromosome, start, and end:

```text
chr1    100    200    peak_a
```

Use `--columns 2,4,5` when the BED coordinates are in different zero-based
columns. Use `--input-format bed` to force BED parsing.

Region-text input defaults to column 0 and accepts delimiters including `:`,
`-`, `*`, `=`, `/`, `^`, `;`, `_`, `%`, `$`, and `,`:

```text
region
chr1:100-200
chr2^300=450
```

Use `--input-format txt`, `--input-format txtnohead`, `--region-column 1`,
and `--header yes` when needed. Every command accepts either a positional
input or the shorter `-i/--input` form.
`--input-format auto` is the default and detects the form from the first row.

## Output formats

Use `--output-format` with `auto`, `bed`, `txt`, or `txtnohead`:

- `auto` (default): follows the detected input format: BED becomes BED,
  headered region text becomes `txt`, and headerless region text becomes
  `txtnohead`.
- `txt`: tab-delimited output with a header.
- `txtnohead`: tab-delimited output without a header.
- `bed`: BED-style chromosome/start/end output without a header. For loop
  commands, use `bedpe` as the equivalent BEDPE output mode.

Examples:

```bash
peak2anno peak2gene tests_data/peaks.bed --tss-bed tests_data/tss.bed \
  --output result.tsv --output-format txt
peak2anno peak2gene tests_data/peaks.bed --tss-bed tests_data/tss.bed \
  --output result.bed --output-format bed
peak2anno peak2gene tests_data/peaks.bed --tss-bed tests_data/tss.bed \
  --output result.tsv --output-format txtnohead
```

## Combining annotations

Multiple single-region annotations can be run once and merged into one table.
The concise syntax accepts consecutive subcommands:

```bash
peak2anno peak2gene narrow2context peaks.bed \
  --tss-bed tests_data/tss.bed \
  --context-dir annotations/hg38 \
  --output combined.tsv --workers 2
```

The equivalent explicit syntax is:

```bash
peak2anno combined \
  --commands peak2gene \
  --commands narrow2context \
  peaks.bed --tss-bed tests_data/tss.bed \
  --context-dir annotations/hg38 --output combined.tsv
```

The combined implementation reuses the same input and reference settings for
each annotation. `--workers N` runs independent annotation steps in up to N
processes; the default is 1. The final merge preserves the input row order.

## Loop/BEDPE commands

`loop2anno`, `loop2context`, and `loop2state` accept BEDPE. By default, the
first anchor uses columns 0-2 and the second anchor uses columns 3-5:

```text
chr1    100    200    chr1    900    1000    loop_a
```

Use `--loop-columns 0,1,2,6,7,8` for another six-column layout. Each anchor
is annotated independently and the result contains `anchor1_...` and
`anchor2_...` columns.

Examples:

```bash
peak2anno loop2anno loops.bedpe --tss-bed tests_data/tss.bed \
  --output loops.annotated.tsv
peak2anno loop2context loops.bedpe --context-dir annotations/hg38 \
  --output loops.context.tsv --context-mode broad
peak2anno loop2state loops.bedpe --states states.bed \
  --output loops.states.bedpe --output-format bedpe
```

## Gene-type filters

`peak2gene --gene-type` filters on BED column 9:

- `all` (default): include every gene type; column 9 is not required.
- `protein_coding`: include only protein-coding genes.
- `lincRNA`: include only long intergenic non-coding RNA genes.
- `nomicro`: use the built-in St. Jude CAB non-micro gene set, including
  protein-coding, lincRNA, macro lncRNA, processed/transcribed pseudogenes,
  processed transcripts, and bidirectional promoter lncRNAs.
- `type1,type2`: provide a comma-separated custom set of BED column 9 values.

When a filter other than `all` is used, the selected gene BED must contain a
gene type in column 9.

## Context database lookup

`narrow2context` and `broad2context` automatically search the database root
from `--db-path`, `$SJCAB_PEAK2ANNO_DB_PATH`, or
`~/.sjcab_peak2anno_db` when `--context-dir` is omitted. They look for the
standard context files such as `2kb.promoter.up.bed` and `2kb.exon.bed` under
`<db>/<species>/context`, `<db>/<species>/features`, or a versioned species
directory. Use `-c/--context-dir` to override this lookup.

## Missing database references

If the selected `--species`/`--ver` reference is not available, the command
prints the exact `sjcab-peak2anno-db` command it proposes and asks:
`Install database files now? [y/N]:` Type `y` or `yes` to run it and retry the
annotation; any other answer leaves the run unchanged. For the default version,
the proposal is `sjcab-peak2anno-db install-bed`; for an explicit version it is
`sjcab-peak2anno-db download-bed SPECIES VERSION`. Context commands propose
`sjcab-peak2anno-db install gencode-feature`. An explicit `--tss-bed`,
`--gene-bed`, or `--context-dir` is never replaced automatically. In a pipe or
other non-interactive session installation is declined safely.
