# Commands

## Single-region commands

```bash
peak2anno peak2gene peaks.bed --species hg38
peak2anno narrow2feature peaks.bed --feature-dir annotations/hg38
peak2anno broad2feature peaks.bed --feature-dir annotations/hg38
peak2anno peak2state peaks.bed --states dense_states.bed
```

- `peak2gene` annotates peaks with nearby/closest genes and promoter or
  enhancer relationships.
- `narrow2feature` assigns one priority-ordered feature label.
- `broad2feature` reports overlap percentages for every feature.
- `peak2state` reports chromatin-state overlap.

Default feature files include promoter, 5'UTR, 3'UTR, exon, intron, TES, and
distal/intergenic regions. Feature priority can be selected for
`narrow2feature`, `broad2feature`, `loop2feature`, and corresponding combined
feature steps:

- `-O/--order-lst def|default|none`: use `DB_PATH/order.lst`.
- `-O/--order-lst utr`: use `DB_PATH/order.utr.lst`.
- `-O/--order-lst PATH`: use a custom order list.

```bash
peak2anno narrow2feature peaks.bed --db-path "$HOME/.sjcab_peak2anno_db" \
  --order-lst order.lst -o features.tsv
peak2anno broad2feature peaks.bed --order-lst utr -o features.utr.tsv
```

The default values `def`, `default`, and `none` use `order.lst` under the
database root. `utr` uses `order.utr.lst` under the database root. A custom
path may also be supplied. If the database-root default is absent, a
neighboring `order.lst` in the feature directory is used for compatibility.
Inspect the selected order list to confirm the meaning of each feature; a
two-column list uses column 2 as the output feature name.

All commands show resolved defaults with `-h` and accept positional input or
`-i/--input`. Omit `-o/--output` to write the main table to stdout.

Common options include:

- `-i/--input`: input BED, BEDPE, or text file.
- `-o/--output`: output file; omit it to write the main table to stdout.
- `-f/--output-format`: `auto`, `bed`, `bedpe`, `txt`, or `txtnohead`.
- `-I/--input-format`: `auto`, `bed`, `txt`, or `txtnohead`.
- `-H/--header`: `auto`, `yes`, or `no`.
- `-C/--columns`: zero-based BED coordinate columns.
- `-R/--region-column`: zero-based text-region column.
- `-b/--backend`: `auto`, `bedtools`, `pybedtools`, or `python`.
- `-x/--overlap-cutoff`: minimum overlap threshold.
- `-m/--summary`: summary output path.
- `-p/--plot`: write summary plots.
- `-a/--auto-install-db`: offer to install missing database files.
- `-A/--no-auto-install-db`: disable database installation.

Gene-reference options include:

- `-s/--species`: species key.
- `--ver/--species-version`: annotation version.
- `--iso/--isoform-set`: `all` or `deduplong`.
- `-g/--gene-bed`: explicit gene BED.
- `-t/--tss-bed`: explicit TSS BED.
- `-G/--gene-type`: gene-type filter.

## Combining commands

Run multiple annotations once and merge their columns:

```bash
peak2anno peak2gene narrow2feature peaks.bed \
  --tss-bed tests_data/tss.bed \
  --feature-dir annotations/hg38 \
  --output combined.tsv --workers 2
```

The explicit form is:

```bash
peak2anno combined \
  --commands peak2gene \
  --commands narrow2feature \
  peaks.bed --tss-bed tests_data/tss.bed \
  --feature-dir annotations/hg38 --output combined.tsv
```

## Loop commands

`loop2gene`, `loop2feature`, and `loop2state` accept BEDPE. Columns 0–2 are
the first anchor and columns 3–5 are the second by default. Each anchor is
annotated independently and output columns are prefixed with `anchor1_` and
`anchor2_`.

```bash
peak2anno loop2gene loops.bedpe --tss-bed tests_data/tss.bed -o loops.tsv
peak2anno loop2feature loops.bedpe --feature-dir annotations/hg38 \
  --feature-mode broad -o loops.features.tsv
peak2anno loop2state loops.bedpe --states states.bed \
  --output-format bedpe -o loops.states.bedpe
```

## Gene filters

`--gene-type` filters BED column 9. Supported values include:

- `all`
- `protein_coding`
- `lincRNA`
- `nomicro`
- A comma-separated custom list.

## Run logs

Every run appends its command line and resolved input/reference files to
`.run.log`. Use `--summary` and `--plot` when a feature/state summary is also
needed.
