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

All commands show resolved defaults with `-h` and accept positional input or
`-i/--input`. Omit `-o/--output` to write the main table to stdout.

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

`--gene-type` filters BED column 9. Supported shortcuts include `all`,
`protein_coding`, `lincRNA`, `nomicro`, or a comma-separated custom list.

## Run logs

Every run appends its command line and resolved input/reference files to
`.run.log`. Use `--summary` and `--plot` when a feature/state summary is also
needed.
