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

- `-O/--order-lst def|default|none` (default: def): use `DB_PATH/order.lst`.
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

- `-i/--input` (default: positional input; required): input BED, BEDPE, or text file.
- `-o/--output` (default: stdout): output file.
- `-f/--output-format` (default: auto): `auto`, `bed`, `bedpe`, `txt`, or `txtnohead`.
- `-I/--input-format` (default: auto): `auto`, `bed`, `txt`, or `txtnohead`.
- `-H/--header` (default: auto): `auto`, `yes`, or `no`.
- `-C/--columns` (default: 0,1,2): zero-based BED coordinate columns.
- `-R/--region-column` (default: 0): zero-based text-region column.
- `-n/--workers` (default: 1): worker processes; loop anchors and combined steps can run in parallel.
- `-x/--overlap-cutoff` (default: 1bp): minimum overlap threshold.
- `-m/--summary` (default: none): summary output path.
- `-p/--plot` (default: false): write summary plots.
- `-a/--auto-install-db` (default: true): automatically install missing database files.
- `-A/--no-auto-install-db` (default: disabled): disable automatic database installation.

The default is also configurable with `SJCAB_PEAK2ANNO_AUTO_INSTALL_DB` in
the rc file or environment (`true`/`false`). CLI flags override that setting.

Gene-reference options include:

- `-s/--species` (default: configured species, usually hg38): species key.
- `--ver/--species-version` (default: configured version, usually v31): annotation version.
- `--iso/--isoform-set` (default: all): `all` or `deduplong`.
- `-g/--gene-bed` (default: database reference): explicit gene BED.
- `-t/--tss-bed` (default: computed from gene BED): explicit TSS BED.
- `-G/--gene-type` (default: nomicro): gene-type filter for `peak2gene` and
  `loop2gene`. `narrow2feature`, `broad2feature`, and `loop2feature` do not
  accept this option.

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
- `nomicro`: protein_coding,processed_transcript,processed_pseudogene,
  transcribed_processed_pseudogene,transcribed_unprocessed_pseudogene,
  translated_unprocessed_pseudogene,transcribed_unitary_pseudogene,
  lincRNA,macro_lncRNA,bidirectional_promoter_lncrna
- A comma-separated custom list from: protein_coding,lncRNA,processed_pseudogene,
  unprocessed_pseudogene,miRNA,snRNA,misc_RNA,TEC,transcribed_unprocessed_pseudogene,
  snoRNA,rRNA_pseudogene,transcribed_processed_pseudogene,IG_V_pseudogene,IG_V_gene,
  transcribed_unitary_pseudogene,TR_V_gene,unitary_pseudogene,TR_J_gene,rRNA,
  polymorphic_pseudogene,IG_D_gene,TR_V_pseudogene,scaRNA,Mt_tRNA,pseudogene,
  IG_J_gene,IG_C_gene,IG_C_pseudogene,ribozyme,TR_C_gene,TR_J_pseudogene,TR_D_gene
  ,IG_J_pseudogene,translated_unprocessed_pseudogene,translated_processed_pseudogene,
  sRNA,Mt_rRNA,vaultRNA,scRNA,IG_pseudogene

## Run logs

Every run appends its command line and resolved input/reference files to
`.run.log`. Use `--summary` and `--plot` when a feature/state summary is also
needed.

The interval backend is configured through `SJCAB_PEAK2ANNO_BACKEND` in the
RC file or environment (`auto`, `bedtools`, or `python`); it is not a CLI
option. When `bedtools` is selected, the exact generated commands are also
appended to `bedtools-peak2anno.sh` for review.

The generated script accepts raw BED input and has both matching strategies:

```bash
./bedtools-peak2anno.sh 1 one-window.matches.tsv  # one wide bedtools window
./bedtools-peak2anno.sh 2 two-stage.matches.tsv  # promoter, then enhancer
```

It uses `bedtools`, `awk`, and `sort` and emits a package-compatible
peak-to-gene table. Mode 2 is the lower-memory strategy; mode 1 uses one
bedtools window and can be faster when process startup dominates.
