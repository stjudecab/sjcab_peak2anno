# peak2anno package

`sjcab_peak2anno` is a small Python command-line package for annotating BED-like peak files.
It exposes the same CLI through `peak2anno` and `sjcab-peak2anno`.

The pip package depends on `sjcab_peak2anno_db==0.1.8`. It also requires the
`bedtools` executable to be installed separately and available in `PATH`.
Python 3.7 or newer is supported; Python 3.6 is not.

## Database

The runtime database root is resolved in this order:

1. `--db-path`
2. `$SJCAB_PEAK2ANNO_DB_PATH`
3. `~/.sjcab_peak2anno_db`

Install the database package from the St. Jude CAB channel when building an environment:

```bash
conda install -c stjudecab sjcab_peak2anno sjcab_peak2anno_db=0.1.8 bedtools pybedtools
```

The local DB manifest currently marks these defaults:

- `hg38`: `v31`
- `mm10`: `vM22`
- `mm9`: `vM17`
- `hg19`: `v31lift37`
- `sacCer3`: `R64-1-1`

For `peak2gene`, the database is searched under the selected species and
version for `all.gene.bed` by default, or `deduplong.gene.bed` with
`--iso deduplong`. The default configuration selects `hg38:v31`; use
`--ver def` to select the database manifest default, or `--ver` to select a
specific version. `--species-version` and `--isoform-version` remain accepted
as longer aliases. Explicit `--tss-bed` or `--gene-bed` overrides discovery.

Gene annotation uses `--prom-enha-cutoffs promoterup,enhancer,promoterdown`.
The third value may be omitted and then reuses promoter-up. `gene0.5` means
50% of gene length, and `transcript2` means twice BED column 5 transcript
length. Species/version-specific environment keys use case-sensitive names such
as `SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS_hg38_v31`.

Defaults may be placed in an rc file selected by `$SJCAB_PEAK2ANNO_CONFIG`,
`$XDG_CONFIG_HOME/sjcab_peak2anno/.sjcab_peak2anno.rc`, or
`~/.sjcab_peak2anno.rc`. Supported keys are
`SJCAB_PEAK2ANNO_DB_PATH`, `SJCAB_PEAK2ANNO_SPECIES_VERSIONS`,
`SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS`, `SJCAB_PEAK2ANNO_GENE_TYPE`, and
`SJCAB_PEAK2ANNO_ISO_SET`, `SJCAB_PEAK2ANNO_2FEATURE_OUT`, and
`SJCAB_PEAK2ANNO_2STATE_OUT`. Environment values override rc values.

`SJCAB_PEAK2ANNO_2FEATURE_OUT` and `SJCAB_PEAK2ANNO_2STATE_OUT` accept
`max`, `percent`, `max,percent`, or `percent,max`. `max` emits one assignment;
`percent` emits ordered percentage columns. Feature percentages allocate each
overlapping base to the first feature in `order.lst`; output columns follow that
same order. State max uses the largest overlap, with `order.lst` as the tie-break
order when one is present.

## Commands

```bash
peak2anno peak2gene peaks.bed -o peaks.peak2gene.tsv --species hg38
peak2anno narrow2feature peaks.bed -o peaks.feature.tsv --species hg38 --feature-dir annotations/hg38
peak2anno broad2feature peaks.bed -o peaks.feature_fractions.tsv --species hg38 --feature-dir annotations/hg38
peak2anno peak2state peaks.bed --states dense_states.bed -o peaks.states.tsv
```

`peak2gene` accepts `--species-version default` (or `def`), `--isoform-version all|deduplong`,
`--prom-enha-cutoffs 2kb,50kb,2kb` and `--gene-type`. Promoter distances are
relative to transcription direction: for a `+` gene, upstream is genomic left,
and for a `-` gene, upstream is genomic right.
For `--gene-type` values other than `all`, the gene/TSS BED must have BED column 9
with gene biotype values. `nomicro` expands to the St. Jude CAB non-micro gene-type set.

`narrow2feature` assigns one feature label by the legacy priority order:
promoter up, promoter down, exon, intron, TES, distal 5', distal 3', intergenic.

`broad2feature` does not collapse by priority. It reports overlap base pairs and
fractions for every feature, plus the primary feature by largest passing fraction.

`peak2state` applies the same broad overlap calculation to a ChromHMM or Segway
dense/segments BED. Use `--state2name` for a two-column state ID to display-name map.

All feature/state commands support `--overlap-cutoff`, where values like `10bp`
mean base pairs and values like `0.1` or `10%` mean a fraction of each input peak.
Use `--summary` to choose the summary TSV path and `--plot` to write PNG/PDF plots
when `matplotlib` is installed.
