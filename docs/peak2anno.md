# peak2anno package

`peak2anno` is a small Python command-line package for annotating BED-like peak files.
It exposes the same CLI through both `peak2anno` and `sjcab_peak2anno`.

## Database

The runtime database root is resolved in this order:

1. `--db-path`
2. `$SJCAB_PEAK2ANNO_DB_PATH`
3. `~/.sjcab_peak2anno_db`

Install the database package from the St. Jude CAB channel when building an environment:

```bash
conda install -c stjudecab sjcab_peak2anno_db peak2anno
```

The local DB manifest currently marks these defaults:

- `hg38`: `v31`
- `mm10`: `vM22`
- `mm9`: `vM17`
- `hg19`: `v31lift37`
- `sacCer3`: `R64-1-1`

## Commands

```bash
peak2anno peak2gene peaks.bed -o peaks.peak2gene.tsv --species hg38
peak2anno narrow2context peaks.bed -o peaks.context.tsv --species hg38 --context-dir annotations/hg38
peak2anno broad2context peaks.bed -o peaks.context_fractions.tsv --species hg38 --context-dir annotations/hg38
peak2anno peak2state peaks.bed --states dense_states.bed -o peaks.states.tsv
```

`peak2gene` accepts `--species-version default` (or `def`), `--isoform-version all|deduplong`,
`--promoter-cutoff 2kb`, `--enhancer-cutoff 50kb`, and `--gene-type`.
For `--gene-type` values other than `all`, the gene/TSS BED must have BED column 9
with gene biotype values. `nomicro` expands to the St. Jude CAB non-micro gene-type set.

`narrow2context` assigns one feature label by the legacy priority order:
promoter up, promoter down, exon, intron, TES, distal 5', distal 3', intergenic.

`broad2context` does not collapse by priority. It reports overlap base pairs and
fractions for every feature, plus the primary feature by largest passing fraction.

`peak2state` applies the same broad overlap calculation to a ChromHMM or Segway
dense/segments BED. Use `--state2name` for a two-column state ID to display-name map.

All context/state commands support `--overlap-cutoff`, where values like `10bp`
mean base pairs and values like `0.1` or `10%` mean a fraction of each input peak.
Use `--summary` to choose the summary TSV path and `--plot` to write PNG/PDF plots
when `matplotlib` is installed.
