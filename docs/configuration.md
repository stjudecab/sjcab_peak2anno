# Configuration

Settings are resolved in this order:

1. Command-line options
2. Environment variables
3. The first existing RC file
4. Built-in defaults

The RC file is selected from `$SJCAB_PEAK2ANNO_CONFIG`,
`$XDG_CONFIG_HOME/sjcab_peak2anno/.sjcab_peak2anno.rc`, or
`~/.sjcab_peak2anno.rc`. When no file exists, peak2anno creates the selected
path and writes supported settings as commented examples.

```text
# SJCAB_PEAK2ANNO_DB_PATH=~/.sjcab_peak2anno_db
# SJCAB_PEAK2ANNO_SPECIES_VERSIONS=hg38:v31
# SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS=2kb,50kb,2kb
# SJCAB_PEAK2ANNO_GENE_TYPE=all
# SJCAB_PEAK2ANNO_ISO_SET=all
# SJCAB_PEAK2ANNO_2FEATURE_OUT=max
# SJCAB_PEAK2ANNO_2STATE_OUT=max,percent
```

Environment variables with the same names override RC values. Species/version
cutoff overrides use case-sensitive names such as
`SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS_hg38_v31`.

## Main settings

| Setting | Default | Purpose |
| --- | --- | --- |
| `SJCAB_PEAK2ANNO_DB_PATH` | `~/.sjcab_peak2anno_db` | Database root |
| `SJCAB_PEAK2ANNO_SPECIES_VERSIONS` | `hg38:v31` | Species/version pairs |
| `SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS` | `2kb,50kb,2kb` | Promoter-up, enhancer, promoter-down cutoffs |
| `SJCAB_PEAK2ANNO_GENE_TYPE` | `all` | Gene-type filter |
| `SJCAB_PEAK2ANNO_ISO_SET` | `all` | `all` or `deduplong` gene BED |
| `SJCAB_PEAK2ANNO_2FEATURE_OUT` | `max` | Feature output mode |
| `SJCAB_PEAK2ANNO_2STATE_OUT` | `max,percent` | State output mode |

Output modes are `max`, `percent`, `max,percent`, or `percent,max`.

## Gene cutoffs

Use `--prom-enha-cutoffs promoterup,enhancer,promoterdown`. The third value
may be omitted and then reuses promoter-up. `gene0.5` means 50% of gene
length; `transcript2` means twice the transcript length in BED column 5.

## Automatic database installation

When an implicit gene or feature BED is missing, the CLI proposes the
appropriate `sjcab-peak2anno-db` command. Confirm interactively, or use:

```bash
peak2anno peak2gene peaks.bed --auto-install-db
```

Use `--no-auto-install-db` to disable installation. If the database CLI is
not installed, the package is installed from
[PyPI](https://pypi.org/project/sjcab-peak2anno-db/) first. Installation is
locked at `<db-path>/.locks/<species>-<version>.lock`, so concurrent runs for
the same reference share one installation while different references proceed
independently.
