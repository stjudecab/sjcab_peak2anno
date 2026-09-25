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
# sjcab_peak2anno_db configuration.
# Uncomment or edit values as needed.
#SJCAB_PEAK2ANNO_DB_PATH=~/.sjcab_peak2anno_db           # Database root
#SJCAB_PEAK2ANNO_SPECIES_VERSIONS=hg38:v31               # Species:version pairs
#SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS=2kb,50kb,2kb          # Promoter-up, enhancer,promoter-down cutoffs
#SJCAB_PEAK2ANNO_GENE_TYPE=nomicro                       # Gene-type filter: nomicro|all|protein_coding|lincRNA
#SJCAB_PEAK2ANNO_ISO_SET=all                             # all or deduplong GeneBEDs
#SJCAB_PEAK2ANNO_2FEATURE_OUT=max                        # Feature output mode: max|percent|max,percent|percent,max
#SJCAB_PEAK2ANNO_2STATE_OUT=max,percent                  # State output mode like Feature output mode
#SJCAB_PEAK2ANNO_TXT_DELIMITER=auto                      # Delimiter inside text regions such as ":-" for chr1:100-200
#SJCAB_PEAK2ANNO_BACKEND=auto                            # default python: python|auto|bedtools
#SJCAB_PEAK2ANNO_AUTO_INSTALL_DB=true                    # Automatically install missing DB files
#SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS_<species>_<version>=2kb,50kb,2kb
```

Environment variables with the same names override RC values. Species/version
cutoff overrides use case-sensitive names such as
`SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS_hg38_v31`.

The backend is selected through the RC file or environment; there is no CLI
backend option. `auto` always uses the indexed Python implementation. Set
`SJCAB_PEAK2ANNO_BACKEND=bedtools` to use native bedtools that must be on
`PATH`, and generated commands are appended to `bedtools-peak2anno.sh` for
review. The script supports mode `1` (one wide window) and mode `2`
(promoter-first, enhancer-second);

`2FEATURE_OUT`,`2STATE_OUT` output modes could be:
- `max`: the assigned feature/state
- `percent`: all feature/states got a percentage and output as table
- `max,percent`: max first, then percent table
- `percent,max`: percent table first, then max

## Gene cutoffs

Use `--prom-enha-cutoffs promoterup,enhancer,promoterdown`. The third value
may be omitted and then reuses promoter-up. `gene0.5` means 50% of gene
length; `transcript2` means twice the transcript length in BED column 5.

## Gene filters

`--gene-type` filters BED column 9 of all.gene.bed. Supported values include:

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



## Automatic database installation

When no explicit `--tss-bed`, `--gene-bed`, `--feature-dir`, or feature list is
provided and the selected reference is missing, the CLI automatically runs the
appropriate `sjcab-peak2anno-db` command for the selected species/version:

The database CLI uses `install-genebed` for gene references and
`install-feature` for feature references.

```bash
peak2anno peak2gene peaks.bed
```

Use `--no-auto-install-db` to disable automatic installation. The explicit
`--auto-install-db` option is retained as an alias for the default behavior.
If the database CLI is not installed, the package is installed from
[PyPI](https://pypi.org/project/sjcab-peak2anno-db/) first. Installation is
locked at `<db-path>/.locks/<species>-<version>.lock`, so concurrent runs for
the same reference share one installation while different references proceed
independently.
