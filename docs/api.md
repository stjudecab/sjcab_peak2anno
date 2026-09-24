# Python API

The canonical Python package name is `sjcab_peak2anno`. The shorter
`peak2anno` name is provided only for the command-line executable.

## Peak-to-gene annotation

```python
from pathlib import Path

from sjcab_peak2anno.peak2gene import PeakGeneConfig, annotate_peak2gene

config = PeakGeneConfig(
    input_path=Path("peaks.bed"),
    output_path=Path("peaks.genes.tsv"),
    species="hg38",
    species_version="v31",
    prom_enha_cutoffs="2kb,50kb,2kb",
    db_path="~/.sjcab_peak2anno_db",
    output_format="txt",
    backend="auto",
)
annotate_peak2gene(config)
```

Use `gene_bed` or `tss_bed` to provide a reference explicitly. Do not set
both. `output_path=None` writes the table to stdout.

## Feature and state annotation

```python
from pathlib import Path

from sjcab_peak2anno.features import (
    FeatureConfig,
    StateConfig,
    annotate_broad_feature,
    annotate_narrow_feature,
    annotate_peak_state,
)

feature_config = FeatureConfig(
    input_path=Path("peaks.bed"),
    output_path=Path("peaks.features.tsv"),
    species="hg38",
    db_path="~/.sjcab_peak2anno_db",
    order_lst="def",       # def/default/none, utr, or a file path
    output_format="txt",
    output_mode="max,percent",
    backend="auto",
)
annotate_narrow_feature(feature_config)

# Use the same FeatureConfig with annotate_broad_feature for per-feature
# overlap output.
annotate_broad_feature(feature_config)

state_config = StateConfig(
    input_path=Path("peaks.bed"),
    states_path=Path("states.bed"),
    output_path=Path("peaks.states.tsv"),
    state2name=Path("state2name.tsv"),
    output_format="txt",
    output_mode="max,percent",
    backend="auto",
)
annotate_peak_state(state_config)
```

Annotation functions return the configured output path (or `None` when
writing to stdout); feature and state functions additionally return their
summary path.

## Loop annotation

`loop2gene`, `loop2feature`, and `loop2state` are available through the
`annotate_loop` function. Each BEDPE anchor is annotated independently.

```python
from pathlib import Path

from sjcab_peak2anno.loops import annotate_loop

annotate_loop(
    "loop2gene",
    Path("loops.bedpe"),
    Path("loops.genes.tsv"),
    "txt",
    species="hg38",
    species_version="v31",
    isoform_version="all",
    db_path="~/.sjcab_peak2anno_db",
    backend="auto",
)
```

## Configuration and backend

The CLI loads RC-file and environment settings with `load_settings()`. API
configuration objects accept explicit values and do not parse command-line
arguments.

```python
from sjcab_peak2anno.config import load_settings

settings = load_settings()
print(settings.default_species, settings.version_for(settings.default_species))
```

The `backend` value may be `auto`, `bedtools`, or `python`.
`auto` prefers the native `bedtools` executable when available.

Feature BED files are discovered by semantic suffixes such as
`.promoter.up.bed`, `.5utr.bed`, and `.exon.bed`, so their numeric cutoff
prefix need not be `2kb`. Check the selected `order.lst` or `order.utr.lst` to confirm feature
meanings. In a two-column order list, column 2 is used as the output feature
name; unmatched entries retain the standard feature names.
