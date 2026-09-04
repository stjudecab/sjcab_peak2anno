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
