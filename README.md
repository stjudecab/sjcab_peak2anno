# sjcab_peak2anno

`sjcab_peak2anno` annotates genomic peaks to genes, genomic features, and
chromatin states. It provides the `peak2anno` and `sjcab-peak2anno` command
line interfaces.

![peak2anno subcommand overview](docs/peak2anno-subcommands.png)

## Install

Pip does not require `bedtools` or `pybedtools`; both are detected when
available. Without them, a slower Python interval fallback is used.

```bash
python -m pip install sjcab_peak2anno
```

The conda package requires both tools:

```bash
conda install -c stjudecab sjcab_peak2anno sjcab_peak2anno_db=0.1.8 bedtools pybedtools
```

## Documentation

Read the full documentation at
[sjcab-peak2anno.readthedocs.io](https://sjcab-peak2anno.readthedocs.io/).

- [Installation](https://github.com/stjudecab/sjcab_peak2anno/blob/main/docs/install.md)
- [Configuration](https://github.com/stjudecab/sjcab_peak2anno/blob/main/docs/configuration.md)
- [Commands](https://github.com/stjudecab/sjcab_peak2anno/blob/main/docs/commands.md)
- [Input and output](https://github.com/stjudecab/sjcab_peak2anno/blob/main/docs/input-output.md)
- [Changelog](https://github.com/stjudecab/sjcab_peak2anno/blob/main/docs/changelog.md)

Package versions are derived from Git tags by `setuptools-scm`. Use a release
tag such as `0.1.8` when building pip or conda packages.
