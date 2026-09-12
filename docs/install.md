# Installation

## Pip

```bash
python -m pip install sjcab_peak2anno
```

Pip does not require `bedtools` or `pybedtools`. Both are detected and used
when available; without them, peak2anno uses a slower Python interval
fallback. For optional acceleration:

```bash
python -m pip install pybedtools
bedtools --version
```

## Conda

The conda package requires both tools:

```bash
conda install -c stjudecab sjcab_peak2anno sjcab_peak2anno_db=0.1.8 bedtools pybedtools
```

## Test a source checkout

```bash
python -m pip install -e ".[test]"
python -m pytest -q
```

Or run without installing the package:

```bash
PYTHONPATH=src python -m pytest -q
```

## Smoke test after installation

```bash
peak2anno peak2gene tests_data/peaks.bed \
  --tss-bed tests_data/tss.bed \
  --prom-enha-cutoffs 100bp,3kb \
  --output peak2anno-smoke.tsv
cat peak2anno-smoke.tsv
```

The same command is available as `sjcab-peak2anno`.
