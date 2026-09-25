# Installation

## Pip

```bash
python -m pip install sjcab_peak2anno
```

Pip does not require `bedtools`; `auto` uses the Python interval backend. To
explicitly use bedtools and generate a reviewable bash script:

```bash
SJCAB_PEAK2ANNO_BACKEND=bedtools peak2anno peak2gene peaks.bed
```

This requires `bedtools` in `PATH`.

## Conda

The conda package does not require external interval tools:

```bash
conda install -c stjudecab sjcab_peak2anno sjcab_peak2anno_db
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
