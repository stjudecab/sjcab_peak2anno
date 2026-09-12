# Input and output

## Input formats

Single-region commands accept BED-like files or region text. BED coordinates
default to columns 0, 1, and 2; use `--columns` to change them.

```text
chr1    100    200    peak_a
```

Region text defaults to column 0 and accepts delimiters such as `:`, `-`,
`*`, `=`, `/`, `^`, `;`, `_`, `%`, `$`, and `,`:

```text
region
chr1:100-200
chr2^300=450
```

`--input-format auto` detects BED, headered text, or headerless text. Force a
format with `bed`, `txt`, or `txtnohead`; use `--region-column` for another
text column.

## Output formats

`--output-format` accepts `auto`, `bed`, `bedpe`, `txt`, and `txtnohead`.

- `auto` follows the detected input format.
- `txt` is tab-delimited with a header.
- `txtnohead` is tab-delimited without a header.
- `bed` writes chromosome/start/end output without a header.
- `bedpe` is the loop equivalent of BED output.

```bash
peak2anno peak2gene peaks.bed --tss-bed tests_data/tss.bed \
  -o result.tsv --output-format txt
peak2anno peak2gene peaks.bed --tss-bed tests_data/tss.bed \
  -o result.bed --output-format bed
```

Feature and state output modes are controlled by `--output-mode` and can emit
the maximum assignment, ordered percentages, or both. Feature percentages
allocate each overlapping base to the first feature in `order.lst`.
