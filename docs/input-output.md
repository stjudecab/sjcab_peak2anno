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

For text input, the region coordinate delimiter can be selected with the
`SJCAB_PEAK2ANNO_TXT_DELIMITER` environment variable or the same setting in
the RC file. The default is `auto`, which recognizes the delimiters shown
above. Set it to a literal delimiter such as `|` when regions use a custom
form:

```bash
SJCAB_PEAK2ANNO_TXT_DELIMITER='|' peak2anno peak2gene regions.txt \
  --input-format txt --tss-bed tests_data/tss.bed
```

With that setting, the text region is written as `chr1|100|200`. The setting
only changes parsing of the chromosome/start/end expression; tabs or spaces
still separate columns in a text table.

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
Use `--order-lst` to choose a custom feature priority list; `def`, `default`,
and `none` select the database-root `order.lst`, while `utr` selects
`order.utr.lst`.
Check the selected order list for feature meanings. For two-column order
lists, the second column is used as the output feature name; standard feature
names are used when a mapping is not found.
