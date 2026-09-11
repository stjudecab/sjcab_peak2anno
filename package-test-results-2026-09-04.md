

no tests ran in 0.00s
ERROR: file or directory not found: tests



==================================== ERRORS ====================================
__________________ ERROR collecting test_peak2anno_package.py __________________
/home/bxu2/.local/lib/python3.6/site-packages/_pytest/python.py:599: in _importtestmodule
    mod = import_path(self.path, mode=importmode, root=self.config.rootpath)
/home/bxu2/.local/lib/python3.6/site-packages/_pytest/pathlib.py:533: in import_path
    importlib.import_module(module_name)
/usr/lib64/python3.6/importlib/__init__.py:126: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
<frozen importlib._bootstrap>:994: in _gcd_import
    ???
<frozen importlib._bootstrap>:971: in _find_and_load
    ???
<frozen importlib._bootstrap>:955: in _find_and_load_unlocked
    ???
<frozen importlib._bootstrap>:665: in _load_unlocked
    ???
/home/bxu2/.local/lib/python3.6/site-packages/_pytest/assertion/rewrite.py:162: in exec_module
    source_stat, co = _rewrite_test(fn, self.config)
/home/bxu2/.local/lib/python3.6/site-packages/_pytest/assertion/rewrite.py:366: in _rewrite_test
    co = compile(tree, strfn, "exec", dont_inherit=True)
E     File "/research/rgs01/home/clusterHome/bxu2/tgz/pypi/sjcab_peak2anno/test_peak2anno_package.py", line 3
E       from __future__ import annotations
E                                        ^
E   SyntaxError: future feature annotations is not defined
=========================== short test summary info ============================
ERROR test_peak2anno_package.py
!!!!!!!!!!!!!!!!!!!! Interrupted: 1 error during collection !!!!!!!!!!!!!!!!!!!!
1 error in 0.23s
Using CPython 3.13.3
Creating virtual environment at: .venv
error: Failed to fetch: `https://pypi.org/simple/matplotlib/`
  Caused by: Could not connect, are you offline?
  Caused by: Request failed after 3 retries
  Caused by: error sending request for url (https://pypi.org/simple/matplotlib/)
  Caused by: client error (Connect)
  Caused by: dns error: failed to lookup address information: Name or service not known
  Caused by: failed to lookup address information: Name or service not known

....                                                                     [100%]
4 passed in 0.08s

....                                                                     [100%]
4 passed in 0.02s

....                                                                     [100%]
4 passed in 0.02s

....                                                                     [100%]
4 passed in 0.03s

....                                                                     [100%]
4 passed in 0.02s

FF....                                                                   [100%]
=================================== FAILURES ===================================
__________________________ test_peak2gene_default_tss __________________________

tmp_path = PosixPath('/tmp/pytest-of-bxu2/pytest-6/test_peak2gene_default_tss0')
toy_db = PosixPath('/tmp/pytest-of-bxu2/pytest-6/test_peak2gene_default_tss0/db')

    def test_peak2gene_default_tss(tmp_path: Path, toy_db: Path) -> None:
        """peak2gene should emit promoter, distal, and closest gene columns."""
        peaks = write(
            tmp_path / "peaks.bed",
            "\n".join(
                [
                    "chr1\t50\t150\tp1",
                    "chr1\t3000\t3050\tp2",
                    "chr1\t9000\t9050\tp3",
                ]
            )
            + "\n",
        )
        out = tmp_path / "peak2gene.tsv"
>       annotate_peak2gene(
            PeakGeneConfig(
                input_path=peaks,
                output_path=out,
                species="toy",
                db_path=str(toy_db),
                prom_enha_cutoffs="100bp,3000bp",
            )
        )

tests/test_package.py:125: 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
src/peak2anno/peak2gene.py:118: in annotate_peak2gene
    tss_records = resolve_tss_records(config)
src/peak2anno/peak2gene.py:50: in resolve_tss_records
    path = gene_annotation_path(
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

species = 'toy', version = 'default', isoform_set = 'all'
root_path = '/tmp/pytest-of-bxu2/pytest-6/test_peak2gene_default_tss0/db'

    def gene_annotation_path(
        species: str,
        version: str = "def",
        isoform_set: str = "all",
        root_path: Optional[str] = None,
    ) -> Path:
        """Find the gene BED for a species, version, and isoform set.
    
        The preferred database layout is ``<root>/<species>/<version>/`` with
        either ``all.gene.bed`` or ``deduplong.gene.bed``. Manifest paths and the
        older ``<root>/<species>/{tss,deduplong>/<version>.bed`` layout are also
        supported.
        """
        if isoform_set not in {"all", "deduplong"}:
            raise ValueError("isoform set must be one of: all, deduplong")
        root = db_root(root_path)
        filename = f"{isoform_set}.gene.bed"
        manifest = load_manifest(root)
        resources = manifest_resources(manifest)
    
        if version in {"default", "def"}:
            for item in resources:
                item_path = str(item.get("path", ""))
                annotation = str(item.get("annotation", ""))
                if (
                    item.get("species") == species
                    and item.get("default") is True
                    and (annotation in {isoform_set, f"{isoform_set}.gene", "gene"}
                         or item_path.endswith(filename))
                ):
                    version = str(item.get("version", "default"))
                    break
            else:
                version_candidates = sorted(
                    path.parent.name
                    for path in (root / species).glob(f"*/{filename}")
                )
                if len(version_candidates) == 1:
                    version = version_candidates[0]
                elif not version_candidates:
                    version = "default"
                else:
                    raise FileNotFoundError(
                        f"No default gene BED version for {species}; choose one with --ver"
                    )
    
        for item in resources:
            if item.get("species") != species or str(item.get("version")) != version:
                continue
            item_path = item.get("path")
            if item_path and (
                str(item.get("annotation", "")) in {isoform_set, f"{isoform_set}.gene", "gene"}
                or str(item_path).endswith(filename)
            ):
                path = root / str(item_path)
                if path.is_file():
                    return path
    
        candidates = [
            root / species / version / filename,
            root / species / f"{version}.{isoform_set}.gene.bed",
            root / species / isoform_set / f"{version}.gene.bed",
        ]
        # Compatibility with the original database package layout.
        candidates.append(root / species / ("tss" if isoform_set == "all" else "deduplong") / f"{version}.bed")
        for path in candidates:
            if path.is_file():
                return path
        searched = ", ".join(str(path) for path in candidates)
>       raise FileNotFoundError(
            f"Missing {filename} for {species} version {version} under {root}; searched: {searched}"
        )
E       FileNotFoundError: Missing all.gene.bed for toy version default under /tmp/pytest-of-bxu2/pytest-6/test_peak2gene_default_tss0/db; searched: /tmp/pytest-of-bxu2/pytest-6/test_peak2gene_default_tss0/db/toy/default/all.gene.bed, /tmp/pytest-of-bxu2/pytest-6/test_peak2gene_default_tss0/db/toy/default.all.gene.bed, /tmp/pytest-of-bxu2/pytest-6/test_peak2gene_default_tss0/db/toy/all/default.gene.bed, /tmp/pytest-of-bxu2/pytest-6/test_peak2gene_default_tss0/db/toy/tss/default.bed

src/peak2anno/db.py:182: FileNotFoundError
____________________ test_peak2gene_finds_default_gene_bed _____________________

tmp_path = PosixPath('/tmp/pytest-of-bxu2/pytest-6/test_peak2gene_finds_default_g0')

    def test_peak2gene_finds_default_gene_bed(tmp_path: Path) -> None:
        """peak2gene should find the default all.gene.bed database file."""
        db = tmp_path / "db"
        write(db / "toy" / "v2" / "all.gene.bed", "chr1\t90\t110\tGeneA\t.\t+\tENSGA\tTXA\n")
        peaks = write(tmp_path / "peaks.bed", "chr1\t100\t101\tpeak1\n")
        output = tmp_path / "output.tsv"
    
        annotate_peak2gene(
            PeakGeneConfig(
                input_path=peaks,
                output_path=output,
                species="toy",
                db_path=str(db),
            )
        )
    
        rows = read_tsv(output)
>       assert rows[1][-3:] == ["GeneA", "ENSGA", "0"]
E       AssertionError: assert ['GeneA', 'ENSGA', '9'] == ['GeneA', 'ENSGA', '0']
E         At index 2 diff: '9' != '0'
E         Use -v to get more diff

tests/test_package.py:167: AssertionError
=========================== short test summary info ============================
FAILED tests/test_package.py::test_peak2gene_default_tss - FileNotFoundError:...
FAILED tests/test_package.py::test_peak2gene_finds_default_gene_bed - Asserti...
2 failed, 4 passed in 0.10s

......                                                                   [100%]
6 passed in 0.04s

.......                                                                  [100%]
7 passed in 0.04s


==================================== ERRORS ====================================
____________________ ERROR collecting tests/test_package.py ____________________
../../../.micromamba/lib/python3.9/site-packages/_pytest/python.py:608: in _importtestmodule
    mod = import_path(self.path, mode=importmode, root=self.config.rootpath)
../../../.micromamba/lib/python3.9/site-packages/_pytest/pathlib.py:533: in import_path
    importlib.import_module(module_name)
../../../.micromamba/lib/python3.9/importlib/__init__.py:127: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
<frozen importlib._bootstrap>:1030: in _gcd_import
    ???
<frozen importlib._bootstrap>:1007: in _find_and_load
    ???
<frozen importlib._bootstrap>:986: in _find_and_load_unlocked
    ???
<frozen importlib._bootstrap>:680: in _load_unlocked
    ???
../../../.micromamba/lib/python3.9/site-packages/_pytest/assertion/rewrite.py:168: in exec_module
    exec(co, module.__dict__)
tests/test_package.py:19: in <module>
    from peak2anno.cli import build_parser, main
E     File "/research/rgs01/home/clusterHome/bxu2/tgz/pypi/sjcab_peak2anno/src/peak2anno/cli.py", line 223
E       references.append(f"gene/TSS: {resolve_tss_path(PeakGeneConfig(
E                                                                      ^
E   SyntaxError: EOL while scanning string literal
=========================== short test summary info ============================
ERROR tests/test_package.py
!!!!!!!!!!!!!!!!!!!! Interrupted: 1 error during collection !!!!!!!!!!!!!!!!!!!!
1 error in 0.25s

....F.FF                                                                 [100%]
=================================== FAILURES ===================================
_________________________ test_narrow_context_priority _________________________

tmp_path = PosixPath('/tmp/pytest-of-bxu2/pytest-9/test_narrow_context_priority0')
context_dir = PosixPath('/tmp/pytest-of-bxu2/pytest-9/test_narrow_context_priority0/context')

    def test_narrow_context_priority(tmp_path: Path, context_dir: Path) -> None:
        """narrow2feature should assign the first priority feature that overlaps."""
        peaks = write(
            tmp_path / "peaks.bed",
            "chr1\t0\t100\tp1\nchr1\t100\t200\tp2\nchr1\t200\t300\tp3\n",
        )
        out = tmp_path / "narrow.tsv"
        output, summary = annotate_narrow_context(
            ContextConfig(
                input_path=peaks,
                output_path=out,
                species="toy",
                context_dir=context_dir,
                overlap_cutoff="0.5",
            )
        )
        rows = read_tsv(output)
        assert [row[-1] for row in rows[1:]] == [
            "Promoter.Up",
            "Exon",
            "Exon",
        ]
>       summary_rows = read_tsv(summary)

tests/test_package.py:219: 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

path = None

    def read_tsv(path: Path) -> list[list[str]]:
        """Read a tab-separated file into rows."""
>       with path.open(newline="", encoding="utf-8") as handle:
E       AttributeError: 'NoneType' object has no attribute 'open'

tests/test_package.py:32: AttributeError
_____________________ test_broad_context_reports_fractions _____________________

tmp_path = PosixPath('/tmp/pytest-of-bxu2/pytest-9/test_broad_context_reports_fra0')
context_dir = PosixPath('/tmp/pytest-of-bxu2/pytest-9/test_broad_context_reports_fra0/context')

    def test_broad_context_reports_fractions(tmp_path: Path, context_dir: Path) -> None:
        """broad2feature should report per-feature fractions instead of priority-only labels."""
        peaks = write(tmp_path / "peaks.bed", "chr1\t0\t100\tp1\n")
        out = tmp_path / "broad.tsv"
        output, summary = annotate_broad_context(
            ContextConfig(
                input_path=peaks,
                output_path=out,
                species="toy",
                context_dir=context_dir,
                overlap_cutoff="1bp",
            )
        )
        rows = read_tsv(output)
        assert "Promoter.Up_bp" in rows[0]
        assert rows[1][rows[0].index("Promoter.Up_bp")] == "100"
        assert rows[1][rows[0].index("Exon_bp")] == "50"
        assert rows[1][-2:] == ["Promoter.Up", "1.000000"]
>       assert read_tsv(summary)[0] == ["Feature", "PrimaryRegions", "OverlapBp", "OverlapFractionOfInputBp"]

tests/test_package.py:263: 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

path = None

    def read_tsv(path: Path) -> list[list[str]]:
        """Read a tab-separated file into rows."""
>       with path.open(newline="", encoding="utf-8") as handle:
E       AttributeError: 'NoneType' object has no attribute 'open'

tests/test_package.py:32: AttributeError
_____________________ test_peak2state_reports_named_states _____________________

tmp_path = PosixPath('/tmp/pytest-of-bxu2/pytest-9/test_peak2state_reports_named_0')

    def test_peak2state_reports_named_states(tmp_path: Path) -> None:
        """peak2state should map state IDs and compute overlap fractions."""
        peaks = write(tmp_path / "peaks.bed", "chr1\t0\t100\tp1\n")
        states = write(tmp_path / "states.bed", "chr1\t0\t60\t1\nchr1\t60\t100\t2\n")
        state2name = write(tmp_path / "state2name.tsv", "1\tActive\n2\tRepressed\n")
        out = tmp_path / "states.tsv"
        output, summary = annotate_peak_state(
            StateConfig(
                input_path=peaks,
                states_path=states,
                state2name=state2name,
                output_path=out,
                overlap_cutoff="1bp",
            )
        )
        rows = read_tsv(output)
        assert rows[1][rows[0].index("Active_bp")] == "60"
        assert rows[1][rows[0].index("Repressed_fraction")] == "0.400000"
        assert rows[1][-2:] == ["Active", "0.600000"]
>       assert read_tsv(summary)[1] == ["Active", "1", "60", "0.600000"]

tests/test_package.py:285: 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 

path = None

    def read_tsv(path: Path) -> list[list[str]]:
        """Read a tab-separated file into rows."""
>       with path.open(newline="", encoding="utf-8") as handle:
E       AttributeError: 'NoneType' object has no attribute 'open'

tests/test_package.py:32: AttributeError
=========================== short test summary info ============================
FAILED tests/test_package.py::test_narrow_context_priority - AttributeError: ...
FAILED tests/test_package.py::test_broad_context_reports_fractions - Attribut...
FAILED tests/test_package.py::test_peak2state_reports_named_states - Attribut...
3 failed, 5 passed in 0.10s

........                                                                 [100%]
8 passed in 0.04s

...........                                                              [100%]
11 passed in 0.06s

...........                                                              [100%]
11 passed in 0.11s

...........                                                              [100%]
11 passed in 0.09s

## 2026-09-07 — missing-reference installation prompt

Command:

```bash
PYTHONPATH=src /research/rgs01/home/clusterHome/bxu2/.micromamba/bin/python -m pytest -q
```

Result:

```text
.............                                                            [100%]
13 passed in 0.10s

Verification rerun after Python 3.7-compatible command quoting:

```text
.............                                                            [100%]
13 passed in 0.12s

## 2026-09-10 — promoter direction and rc configuration

Validation commands:

```bash
ruff check src tests
PYTHONPATH=src .venv/bin/python -m compileall -q src tests
```

Result:

```text
All checks passed!
```

The existing pytest suite was not rerun because the available environments did
not provide a compatible pytest installation.

Additional configuration check on 2026-09-10:

```text
SJCAB_PEAK2ANNO_GENE_TYPE=nomicro SJCAB_PEAK2ANNO_ISO_SET=deduplong
parser defaults: nomicro deduplong
ruff: All checks passed!

2026-09-11 validation:

```text
cutoff/CLI checks passed
ruff: All checks passed!
```

Feature/state output-mode validation on 2026-09-10:

```text
max,percent -> [max, percent]
percent,max -> [percent, max]
ruff: All checks passed!
```

Output-mode configuration validation on 2026-09-10:

```text
SJCAB_PEAK2ANNO_2FEATURE_OUT=percent,max -> percent,max
SJCAB_PEAK2ANNO_2STATE_OUT=max -> max
ruff: All checks passed!
```
```
```
```
