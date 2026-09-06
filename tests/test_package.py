"""Tests for the peak2anno package commands."""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pytest

from peak2anno.context import (
    ContextConfig,
    StateConfig,
    annotate_broad_context,
    annotate_narrow_context,
    annotate_peak_state,
)
from peak2anno.cli import build_parser, main
from peak2anno.intervals import read_regions
from peak2anno.peak2gene import PeakGeneConfig, annotate_peak2gene


def write(path: Path, text: str) -> Path:
    """Write text to a path for tests."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def read_tsv(path: Path) -> list[list[str]]:
    """Read a tab-separated file into rows."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def make_toy_db(tmp_path: Path) -> Path:
    """Create a minimal sjcab_peak2anno_db-like directory."""
    root = tmp_path / "db"
    manifest = {
        "files": [
            {
                "species": "toy",
                "annotation": "tss",
                "version": "v1",
                "default": True,
                "path": "toy/tss/v1.bed",
            },
            {
                "species": "toy",
                "annotation": "deduplong",
                "version": "v1",
                "default": True,
                "path": "toy/deduplong/v1.bed",
            },
        ]
    }
    write(root / "manifest.json", json.dumps(manifest))
    write(
        root / "toy" / "tss" / "v1.bed",
        "\n".join(
            [
                "chr1\t99\t100\tGeneA\t.\t+\tENSGA\tTXA",
                "chr1\t5000\t5001\tGeneB\t.\t+\tENSGB\tTXB",
                "chr1\t10000\t10001\tGeneC\t.\t-\tENSGC\tTXC",
            ]
        )
        + "\n",
    )
    write(
        root / "toy" / "deduplong" / "v1.bed",
        "\n".join(
            [
                "chr1\t99\t500\tGeneA\t.\t+\tENSGA\tTXA",
                "chr1\t4900\t5001\tGeneB\t.\t-\tENSGB\tTXB",
            ]
        )
        + "\n",
    )
    return root


def make_context_dir(tmp_path: Path) -> Path:
    """Create a minimal default context annotation directory."""
    context = tmp_path / "context"
    files = {
        "2kb.promoter.up.bed": "chr1\t0\t100\n",
        "2kb.promoter.down.bed": "",
        "2kb.exon.bed": "chr1\t50\t250\n",
        "2kb.intron.bed": "",
        "2kb.tes.bed": "",
        "2kb.dis5.bed": "",
        "2kb.dis3.bed": "",
        "2kb.intergenic.bed": "chr1\t200\t300\n",
    }
    for filename, text in files.items():
        write(context / filename, text)
    return context


@pytest.fixture
def toy_db(tmp_path: Path) -> Path:
    """Provide a toy annotation database for a test."""
    return make_toy_db(tmp_path)


@pytest.fixture
def context_dir(tmp_path: Path) -> Path:
    """Provide a toy context annotation directory for a test."""
    return make_context_dir(tmp_path)


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
    annotate_peak2gene(
        PeakGeneConfig(
            input_path=peaks,
            output_path=out,
            species="toy",
            db_path=str(toy_db),
            promoter_cutoff="100bp",
            enhancer_cutoff="3000bp",
        )
    )
    rows = read_tsv(out)
    assert rows[0][-7:] == [
        "Gene_100bp",
        "Gencode_ids",
        "Gene_100bp-3kb",
        "Gencode_ids",
        "Closest_Gene",
        "Gencode_id",
        "Distance",
    ]
    assert rows[1][-7:] == ["GeneA", "ENSGA", ".", ".", "GeneA", "ENSGA", "0"]
    assert rows[2][-7:] == [".", ".", "GeneB,GeneA", "ENSGB,ENSGA", "GeneB", "ENSGB", "1950"]
    assert rows[3][-3:] == ["GeneC", "ENSGC", "950"]


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
    assert rows[1][-3:] == ["GeneA", "ENSGA", "9"]


def test_peak2gene_help_shows_short_options_and_defaults(capsys: pytest.CaptureFixture[str]) -> None:
    """The peak2gene help should show concise options and their defaults."""
    with pytest.raises(SystemExit):
        build_parser().parse_args(["peak2gene", "--help"])
    help_text = capsys.readouterr().out
    assert "--ver" in help_text
    assert "--iso" in help_text
    assert "(default: hg38)" in help_text
    assert "(default: def)" in help_text


def test_cli_writes_stdout_and_run_log(tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    """The CLI should stream output and record resolved references without -o."""
    peaks = write(tmp_path / "peaks.bed", "chr1\t100\t101\tpeak1\n")
    tss = write(tmp_path / "tss.bed", "chr1\t99\t100\tGeneA\t.\t+\tENSGA\tTXA\n")
    monkeypatch.chdir(tmp_path)

    assert main(["peak2gene", str(peaks), "--tss-bed", str(tss), "--promoter-cutoff", "100bp"]) == 0

    assert "GeneA" in capsys.readouterr().out
    run_log = (tmp_path / ".run.log").read_text(encoding="utf-8")
    assert "command:" in run_log
    assert "output: stdout" in run_log
    assert str(tss.resolve()) in run_log


def test_region_text_accepts_header_and_common_delimiters(tmp_path: Path) -> None:
    """Text input should parse a header and non-colon region delimiters."""
    path = write(tmp_path / "regions.txt", "region\nchr1^100=200\n")
    header, regions = read_regions(path, input_format="txt")
    assert header == ["region"]
    assert (regions[0].chrom, regions[0].start, regions[0].end) == ("chr1", 100, 200)


def test_loop2anno_merges_two_anchor_annotations(tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    """loop2anno should annotate both BEDPE anchors in one table."""
    loops = write(tmp_path / "loops.bedpe", "chr1\t50\t150\tchr1\t300\t400\n")
    tss = write(tmp_path / "tss.bed", "chr1\t99\t100\tGeneA\t.\t+\tENSGA\tTXA\n")
    monkeypatch.chdir(tmp_path)
    assert main(["loop2anno", str(loops), "--tss-bed", str(tss), "--output-format", "txt"]) == 0
    output = capsys.readouterr().out
    assert "anchor1_Closest_Gene" in output
    assert "anchor2_Closest_Gene" in output


def test_combined_annotations_merge_columns(tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    """The combined command should merge multiple annotation result blocks."""
    peaks = write(tmp_path / "peaks.bed", "chr1\t50\t150\tpeak1\n")
    tss = write(tmp_path / "tss.bed", "chr1\t99\t100\tGeneA\t.\t+\tENSGA\tTXA\n")
    context = make_context_dir(tmp_path)
    monkeypatch.chdir(tmp_path)
    assert main([
        "peak2gene", "narrow2context", str(peaks), "--tss-bed", str(tss),
        "--context-dir", str(context), "--workers", "2", "--output-format", "txt",
    ]) == 0
    output = capsys.readouterr().out
    assert "Closest_Gene" in output
    assert "FeatureAssignment" in output


def test_narrow_context_priority(tmp_path: Path, context_dir: Path) -> None:
    """narrow2context should assign the first priority feature that overlaps."""
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
    summary_rows = read_tsv(summary)
    assert summary_rows[0][4:7] == ["Promoter.Up", "Promoter.Down", "Exon"]
    assert summary_rows[1][4:7] == ["1", "0", "2"]


def test_narrow_context_finds_context_dir_under_db_path(tmp_path: Path) -> None:
    """Context commands should search the database root when no directory is given."""
    context = make_context_dir(tmp_path)
    db_context = tmp_path / "db" / "toy" / "context"
    shutil.copytree(context, db_context)
    peaks = write(tmp_path / "peaks.bed", "chr1\t0\t100\tp1\n")
    output = tmp_path / "narrow.tsv"

    annotate_narrow_context(
        ContextConfig(
            input_path=peaks,
            output_path=output,
            species="toy",
            db_path=str(tmp_path / "db"),
            overlap_cutoff="0.5",
        )
    )

    assert read_tsv(output)[1][-1] == "Promoter.Up"


def test_broad_context_reports_fractions(tmp_path: Path, context_dir: Path) -> None:
    """broad2context should report per-feature fractions instead of priority-only labels."""
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
    assert read_tsv(summary)[0] == ["Feature", "PrimaryRegions", "OverlapBp", "OverlapFractionOfInputBp"]


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
    assert read_tsv(summary)[1] == ["Active", "1", "60", "0.600000"]
