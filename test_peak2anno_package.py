"""Tests for the peak2anno package commands."""

from __future__ import annotations

import json
from pathlib import Path

from peak2anno.context import (
    ContextConfig,
    StateConfig,
    annotate_broad_context,
    annotate_narrow_context,
    annotate_peak_state,
)
from peak2anno.peak2gene import PeakGeneConfig, annotate_peak2gene


def write(path: Path, text: str) -> Path:
    """Write text to a path for tests."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def read_tsv(path: Path) -> list[list[str]]:
    """Read a tab-separated file into rows."""
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()]


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


def test_peak2gene_default_tss(tmp_path: Path) -> None:
    """peak2gene should emit promoter, distal, and closest gene columns."""
    db = make_toy_db(tmp_path)
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
            db_path=str(db),
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


def test_narrow_context_priority(tmp_path: Path) -> None:
    """narrow2context should assign the first priority feature that overlaps."""
    context = make_context_dir(tmp_path)
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
            context_dir=context,
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


def test_broad_context_reports_fractions(tmp_path: Path) -> None:
    """broad2context should report per-feature fractions instead of priority-only labels."""
    context = make_context_dir(tmp_path)
    peaks = write(tmp_path / "peaks.bed", "chr1\t0\t100\tp1\n")
    out = tmp_path / "broad.tsv"
    output, summary = annotate_broad_context(
        ContextConfig(
            input_path=peaks,
            output_path=out,
            species="toy",
            context_dir=context,
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
