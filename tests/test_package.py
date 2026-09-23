"""Tests for the peak2anno package commands."""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from peak2anno.feature import (
    FeatureConfig,
    StateConfig,
    annotate_broad_feature,
    annotate_narrow_feature,
    annotate_peak_state,
)
from peak2anno.cli import build_parser, main
from peak2anno.config import load_settings
from peak2anno.intervals import read_regions
from peak2anno.peak2gene import PeakGeneConfig, annotate_peak2gene, promoter_records
from peak2anno.runtime import ToolStatus, resolve_backend


def write(path: Path, text: str) -> Path:
    """Write text to a path for tests."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def read_tsv(path: Path) -> list[list[str]]:
    """Read a tab-separated file into rows."""
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def test_backend_auto_prefers_bedtools() -> None:
    """auto should prefer the native executable over the Python wrapper."""
    status = ToolStatus(bedtools="/usr/bin/bedtools", pybedtools=True)
    assert resolve_backend("auto", status) == "bedtools"
    assert resolve_backend("pybedtools", status) == "pybedtools"
    assert resolve_backend("auto", ToolStatus(bedtools=None, pybedtools=False)) == "python"


def test_backend_rejects_unavailable_tool() -> None:
    """Explicit backend requests should fail clearly when unavailable."""
    with pytest.raises(ValueError, match="bedtools was not found"):
        resolve_backend("bedtools", ToolStatus(bedtools=None, pybedtools=False))


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


def make_feature_dir(tmp_path: Path) -> Path:
    """Create a minimal default feature annotation directory."""
    feature = tmp_path / "feature"
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
        write(feature / filename, text)
    return feature


@pytest.fixture
def toy_db(tmp_path: Path) -> Path:
    """Provide a toy annotation database for a test."""
    return make_toy_db(tmp_path)


@pytest.fixture
def feature_dir(tmp_path: Path) -> Path:
    """Provide a toy feature annotation directory for a test."""
    return make_feature_dir(tmp_path)


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
            prom_enha_cutoffs="100bp,3000bp",
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
    assert rows[2][-7:] == [".", ".", "GeneB,GeneA", "ENSGB,ENSGA", "GeneB", "ENSGB", "1951"]
    assert rows[3][-3:] == ["GeneC", "ENSGC", "950"]


def test_peak2gene_keeps_promoter_and_distal_columns_exclusive(tmp_path: Path, toy_db: Path) -> None:
    """voom-compatible output reports distal genes only without promoters."""
    peaks = write(tmp_path / "peaks.bed", "chr1\t3000\t3050\np1\n")
    output = tmp_path / "output.tsv"
    annotate_peak2gene(
        PeakGeneConfig(
            input_path=peaks,
            output_path=output,
            species="toy",
            db_path=str(toy_db),
            prom_enha_cutoffs="2500bp,6000bp,2500bp",
        )
    )
    row = read_tsv(output)[1]
    assert row[3:7] == ["GeneB", "ENSGB", ".", "."]


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
    assert "(default: v31)" in help_text
    assert "--prom-enha-cutoffs" in help_text
    assert "(default: 2kb,50kb,2kb)" in help_text


def test_common_short_options_parse() -> None:
    """Frequently used long options should have concise aliases."""
    args = build_parser().parse_args(
        [
            "peak2gene",
            "-i", "peaks.txt",
            "-f", "txt",
            "-I", "txt",
            "-H", "yes",
            "-C", "0,1,2",
            "-R", "0",
            "-g", "genes.bed",
            "-t", "tss.bed",
        ]
    )
    assert args.input_option == Path("peaks.txt")
    assert args.output_format == "txt"
    assert args.input_format == "txt"
    assert args.header == "yes"
    assert args.columns == "0,1,2"
    assert args.region_column == 0
    assert args.gene_bed == Path("genes.bed")
    assert args.tss_bed == Path("tss.bed")


def test_rc_settings_are_overridden_by_environment(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The selected rc file supplies defaults and same-name environment values win."""
    rc = write(
        tmp_path / "settings.rc",
        "\n".join(
            [
                "SJCAB_PEAK2ANNO_SPECIES_VERSIONS=mm10:vM22",
                "SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS=111bp,333bp,222bp",
                "SJCAB_PEAK2ANNO_GENE_TYPE=nomicro",
                "SJCAB_PEAK2ANNO_ISO_SET=deduplong",
                "SJCAB_PEAK2ANNO_2FEATURE_OUT=percent,max",
                "SJCAB_PEAK2ANNO_2STATE_OUT=max",
            ]
        )
        + "\n",
    )
    monkeypatch.setenv("SJCAB_PEAK2ANNO_CONFIG", str(rc))
    monkeypatch.setenv("SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS_hg38_v31", "gene0.5,transcript2")
    settings = load_settings()
    assert settings.default_species == "mm10"
    assert settings.version_for("mm10") == "vM22"
    assert settings.prom_enha_cutoffs == "111bp,333bp,222bp"
    assert settings.prom_enha_cutoffs_for("hg38", "v31") == "gene0.5,transcript2"
    assert settings.gene_type == "nomicro"
    assert settings.iso_set == "deduplong"
    assert settings.feature_out == "percent,max"
    assert settings.state_out == "max"


def test_feature_percent_mode_uses_order_list(tmp_path: Path) -> None:
    """Feature percentages should partition overlap in order-list order."""
    feature = make_feature_dir(tmp_path)
    write(feature / "order.lst", "exon\npromoter.up\npromoter.down\nintron\ntes\ndis5\ndis3\nintergenic\n")
    peaks = write(tmp_path / "peaks.bed", "chr1\t0\t200\tp1\n")
    output = tmp_path / "features.tsv"
    annotate_narrow_feature(
        FeatureConfig(
            input_path=peaks,
            output_path=output,
            species="toy",
            feature_dir=feature,
            output_format="txt",
            output_mode="percent",
        )
    )
    rows = read_tsv(output)
    assert rows[0][-8:-5] == ["Exon_percent", "Promoter.Up_percent", "Promoter.Down_percent"]
    assert rows[1][-8:-5] == ["75.000000", "25.000000", "0.000000"]


def test_promoter_cutoffs_are_strand_aware() -> None:
    """Upstream/downstream promoter windows should follow the gene strand."""
    from peak2anno.intervals import BedRecord, InputRegion, IntervalIndex

    records = [
        BedRecord("chr1", 100, 101, ("chr1", "100", "101", "plus", ".", "+")),
        BedRecord("chr1", 300, 301, ("chr1", "300", "301", "minus", ".", "-")),
    ]
    index = IntervalIndex(records)
    assert [record.name for record in promoter_records(index, InputRegion("chr1", 0, 50, ()), 60, 10)] == ["plus"]
    assert promoter_records(index, InputRegion("chr1", 150, 200, ()), 60, 10) == []
    assert [record.name for record in promoter_records(index, InputRegion("chr1", 350, 400, ()), 60, 10)] == ["minus"]


def test_missing_database_reference_can_be_declined(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A missing selected reference should ask before running the installer."""
    peaks = write(tmp_path / "peaks.bed", "chr1\t100\t101\tpeak1\n")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", type("TTY", (), {"isatty": lambda self: True})())
    monkeypatch.setattr("builtins.input", lambda _prompt: "no")

    with pytest.raises(SystemExit):
        main(["peak2gene", str(peaks), "-s", "missing", "--ver", "v9", "-d", str(tmp_path / "db")])

    error = capsys.readouterr().err
    assert "sjcab-peak2anno-db download-gencode-bed missing v9" in error
    assert "installation declined" in error


def test_missing_database_reference_can_be_installed(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A yes response should run the proposed installer and retry the command."""
    peaks = write(tmp_path / "peaks.bed", "chr1\t100\t101\tpeak1\n")
    db = tmp_path / "db"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "stdin", type("TTY", (), {"isatty": lambda self: True})())
    monkeypatch.setattr("builtins.input", lambda _prompt: "yes")

    def fake_install(command: list[str], check: bool) -> subprocess.CompletedProcess[str]:
        assert check is True
        assert command[1:5] == ["sjcab-peak2anno-db", "download-gencode-bed", "missing", "v9"]
        write(db / "missing" / "v9" / "all.gene.bed", "chr1\t99\t100\tGeneA\t.\t+\tENSGA\tTXA\n")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("peak2anno.cli.shutil.which", lambda _name: "/fake/sjcab-peak2anno-db")
    monkeypatch.setattr("peak2anno.cli.subprocess.run", fake_install)
    assert main(["peak2gene", str(peaks), "-s", "missing", "--ver", "v9", "-d", str(db)]) == 0
    assert "GeneA" in capsys.readouterr().out


def test_cli_writes_stdout_and_run_log(tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    """The CLI should stream output and record resolved references without -o."""
    peaks = write(tmp_path / "peaks.bed", "chr1\t100\t101\tpeak1\n")
    tss = write(tmp_path / "tss.bed", "chr1\t99\t100\tGeneA\t.\t+\tENSGA\tTXA\n")
    monkeypatch.chdir(tmp_path)

    assert main(["peak2gene", str(peaks), "--tss-bed", str(tss), "--prom-enha-cutoffs", "100bp,3kb"]) == 0

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


def test_region_text_accepts_custom_delimiter(tmp_path: Path) -> None:
    """Text input should parse the configured custom region delimiter."""
    path = write(tmp_path / "regions.txt", "region\nchr1|100|200\n")
    header, regions = read_regions(path, input_format="txt", txt_delimiter="|")
    assert header == ["region"]
    assert regions[0].region_name == "chr1:100-200"


def test_loop2gene_merges_two_anchor_annotations(tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    """loop2gene should annotate both BEDPE anchors in one table."""
    loops = write(tmp_path / "loops.bedpe", "chr1\t50\t150\tchr1\t300\t400\n")
    tss = write(tmp_path / "tss.bed", "chr1\t99\t100\tGeneA\t.\t+\tENSGA\tTXA\n")
    monkeypatch.chdir(tmp_path)
    assert main(["loop2gene", str(loops), "--tss-bed", str(tss), "--output-format", "txt"]) == 0
    output = capsys.readouterr().out
    assert "anchor1_Closest_Gene" in output
    assert "anchor2_Closest_Gene" in output


def test_combined_annotations_merge_columns(tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    """The combined command should merge multiple annotation result blocks."""
    peaks = write(tmp_path / "peaks.bed", "chr1\t50\t150\tpeak1\n")
    tss = write(tmp_path / "tss.bed", "chr1\t99\t100\tGeneA\t.\t+\tENSGA\tTXA\n")
    feature = make_feature_dir(tmp_path)
    monkeypatch.chdir(tmp_path)
    assert main([
        "peak2gene", "narrow2feature", str(peaks), "--tss-bed", str(tss),
        "--feature-dir", str(feature), "--workers", "2", "--output-format", "txt",
    ]) == 0
    output = capsys.readouterr().out
    assert "Closest_Gene" in output
    assert "FeatureAssignment" in output


def test_narrow_feature_priority(tmp_path: Path, feature_dir: Path) -> None:
    """narrow2feature should assign the first priority feature that overlaps."""
    peaks = write(
        tmp_path / "peaks.bed",
        "chr1\t0\t100\tp1\nchr1\t100\t200\tp2\nchr1\t200\t300\tp3\n",
    )
    out = tmp_path / "narrow.tsv"
    output, summary = annotate_narrow_feature(
        FeatureConfig(
            input_path=peaks,
            output_path=out,
            species="toy",
            feature_dir=feature_dir,
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


def test_narrow_feature_finds_feature_dir_under_db_path(tmp_path: Path) -> None:
    """Feature commands should search the database root when no directory is given."""
    feature = make_feature_dir(tmp_path)
    db_feature = tmp_path / "db" / "toy" / "features"
    shutil.copytree(feature, db_feature)
    peaks = write(tmp_path / "peaks.bed", "chr1\t0\t100\tp1\n")
    output = tmp_path / "narrow.tsv"

    annotate_narrow_feature(
        FeatureConfig(
            input_path=peaks,
            output_path=output,
            species="toy",
            db_path=str(tmp_path / "db"),
            overlap_cutoff="0.5",
        )
    )

    assert read_tsv(output)[1][-1] == "Promoter.Up"


def test_broad_feature_reports_fractions(tmp_path: Path, feature_dir: Path) -> None:
    """broad2feature should report per-feature fractions instead of priority-only labels."""
    peaks = write(tmp_path / "peaks.bed", "chr1\t0\t100\tp1\n")
    out = tmp_path / "broad.tsv"
    output, summary = annotate_broad_feature(
        FeatureConfig(
            input_path=peaks,
            output_path=out,
            species="toy",
            feature_dir=feature_dir,
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
