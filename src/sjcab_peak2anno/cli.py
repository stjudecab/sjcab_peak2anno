"""Console entry point for peak2anno."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from contextlib import contextmanager
import fcntl
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, Sequence

from . import __version__
from .features import (
    FeatureConfig,
    StateConfig,
    annotate_broad_feature,
    annotate_narrow_feature,
    annotate_peak_state,
    resolve_order_path,
    resolve_features,
)
from .config import Settings, load_settings
from .db import available_versions, db_root, gene_annotation_path
from .peak2gene import PeakGeneConfig, annotate_peak2gene, resolve_tss_path
from .loops import annotate_loop
from .intervals import detect_output_format, read_regions, write_table
from .runtime import BACKENDS, detect_tools, resolve_backend, warn_if_slow


DB_PACKAGE = "sjcab_peak2anno_db==0.1.8"


def add_common_feature_args(parser: argparse.ArgumentParser, settings: Settings) -> None:
    """Add arguments shared by narrow2feature and broad2feature."""
    add_input_args(parser, "Input BED/TSV or region-text file.", settings.txt_delimiter)
    parser.add_argument("-o", "--output", type=Path, help="Output TSV path; defaults to stdout.")
    parser.add_argument("-s", "--species", default=settings.default_species, help="Species key, for example hg38 or mm10.")
    parser.add_argument("-d", "--db-path", default=settings.db_path, help="Database root; defaults to rc/env or ~/.sjcab_peak2anno_db.")
    parser.add_argument(
        "-c", "--feature-dir", type=Path,
        help="Feature BED directory. If omitted, search --db-path for the species feature files.",
    )
    parser.add_argument("-I", "--input-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Input format.")
    parser.add_argument("-C", "--columns", help="BED columns as comma-separated zero-based indexes, for example 0,1,2.")
    parser.add_argument("-R", "--region-column", type=int, default=0, help="Zero-based region column for txt input.")
    parser.add_argument("-f", "--output-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Output format; auto follows input format.")
    parser.add_argument("-F", "--features", help="Comma-separated feature BEDs or .lst file.")
    parser.add_argument("-L", "--feature-labels", help="Comma-separated feature labels or .labels.lst file.")
    parser.add_argument(
        "-O", "--order-lst", default="def",
        help="Feature priority list; inspect order.lst for feature meanings. def/default/none uses DB_PATH/order.lst; utr uses DB_PATH/order.utr.lst; or provide a path.",
    )
    parser.add_argument(
        "-x", "--overlap-cutoff",
        default="1bp",
        help="Minimum overlap. Examples: 1bp, 10bp, 0.1 for 10%% of peak, 10%%.",
    )
    parser.add_argument("-H", "--header", choices=["auto", "yes", "no"], default="auto", help="Input header handling.")
    parser.add_argument("-m", "--summary", type=Path, help="Summary TSV path.")
    parser.add_argument("-p", "--plot", action="store_true", help="Write PNG/PDF bar and pie plots for summary counts.")


def add_input_args(parser: argparse.ArgumentParser, help_text: str, txt_delimiter: str = "auto") -> None:
    """Add positional and short-option input forms."""
    parser.add_argument("input", nargs="?", type=Path, help=help_text)
    parser.add_argument("-i", "--input", dest="input_option", type=Path, metavar="INPUT", help="Input file.")
    parser.set_defaults(txt_delimiter=txt_delimiter)


def add_backend_arg(parser: argparse.ArgumentParser, settings: Settings) -> None:
    """Add the interval backend selector shared by annotation commands."""
    parser.add_argument(
        "-b", "--backend", choices=BACKENDS, default=settings.backend,
        help="Interval backend. auto prefers bedtools, then pybedtools, then Python.",
    )


def add_db_install_args(parser: argparse.ArgumentParser) -> None:
    """Add automatic database installation controls."""
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-a", "--auto-install-db",
        dest="auto_install_db",
        action="store_true",
        default=False,
        help="Install missing gene/feature BEDs automatically.",
    )
    group.add_argument(
        "-A", "--no-auto-install-db",
        dest="auto_install_db",
        action="store_false",
        help="Do not install missing database files; fail with the missing path.",
    )


def add_gene_cutoff_args(parser: argparse.ArgumentParser, settings: Settings) -> None:
    """Add the three-part promoter/enhancer cutoff option."""
    parser.add_argument(
        "-P", "--prom-enha-cutoffs",
        default=settings.prom_enha_cutoffs,
        help="promoter-up,enhancer,promoter-down; supports geneN and transcriptN.",
    )


def add_output_mode_arg(parser: argparse.ArgumentParser, default: str, help_text: str) -> None:
    """Add max/percent output selection with its resolved default."""
    parser.add_argument(
        "-M", "--output-mode",
        choices=["max", "percent", "max,percent", "percent,max"],
        default=default,
        help=help_text,
    )


def parse_columns(value: Optional[str]) -> Optional[tuple[int, int, int]]:
    """Parse three zero-based BED coordinate column indexes."""
    if value is None:
        return None
    values = tuple(int(item.strip()) for item in value.split(","))
    if len(values) != 3 or min(values) < 0:
        raise ValueError("--columns must contain three non-negative indexes, for example 0,1,2")
    return values  # type: ignore[return-value]


def parse_loop_columns(value: str) -> tuple[int, int, int, int, int, int]:
    """Parse six zero-based BEDPE coordinate indexes."""
    values = tuple(int(item.strip()) for item in value.split(","))
    if len(values) != 6 or min(values) < 0:
        raise ValueError("--loop-columns must contain six non-negative indexes")
    return values  # type: ignore[return-value]


def build_parser(settings: Optional[Settings] = None) -> argparse.ArgumentParser:
    """Build the command-line parser."""
    settings = settings or load_settings()
    parser = argparse.ArgumentParser(
        prog="peak2anno",
        description="Annotate genomic peaks to nearby genes, genomic features, and chromatin states.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"peak2anno {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    peak2gene = subparsers.add_parser(
        "peak2gene",
        help="Annotate peaks to nearby and closest genes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_input_args(peak2gene, "Input BED/TSV or region-text file.", settings.txt_delimiter)
    add_backend_arg(peak2gene, settings)
    add_db_install_args(peak2gene)
    peak2gene.add_argument("-o", "--output", type=Path, help="Output TSV path; defaults to stdout.")
    peak2gene.add_argument("-s", "--species", default=settings.default_species, help="Species key, for example hg38 or mm10.")
    peak2gene.add_argument(
        "--ver", "--species-version",
        dest="species_version",
        default=settings.version_for(settings.default_species),
        help="Annotation version; def selects the database default.",
    )
    peak2gene.add_argument(
        "--iso", "--isoform-set", "--isoform-version",
        dest="isoform_version",
        choices=["all", "deduplong"],
        default=settings.iso_set,
        help="Gene BED to use: all isoforms or one longest isoform per gene.",
    )
    add_gene_cutoff_args(peak2gene, settings)
    peak2gene.add_argument(
        "-G", "--gene-type",
        default=settings.gene_type,
        help="Gene type selector from BED column 9: all, protein_coding, lincRNA, nomicro, or comma list.",
    )
    peak2gene.add_argument("-d", "--db-path", default=settings.db_path, help="Database root; defaults to rc/env or ~/.sjcab_peak2anno_db.")
    peak2gene.add_argument("-g", "--gene-bed", type=Path, help="Override gene BED; TSS is computed from strand.")
    peak2gene.add_argument("-t", "--tss-bed", type=Path, help="Override TSS BED.")
    peak2gene.add_argument("-H", "--header", choices=["auto", "yes", "no"], default="auto", help="Input header handling.")
    peak2gene.add_argument("-I", "--input-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Input format.")
    peak2gene.add_argument("-C", "--columns", help="BED columns as comma-separated zero-based indexes, for example 0,1,2.")
    peak2gene.add_argument("-R", "--region-column", type=int, default=0, help="Zero-based region column for txt input.")
    peak2gene.add_argument("-f", "--output-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Output format; auto follows input format.")

    narrow = subparsers.add_parser(
        "narrow2feature",
        help="Assign one priority-ordered genomic feature label.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_common_feature_args(narrow, settings)
    add_db_install_args(narrow)
    add_backend_arg(narrow, settings)
    narrow.add_argument("--column-name", default="FeatureAssignment", help="Output annotation column name.")
    add_output_mode_arg(narrow, settings.feature_out, "Output max assignment, ordered percentages, or both.")

    broad = subparsers.add_parser(
        "broad2feature",
        help="Report per-feature genomic feature overlap fractions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_common_feature_args(broad, settings)
    add_db_install_args(broad)
    add_backend_arg(broad, settings)
    add_output_mode_arg(broad, settings.feature_out, "Output max assignment, ordered percentages, or both.")

    state = subparsers.add_parser(
        "peak2state",
        help="Report per-state chromatin overlap fractions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_input_args(state, "Input BED/TSV or region-text file.", settings.txt_delimiter)
    add_backend_arg(state, settings)
    add_db_install_args(state)
    state.add_argument("-s", "-S", "--states", type=Path, required=True, help="Chromatin state dense/segments BED.")
    state.add_argument("-o", "--output", type=Path, help="Output TSV path; defaults to stdout.")
    state.add_argument("-N", "--state2name", type=Path, help="Optional two-column state ID to label mapping.")
    state.add_argument(
        "-x", "--overlap-cutoff",
        default="1bp",
        help="Minimum overlap. Examples: 1bp, 10bp, 0.1 for 10%% of peak, 10%%.",
    )
    state.add_argument("-H", "--header", choices=["auto", "yes", "no"], default="auto", help="Input header handling.")
    state.add_argument("-I", "--input-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Input format.")
    state.add_argument("-C", "--columns", help="BED columns as comma-separated zero-based indexes, for example 0,1,2.")
    state.add_argument("-R", "--region-column", type=int, default=0, help="Zero-based region column for txt input.")
    state.add_argument("-f", "--output-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Output format; auto follows input format.")
    add_output_mode_arg(state, settings.state_out, "Output max state, percentages, or both.")

    for name, help_text in (
        ("loop2gene", "Annotate both anchors of BEDPE loops to genes."),
        ("loop2feature", "Annotate both anchors of BEDPE loops to genomic features."),
        ("loop2state", "Annotate both anchors of BEDPE loops to chromatin states."),
    ):
        loop = subparsers.add_parser(name, help=help_text, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        add_input_args(loop, "BEDPE input.", settings.txt_delimiter)
        add_backend_arg(loop, settings)
        add_db_install_args(loop)
        loop.add_argument("-o", "--output", type=Path, help="Output path; defaults to stdout.")
        loop.add_argument("-H", "--header", choices=["auto", "yes", "no"], default="auto", help="Input header handling.")
        loop.add_argument("-C", "--loop-columns", default="0,1,2,3,4,5", help="BEDPE coordinate columns.")
        loop.add_argument("-f", "--output-format", choices=["auto", "bedpe", "bed", "txt", "txtnohead"], default="auto", help="Output format; auto follows BEDPE input.")
        loop.add_argument("-s", "--species", default=settings.default_species, help="Species key.")
        loop.add_argument("-d", "--db-path", default=settings.db_path, help="Database root.")
        loop.add_argument("--ver", "--species-version", dest="species_version", default=settings.version_for(settings.default_species), help="Annotation version.")
        loop.add_argument("--iso", "--isoform-set", "--isoform-version", dest="isoform_version", choices=["all", "deduplong"], default=settings.iso_set, help="Isoform set.")
        loop.add_argument("--gene-type", default=settings.gene_type, help="Gene type selector for loop2gene.")
        loop.add_argument("-t", "--tss-bed", type=Path, help="Override TSS BED.")
        loop.add_argument("-g", "--gene-bed", type=Path, help="Override gene BED.")
        loop.add_argument("-c", "--feature-dir", type=Path, help="Feature BED directory.")
        loop.add_argument("-O", "--order-lst", default="def", help="Feature priority list; inspect order.lst for feature meanings. def/default/none uses DB_PATH/order.lst, utr uses DB_PATH/order.utr.lst.")
        loop.add_argument("-B", "--feature-mode", choices=["narrow", "broad"], default="narrow", help="Feature mode.")
        loop.add_argument("-x", "--overlap-cutoff", default="1bp", help="Minimum overlap.")
        loop.add_argument("-S", "--states", type=Path, help="Chromatin state BED.")
        loop.add_argument("-N", "--state2name", type=Path, help="State ID/name map.")
        add_gene_cutoff_args(loop, settings)
        add_output_mode_arg(loop, settings.state_out if name == "loop2state" else settings.feature_out, "Output max assignment, ordered percentages, or both.")

    combined = subparsers.add_parser("combined", help="Run multiple annotations and merge their columns.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    combined.add_argument("--commands", action="append", choices=["peak2gene", "narrow2feature", "broad2feature", "peak2state"], required=True, help="Annotation step; repeat for multiple steps.")
    add_input_args(combined, "Input BED/TSV or region-text file.", settings.txt_delimiter)
    add_backend_arg(combined, settings)
    add_db_install_args(combined)
    combined.add_argument("-o", "--output", type=Path, help="Output path; defaults to stdout.")
    combined.add_argument("-s", "--species", default=settings.default_species, help="Species key.")
    combined.add_argument("-d", "--db-path", default=settings.db_path, help="Database root.")
    combined.add_argument("--ver", "--species-version", dest="species_version", default=settings.version_for(settings.default_species), help="Annotation version.")
    combined.add_argument("--iso", "--isoform-set", "--isoform-version", dest="isoform_version", choices=["all", "deduplong"], default=settings.iso_set, help="Isoform set.")
    combined.add_argument("-t", "--tss-bed", type=Path, help="Override TSS BED.")
    combined.add_argument("-g", "--gene-bed", type=Path, help="Override gene BED.")
    combined.add_argument("-c", "--feature-dir", type=Path, help="Feature BED directory.")
    combined.add_argument("-O", "--order-lst", default="def", help="Feature priority list; inspect order.lst for feature meanings. def/default/none uses DB_PATH/order.lst, utr uses DB_PATH/order.utr.lst.")
    combined.add_argument("--features", help="Feature BEDs.")
    combined.add_argument("--feature-labels", help="Feature labels.")
    combined.add_argument("--gene-type", default=settings.gene_type, help="Gene type filter.")
    combined.add_argument("-S", "--states", type=Path, help="Chromatin state BED.")
    combined.add_argument("-N", "--state2name", type=Path, help="State ID/name map.")
    add_gene_cutoff_args(combined, settings)
    combined.add_argument("--column-name", default="FeatureAssignment", help="Feature output column name.")
    combined.add_argument("-x", "--overlap-cutoff", default="1bp", help="Minimum overlap.")
    combined.add_argument("-H", "--header", choices=["auto", "yes", "no"], default="auto", help="Input header handling.")
    combined.add_argument("-I", "--input-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto")
    combined.add_argument("-C", "--columns")
    combined.add_argument("-R", "--region-column", type=int, default=0)
    combined.add_argument("-f", "--output-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Output format; auto follows input format.")
    combined.add_argument("--workers", type=int, default=1, help="Processes for independent annotations.")
    combined.add_argument("--feature-output-mode", choices=["max", "percent", "max,percent", "percent,max"], default=settings.feature_out, help="Output mode for feature steps.")
    combined.add_argument("--state-output-mode", choices=["max", "percent", "max,percent", "percent,max"], default=settings.state_out, help="Output mode for state steps.")
    state.add_argument("-m", "--summary", type=Path, help="Summary TSV path.")
    state.add_argument("-p", "--plot", action="store_true", help="Write PNG/PDF bar and pie plots for summary counts.")

    versions = subparsers.add_parser(
        "list-db",
        help="List available sjcab_peak2anno_db annotations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    versions.add_argument("-d", "--db-path", default=settings.db_path, help="sjcab_peak2anno_db root.")
    return parser


def run(args: argparse.Namespace) -> Optional[Path]:
    """Dispatch parsed CLI arguments."""
    if args.command == "peak2gene":
        return annotate_peak2gene(
            PeakGeneConfig(
                input_path=args.input,
                output_path=args.output,
                species=args.species,
                species_version=args.species_version,
                isoform_version=args.isoform_version,
                prom_enha_cutoffs=args.prom_enha_cutoffs,
                gene_type=args.gene_type,
                db_path=args.db_path,
                gene_bed=args.gene_bed,
                tss_bed=args.tss_bed,
                header=args.header,
                input_format=args.input_format,
                columns=parse_columns(args.columns),
                region_column=args.region_column,
                output_format=args.output_format,
                txt_delimiter=args.txt_delimiter,
                backend=args.backend,
            )
        )
    if args.command == "narrow2feature":
        output, _summary = annotate_narrow_feature(
            FeatureConfig(
                input_path=args.input,
                output_path=args.output,
                species=args.species,
                db_path=args.db_path,
                feature_dir=args.feature_dir,
                features=args.features,
                feature_labels=args.feature_labels,
                overlap_cutoff=args.overlap_cutoff,
                header=args.header,
                column_name=args.column_name,
                summary_path=args.summary,
                plot=args.plot,
                input_format=args.input_format,
                columns=parse_columns(args.columns),
                region_column=args.region_column,
                output_format=args.output_format,
                output_mode=args.output_mode,
                txt_delimiter=args.txt_delimiter,
                backend=args.backend,
                order_lst=args.order_lst,
            )
        )
        return output
    if args.command == "broad2feature":
        output, _summary = annotate_broad_feature(
            FeatureConfig(
                input_path=args.input,
                output_path=args.output,
                species=args.species,
                db_path=args.db_path,
                feature_dir=args.feature_dir,
                features=args.features,
                feature_labels=args.feature_labels,
                overlap_cutoff=args.overlap_cutoff,
                header=args.header,
                summary_path=args.summary,
                plot=args.plot,
                input_format=args.input_format,
                columns=parse_columns(args.columns),
                region_column=args.region_column,
                output_format=args.output_format,
                output_mode=args.output_mode,
                txt_delimiter=args.txt_delimiter,
                backend=args.backend,
                order_lst=args.order_lst,
            )
        )
        return output
    if args.command == "peak2state":
        output, _summary = annotate_peak_state(
            StateConfig(
                input_path=args.input,
                states_path=args.states,
                output_path=args.output,
                state2name=args.state2name,
                overlap_cutoff=args.overlap_cutoff,
                header=args.header,
                summary_path=args.summary,
                plot=args.plot,
                input_format=args.input_format,
                columns=parse_columns(args.columns),
                region_column=args.region_column,
                output_format=args.output_format,
                output_mode=args.output_mode,
                txt_delimiter=args.txt_delimiter,
                backend=args.backend,
            )
        )
        return output
    if args.command in {"loop2gene", "loop2feature", "loop2state"}:
        if args.command == "loop2state" and args.states is None:
            raise ValueError("loop2state requires --states")
        return annotate_loop(
            args.command,
            args.input,
            args.output,
            args.output_format,
            header=args.header,
            loop_columns=parse_loop_columns(args.loop_columns),
            species=args.species,
            species_version=args.species_version,
            isoform_version=args.isoform_version,
            db_path=args.db_path,
            tss_bed=args.tss_bed,
            gene_bed=args.gene_bed,
            gene_type=args.gene_type,
            feature_dir=args.feature_dir,
            feature_mode=args.feature_mode,
            overlap_cutoff=args.overlap_cutoff,
            prom_enha_cutoffs=args.prom_enha_cutoffs,
            output_mode=args.output_mode,
            states=args.states,
            state2name=args.state2name,
            backend=args.backend,
            order_lst=args.order_lst,
        )
    if args.command == "combined":
        return run_combined(args)
    if args.command == "list-db":
        rows = available_versions(args.db_path)
        print("species\tannotation\tversion\tdefault\tpath")
        for row in rows:
            print(
                "\t".join(
                    [
                        str(row.get("species", ".")),
                        str(row.get("annotation", ".")),
                        str(row.get("version", ".")),
                        str(row.get("default", ".")),
                        str(row.get("path", ".")),
                    ]
                )
            )
        return None
    raise ValueError(f"Unknown command: {args.command}")


def run_combined(args: argparse.Namespace) -> Optional[Path]:
    """Run multiple single-region annotations and merge their columns."""
    input_header, _regions = read_regions(
        args.input,
        header=args.header,
        input_format=args.input_format,
        columns=parse_columns(args.columns),
        region_column=args.region_column,
        txt_delimiter=args.txt_delimiter,
    )
    with tempfile.TemporaryDirectory(prefix="peak2anno-combined-") as temp:
        temp_root = Path(temp)
        step_args = []
        for index, command in enumerate(args.commands):
            values = vars(args).copy()
            values.update(
                command=command,
                output=temp_root / f"{index}.tsv",
                output_format="txt",
                summary=None,
                plot=False,
                features=None,
                feature_labels=None,
                column_name="FeatureAssignment",
                output_mode=(args.state_output_mode if command == "peak2state" else args.feature_output_mode),
            )
            run_args = argparse.Namespace(**values)
            if command == "peak2state" and args.states is None:
                raise ValueError("combined peak2state requires --states")
            step_args.append(run_args)
        if args.workers > 1 and len(step_args) > 1:
            with ProcessPoolExecutor(max_workers=args.workers) as executor:
                tables = list(executor.map(run_combined_step, step_args))
        else:
            tables = [run_combined_step(step_args_item) for step_args_item in step_args]
        output_header = list(input_header)
        for header, _rows in tables:
            output_header.extend(header[len(input_header):])
        output_rows = []
        for row_index in range(len(tables[0][1])):
            row = list(tables[0][1][row_index])
            for _header, rows in tables[1:]:
                row.extend(rows[row_index][len(input_header):])
            output_rows.append(row)
        output_format = args.output_format if args.output_format != "auto" else detect_output_format(args.input, args.header, args.input_format, parse_columns(args.columns), args.region_column, args.txt_delimiter)
        write_table(args.output, output_header, output_rows, include_header=output_format == "txt")
    return args.output


def run_combined_step(args: argparse.Namespace) -> tuple[list[str], list[list[str]]]:
    """Run one combined annotation step in the current or a worker process."""
    run(args)
    with args.output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return rows[0], rows[1:]


def normalize_argv(argv: Optional[Sequence[str]]) -> Optional[List[str]]:
    """Accept ``peak2anno peak2gene narrow2feature ...`` as combined syntax."""
    values = list(sys.argv[1:] if argv is None else argv)
    commands = {"peak2gene", "narrow2feature", "broad2feature", "peak2state"}
    if len(values) >= 2 and values[0] in commands and values[1] in commands:
        index = 0
        selected = []
        while index < len(values) and values[index] in commands:
            selected.append(values[index])
            index += 1
        normalized = ["combined"]
        for command in selected:
            normalized.extend(["--commands", command])
        return normalized + values[index:]
    return values


def _has_option(argv: Sequence[str], *options: str) -> bool:
    """Return whether one of the options was explicitly supplied."""
    return any(value == option or value.startswith(option + "=") for value in argv for option in options)


def apply_settings(args: argparse.Namespace, argv: Sequence[str], settings: Settings) -> None:
    """Apply species-aware config defaults while preserving explicit CLI values."""
    if hasattr(args, "species") and not _has_option(argv, "-s", "--species"):
        args.species = settings.default_species
    if hasattr(args, "species_version") and not _has_option(argv, "--ver", "--species-version"):
        args.species_version = settings.version_for(args.species)
    if hasattr(args, "prom_enha_cutoffs") and not _has_option(argv, "--prom-enha-cutoffs"):
        args.prom_enha_cutoffs = settings.prom_enha_cutoffs_for(args.species, args.species_version)


def resolved_references(args: argparse.Namespace) -> list[str]:
    """Return input and reference files resolved for a parsed command."""
    references = [f"input: {args.input.expanduser().resolve()}"]
    if args.command == "peak2gene":
        config = PeakGeneConfig(
            input_path=args.input,
            output_path=args.output,
            species=args.species,
            species_version=args.species_version,
            isoform_version=args.isoform_version,
            db_path=args.db_path,
            gene_bed=args.gene_bed,
            tss_bed=args.tss_bed,
        )
        references.append(f"gene/TSS: {resolve_tss_path(config).expanduser().resolve()}")
    elif args.command in {"narrow2feature", "broad2feature"}:
        config = FeatureConfig(
            input_path=args.input,
            output_path=args.output,
            species=args.species,
            db_path=args.db_path,
            feature_dir=args.feature_dir,
            features=args.features,
            feature_labels=args.feature_labels,
            order_lst=args.order_lst,
        )
        specs = resolve_features(config)
        for spec in specs:
            references.append(f"feature {spec.label}: {spec.path.expanduser().resolve()}")
        feature_dir = specs[0].path.parent if specs else (config.feature_dir or config.input_path.parent)
        order_path = resolve_order_path(config, feature_dir)
        if order_path is not None:
            references.append(f"feature order: {order_path.expanduser().resolve()}")
    elif args.command == "peak2state":
        references.append(f"states: {args.states.expanduser().resolve()}")
        if args.state2name is not None:
            references.append(f"state names: {args.state2name.expanduser().resolve()}")
    elif args.command == "list-db":
        references.append(f"database root: {db_root(args.db_path).resolve()}")
        references.append(f"manifest: {(db_root(args.db_path) / 'manifest.json').resolve()}")
    return references


def record_run(args: argparse.Namespace, command_args: Optional[Sequence[str]] = None) -> None:
    """Append the command line and resolved reference files to ``.run.log``."""
    actual_args = sys.argv[1:] if command_args is None else command_args
    command = " ".join(shlex.quote(value) for value in [sys.argv[0], *actual_args])
    references = []
    try:
        references = resolved_references(args)
    except Exception as exc:  # noqa: BLE001
        references.append(f"resolution error: {exc}")
    output = getattr(args, "output", None)
    with Path(".run.log").open("a", encoding="utf-8") as handle:
        handle.write(f"[{datetime.now().isoformat(timespec='seconds')}]\n")
        handle.write(f"command: {command}\n")
        if hasattr(args, "backend"):
            handle.write(f"backend: {args.backend}\n")
        handle.write(f"output: {output.expanduser().resolve() if output else 'stdout'}\n")
        for reference in references:
            handle.write(f"reference: {reference}\n")
        handle.write("\n")


def database_install_command(args: argparse.Namespace) -> Optional[list[str]]:
    """Return the database-package command that can install missing references."""
    species = getattr(args, "species", "hg38")
    version = getattr(args, "species_version", "def")
    db_path = getattr(args, "db_path", None)
    data_dir = str(Path(db_path).expanduser()) if db_path else None

    commands = getattr(args, "commands", [getattr(args, "command", "")])
    if isinstance(commands, str):
        commands = [commands]
    needs_gene = "peak2gene" in commands or "loop2gene" in commands
    needs_feature = any(command in {"narrow2feature", "broad2feature", "loop2feature"} for command in commands)
    has_gene_override = getattr(args, "gene_bed", None) is not None or getattr(args, "tss_bed", None) is not None
    has_feature_override = getattr(args, "feature_dir", None) is not None or getattr(args, "features", None) is not None

    if needs_gene and not has_gene_override:
        if version in {"def", "default"}:
            command = ["sjcab-peak2anno-db", "install-gencode-bed"]
            if data_dir:
                command.extend(["-d", data_dir])
        else:
            command = ["sjcab-peak2anno-db", "download-gencode-bed", str(species), str(version)]
            if data_dir:
                command.extend(["-o", data_dir])
        return command
    if needs_feature and not has_feature_override:
        # Feature files are generated by the database package's feature command.
        if version in {"def", "default"}:
            command = ["sjcab-peak2anno-db", "install-gencode-feature"]
            if data_dir:
                command.extend(["-d", data_dir])
        else:
            command = ["sjcab-peak2anno-db", "download-gencode-feature", str(species), str(version)]
            if data_dir:
                command.extend(["-o", data_dir])
        return command
    return None


def database_references_present(args: argparse.Namespace) -> bool:
    """Return whether all implicit gene and feature references now exist."""
    commands = getattr(args, "commands", [getattr(args, "command", "")])
    if isinstance(commands, str):
        commands = [commands]
    needs_gene = "peak2gene" in commands or "loop2gene" in commands
    needs_feature = any(command in {"narrow2feature", "broad2feature", "loop2feature"} for command in commands)
    if needs_gene and getattr(args, "gene_bed", None) is None and getattr(args, "tss_bed", None) is None:
        gene_annotation_path(
            str(args.species),
            version=str(args.species_version),
            isoform_set=str(args.isoform_version),
            root_path=getattr(args, "db_path", None),
        )
    if needs_feature and getattr(args, "feature_dir", None) is None and getattr(args, "features", None) is None:
        resolve_features(
            FeatureConfig(
                input_path=args.input,
                output_path=getattr(args, "output", None),
                species=str(args.species),
                db_path=getattr(args, "db_path", None),
            )
        )
    return True


def database_lock_path(args: argparse.Namespace) -> Path:
    """Return the per-species/version installation lock path."""
    root = db_root(getattr(args, "db_path", None))
    species = str(getattr(args, "species", "unknown"))
    version = str(getattr(args, "species_version", "def"))
    safe_key = "".join(char if char.isalnum() or char in "._-" else "_" for char in f"{species}-{version}")
    return root / ".locks" / f"{safe_key}.lock"


@contextmanager
def database_install_lock(args: argparse.Namespace):
    """Serialize automatic installation for one species/version pair."""
    path = database_lock_path(args)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def ensure_database_cli() -> str:
    """Return the DB CLI, installing the PyPI package when necessary."""
    command = shutil.which("sjcab-peak2anno-db")
    if command:
        return command
    subprocess.run([sys.executable, "-m", "pip", "install", DB_PACKAGE], check=True)
    command = shutil.which("sjcab-peak2anno-db")
    if not command:
        raise FileNotFoundError("sjcab-peak2anno-db was not found after installing sjcab_peak2anno_db")
    return command


def offer_database_install(args: argparse.Namespace, missing: Exception) -> bool:
    """Offer to install missing database files, returning whether to retry."""
    command = database_install_command(args)
    if command is None:
        return False
    # ``shlex.join`` is only available from Python 3.8; keep the CLI usable on
    # the package's Python 3.7 baseline.
    command_text = " ".join(shlex.quote(value) for value in command)
    print(f"peak2anno: selected database reference was not found: {missing}", file=sys.stderr)
    print("You can install it with:", file=sys.stderr)
    print(f"  {command_text}", file=sys.stderr)
    automatic = bool(getattr(args, "auto_install_db", False))
    if not automatic and not sys.stdin.isatty():
        print("peak2anno: non-interactive input; installation declined", file=sys.stderr)
        return False
    if not automatic:
        try:
            answer = input("Install database files now? [y/N]: ").strip().lower()
        except EOFError:
            return False
        if answer not in {"y", "yes"}:
            print("peak2anno: installation declined", file=sys.stderr)
            return False
    try:
        with database_install_lock(args):
            try:
                database_references_present(args)
                print("peak2anno: database references became available while waiting for the installer", file=sys.stderr)
                return True
            except FileNotFoundError:
                db_cli = ensure_database_cli()
                subprocess.run([db_cli, *command[1:]], check=True)
    except (OSError, subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"peak2anno: database installation failed: {exc}", file=sys.stderr)
        return False
    return True


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the peak2anno command-line interface."""
    settings = load_settings()
    parser = build_parser(settings)
    normalized = normalize_argv(argv)
    args = parser.parse_args(normalized)
    apply_settings(args, normalized or [], settings)
    if hasattr(args, "backend"):
        tool_status = detect_tools()
        try:
            args.backend = resolve_backend(args.backend, tool_status)
        except ValueError as exc:
            parser.error(str(exc))
        warn_if_slow(tool_status, args.backend)
    if args.command != "list-db":
        if getattr(args, "input", None) is None:
            args.input = getattr(args, "input_option", None)
        elif getattr(args, "input_option", None) is not None:
            parser.error("provide input either positionally or with -i/--input, not both")
        if args.input is None:
            parser.error("an input file is required; provide it positionally or with -i/--input")
    record_run(args, normalized)
    try:
        run(args)
    except FileNotFoundError as exc:
        if offer_database_install(args, exc):
            record_run(args, normalized)
            try:
                run(args)
            except Exception as retry_exc:  # noqa: BLE001
                parser.exit(1, f"peak2anno: error after database installation: {retry_exc}\n")
        else:
            parser.exit(1, f"peak2anno: error: {exc}\n")
    except Exception as exc:  # noqa: BLE001
        parser.exit(1, f"peak2anno: error: {exc}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
