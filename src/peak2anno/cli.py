"""Console entry point for peak2anno."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
import shlex
import sys
import tempfile
from pathlib import Path
from typing import Optional, Sequence

from . import __version__
from .context import (
    ContextConfig,
    StateConfig,
    annotate_broad_context,
    annotate_narrow_context,
    annotate_peak_state,
    resolve_features,
)
from .db import available_versions, db_root
from .peak2gene import PeakGeneConfig, annotate_peak2gene, resolve_tss_path
from .loops import annotate_loop
from .intervals import detect_output_format, read_regions, write_table


def add_common_context_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments shared by narrow2context and broad2context."""
    add_input_args(parser, "Input BED/TSV or region-text file.")
    parser.add_argument("-o", "--output", type=Path, help="Output TSV path; defaults to stdout.")
    parser.add_argument("-s", "--species", default="hg38", help="Species key, for example hg38 or mm10.")
    parser.add_argument("-d", "--db-path", help="Database root; defaults to $SJCAB_PEAK2ANNO_DB_PATH or ~/.sjcab_peak2anno_db.")
    parser.add_argument(
        "-c", "--context-dir", type=Path,
        help="Context BED directory. If omitted, search --db-path for the species context files.",
    )
    parser.add_argument("--input-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Input format.")
    parser.add_argument("--columns", help="BED columns as comma-separated zero-based indexes, for example 0,1,2.")
    parser.add_argument("--region-column", type=int, default=0, help="Zero-based region column for txt input.")
    parser.add_argument("--output-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Output format; auto follows input format.")
    parser.add_argument("--features", help="Comma-separated feature BEDs or .lst file.")
    parser.add_argument("--feature-labels", help="Comma-separated feature labels or .labels.lst file.")
    parser.add_argument(
        "--overlap-cutoff",
        default="1bp",
        help="Minimum overlap. Examples: 1bp, 10bp, 0.1 for 10%% of peak, 10%%.",
    )
    parser.add_argument("--header", choices=["auto", "yes", "no"], default="auto", help="Input header handling.")
    parser.add_argument("--summary", type=Path, help="Summary TSV path.")
    parser.add_argument("--plot", action="store_true", help="Write PNG/PDF bar and pie plots for summary counts.")


def add_input_args(parser: argparse.ArgumentParser, help_text: str) -> None:
    """Add positional and short-option input forms."""
    parser.add_argument("input", nargs="?", type=Path, help=help_text)
    parser.add_argument("-i", "--input", dest="input_option", type=Path, metavar="INPUT", help="Input file.")


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


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="peak2anno",
        description="Annotate genomic peaks to nearby genes, genomic contexts, and chromatin states.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"peak2anno {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    peak2gene = subparsers.add_parser(
        "peak2gene",
        help="Annotate peaks to nearby and closest genes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_input_args(peak2gene, "Input BED/TSV or region-text file.")
    peak2gene.add_argument("-o", "--output", type=Path, help="Output TSV path; defaults to stdout.")
    peak2gene.add_argument("-s", "--species", default="hg38", help="Species key, for example hg38 or mm10.")
    peak2gene.add_argument(
        "--ver", "--species-version",
        dest="species_version",
        default="def",
        help="Annotation version; def selects the database default.",
    )
    peak2gene.add_argument(
        "--iso", "--isoform-set", "--isoform-version",
        dest="isoform_version",
        choices=["all", "deduplong"],
        default="all",
        help="Gene BED to use: all isoforms or one longest isoform per gene.",
    )
    peak2gene.add_argument("--promoter-cutoff", default="2kb", help="Promoter TSS distance cutoff.")
    peak2gene.add_argument("--enhancer-cutoff", default="50kb", help="Putative enhancer outer TSS distance cutoff.")
    peak2gene.add_argument(
        "--gene-type",
        default="all",
        help="Gene type selector from BED column 9: all, protein_coding, lincRNA, nomicro, or comma list.",
    )
    peak2gene.add_argument("-d", "--db-path", help="Database root; defaults to $SJCAB_PEAK2ANNO_DB_PATH or ~/.sjcab_peak2anno_db.")
    peak2gene.add_argument("--gene-bed", type=Path, help="Override gene BED; TSS is computed from strand.")
    peak2gene.add_argument("--tss-bed", type=Path, help="Override TSS BED.")
    peak2gene.add_argument("--header", choices=["auto", "yes", "no"], default="auto", help="Input header handling.")
    peak2gene.add_argument("--input-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Input format.")
    peak2gene.add_argument("--columns", help="BED columns as comma-separated zero-based indexes, for example 0,1,2.")
    peak2gene.add_argument("--region-column", type=int, default=0, help="Zero-based region column for txt input.")
    peak2gene.add_argument("--output-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Output format; auto follows input format.")

    narrow = subparsers.add_parser(
        "narrow2context",
        help="Assign one priority-ordered genomic context label.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_common_context_args(narrow)
    narrow.add_argument("--column-name", default="FeatureAssignment", help="Output annotation column name.")

    broad = subparsers.add_parser(
        "broad2context",
        help="Report per-feature genomic context overlap fractions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_common_context_args(broad)

    state = subparsers.add_parser(
        "peak2state",
        help="Report per-state chromatin overlap fractions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_input_args(state, "Input BED/TSV or region-text file.")
    state.add_argument("-s", "--states", type=Path, required=True, help="Chromatin state dense/segments BED.")
    state.add_argument("-o", "--output", type=Path, help="Output TSV path; defaults to stdout.")
    state.add_argument("--state2name", type=Path, help="Optional two-column state ID to label mapping.")
    state.add_argument(
        "--overlap-cutoff",
        default="1bp",
        help="Minimum overlap. Examples: 1bp, 10bp, 0.1 for 10%% of peak, 10%%.",
    )
    state.add_argument("--header", choices=["auto", "yes", "no"], default="auto", help="Input header handling.")
    state.add_argument("--input-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Input format.")
    state.add_argument("--columns", help="BED columns as comma-separated zero-based indexes, for example 0,1,2.")
    state.add_argument("--region-column", type=int, default=0, help="Zero-based region column for txt input.")
    state.add_argument("--output-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Output format; auto follows input format.")

    for name, help_text in (
        ("loop2anno", "Annotate both anchors of BEDPE loops to genes."),
        ("loop2context", "Annotate both anchors of BEDPE loops to genomic contexts."),
        ("loop2state", "Annotate both anchors of BEDPE loops to chromatin states."),
    ):
        loop = subparsers.add_parser(name, help=help_text, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        add_input_args(loop, "BEDPE input.")
        loop.add_argument("-o", "--output", type=Path, help="Output path; defaults to stdout.")
        loop.add_argument("--header", choices=["auto", "yes", "no"], default="auto", help="Input header handling.")
        loop.add_argument("--loop-columns", default="0,1,2,3,4,5", help="BEDPE coordinate columns.")
        loop.add_argument("--output-format", choices=["auto", "bedpe", "bed", "txt", "txtnohead"], default="auto", help="Output format; auto follows BEDPE input.")
        loop.add_argument("-s", "--species", default="hg38", help="Species key.")
        loop.add_argument("-d", "--db-path", help="Database root.")
        loop.add_argument("--ver", "--species-version", dest="species_version", default="def", help="Annotation version.")
        loop.add_argument("--iso", "--isoform-set", "--isoform-version", dest="isoform_version", choices=["all", "deduplong"], default="all", help="Isoform set.")
        loop.add_argument("--tss-bed", type=Path, help="Override TSS BED.")
        loop.add_argument("--gene-bed", type=Path, help="Override gene BED.")
        loop.add_argument("-c", "--context-dir", type=Path, help="Context BED directory.")
        loop.add_argument("--context-mode", choices=["narrow", "broad"], default="narrow", help="Context mode.")
        loop.add_argument("--overlap-cutoff", default="1bp", help="Minimum overlap.")
        loop.add_argument("--states", type=Path, help="Chromatin state BED.")
        loop.add_argument("--state2name", type=Path, help="State ID/name map.")

    combined = subparsers.add_parser("combined", help="Run multiple annotations and merge their columns.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    combined.add_argument("--commands", action="append", choices=["peak2gene", "narrow2context", "broad2context", "peak2state"], required=True, help="Annotation step; repeat for multiple steps.")
    add_input_args(combined, "Input BED/TSV or region-text file.")
    combined.add_argument("-o", "--output", type=Path, help="Output path; defaults to stdout.")
    combined.add_argument("-s", "--species", default="hg38", help="Species key.")
    combined.add_argument("-d", "--db-path", help="Database root.")
    combined.add_argument("--ver", "--species-version", dest="species_version", default="def", help="Annotation version.")
    combined.add_argument("--iso", "--isoform-set", "--isoform-version", dest="isoform_version", choices=["all", "deduplong"], default="all", help="Isoform set.")
    combined.add_argument("--tss-bed", type=Path, help="Override TSS BED.")
    combined.add_argument("--gene-bed", type=Path, help="Override gene BED.")
    combined.add_argument("-c", "--context-dir", type=Path, help="Context BED directory.")
    combined.add_argument("--features", help="Context feature BEDs.")
    combined.add_argument("--feature-labels", help="Context feature labels.")
    combined.add_argument("--gene-type", default="all", help="Gene type filter.")
    combined.add_argument("--states", type=Path, help="Chromatin state BED.")
    combined.add_argument("--state2name", type=Path, help="State ID/name map.")
    combined.add_argument("--promoter-cutoff", default="2kb", help="Promoter TSS cutoff.")
    combined.add_argument("--enhancer-cutoff", default="50kb", help="Enhancer outer cutoff.")
    combined.add_argument("--column-name", default="FeatureAssignment", help="Context output column name.")
    combined.add_argument("--overlap-cutoff", default="1bp", help="Minimum overlap.")
    combined.add_argument("--header", choices=["auto", "yes", "no"], default="auto", help="Input header handling.")
    combined.add_argument("--input-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto")
    combined.add_argument("--columns")
    combined.add_argument("--region-column", type=int, default=0)
    combined.add_argument("--output-format", choices=["auto", "bed", "txt", "txtnohead"], default="auto", help="Output format; auto follows input format.")
    combined.add_argument("--workers", type=int, default=1, help="Processes for independent annotations.")
    state.add_argument("--summary", type=Path, help="Summary TSV path.")
    state.add_argument("--plot", action="store_true", help="Write PNG/PDF bar and pie plots for summary counts.")

    versions = subparsers.add_parser(
        "list-db",
        help="List available sjcab_peak2anno_db annotations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    versions.add_argument("--db-path", help="sjcab_peak2anno_db root.")
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
                promoter_cutoff=args.promoter_cutoff,
                enhancer_cutoff=args.enhancer_cutoff,
                gene_type=args.gene_type,
                db_path=args.db_path,
                gene_bed=args.gene_bed,
                tss_bed=args.tss_bed,
                header=args.header,
                input_format=args.input_format,
                columns=parse_columns(args.columns),
                region_column=args.region_column,
                output_format=args.output_format,
            )
        )
    if args.command == "narrow2context":
        output, _summary = annotate_narrow_context(
            ContextConfig(
                input_path=args.input,
                output_path=args.output,
                species=args.species,
                db_path=args.db_path,
                context_dir=args.context_dir,
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
            )
        )
        return output
    if args.command == "broad2context":
        output, _summary = annotate_broad_context(
            ContextConfig(
                input_path=args.input,
                output_path=args.output,
                species=args.species,
                db_path=args.db_path,
                context_dir=args.context_dir,
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
            )
        )
        return output
    if args.command in {"loop2anno", "loop2context", "loop2state"}:
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
            context_dir=args.context_dir,
            context_mode=args.context_mode,
            overlap_cutoff=args.overlap_cutoff,
            states=args.states,
            state2name=args.state2name,
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
        output_format = args.output_format if args.output_format != "auto" else detect_output_format(args.input, args.header, args.input_format, parse_columns(args.columns), args.region_column)
        write_table(args.output, output_header, output_rows, include_header=output_format == "txt")
    return args.output


def run_combined_step(args: argparse.Namespace) -> tuple[list[str], list[list[str]]]:
    """Run one combined annotation step in the current or a worker process."""
    run(args)
    with args.output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return rows[0], rows[1:]


def normalize_argv(argv: Optional[Sequence[str]]) -> Optional[List[str]]:
    """Accept ``peak2anno peak2gene narrow2context ...`` as combined syntax."""
    values = list(sys.argv[1:] if argv is None else argv)
    commands = {"peak2gene", "narrow2context", "broad2context", "peak2state"}
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
    elif args.command in {"narrow2context", "broad2context"}:
        config = ContextConfig(
            input_path=args.input,
            output_path=args.output,
            species=args.species,
            db_path=args.db_path,
            context_dir=args.context_dir,
            features=args.features,
            feature_labels=args.feature_labels,
        )
        for spec in resolve_features(config):
            references.append(f"context {spec.label}: {spec.path.expanduser().resolve()}")
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
        handle.write(f"output: {output.expanduser().resolve() if output else 'stdout'}\n")
        for reference in references:
            handle.write(f"reference: {reference}\n")
        handle.write("\n")


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the peak2anno command-line interface."""
    parser = build_parser()
    normalized = normalize_argv(argv)
    args = parser.parse_args(normalized)
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
    except Exception as exc:  # noqa: BLE001
        parser.exit(1, f"peak2anno: error: {exc}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
