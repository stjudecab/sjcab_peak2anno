"""Console entry point for peak2anno."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from . import __version__
from .context import (
    ContextConfig,
    StateConfig,
    annotate_broad_context,
    annotate_narrow_context,
    annotate_peak_state,
)
from .db import available_versions
from .peak2gene import PeakGeneConfig, annotate_peak2gene


def add_common_context_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments shared by narrow2context and broad2context."""
    parser.add_argument("input", type=Path, help="Input BED/TSV with BED columns or a Region column.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Output TSV path.")
    parser.add_argument("-s", "--species", default="hg38", help="Species key, for example hg38 or mm10.")
    parser.add_argument("--db-path", help="sjcab_peak2anno_db root. Defaults to $SJCAB_PEAK2ANNO_DB_PATH or ~/.sjcab_peak2anno_db.")
    parser.add_argument("--context-dir", type=Path, help="Directory containing context BEDs such as 2kb.exon.bed.")
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


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="peak2anno",
        description="Annotate genomic peaks to nearby genes, genomic contexts, and chromatin states.",
    )
    parser.add_argument("--version", action="version", version=f"peak2anno {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    peak2gene = subparsers.add_parser("peak2gene", help="Annotate peaks to nearby and closest genes.")
    peak2gene.add_argument("input", type=Path, help="Input BED/TSV with BED columns or a Region column.")
    peak2gene.add_argument("-o", "--output", type=Path, required=True, help="Output TSV path.")
    peak2gene.add_argument("-s", "--species", default="hg38", help="Species key, for example hg38 or mm10.")
    peak2gene.add_argument(
        "--species-version",
        default="default",
        help="Annotation version, for example v31 or vM22. Use default/def to read manifest default.",
    )
    peak2gene.add_argument(
        "--isoform-version",
        choices=["all", "deduplong"],
        default="all",
        help="Use all TSS isoforms or one longest isoform per gene.",
    )
    peak2gene.add_argument("--promoter-cutoff", default="2kb", help="Promoter TSS distance cutoff.")
    peak2gene.add_argument("--enhancer-cutoff", default="50kb", help="Putative enhancer outer TSS distance cutoff.")
    peak2gene.add_argument(
        "--gene-type",
        default="all",
        help="Gene type selector from BED column 9: all, protein_coding, lincRNA, nomicro, or comma list.",
    )
    peak2gene.add_argument("--db-path", help="sjcab_peak2anno_db root.")
    peak2gene.add_argument("--gene-bed", type=Path, help="Override gene BED; TSS is computed from strand.")
    peak2gene.add_argument("--tss-bed", type=Path, help="Override TSS BED.")
    peak2gene.add_argument("--header", choices=["auto", "yes", "no"], default="auto", help="Input header handling.")

    narrow = subparsers.add_parser("narrow2context", help="Assign one priority-ordered genomic context label.")
    add_common_context_args(narrow)
    narrow.add_argument("--column-name", default="FeatureAssignment", help="Output annotation column name.")

    broad = subparsers.add_parser("broad2context", help="Report per-feature genomic context overlap fractions.")
    add_common_context_args(broad)

    state = subparsers.add_parser("peak2state", help="Report per-state chromatin overlap fractions.")
    state.add_argument("input", type=Path, help="Input BED/TSV with BED columns or a Region column.")
    state.add_argument("-s", "--states", type=Path, required=True, help="Chromatin state dense/segments BED.")
    state.add_argument("-o", "--output", type=Path, required=True, help="Output TSV path.")
    state.add_argument("--state2name", type=Path, help="Optional two-column state ID to label mapping.")
    state.add_argument(
        "--overlap-cutoff",
        default="1bp",
        help="Minimum overlap. Examples: 1bp, 10bp, 0.1 for 10%% of peak, 10%%.",
    )
    state.add_argument("--header", choices=["auto", "yes", "no"], default="auto", help="Input header handling.")
    state.add_argument("--summary", type=Path, help="Summary TSV path.")
    state.add_argument("--plot", action="store_true", help="Write PNG/PDF bar and pie plots for summary counts.")

    versions = subparsers.add_parser("list-db", help="List available sjcab_peak2anno_db annotations.")
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
            )
        )
        return output
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


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the peak2anno command-line interface."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        run(args)
    except Exception as exc:  # noqa: BLE001
        parser.exit(1, f"peak2anno: error: {exc}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
