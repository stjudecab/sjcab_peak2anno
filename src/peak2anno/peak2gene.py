"""Peak-to-gene annotation command implementation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

from .db import gene_annotation_path
from .intervals import (
    BedRecord,
    InputRegion,
    IntervalIndex,
    distance_bp,
    format_distance,
    parse_distance,
    output_region_values,
    output_region_header,
    detect_output_format,
    read_bed_records,
    read_regions,
    unique_join,
    write_table,
)


@dataclass(frozen=True)
class PeakGeneConfig:
    """Configuration for peak-to-gene annotation."""

    input_path: Path
    output_path: Optional[Path]
    species: str
    species_version: str = "default"
    isoform_version: str = "all"
    promoter_cutoff: str = "2kb"
    enhancer_cutoff: str = "50kb"
    gene_type: str = "all"
    db_path: Optional[str] = None
    gene_bed: Optional[Path] = None
    tss_bed: Optional[Path] = None
    header: str = "auto"
    input_format: str = "auto"
    columns: Optional[Tuple[int, int, int]] = None
    region_column: int = 0
    output_format: str = "txt"


def resolve_tss_records(config: PeakGeneConfig) -> List[BedRecord]:
    """Load TSS records for peak-to-gene annotation."""
    if config.tss_bed is not None and config.gene_bed is not None:
        raise ValueError("Use either --tss-bed or --gene-bed, not both")
    path = resolve_tss_path(config)
    return read_bed_records(path, gene_type=config.gene_type, as_tss=config.tss_bed is None)


def resolve_tss_path(config: PeakGeneConfig) -> Path:
    """Resolve the gene/TSS reference path used by a peak2gene run."""
    if config.tss_bed is not None and config.gene_bed is not None:
        raise ValueError("Use either --tss-bed or --gene-bed, not both")
    if config.tss_bed is not None:
        return config.tss_bed
    if config.gene_bed is not None:
        return config.gene_bed
    return gene_annotation_path(
        config.species,
        version=config.species_version,
        isoform_set=config.isoform_version,
        root_path=config.db_path,
    )


def nearby_records(
    index: IntervalIndex,
    region: InputRegion,
    cutoff_bp: int,
) -> List[Tuple[BedRecord, int]]:
    """Return TSS records within a distance cutoff of a peak."""
    query_start = max(0, region.start - cutoff_bp)
    query_end = region.end + cutoff_bp
    matches: List[Tuple[BedRecord, int]] = []
    for record in index.query(region.chrom, query_start, query_end):
        dist = distance_bp(region.start, region.end, record.start, record.end)
        if dist <= cutoff_bp:
            matches.append((record, dist))
    matches.sort(key=lambda item: (item[1], item[0].start, item[0].name, item[0].gene_id))
    return matches


def closest_record(index: IntervalIndex, region: InputRegion) -> Tuple[Optional[BedRecord], Optional[int]]:
    """Return the closest TSS record and distance for a peak."""
    records = index.records_for_chrom(region.chrom)
    if not records:
        return None, None
    best: Optional[Tuple[BedRecord, int]] = None
    for record in records:
        if best is not None and record.start > region.end + best[1]:
            break
        dist = distance_bp(region.start, region.end, record.start, record.end)
        if best is None or (dist, record.start, record.name, record.gene_id) < (
            best[1],
            best[0].start,
            best[0].name,
            best[0].gene_id,
        ):
            best = (record, dist)
    if best is None:
        return None, None
    return best


def names_and_ids(records: Iterable[BedRecord]) -> Tuple[str, str]:
    """Return comma-separated unique gene symbols and IDs."""
    rows = list(records)
    return unique_join(record.name for record in rows), unique_join(record.gene_id for record in rows)


def annotate_peak2gene(config: PeakGeneConfig) -> Path:
    """Annotate peaks with nearby and closest genes.

    Args:
        config (PeakGeneConfig): Annotation configuration.

    Returns:
        Path: Output TSV path.
    """
    promoter_bp = parse_distance(config.promoter_cutoff)
    enhancer_bp = parse_distance(config.enhancer_cutoff)
    if enhancer_bp < promoter_bp:
        raise ValueError("--enhancer-cutoff must be greater than or equal to --promoter-cutoff")
    output_format = config.output_format if config.output_format != "auto" else detect_output_format(
        config.input_path, config.header, config.input_format, config.columns, config.region_column
    )
    header, regions = read_regions(
        config.input_path,
        header=config.header,
        input_format=config.input_format,
        columns=config.columns,
        region_column=config.region_column,
    )
    tss_records = resolve_tss_records(config)
    index = IntervalIndex(tss_records)

    promoter_label = format_distance(promoter_bp)
    enhancer_label = format_distance(enhancer_bp)
    distal_label = f"{promoter_label}-{enhancer_label}"
    if promoter_label.endswith("kb") and enhancer_label.endswith("kb"):
        distal_label = f"{promoter_label[:-2]}-{enhancer_label}"
    elif promoter_label.endswith("Mb") and enhancer_label.endswith("Mb"):
        distal_label = f"{promoter_label[:-2]}-{enhancer_label}"

    out_header = output_region_header(header, output_format) + [
        f"Gene_{promoter_label}",
        "Gencode_ids",
        f"Gene_{distal_label}",
        "Gencode_ids",
        "Closest_Gene",
        "Gencode_id",
        "Distance",
    ]
    rows: List[Sequence[object]] = []
    for region in regions:
        promoter = nearby_records(index, region, promoter_bp)
        outer = nearby_records(index, region, enhancer_bp)
        distal = [(record, dist) for record, dist in outer if dist > promoter_bp]
        promoter_names, promoter_ids = names_and_ids(record for record, _ in promoter)
        distal_names, distal_ids = names_and_ids(record for record, _ in distal)
        closest, distance = closest_record(index, region)
        rows.append(
            output_region_values(region, region.values, output_format)
            + [
                promoter_names,
                promoter_ids,
                distal_names,
                distal_ids,
                closest.name if closest else ".",
                closest.gene_id if closest else ".",
                distance if distance is not None else ".",
            ]
        )
    write_table(config.output_path, out_header, rows, include_header=output_format == "txt")
    return config.output_path
