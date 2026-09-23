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
    prom_enha_cutoffs: str = "2kb,50kb,2kb"
    gene_type: str = "all"
    db_path: Optional[str] = None
    gene_bed: Optional[Path] = None
    tss_bed: Optional[Path] = None
    header: str = "auto"
    input_format: str = "auto"
    columns: Optional[Tuple[int, int, int]] = None
    region_column: int = 0
    output_format: str = "txt"
    txt_delimiter: str = "auto"
    backend: str = "python"


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
    cutoff_bp: object,
) -> List[Tuple[BedRecord, int]]:
    """Return TSS records within a distance cutoff of a peak."""
    fixed_cutoff = _fixed_cutoff(cutoff_bp)
    candidates = (
        index.query_window(region.chrom, region.start, region.end, fixed_cutoff)
        if fixed_cutoff is not None
        else index.records_for_chrom(region.chrom)
    )
    matches: List[Tuple[BedRecord, int]] = []
    for record in candidates:
        dist = distance_bp(region.start, region.end, record.start, record.end)
        if dist <= cutoff_value(cutoff_bp, record):
            matches.append((record, dist))
    matches.sort(key=lambda item: (item[1], item[0].start, item[0].name, item[0].gene_id))
    return matches


def promoter_records(
    index: IntervalIndex,
    region: InputRegion,
    upstream_bp: object,
    downstream_bp: object,
) -> List[BedRecord]:
    """Return records in strand-aware upstream/downstream promoter windows."""
    fixed_upstream = _fixed_cutoff(upstream_bp)
    fixed_downstream = _fixed_cutoff(downstream_bp)
    fixed_window = (
        max(fixed_upstream, fixed_downstream)
        if fixed_upstream is not None and fixed_downstream is not None
        else None
    )
    candidates = (
        index.query_window(region.chrom, region.start, region.end, fixed_window)
        if fixed_window is not None
        else index.records_for_chrom(region.chrom)
    )
    matches: List[Tuple[BedRecord, int]] = []
    for record in candidates:
        if region.end <= record.start:
            distance = record.start - region.end
            side = "left"
        elif region.start >= record.end:
            distance = region.start - record.end
            side = "right"
        else:
            distance = 0
            side = "overlap"
        upstream = cutoff_value(upstream_bp, record)
        downstream = cutoff_value(downstream_bp, record)
        if side == "overlap":
            allowed = max(upstream, downstream)
        elif record.strand == "+":
            allowed = upstream if side == "left" else downstream
        elif record.strand == "-":
            allowed = downstream if side == "left" else upstream
        else:
            allowed = max(upstream, downstream)
        if distance <= allowed:
            matches.append((record, distance))
    matches.sort(key=lambda item: (item[1], item[0].start, item[0].name, item[0].gene_id))
    return [record for record, _distance in matches]


def format_promoter_label(upstream_bp: int, downstream_bp: int) -> str:
    """Format the promoter range for output column names."""
    upstream = format_distance(upstream_bp)
    downstream = format_distance(downstream_bp)
    if upstream == downstream:
        return upstream
    return f"{upstream}up_{downstream}down"


def parse_prom_enha_cutoffs(value: str) -> Tuple[str, str, str]:
    """Parse promoter-up, enhancer, and optional promoter-down cutoffs."""
    values = [item.strip() for item in value.split(",")]
    if len(values) == 2:
        values.append(values[0])
    if len(values) != 3 or any(not item for item in values):
        raise ValueError("--prom-enha-cutoffs must be promoterup,enhancer[,promoterdown]")
    return values[0], values[1], values[2]


def cutoff_value(value: object, record: BedRecord) -> int:
    """Resolve a fixed, gene-relative, or transcript-relative cutoff."""
    text = str(value).strip()
    lowered = text.lower()
    if lowered.startswith("gene"):
        fraction = float(text[4:])
        length = (record.source_end - record.source_start) if record.source_start is not None and record.source_end is not None else record.length
        return int(length * fraction)
    if lowered.startswith("transcript"):
        multiplier = float(text[10:])
        if len(record.fields) <= 4:
            raise ValueError(f"{text!r} requires BED column 5 transcript length")
        try:
            transcript_length = int(float(record.fields[4]))
        except ValueError as exc:
            raise ValueError(f"{text!r} requires numeric BED column 5 transcript length") from exc
        return int(transcript_length * multiplier)
    return parse_distance(text)


def _fixed_cutoff(value: object) -> Optional[int]:
    """Return a fixed cutoff in base pairs, or ``None`` for relative cutoffs."""
    text = str(value).strip().lower()
    if text.startswith("gene") or text.startswith("transcript"):
        return None
    try:
        return parse_distance(text)
    except ValueError:
        return None


def cutoff_label(value: str) -> str:
    """Format a cutoff for output column names."""
    try:
        return format_distance(parse_distance(value))
    except ValueError:
        return value.replace(".", "p")


def closest_record(index: IntervalIndex, region: InputRegion) -> Tuple[Optional[BedRecord], Optional[int]]:
    """Return the closest TSS record and distance for a peak."""
    candidates = index.nearest_candidates(region.chrom, region.start, region.end)
    best: Optional[Tuple[BedRecord, int]] = None
    for record in candidates:
        dist = distance_bp(region.start, region.end, record.start, record.end)
        if best is None or dist < best[1]:
            best = (record, dist)
    if best is None:
        return None, None
    return best


def names_and_ids(records: Iterable[BedRecord]) -> Tuple[str, str]:
    """Return comma-separated unique gene symbols and IDs."""
    rows = list(records)
    # winandgroup.sh emits each distinct grouped column in lexical order;
    # sort names and IDs independently to match its two output columns.
    names = unique_join(sorted(record.name for record in rows))
    ids = unique_join(sorted(record.gene_id for record in rows))
    return names, ids


def annotate_peak2gene(config: PeakGeneConfig) -> Path:
    """Annotate peaks with nearby and closest genes.

    Args:
        config (PeakGeneConfig): Annotation configuration.

    Returns:
        Path: Output TSV path.
    """
    upstream_cutoff, enhancer_cutoff, downstream_cutoff = parse_prom_enha_cutoffs(config.prom_enha_cutoffs)
    output_format = config.output_format if config.output_format != "auto" else detect_output_format(
        config.input_path, config.header, config.input_format, config.columns, config.region_column,
        config.txt_delimiter,
    )
    header, regions = read_regions(
        config.input_path,
        header=config.header,
        input_format=config.input_format,
        columns=config.columns,
        region_column=config.region_column,
        txt_delimiter=config.txt_delimiter,
    )
    tss_records = resolve_tss_records(config)
    index = IntervalIndex(tss_records)

    promoter_label = f"{cutoff_label(upstream_cutoff)}up_{cutoff_label(downstream_cutoff)}down"
    if upstream_cutoff == downstream_cutoff:
        promoter_label = cutoff_label(upstream_cutoff)
    enhancer_label = cutoff_label(enhancer_cutoff)
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
    rows: List[Tuple[str, int, int, Sequence[object]]] = []
    for region in regions:
        promoter_records_for_region = promoter_records(index, region, upstream_cutoff, downstream_cutoff)
        # Match voom2anno.sh: distal/enhancer genes are reported only for
        # peaks without a promoter assignment.  This keeps the two columns
        # mutually exclusive and avoids substantially larger output tables.
        distal = [] if promoter_records_for_region else nearby_records(index, region, enhancer_cutoff)
        promoter_names, promoter_ids = names_and_ids(promoter_records_for_region)
        distal_names, distal_ids = names_and_ids(record for record, _ in distal)
        closest, distance = closest_record(index, region)
        rows.append((
            region.chrom,
            region.start,
            region.end,
            output_region_values(region, region.values, output_format) + [
                promoter_names,
                promoter_ids,
                distal_names,
                distal_ids,
                closest.name if closest else ".",
                closest.gene_id if closest else ".",
                distance if distance is not None else ".",
            ],
        ))
    rows.sort(key=lambda row: (row[0], row[1], row[2]))
    write_table(config.output_path, out_header, (row[3] for row in rows), include_header=output_format == "txt")
    return config.output_path
