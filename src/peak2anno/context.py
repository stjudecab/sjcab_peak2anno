"""Context and chromatin-state annotation implementations."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from .db import candidate_context_dirs
from .intervals import (
    BedRecord,
    InputRegion,
    IntervalIndex,
    OverlapCutoff,
    overlap_bp,
    output_region_values,
    output_region_header,
    detect_output_format,
    parse_overlap_cutoff,
    read_bed_records,
    read_regions,
    split_fields,
    write_table,
)
from .plots import write_count_plots


DEFAULT_FEATURES = [
    ("2kb.promoter.up.bed", "Promoter.Up"),
    ("2kb.promoter.down.bed", "Promoter.Down"),
    ("2kb.exon.bed", "Exon"),
    ("2kb.intron.bed", "Intron"),
    ("2kb.tes.bed", "TES (transcription end sites)"),
    ("2kb.dis5.bed", "Dis5 (5' distal regions)"),
    ("2kb.dis3.bed", "Dis3 (3' distal regions)"),
    ("2kb.intergenic.bed", "Intergenic"),
]


@dataclass(frozen=True)
class FeatureSpec:
    """One annotation feature BED and output label."""

    path: Path
    label: str


@dataclass(frozen=True)
class ContextConfig:
    """Configuration for narrow or broad genomic context annotation."""

    input_path: Path
    output_path: Optional[Path]
    species: str
    db_path: Optional[str] = None
    context_dir: Optional[Path] = None
    features: Optional[str] = None
    feature_labels: Optional[str] = None
    overlap_cutoff: str = "1bp"
    header: str = "auto"
    column_name: str = "FeatureAssignment"
    summary_path: Optional[Path] = None
    plot: bool = False
    input_format: str = "auto"
    columns: Optional[Tuple[int, int, int]] = None
    region_column: int = 0
    output_format: str = "txt"


@dataclass(frozen=True)
class StateConfig:
    """Configuration for chromatin-state annotation."""

    input_path: Path
    states_path: Path
    output_path: Optional[Path]
    state2name: Optional[Path] = None
    overlap_cutoff: str = "1bp"
    header: str = "auto"
    summary_path: Optional[Path] = None
    plot: bool = False
    input_format: str = "auto"
    columns: Optional[Tuple[int, int, int]] = None
    region_column: int = 0
    output_format: str = "txt"


def read_list_or_csv(value: str, base_dir: Optional[Path] = None) -> List[str]:
    """Read values from a .lst file or comma-separated text."""
    path = Path(value)
    if base_dir is not None and not path.is_absolute():
        candidate = base_dir / path
        if candidate.is_file():
            path = candidate
    if path.is_file():
        return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return [item.strip() for item in value.split(",") if item.strip()]


def resolve_context_dir(species: str, db_path: Optional[str], explicit: Optional[Path]) -> Path:
    """Resolve the feature BED directory for genomic context commands."""
    if explicit is not None:
        if not explicit.is_dir():
            raise FileNotFoundError(f"--context-dir does not exist: {explicit}")
        return explicit
    for candidate in candidate_context_dirs(species, root_path=db_path):
        if (candidate / DEFAULT_FEATURES[0][0]).is_file():
            return candidate
    searched = "\n".join(str(path) for path in candidate_context_dirs(species, root_path=db_path))
    raise FileNotFoundError(
        "Could not find default context feature BEDs. Use --context-dir or --features/--feature-labels.\n"
        f"Searched:\n{searched}"
    )


def resolve_features(config: ContextConfig) -> List[FeatureSpec]:
    """Resolve context feature BED files and labels."""
    context_dir = resolve_context_dir(config.species, config.db_path, config.context_dir)
    if config.features is None and config.feature_labels is None:
        return [FeatureSpec(context_dir / filename, label) for filename, label in DEFAULT_FEATURES]
    if config.features is None or config.feature_labels is None:
        raise ValueError("--features and --feature-labels must be provided together")
    feature_values = read_list_or_csv(config.features, base_dir=context_dir)
    label_values = read_list_or_csv(config.feature_labels, base_dir=context_dir)
    if len(feature_values) != len(label_values):
        raise ValueError("--features and --feature-labels must have the same number of entries")
    specs: List[FeatureSpec] = []
    for value, label in zip(feature_values, label_values):
        path = Path(value)
        if not path.is_absolute():
            path = context_dir / path
        specs.append(FeatureSpec(path, label))
    return specs


def feature_indexes(specs: Sequence[FeatureSpec]) -> OrderedDict[str, IntervalIndex]:
    """Load feature BEDs into interval indexes keyed by label."""
    indexes: OrderedDict[str, IntervalIndex] = OrderedDict()
    for spec in specs:
        if not spec.path.is_file():
            raise FileNotFoundError(f"Missing feature BED for {spec.label}: {spec.path}")
        indexes[spec.label] = IntervalIndex(read_bed_records(spec.path))
    return indexes


def total_overlap(index: IntervalIndex, region: InputRegion) -> int:
    """Return merged overlap bases between a region and an indexed feature."""
    intervals: List[Tuple[int, int]] = []
    for record in index.query(region.chrom, region.start, region.end):
        bp = overlap_bp(region.start, region.end, record.start, record.end)
        if bp > 0:
            intervals.append((max(region.start, record.start), min(region.end, record.end)))
    if not intervals:
        return 0
    intervals.sort()
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return sum(end - start for start, end in merged)


def assign_priority(
    region: InputRegion,
    indexes: Mapping[str, IntervalIndex],
    cutoff: OverlapCutoff,
) -> str:
    """Assign the first feature satisfying the overlap cutoff."""
    for label, index in indexes.items():
        bp = total_overlap(index, region)
        if cutoff.passes(bp, region.length):
            return label
    return "False"


def summarize_counts(labels: Sequence[str], assignments: Iterable[str]) -> OrderedDict[str, int]:
    """Count assignments while preserving requested label order."""
    counts: OrderedDict[str, int] = OrderedDict((label, 0) for label in labels)
    counts["False"] = 0
    for assignment in assignments:
        counts[assignment] = counts.get(assignment, 0) + 1
    return counts


def write_count_summary(
    path: Path,
    source_name: str,
    counts: Mapping[str, int],
    total_regions: int,
    threshold: str,
) -> None:
    """Write a one-row feature count summary."""
    header = ["Annotated file", "Category", "Threshold", "#Regions"] + list(counts.keys())
    row = [source_name, "AllRegions", threshold, total_regions] + [counts[key] for key in counts]
    write_table(path, header, [row])


def annotate_narrow_context(config: ContextConfig) -> Tuple[Path, Path]:
    """Annotate each peak to one prioritized genomic context feature."""
    output_format = config.output_format if config.output_format != "auto" else detect_output_format(config.input_path, config.header, config.input_format, config.columns, config.region_column)
    header, regions = read_regions(config.input_path, header=config.header, input_format=config.input_format, columns=config.columns, region_column=config.region_column)
    specs = resolve_features(config)
    indexes = feature_indexes(specs)
    cutoff = parse_overlap_cutoff(config.overlap_cutoff)
    assignments = [assign_priority(region, indexes, cutoff) for region in regions]
    out_header = output_region_header(header, output_format) + [config.column_name]
    rows = [output_region_values(region, region.values, output_format) + [assignment] for region, assignment in zip(regions, assignments)]
    write_table(config.output_path, out_header, rows, include_header=output_format == "txt")

    counts = summarize_counts([spec.label for spec in specs], assignments)
    summary = config.summary_path or (
        config.output_path.with_suffix(config.output_path.suffix + ".summary.tsv")
        if config.output_path is not None
        else None
    )
    if summary is not None:
        write_count_summary(summary, config.input_path.name, counts, len(regions), config.overlap_cutoff)
    if config.plot:
        if summary is None:
            raise ValueError("--plot requires --summary when --output is omitted")
        write_count_plots(counts, summary.with_suffix(""), f"Genomic context: {config.input_path.name}")
    return config.output_path, summary


def broad_rows(
    regions: Sequence[InputRegion],
    indexes: Mapping[str, IntervalIndex],
    cutoff: OverlapCutoff,
) -> Tuple[List[str], List[List[object]], OrderedDict[str, int], OrderedDict[str, int]]:
    """Compute broad overlap rows and summary counts."""
    labels = list(indexes.keys())
    primary_counts: OrderedDict[str, int] = OrderedDict((label, 0) for label in labels)
    primary_counts["False"] = 0
    bp_totals: OrderedDict[str, int] = OrderedDict((label, 0) for label in labels)
    out_header: List[str] = []
    for label in labels:
        safe_label = label.replace(" ", "_")
        out_header.extend([f"{safe_label}_bp", f"{safe_label}_fraction"])
    out_header.extend(["PrimaryFeature", "PrimaryFeatureFraction"])
    rows: List[List[object]] = []
    for region in regions:
        overlaps = OrderedDict((label, total_overlap(index, region)) for label, index in indexes.items())
        fractions = OrderedDict(
            (label, (bp / region.length if region.length > 0 else 0.0)) for label, bp in overlaps.items()
        )
        for label, bp in overlaps.items():
            bp_totals[label] += bp
        passing = [
            (label, fractions[label], overlaps[label])
            for label in labels
            if cutoff.passes(overlaps[label], region.length)
        ]
        if passing:
            primary, primary_fraction, _ = max(passing, key=lambda item: (item[1], item[2], -labels.index(item[0])))
        else:
            primary, primary_fraction = "False", 0.0
        primary_counts[primary] = primary_counts.get(primary, 0) + 1
        row: List[object] = []
        for label in labels:
            row.extend([overlaps[label], f"{fractions[label]:.6f}"])
        row.extend([primary, f"{primary_fraction:.6f}"])
        rows.append(row)
    return out_header, rows, primary_counts, bp_totals


def annotate_broad_context(config: ContextConfig) -> Tuple[Path, Path]:
    """Annotate each peak with per-feature context overlap fractions."""
    output_format = config.output_format if config.output_format != "auto" else detect_output_format(config.input_path, config.header, config.input_format, config.columns, config.region_column)
    header, regions = read_regions(config.input_path, header=config.header, input_format=config.input_format, columns=config.columns, region_column=config.region_column)
    specs = resolve_features(config)
    indexes = feature_indexes(specs)
    cutoff = parse_overlap_cutoff(config.overlap_cutoff)
    extra_header, extra_rows, primary_counts, bp_totals = broad_rows(regions, indexes, cutoff)
    rows = [output_region_values(region, region.values, output_format) + extra for region, extra in zip(regions, extra_rows)]
    write_table(config.output_path, output_region_header(header, output_format) + extra_header, rows, include_header=output_format == "txt")

    summary = config.summary_path or (
        config.output_path.with_suffix(config.output_path.suffix + ".summary.tsv")
        if config.output_path is not None
        else None
    )
    total_bp = sum(region.length for region in regions)
    summary_header = ["Feature", "PrimaryRegions", "OverlapBp", "OverlapFractionOfInputBp"]
    summary_rows = [
        [
            label,
            primary_counts.get(label, 0),
            bp_totals.get(label, 0),
            f"{(bp_totals.get(label, 0) / total_bp if total_bp else 0.0):.6f}",
        ]
        for label in indexes.keys()
    ]
    summary_rows.append(["False", primary_counts.get("False", 0), 0, "0.000000"])
    if summary is not None:
        write_table(summary, summary_header, summary_rows)
    if config.plot:
        if summary is None:
            raise ValueError("--plot requires --summary when --output is omitted")
        write_count_plots(primary_counts, summary.with_suffix(""), f"Broad genomic context: {config.input_path.name}")
    return config.output_path, summary


def read_state_names(path: Optional[Path]) -> Dict[str, str]:
    """Read state ID to state label mappings."""
    if path is None:
        return {}
    names: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = split_fields(raw)
            if len(fields) < 2:
                continue
            if fields[0].lower() in {"state", "id"}:
                continue
            names[fields[0]] = fields[1]
    return names


def load_state_index(states_path: Path, state_names: Mapping[str, str]) -> OrderedDict[str, IntervalIndex]:
    """Load a chromatin-state BED into one interval index per state label."""
    grouped: OrderedDict[str, List[BedRecord]] = OrderedDict()
    for record in read_bed_records(states_path):
        state_id = record.name
        label = state_names.get(state_id, state_id)
        grouped.setdefault(label, []).append(record)
    return OrderedDict((label, IntervalIndex(records)) for label, records in grouped.items())


def annotate_peak_state(config: StateConfig) -> Tuple[Path, Path]:
    """Annotate peaks with chromatin-state overlap fractions."""
    output_format = config.output_format if config.output_format != "auto" else detect_output_format(config.input_path, config.header, config.input_format, config.columns, config.region_column)
    header, regions = read_regions(config.input_path, header=config.header, input_format=config.input_format, columns=config.columns, region_column=config.region_column)
    state_names = read_state_names(config.state2name)
    indexes = load_state_index(config.states_path, state_names)
    cutoff = parse_overlap_cutoff(config.overlap_cutoff)
    extra_header, extra_rows, primary_counts, bp_totals = broad_rows(regions, indexes, cutoff)
    rows = [output_region_values(region, region.values, output_format) + extra for region, extra in zip(regions, extra_rows)]
    write_table(config.output_path, output_region_header(header, output_format) + extra_header, rows, include_header=output_format == "txt")

    summary = config.summary_path or (
        config.output_path.with_suffix(config.output_path.suffix + ".summary.tsv")
        if config.output_path is not None
        else None
    )
    total_bp = sum(region.length for region in regions)
    summary_header = ["State", "PrimaryRegions", "OverlapBp", "OverlapFractionOfInputBp"]
    summary_rows = [
        [
            label,
            primary_counts.get(label, 0),
            bp_totals.get(label, 0),
            f"{(bp_totals.get(label, 0) / total_bp if total_bp else 0.0):.6f}",
        ]
        for label in indexes.keys()
    ]
    summary_rows.append(["False", primary_counts.get("False", 0), 0, "0.000000"])
    if summary is not None:
        write_table(summary, summary_header, summary_rows)
    if config.plot:
        if summary is None:
            raise ValueError("--plot requires --summary when --output is omitted")
        write_count_plots(primary_counts, summary.with_suffix(""), f"Peak state context: {config.input_path.name}")
    return config.output_path, summary
