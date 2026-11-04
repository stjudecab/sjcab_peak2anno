"""Fast Python implementations of common paired-locus operations."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

from .intervals import BedRecord, InputRegion, IntervalIndex, read_bed_records, write_table
from .loops import read_bedpe


LoopRow = Tuple[List[str], InputRegion, InputRegion]


def _overlap(a: InputRegion, b: BedRecord, distance: int = 0) -> bool:
    """Return whether two half-open intervals overlap within ``distance``."""
    return a.chrom == b.chrom and a.start - distance < b.end and b.start - distance < a.end


def _loop_overlap(a: InputRegion, b: InputRegion, distance: int = 0) -> bool:
    """Return whether two anchors overlap within ``distance``."""
    return a.chrom == b.chrom and a.start - distance < b.end and b.start - distance < a.end


def _loop_output(
    output_path: Optional[Path],
    output_format: str,
    header: Sequence[str],
    rows: Iterable[Sequence[str]],
) -> Optional[Path]:
    """Write loop rows in BEDPE or tabular form."""
    write_table(output_path, header, rows, include_header=output_format == "txt")
    return output_path


def loop_inter1d(
    input_path: Path,
    bed_path: Path,
    output_path: Optional[Path],
    output_format: str = "txt",
    header: str = "auto",
    loop_columns: Tuple[int, int, int, int, int, int] = (0, 1, 2, 3, 4, 5),
    distance: int = 0,
    invert: bool = False,
    unique: bool = False,
    direction: str = "both",
) -> Optional[Path]:
    """Intersect each loop with BED and annotate the overlapping anchor.

    ``direction`` selects both anchors, only the source anchor, or only the
    target anchor. ``-u`` emits one row per loop and ``-v`` emits loops without
    an overlap.
    """
    if direction not in {"both", "source", "target"}:
        raise ValueError("direction must be one of: both, source, target")
    input_header, loops = read_bedpe(input_path, header=header, columns=loop_columns)
    bed_records = read_bed_records(bed_path, gene_type="all")
    index = IntervalIndex(bed_records)
    output_header = list(input_header) + ["Intersected_Anchor"]
    output_rows: List[List[str]] = []
    for values, left, right in loops:
        left_hits = list(index.query_window(left.chrom, left.start, left.end, distance)) if direction in {"both", "source"} else []
        right_hits = list(index.query_window(right.chrom, right.start, right.end, distance)) if direction in {"both", "target"} else []
        left_hit = bool(left_hits)
        right_hit = bool(right_hits)
        hit = left_hit or right_hit
        if invert:
            if hit:
                continue
            label = "."
        else:
            if not hit:
                continue
            label = "AB" if left_hit and right_hit else ("A" if left_hit else "B")
        if unique:
            output_rows.append(values + [label])
            continue
        if invert:
            output_rows.append(values + [label])
            continue
        # Match pgltools' default 1D behavior: emit one row per intersecting
        # anchor and replace that anchor with its actual overlap interval.
        for record in left_hits:
            output_rows.append(
                [
                    values[0], str(max(left.start, record.start)), str(min(left.end, record.end)),
                    *values[3:6], *values[6:], "A",
                ]
            )
        for record in right_hits:
            output_rows.append(
                [
                    *values[0:3], values[3], str(max(right.start, record.start)),
                    str(min(right.end, record.end)), *values[6:], "B",
                ]
            )
    return _loop_output(output_path, output_format, output_header, output_rows)


def loop_inter(
    input_path: Path,
    other_path: Path,
    output_path: Optional[Path],
    output_format: str = "txt",
    header: str = "auto",
    loop_columns: Tuple[int, int, int, int, int, int] = (0, 1, 2, 3, 4, 5),
    distance: int = 0,
    union: bool = False,
    keep: str = "intersection",
    direction: str = "direction",
) -> Optional[Path]:
    """Intersect paired loci using same-orientation or reverse anchors."""
    if direction not in {"direction", "reverse"}:
        raise ValueError("direction must be one of: direction, reverse")
    input_header, loops_a = read_bedpe(input_path, header=header, columns=loop_columns)
    other_header, loops_b = read_bedpe(other_path, header=header, columns=loop_columns)
    first_anchor_index = IntervalIndex(
        BedRecord(
            (left if direction == "direction" else right).chrom,
            (left if direction == "direction" else right).start,
            (left if direction == "direction" else right).end,
            ((left if direction == "direction" else right).chrom,
             str((left if direction == "direction" else right).start),
             str((left if direction == "direction" else right).end), str(index)),
        )
        for index, (_values, left, right) in enumerate(loops_b)
    )
    output_header = list(input_header) + [f"B_{name}" for name in other_header[6:]]
    output_rows: List[List[str]] = []
    for values_a, a1, a2 in loops_a:
        candidates = first_anchor_index.query_window(a1.chrom, a1.start, a1.end, distance)
        for candidate in candidates:
            values_b, b1, b2 = loops_b[int(candidate.fields[3])]
            b_source, b_target = (b1, b2) if direction == "direction" else (b2, b1)
            if not _loop_overlap(a1, b_source, distance) or not _loop_overlap(a2, b_target, distance):
                continue
            if keep == "a":
                values = values_a[:]
            elif keep == "b":
                values = values_b[:]
            elif union:
                values = [
                    a1.chrom, str(min(a1.start, b_source.start)), str(max(a1.end, b_source.end)),
                    a2.chrom, str(min(a2.start, b_target.start)), str(max(a2.end, b_target.end)),
                ] + values_a[6:]
            else:
                values = [
                    a1.chrom, str(max(a1.start, b_source.start)), str(min(a1.end, b_source.end)),
                    a2.chrom, str(max(a2.start, b_target.start)), str(min(a2.end, b_target.end)),
                ] + values_a[6:]
            if keep not in {"a", "b"}:
                values.extend(values_b[6:])
            output_rows.append(values)
    return _loop_output(output_path, output_format, output_header, output_rows)


def _merge_anchor_groups(loops: Sequence[LoopRow], side: int, distance: int) -> List[LoopRow]:
    """Merge one anchor while retaining every loop's opposite anchor."""
    ordered = sorted(loops, key=lambda row: ((row[1] if side == 1 else row[2]).chrom,
                                             (row[1] if side == 1 else row[2]).start,
                                             (row[1] if side == 1 else row[2]).end))
    groups: List[List[LoopRow]] = []
    current: List[LoopRow] = []
    current_chrom = ""
    current_end = -1
    for row in ordered:
        anchor = row[1] if side == 1 else row[2]
        if current and (anchor.chrom != current_chrom or anchor.start - distance > current_end):
            groups.append(current)
            current = []
        current.append(row)
        current_chrom = anchor.chrom
        current_end = max(current_end, anchor.end)
    if current:
        groups.append(current)
    output: List[LoopRow] = []
    for group in groups:
        merged_anchor = group[0][1] if side == 1 else group[0][2]
        for _values, left, right in group[1:]:
            anchor = left if side == 1 else right
            merged_anchor = InputRegion(merged_anchor.chrom, min(merged_anchor.start, anchor.start), max(merged_anchor.end, anchor.end), ())
        for values, left, right in group:
            output.append((values[:], merged_anchor if side == 1 else left, right if side == 1 else merged_anchor))
    return output


def loop_merge(
    input_path: Path,
    output_path: Optional[Path],
    output_format: str = "txt",
    header: str = "auto",
    loop_columns: Tuple[int, int, int, int, int, int] = (0, 1, 2, 3, 4, 5),
    distance: int = 0,
    direction: str = "both",
) -> Optional[Path]:
    """Merge both anchors, or one selected anchor while retaining the other."""
    if direction not in {"both", "source", "target"}:
        raise ValueError("direction must be one of: both, source, target")
    input_header, loops = read_bedpe(input_path, header=header, columns=loop_columns)
    if direction == "source":
        merged = _merge_anchor_groups(loops, 1, distance)
    elif direction == "target":
        merged = _merge_anchor_groups(loops, 2, distance)
    else:
        ordered = sorted(loops, key=lambda row: (row[1].chrom, row[1].start, row[1].end, row[2].chrom, row[2].start, row[2].end))
        merged = []
        for values, left, right in ordered:
            if merged:
                old_values, old_left, old_right = merged[-1]
                if _loop_overlap(old_left, left, distance) and _loop_overlap(old_right, right, distance):
                    merged[-1] = (old_values + [value for value in values[6:] if value not in old_values[6:]],
                                  InputRegion(old_left.chrom, min(old_left.start, left.start), max(old_left.end, left.end), ()),
                                  InputRegion(old_right.chrom, min(old_right.start, right.start), max(old_right.end, right.end), ()))
                    continue
            merged.append((values[:], left, right))
    output_rows = [
        [left.chrom, str(left.start), str(left.end), right.chrom, str(right.start), str(right.end)] + values[6:]
        for values, left, right in merged
    ]
    return _loop_output(output_path, output_format, input_header, output_rows)
