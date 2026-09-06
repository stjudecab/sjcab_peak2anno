"""BEDPE loop annotation helpers."""

from __future__ import annotations

import csv
import tempfile
from pathlib import Path
from typing import Callable, List, Optional, Sequence, Tuple

from .context import ContextConfig, StateConfig, annotate_broad_context, annotate_narrow_context, annotate_peak_state
from .intervals import InputRegion, read_regions, split_fields, write_table
from .peak2gene import PeakGeneConfig, annotate_peak2gene


def read_bedpe(path: Path, header: str = "auto", columns: Tuple[int, int, int, int, int, int] = (0, 1, 2, 3, 4, 5)) -> Tuple[List[str], List[Tuple[List[str], InputRegion, InputRegion]]]:
    """Read BEDPE rows and return original values plus two anchor regions."""
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if raw.strip() and not raw.startswith(("#", "track", "browser")):
                rows.append(split_fields(raw))
    if not rows:
        raise ValueError(f"No BEDPE rows found in {path}")
    first = rows[0]
    has_header = header == "yes" or (
        header == "auto"
        and (max(columns) >= len(first) or any(not value.lstrip("-").isdigit() for value in (first[columns[1]], first[columns[2]], first[columns[4]], first[columns[5]])))
    )
    if has_header:
        out_header = first
        data = rows[1:]
    else:
        out_header = ["chr1", "start1", "end1", "chr2", "start2", "end2"]
        if len(first) > 6:
            out_header.extend(f"field{i}" for i in range(7, len(first) + 1))
        data = rows
    result = []
    for number, row in enumerate(data, start=2 if has_header else 1):
        if max(columns) >= len(row):
            raise ValueError(f"BEDPE columns {columns} are not present at {path}:{number}")
        values = row + ["."] * (len(out_header) - len(row))
        c1, s1, e1, c2, s2, e2 = columns
        try:
            anchor1 = InputRegion(values[c1], int(values[s1]), int(values[e1]), (values[c1], values[s1], values[e1], f"loop{number}_1"))
            anchor2 = InputRegion(values[c2], int(values[s2]), int(values[e2]), (values[c2], values[s2], values[e2], f"loop{number}_2"))
        except ValueError as exc:
            raise ValueError(f"Invalid BEDPE coordinates at {path}:{number}: {row!r}") from exc
        result.append((values, anchor1, anchor2))
    return out_header, result


def _write_anchor(path: Path, rows: Sequence[Tuple[List[str], InputRegion, InputRegion]], side: int) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for number, (_values, left, right) in enumerate(rows, start=1):
            region = left if side == 1 else right
            handle.write("\t".join([region.chrom, str(region.start), str(region.end), f"loop{number}_{side}"]) + "\n")


def _read_table(path: Path) -> Tuple[List[str], List[List[str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    return rows[0], rows[1:]


def annotate_loop(
    command: str,
    input_path: Path,
    output_path: Optional[Path],
    output_format: str,
    header: str = "auto",
    loop_columns: Tuple[int, int, int, int, int, int] = (0, 1, 2, 3, 4, 5),
    **kwargs: object,
) -> Optional[Path]:
    """Annotate both anchors of a BEDPE file and merge the result."""
    input_header, rows = read_bedpe(input_path, header=header, columns=loop_columns)
    if output_format == "auto":
        if header == "yes":
            output_format = "txt"
        elif header == "no":
            output_format = "bedpe"
        else:
            first = next(line for line in input_path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#"))
            fields = split_fields(first)
            output_format = "bedpe" if all(fields[index].lstrip("-").isdigit() for index in (loop_columns[1], loop_columns[2], loop_columns[4], loop_columns[5])) else "txt"
    with tempfile.TemporaryDirectory(prefix="peak2anno-loop-") as temp:
        temp_root = Path(temp)
        anchor_outputs = []
        for side in (1, 2):
            anchor_input = temp_root / f"anchor{side}.bed"
            anchor_output = temp_root / f"anchor{side}.tsv"
            _write_anchor(anchor_input, rows, side)
            if command == "loop2anno":
                annotate_peak2gene(PeakGeneConfig(input_path=anchor_input, output_path=anchor_output, species=str(kwargs["species"]), species_version=str(kwargs["species_version"]), isoform_version=str(kwargs["isoform_version"]), db_path=kwargs.get("db_path"), tss_bed=kwargs.get("tss_bed"), gene_bed=kwargs.get("gene_bed"), output_format="txt"))
            elif command == "loop2context":
                config = ContextConfig(input_path=anchor_input, output_path=anchor_output, species=str(kwargs["species"]), db_path=kwargs.get("db_path"), context_dir=kwargs.get("context_dir"), overlap_cutoff=str(kwargs["overlap_cutoff"]), output_format="txt")
                (annotate_broad_context if kwargs.get("context_mode") == "broad" else annotate_narrow_context)(config)
            elif command == "loop2state":
                annotate_peak_state(StateConfig(input_path=anchor_input, output_path=anchor_output, states_path=Path(str(kwargs["states"])), state2name=kwargs.get("state2name"), overlap_cutoff=str(kwargs["overlap_cutoff"]), output_format="txt"))
            else:
                raise ValueError(f"Unknown loop command: {command}")
            anchor_outputs.append(_read_table(anchor_output))
        left_header, left_rows = anchor_outputs[0]
        right_header, right_rows = anchor_outputs[1]
        output_header = list(input_header) + [f"anchor1_{name}" for name in left_header[4:]] + [f"anchor2_{name}" for name in right_header[4:]]
        output_rows = [values + left[4:] + right[4:] for (values, _a, _b), left, right in zip(rows, left_rows, right_rows)]
        write_table(output_path, output_header, output_rows, include_header=output_format == "txt")
    return output_path
