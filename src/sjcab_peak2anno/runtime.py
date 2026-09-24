"""Optional external tool detection and runtime warnings."""

from __future__ import annotations

import shutil
import sys
from dataclasses import dataclass
from typing import Optional


BACKENDS = ("auto", "bedtools", "python")


@dataclass(frozen=True)
class ToolStatus:
    """Availability of optional interval acceleration tools."""

    bedtools: Optional[str]


def detect_tools() -> ToolStatus:
    """Detect the optional native bedtools executable."""
    return ToolStatus(bedtools=shutil.which("bedtools"))


def resolve_backend(requested: str, status: ToolStatus) -> str:
    """Resolve a requested backend against the tools available on ``PATH``.

    ``auto`` deliberately selects the indexed Python implementation.
    """
    if requested not in BACKENDS:
        raise ValueError(f"backend must be one of: {', '.join(BACKENDS)}")
    if requested == "auto":
        return "python"
    if requested == "bedtools" and status.bedtools is None:
        raise ValueError("SJCAB_PEAK2ANNO_BACKEND=bedtools requested, but bedtools was not found on PATH")
    return requested


def warn_if_slow(status: ToolStatus, backend: str = "auto") -> None:
    """Explain the Python fallback when no external backend is selected."""
    if backend != "python":
        return
    missing = []
    if status.bedtools is None:
        missing.append("bedtools")
    if not missing:
        return
    tools = " and ".join(missing)
    print(
        f"peak2anno: {tools} not found; using the slower Python interval fallback. "
        "Set SJCAB_PEAK2ANNO_BACKEND=bedtools if you want to generate and review a bedtools bash script; bedtools must be on PATH.",
        file=sys.stderr,
    )
