"""Optional external tool detection and runtime warnings."""

from __future__ import annotations

import os
import shutil
import sys
from dataclasses import dataclass
from typing import Optional


BACKENDS = ("auto", "bedtools", "pybedtools", "python")


@dataclass(frozen=True)
class ToolStatus:
    """Availability of optional interval acceleration tools."""

    bedtools: Optional[str]
    pybedtools: bool

    @property
    def accelerated(self) -> bool:
        """Return whether both pybedtools and bedtools are available."""
        return self.bedtools is not None and self.pybedtools


def detect_tools() -> ToolStatus:
    """Detect bedtools and pybedtools without making either a pip requirement."""
    bedtools = shutil.which("bedtools")
    try:
        import pybedtools  # type: ignore[import-not-found]
    except ImportError:
        return ToolStatus(bedtools=bedtools, pybedtools=False)

    if bedtools is not None:
        pybedtools.set_bedtools_path(os.path.dirname(bedtools))
    return ToolStatus(bedtools=bedtools, pybedtools=True)


def resolve_backend(requested: str, status: ToolStatus) -> str:
    """Resolve a requested backend against the tools available on ``PATH``.

    ``bedtools`` is preferred for ``auto`` because it avoids the Python
    wrapper overhead.  The pybedtools backend still requires the bedtools
    executable for interval operations.
    """
    if requested not in BACKENDS:
        raise ValueError(f"backend must be one of: {', '.join(BACKENDS)}")
    if requested == "auto":
        if status.bedtools is not None:
            return "bedtools"
        if status.pybedtools and status.bedtools is not None:
            return "pybedtools"
        return "python"
    if requested == "bedtools" and status.bedtools is None:
        raise ValueError("--backend bedtools requested, but bedtools was not found on PATH")
    if requested == "pybedtools" and not status.pybedtools:
        raise ValueError("--backend pybedtools requested, but pybedtools is not installed")
    if requested == "pybedtools" and status.bedtools is None:
        raise ValueError("--backend pybedtools requested, but its bedtools executable was not found on PATH")
    return requested


def warn_if_slow(status: ToolStatus, backend: str = "auto") -> None:
    """Explain the Python fallback when no external backend is selected."""
    if backend != "python":
        return
    missing = []
    if status.bedtools is None:
        missing.append("bedtools")
    if not status.pybedtools:
        missing.append("pybedtools")
    tools = " and ".join(missing)
    print(
        f"peak2anno: {tools} not found; using the slower Python interval fallback. "
        "Install bedtools and pybedtools (the conda package includes both) for faster runs.",
        file=sys.stderr,
    )
