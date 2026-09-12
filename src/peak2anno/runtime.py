"""Optional external tool detection and runtime warnings."""

from __future__ import annotations

import os
import shutil
import sys
from dataclasses import dataclass
from typing import Optional


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


def warn_if_slow(status: ToolStatus) -> None:
    """Explain the Python fallback when bedtools acceleration is unavailable."""
    if status.accelerated:
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
