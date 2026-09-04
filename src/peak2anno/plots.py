"""Optional plotting helpers for peak2anno summary tables."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping


def write_count_plots(counts: Mapping[str, int], prefix: Path, title: str) -> None:
    """Write bar and pie plots for feature counts.

    Args:
        counts (Mapping[str, int]): Feature counts.
        prefix (Path): Output prefix without extension.
        title (str): Plot title.

    Raises:
        RuntimeError: If matplotlib is not installed.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Plotting requires matplotlib; install matplotlib-base or omit --plot") from exc

    labels = list(counts.keys())
    values = [counts[label] for label in labels]
    prefix.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, max(4, len(labels) * 0.35)))
    ax.barh(labels, values, color="#4c78a8")
    ax.set_xlabel("Regions")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(prefix.with_suffix(".barPlot.png"), dpi=200)
    fig.savefig(prefix.with_suffix(".barPlot.pdf"))
    plt.close(fig)

    if sum(values) > 0:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.pie(values, labels=labels, autopct="%1.1f%%")
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(prefix.with_suffix(".piePlot.png"), dpi=200)
        fig.savefig(prefix.with_suffix(".piePlot.pdf"))
        plt.close(fig)

