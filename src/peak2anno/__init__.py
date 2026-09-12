"""Command-line tools for annotating genomic peaks."""

try:
    from ._version import version as __version__
except ImportError:  # pragma: no cover - source checkout before a build
    __version__ = "0+unknown"
