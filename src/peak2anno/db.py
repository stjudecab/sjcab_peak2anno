"""Resolve sjcab_peak2anno_db annotation files."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional


DEFAULT_DB_PATH = Path("~/.sjcab_peak2anno_db").expanduser()


def db_root(path: Optional[str] = None) -> Path:
    """Return the annotation database root path.

    Args:
        path (Optional[str]): Explicit database path.

    Returns:
        Path: Resolved database root.
    """
    if path:
        return Path(path).expanduser()
    env_path = os.environ.get("SJCAB_PEAK2ANNO_DB_PATH")
    if env_path:
        return Path(env_path).expanduser()
    return DEFAULT_DB_PATH


def load_manifest(root: Path) -> Dict[str, Any]:
    """Load the database manifest when available."""
    manifest = root / "manifest.json"
    if not manifest.is_file():
        return {}
    with manifest.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def manifest_resources(manifest: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Return annotation resources from supported manifest schema versions."""
    if "resources" in manifest:
        return list(manifest["resources"])
    if "files" in manifest:
        return list(manifest["files"])
    return []


def resolve_version(root: Path, species: str, version: str, annotation: str) -> str:
    """Resolve ``default`` into a concrete species annotation version.

    Args:
        root (Path): Database root.
        species (str): Genome species key, for example ``hg38``.
        version (str): Requested version or ``default``.
        annotation (str): Annotation type, for example ``tss``.

    Returns:
        str: Concrete version.

    Raises:
        FileNotFoundError: If no matching version can be found.
    """
    if version not in {"default", "def"}:
        return version
    manifest = load_manifest(root)
    for item in manifest_resources(manifest):
        if (
            item.get("species") == species
            and item.get("annotation") == annotation
            and item.get("default") is True
        ):
            return str(item["version"])
    folder = root / species / annotation
    candidates = sorted(path.stem for path in folder.glob("*.bed")) if folder.is_dir() else []
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(f"No {annotation} BED files found for {species} under {root}")
    raise FileNotFoundError(
        f"No manifest default for {species}/{annotation}; choose one of: {', '.join(candidates)}"
    )


def annotation_path(
    species: str,
    version: str = "default",
    annotation: str = "tss",
    root_path: Optional[str] = None,
) -> Path:
    """Return a database annotation BED path.

    Args:
        species (str): Genome species key.
        version (str): Version or ``default``.
        annotation (str): Annotation type.
        root_path (Optional[str]): Explicit database path.

    Returns:
        Path: BED path.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    root = db_root(root_path)
    concrete = resolve_version(root, species, version, annotation)
    path = root / species / annotation / f"{concrete}.bed"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {annotation} BED: {path}")
    return path


def available_versions(root_path: Optional[str] = None) -> List[Dict[str, Any]]:
    """Return available database annotations from the manifest or filesystem."""
    root = db_root(root_path)
    manifest = load_manifest(root)
    resources = manifest_resources(manifest)
    if resources:
        return resources
    rows: List[Dict[str, Any]] = []
    for species_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        if species_dir.name in {"blacklists", "cgi"}:
            continue
        for annotation_dir in sorted(path for path in species_dir.iterdir() if path.is_dir()):
            for bed in sorted(annotation_dir.glob("*.bed")):
                rows.append(
                    {
                        "species": species_dir.name,
                        "annotation": annotation_dir.name,
                        "version": bed.stem,
                        "path": str(bed.relative_to(root)),
                    }
                )
    return rows


def candidate_context_dirs(species: str, root_path: Optional[str] = None) -> List[Path]:
    """Return context feature directory candidates for a species."""
    root = db_root(root_path)
    package_root = Path(__file__).resolve().parents[2]
    cwd = Path.cwd()
    return [
        root / species / "context",
        root / species / "features",
        root / species,
        cwd / "annotations" / species,
        package_root / "annotations" / species,
    ]
