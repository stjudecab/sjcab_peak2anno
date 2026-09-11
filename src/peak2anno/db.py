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


def gene_annotation_path(
    species: str,
    version: str = "def",
    isoform_set: str = "all",
    root_path: Optional[str] = None,
) -> Path:
    """Find the gene BED for a species, version, and isoform set.

    The preferred database layout is ``<root>/<species>/<version>/`` with
    either ``all.gene.bed`` or ``deduplong.gene.bed``. Manifest paths and the
    older ``<root>/<species>/{tss,deduplong>/<version>.bed`` layout are also
    supported.
    """
    if isoform_set not in {"all", "deduplong"}:
        raise ValueError("isoform set must be one of: all, deduplong")
    root = db_root(root_path)
    filename = f"{isoform_set}.gene.bed"
    manifest = load_manifest(root)
    resources = manifest_resources(manifest)

    if version in {"default", "def"}:
        for item in resources:
            item_path = str(item.get("path", ""))
            annotation = str(item.get("annotation", ""))
            if (
                item.get("species") == species
                and item.get("default") is True
                and (annotation in {isoform_set, f"{isoform_set}.gene", "gene"}
                     or (isoform_set == "all" and annotation == "tss")
                     or item_path.endswith(filename))
            ):
                version = str(item.get("version", "default"))
                break
        else:
            version_candidates = sorted(
                path.parent.name
                for path in (root / species).glob(f"*/{filename}")
            )
            if len(version_candidates) == 1:
                version = version_candidates[0]
            elif not version_candidates:
                version = "default"
            else:
                raise FileNotFoundError(
                    f"No default gene BED version for {species}; choose one with --ver"
                )

    for item in resources:
        if item.get("species") != species or str(item.get("version")) != version:
            continue
        item_path = item.get("path")
        if item_path and (
            str(item.get("annotation", "")) in {isoform_set, f"{isoform_set}.gene", "gene"}
            or (isoform_set == "all" and str(item.get("annotation", "")) == "tss")
            or str(item_path).endswith(filename)
        ):
            path = root / str(item_path)
            if path.is_file():
                return path

    candidates = [
        root / species / version / filename,
        root / species / f"{version}.{isoform_set}.gene.bed",
        root / species / isoform_set / f"{version}.gene.bed",
    ]
    # Compatibility with the original database package layout.
    candidates.append(root / species / ("tss" if isoform_set == "all" else "deduplong") / f"{version}.bed")
    for path in candidates:
        if path.is_file():
            return path
    searched = ", ".join(str(path) for path in candidates)
    raise FileNotFoundError(
        f"Missing {filename} for {species} version {version} under {root}; searched: {searched}"
    )


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


def candidate_feature_dirs(species: str, root_path: Optional[str] = None) -> List[Path]:
    """Return feature directory candidates for a species."""
    root = db_root(root_path)
    package_root = Path(__file__).resolve().parents[2]
    cwd = Path.cwd()
    candidates = [
        root / species / "features",
        root / species,
    ]
    species_root = root / species
    if species_root.is_dir():
        for version_root in sorted(path for path in species_root.iterdir() if path.is_dir()):
            candidates.extend([version_root / "features", version_root])
    candidates.extend([cwd / "annotations" / species, package_root / "annotations" / species])
    return candidates
