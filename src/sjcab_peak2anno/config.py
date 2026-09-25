"""Configuration-file and environment defaults for the command-line tools."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import shlex
from typing import Dict, Optional


DEFAULT_SPECIES_VERSIONS = "hg38:v31"
DEFAULT_PROM_ENHA_CUTOFFS = "2kb,50kb,2kb"
DEFAULT_GENE_TYPE = "nomicro"
DEFAULT_ISO_SET = "all"
DEFAULT_2FEATURE_OUT = "max"
DEFAULT_2STATE_OUT = "max,percent"
DEFAULT_TXT_DELIMITER = "auto"
DEFAULT_BACKEND = "python"

RC_DEFAULTS = {
    "SJCAB_PEAK2ANNO_DB_PATH": "~/.sjcab_peak2anno_db",
    "SJCAB_PEAK2ANNO_SPECIES_VERSIONS": DEFAULT_SPECIES_VERSIONS,
    "SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS": DEFAULT_PROM_ENHA_CUTOFFS,
    "SJCAB_PEAK2ANNO_GENE_TYPE": DEFAULT_GENE_TYPE,
    "SJCAB_PEAK2ANNO_ISO_SET": DEFAULT_ISO_SET,
    "SJCAB_PEAK2ANNO_2FEATURE_OUT": DEFAULT_2FEATURE_OUT,
    "SJCAB_PEAK2ANNO_2STATE_OUT": DEFAULT_2STATE_OUT,
    "SJCAB_PEAK2ANNO_TXT_DELIMITER": DEFAULT_TXT_DELIMITER,
    "SJCAB_PEAK2ANNO_BACKEND": DEFAULT_BACKEND,
}


@dataclass(frozen=True)
class Settings:
    """Resolved package settings after rc-file and environment overrides."""

    db_path: Optional[str]
    species_versions: str
    prom_enha_cutoffs: str
    prom_enha_overrides: Dict[str, str]
    gene_type: str
    iso_set: str
    feature_out: str
    state_out: str
    txt_delimiter: str
    backend: str

    @property
    def default_species(self) -> str:
        """Return the first species in the species/version setting."""
        return next(iter(self.species_version_map), "hg38")

    def version_for(self, species: str) -> str:
        """Return the configured version for ``species`` or ``def``."""
        return self.species_version_map.get(species, "def")

    def prom_enha_cutoffs_for(self, species: str, version: str) -> str:
        """Return species/version-specific cutoffs, preserving case."""
        for key in (f"{species}_{version}", f"{species}:{version}"):
            if key in self.prom_enha_overrides:
                return self.prom_enha_overrides[key]
        return self.prom_enha_cutoffs

    @property
    def species_version_map(self) -> Dict[str, str]:
        """Parse ``species:version`` pairs from the configured value."""
        result: Dict[str, str] = {}
        for item in self.species_versions.replace(";", ",").split(","):
            item = item.strip()
            if not item:
                continue
            if ":" not in item:
                raise ValueError(
                    "SJCAB_PEAK2ANNO_SPECIES_VERSIONS must contain species:version pairs"
                )
            species, version = (part.strip() for part in item.split(":", 1))
            if species and version:
                result[species] = version
        return result


def rc_paths() -> list[Path]:
    """Return rc-file candidates in precedence order."""
    explicit = os.environ.get("SJCAB_PEAK2ANNO_CONFIG")
    if explicit:
        return [Path(explicit).expanduser()]
    xdg_home = Path(os.environ.get("XDG_CONFIG_HOME") or "~/.config").expanduser()
    return [
        xdg_home / "sjcab_peak2anno" / ".sjcab_peak2anno.rc",
        Path("~/.sjcab_peak2anno.rc").expanduser(),
    ]


def _read_rc(path: Path) -> Dict[str, str]:
    """Read simple ``KEY=VALUE`` rc syntax."""
    values: Dict[str, str] = {}
    if not path.is_file():
        return values
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("export "):
                line = line[7:].lstrip()
            if "=" not in line:
                raise ValueError(f"Invalid rc setting at {path}:{line_number}: {raw.rstrip()}")
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()
            try:
                tokens = shlex.split(value, comments=True)
            except ValueError as exc:
                raise ValueError(f"Invalid rc setting at {path}:{line_number}: {exc}") from exc
            values[key] = " ".join(tokens) if tokens else ""
    return values


def _ensure_rc_template(path: Path) -> None:
    """Create or extend an RC file with commented supported settings."""
    try:
        existing = path.read_text(encoding="utf-8") if path.is_file() else ""
        present = set()
        for raw in existing.splitlines():
            line = raw.strip().lstrip("#").strip()
            if "=" in line:
                present.add(line.split("=", 1)[0].strip())
        missing = [
            f"#{key}={value}"
            for key, value in RC_DEFAULTS.items()
            if key not in present
        ]
        if not any(key.startswith("SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS_") for key in present):
            missing.append(
                "#SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS_<species>_<version>=2kb,50kb,2kb"
            )
        if not missing:
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        prefix = "" if not existing or existing.endswith("\n") else "\n"
        with path.open("a", encoding="utf-8") as handle:
            handle.write(prefix)
            if not existing:
                handle.write("# peak2anno configuration; uncomment settings to override defaults.\n")
            handle.write("\n".join(missing) + "\n")
    except OSError:
        # A read-only home should not prevent command execution.
        return


def load_settings() -> Settings:
    """Load defaults, the first available rc file, then environment values."""
    values: Dict[str, str] = {
        "SJCAB_PEAK2ANNO_SPECIES_VERSIONS": DEFAULT_SPECIES_VERSIONS,
        "SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS": DEFAULT_PROM_ENHA_CUTOFFS,
        "SJCAB_PEAK2ANNO_GENE_TYPE": DEFAULT_GENE_TYPE,
        "SJCAB_PEAK2ANNO_ISO_SET": DEFAULT_ISO_SET,
        "SJCAB_PEAK2ANNO_2FEATURE_OUT": DEFAULT_2FEATURE_OUT,
        "SJCAB_PEAK2ANNO_2STATE_OUT": DEFAULT_2STATE_OUT,
        "SJCAB_PEAK2ANNO_TXT_DELIMITER": DEFAULT_TXT_DELIMITER,
        "SJCAB_PEAK2ANNO_BACKEND": DEFAULT_BACKEND,
    }
    candidates = rc_paths()
    rc_file = next((path for path in candidates if path.is_file()), candidates[0])
    _ensure_rc_template(rc_file)
    if rc_file.is_file():
        values.update(_read_rc(rc_file))
    for key in tuple(values) + ("SJCAB_PEAK2ANNO_DB_PATH",):
        if key in os.environ:
            values[key] = os.environ[key]
    cutoff_prefix = "SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS_"
    for key, value in os.environ.items():
        if key.startswith(cutoff_prefix):
            values[key] = value

    settings = Settings(
        db_path=values.get("SJCAB_PEAK2ANNO_DB_PATH") or None,
        species_versions=values["SJCAB_PEAK2ANNO_SPECIES_VERSIONS"],
        prom_enha_cutoffs=values["SJCAB_PEAK2ANNO_PROM_ENHA_CUTOFFS"],
        prom_enha_overrides={
            key[len(cutoff_prefix):]: value
            for key, value in values.items()
            if key.startswith(cutoff_prefix)
        },
        gene_type=values["SJCAB_PEAK2ANNO_GENE_TYPE"],
        iso_set=values["SJCAB_PEAK2ANNO_ISO_SET"],
        feature_out=values["SJCAB_PEAK2ANNO_2FEATURE_OUT"],
        state_out=values["SJCAB_PEAK2ANNO_2STATE_OUT"],
        txt_delimiter=values["SJCAB_PEAK2ANNO_TXT_DELIMITER"],
        backend=values["SJCAB_PEAK2ANNO_BACKEND"],
    )
    if settings.iso_set not in {"all", "deduplong"}:
        raise ValueError("SJCAB_PEAK2ANNO_ISO_SET must be all or deduplong")
    for name, value in (("SJCAB_PEAK2ANNO_2FEATURE_OUT", settings.feature_out), ("SJCAB_PEAK2ANNO_2STATE_OUT", settings.state_out)):
        modes = value.split(",")
        if not modes or any(mode not in {"max", "percent"} for mode in modes) or len(set(modes)) != len(modes):
            raise ValueError(f"{name} must be max, percent, max,percent, or percent,max")
    # Validate the mapping while reporting configuration errors before parsing CLI args.
    settings.species_version_map
    return settings
