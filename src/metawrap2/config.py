"""MetaWrap2 configuration and environment preflight.

Configuration lives in a small TOML file (``metawrap2.toml``): database paths, a few
settings, and (optionally) overrides for any surfaced module constant. Resolution order
for the config file: an explicit ``--config`` path, then ``$METAWRAP2_CONFIG``, then
``~/.metawrap2/config.toml``, then the packaged default. Values a user does not set fall
back to the in-code module defaults.

Also provides the preflight check: report missing tools/databases (and whether a module's
conda env exists) *before* a long run starts, instead of failing hours in.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

try:  # Python 3.11+
    import tomllib as _toml
except ModuleNotFoundError:  # pragma: no cover - exercised on <3.11
    import tomli as _toml  # type: ignore

# One conda env per module (see envs/<module>.yaml), holding only that module's tools.
MODULE_ENVS: Dict[str, str] = {
    "read_qc": "metawrap2-read_qc",
    "assembly": "metawrap2-assembly",
    "binning": "metawrap2-binning",
    "bin_refinement": "metawrap2-bin_refinement",
    "reassemble_bins": "metawrap2-reassemble_bins",
    "quant_bins": "metawrap2-quant_bins",
    "classify_bins": "metawrap2-classify_bins",
    "annotate_bins": "metawrap2-annotate_bins",
    "kraken2": "metawrap2-kraken2",
    "blobology": "metawrap2-blobology",
}

# External tools each module invokes (used by the preflight check).
MODULE_TOOLS: Dict[str, List[str]] = {
    "read_qc": ["trim_galore", "fastqc", "bmtagger.sh", "bowtie2"],
    "assembly": ["metaspades.py", "megahit", "quast", "bwa"],
    "binning": ["metabat2", "run_MaxBin.pl", "concoct", "bwa", "bowtie2", "samtools"],
    "bin_refinement": ["checkm"],
    "reassemble_bins": ["spades.py", "bwa", "samtools", "minimap2"],
    "quant_bins": ["salmon"],
    "classify_bins": ["blastn", "taxator", "binner", "taxknife"],
    "annotate_bins": ["prokka"],
    "kraken2": ["kraken2", "ktImportText"],
    "blobology": ["blastn", "bowtie2", "samtools"],
}

# Database keys and which module(s) use them.
DB_KEYS = {
    "KRAKEN2_DB": "kraken2",
    "BMTAGGER_DB": "read_qc (host removal)",
    "BLASTDB": "blobology",
    "TAXDUMP": "blobology / classify_bins",
}


@dataclass
class Settings:
    databases: Dict[str, str] = field(default_factory=dict)
    threads: int = 1
    use_conda_envs: bool = True
    overrides: Dict[str, Dict[str, object]] = field(default_factory=dict)

    def db(self, key: str) -> str:
        return self.databases.get(key, "")

    def module_override(self, module: str) -> Dict[str, object]:
        return self.overrides.get(module, {})


def find_config(explicit: Optional[str] = None) -> Optional[str]:
    """Locate the metawrap2.toml to use, or None if only in-code defaults apply."""
    candidates = [
        explicit,
        os.environ.get("METAWRAP2_CONFIG"),
        os.path.expanduser("~/.metawrap2/config.toml"),
    ]
    for path in candidates:
        if path and os.path.isfile(path):
            return path
    return None


def load_settings(explicit: Optional[str] = None) -> Settings:
    """Load Settings from the resolved config file (or defaults if none found)."""
    path = find_config(explicit)
    if not path:
        return Settings()
    with open(path, "rb") as fh:
        data = _toml.load(fh)
    return Settings(
        databases=dict(data.get("databases", {})),
        threads=int(data.get("settings", {}).get("threads", 1)),
        use_conda_envs=bool(data.get("settings", {}).get("use_conda_envs", True)),
        overrides={k: dict(v) for k, v in data.get("modules", {}).items()},
    )


def conda_env_exists(env: str) -> bool:
    """True if a conda env named *env* appears to exist."""
    base = os.environ.get("CONDA_PREFIX", "")
    roots = []
    if base:
        # .../envs/<name> or the base prefix itself
        roots.append(os.path.join(os.path.dirname(base), env))
        roots.append(os.path.join(base, "envs", env))
    home_envs = os.path.expanduser("~/miniconda3/envs/%s" % env)
    roots.append(home_envs)
    return any(os.path.isdir(r) for r in roots)


def check_databases(settings: Settings) -> List[Tuple[str, str, bool, str]]:
    """Return (key, path, present, used_by) for each known database path."""
    results = []
    for key, used_by in DB_KEYS.items():
        path = settings.db(key)
        present = bool(path) and os.path.exists(path)
        results.append((key, path, present, used_by))
    return results


def check_tools(modules: List[str]) -> List[Tuple[str, str, bool]]:
    """Return (tool, module, present_on_PATH) for the tools used by the given modules."""
    results = []
    seen = set()
    for module in modules:
        for tool in MODULE_TOOLS.get(module, []):
            if (tool, module) in seen:
                continue
            seen.add((tool, module))
            results.append((tool, module, shutil.which(tool) is not None))
    return results
