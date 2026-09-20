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

import difflib
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Tuple

try:  # Python 3.11+
    import tomllib as _toml
except ModuleNotFoundError:  # pragma: no cover - exercised on <3.11
    import tomli as _toml  # only reachable on Python < 3.11

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

# Opt-in alternative envs, keyed by the module flag that selects them. These hold a modern
# replacement for the module's default tool; scores/labels differ from the default, so they
# are never selected automatically. See envs/<module>-<tool>.yaml.
OPTIONAL_ENVS: Dict[str, str] = {
    "checkm2": "metawrap2-bin_refinement-checkm2",
    "gtdbtk": "metawrap2-classify_bins-gtdbtk",
    "bakta": "metawrap2-annotate_bins-bakta",
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
    "BLASTDB": "blobology / classify_bins",
    "BLASTDB_NAME": "blobology / classify_bins (db name inside BLASTDB; default nt)",
    "TAXDUMP": "blobology / classify_bins",
    "CHECKM_DB": "bin_refinement (CheckM1 reference data)",
    "BAKTA_DB": "annotate_bins (--bakta)",
    "GTDBTK_DATA_PATH": "classify_bins (--gtdbtk)",
}

# Databases only an opt-in code path needs; absent ones are not counted as a failure by
# `metawrap2 check` unless that path was actually requested.
OPTIONAL_DB_KEYS = {"BAKTA_DB", "GTDBTK_DATA_PATH", "CHECKM_DB", "BLASTDB_NAME"}

# Keys whose value is a setting, not a filesystem path - checking them for existence is
# meaningless and reporting them as MISSING is just confusing.
NON_PATH_DB_KEYS = {"BLASTDB_NAME"}


class ConfigError(Exception):
    """The config file cannot be used as written."""


@dataclass
class Settings:
    databases: Dict[str, str] = field(default_factory=dict)
    threads: int = 1
    use_conda_envs: bool = True
    #: Where modules put temporary space. Empty means "inside the run's output directory".
    #: On a cluster this matters: /tmp is often small and node-local, and a metaSPAdes
    #: intermediate can be hundreds of GB. See metawrap2.scratch.
    scratch_dir: str = ""
    #: Set true to skip the pre-flight disk-space estimate.
    skip_space_check: bool = False
    #: Set true to make a failed disk-space estimate stop the run rather than warn. Worth doing
    #: on shared infrastructure, where one run filling the filesystem is everyone's problem.
    strict_space_check: bool = False
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


#: The keys ``[settings]`` understands. A key not in here was silently ignored before, so a
#: typo like ``use_conda_env`` looked like it worked and then did nothing.
SETTINGS_KEYS = {
    "threads": "default thread count for every module",
    "use_conda_envs": "run tools inside metawrap2-<module> envs (true) or off PATH (false)",
    "scratch_dir": "where modules put temporary space",
    "skip_space_check": "skip the pre-flight disk-space estimate",
    "strict_space_check": "let a failed disk-space estimate stop the run, not just warn",
}

#: The top-level tables a config may contain.
CONFIG_TABLES = {
    "settings": "global defaults - see SETTINGS_KEYS",
    "databases": "database paths - see DB_KEYS",
    "modules": "per-module option overrides, e.g. [modules.binning]",
}


def _did_you_mean(name: str, known: Iterable[str]) -> str:
    """`` (did you mean 'x'?)`` for the closest known name, or "" if nothing is close."""
    matches = difflib.get_close_matches(name, list(known), n=1, cutoff=0.6)
    return " (did you mean %r?)" % matches[0] if matches else ""


#: Database keys that are recognised but no longer used, and what to do instead. These are worth
#: naming rather than treating as a typo: they are the leftovers in a config carried forward from
#: an older setup, and they are harmless, so noticing them must not stop a run.
RETIRED_DB_KEYS = {
    "KRAKEN_DB": "the KRAKEN1 module no longer exists - use KRAKEN2_DB with the kraken2 module",
}


def check_config(data: Dict[str, Any], path: str = "") -> Tuple[List[str], List[str]]:
    """Check a parsed config. Returns ``(problems, notes)``.

    *problems* stop a run; *notes* are printed and then ignored. Which is which follows from what
    the mistake costs:

    * An unknown **table**, or an unknown ``[settings]`` key, is a problem. That namespace is
      small, entirely MetaWrap2's own, and a typo there means the setting silently does not apply -
      indistinguishable from never having written it. A did-you-mean makes it a five-second fix.
    * An unknown or retired **database** key is only a note. ``[databases]`` is deliberately open,
      so a locally patched module may look up a key nothing here has heard of; and a leftover key
      from an older configuration is harmless. Refusing to run ``binning`` because the config
      still mentions a database that only a removed module ever read would be absurd.

    Reporting everything at once matters for the same reason it does for inputs: one round trip,
    not one finding per re-run.
    """
    problems: List[str] = []
    notes: List[str] = []
    where = " in %s" % path if path else ""
    for table in data:
        if table not in CONFIG_TABLES:
            problems.append(
                "unknown table [%s]%s%s. Known tables: %s"
                % (table, where, _did_you_mean(table, CONFIG_TABLES), ", ".join(CONFIG_TABLES))
            )
    settings = data.get("settings")
    if isinstance(settings, dict):
        for key in settings:
            if key not in SETTINGS_KEYS:
                problems.append(
                    "unknown key 'settings.%s'%s%s. Known keys: %s"
                    % (key, where, _did_you_mean(key, SETTINGS_KEYS), ", ".join(SETTINGS_KEYS))
                )
    databases = data.get("databases")
    if isinstance(databases, dict):
        known = set(DB_KEYS) | OPTIONAL_DB_KEYS
        for key in databases:
            if key in known:
                continue
            if key in RETIRED_DB_KEYS:
                notes.append(
                    "%r%s is no longer used: %s. You can delete the line."
                    % (key, where, RETIRED_DB_KEYS[key])
                )
                continue
            suggestion = _did_you_mean(key, known)
            if suggestion:  # close to a real key, so probably a typo rather than a custom key
                notes.append(
                    "%r%s is not a database key MetaWrap2 reads%s - it will be ignored."
                    % (key, where, suggestion)
                )
    return problems, notes


def load_settings(explicit: Optional[str] = None, strict: bool = True) -> Settings:
    """Load Settings from the resolved config file (or defaults if none found).

    With *strict* (the default), a config containing unknown tables or ``[settings]`` keys is an
    error rather than being silently half-applied. Pass ``strict=False`` to load it anyway.
    """
    path = find_config(explicit)
    if not path:
        return Settings()
    with open(path, "rb") as fh:
        data = _toml.load(fh)
    problems, notes = check_config(data, path)
    for note in notes:
        # A note is not worth failing over, but it is worth saying once, where the user will see
        # it next to the run it affects.
        print("NOTE: %s" % note, file=sys.stderr)
    if problems and strict:
        raise ConfigError(
            "\n".join(
                ["%d problem(s) in your MetaWrap2 config:" % len(problems)]
                + ["  - " + p for p in problems]
            )
        )
    return Settings(
        databases=dict(data.get("databases", {})),
        threads=int(data.get("settings", {}).get("threads", 1)),
        use_conda_envs=bool(data.get("settings", {}).get("use_conda_envs", True)),
        scratch_dir=str(data.get("settings", {}).get("scratch_dir", "") or ""),
        skip_space_check=bool(data.get("settings", {}).get("skip_space_check", False)),
        strict_space_check=bool(data.get("settings", {}).get("strict_space_check", False)),
        overrides={k: dict(v) for k, v in data.get("modules", {}).items()},
    )


def env_dirs() -> List[str]:
    """Directories conda/mamba keep environments in, most authoritative first.

    Asking the solver (``conda env list``) is the only reliable answer, because envs can
    live anywhere ``envs_dirs`` points - not just under the active prefix. The path guesses
    are a fallback for when neither conda nor mamba is on PATH.
    """
    dirs: List[str] = []
    for prog in ("conda", "mamba", "micromamba"):
        exe = shutil.which(prog)
        if not exe:
            continue
        try:
            out = subprocess.run(
                [exe, "env", "list"],
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
                text=True,
                timeout=60,
                check=False,
            )
        except (OSError, subprocess.SubprocessError):
            continue
        for line in (out.stdout or "").splitlines():
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            path = line.replace("*", " ").split()[-1]
            if os.path.isdir(path):
                dirs.append(
                    os.path.dirname(path)
                    if os.path.basename(os.path.dirname(path)) == "envs"
                    else path
                )
        if dirs:
            break

    prefix = os.environ.get("CONDA_PREFIX", "")
    if prefix:
        dirs.append(os.path.join(prefix, "envs"))
        if os.path.basename(os.path.dirname(prefix)) == "envs":
            dirs.append(os.path.dirname(prefix))
    for guess in (
        "~/miniforge3/envs",
        "~/miniconda3/envs",
        "~/mambaforge/envs",
        "~/anaconda3/envs",
        "~/micromamba/envs",
    ):
        dirs.append(os.path.expanduser(guess))

    seen, unique = set(), []
    for d in dirs:
        if d not in seen:
            seen.add(d)
            unique.append(d)
    return unique


def conda_env_exists(env: str) -> bool:
    """True if a conda env named *env* exists in any of conda's environment directories."""
    return any(os.path.isdir(os.path.join(d, env)) for d in env_dirs())


def check_databases(settings: Settings) -> List[Tuple[str, str, bool, str]]:
    """Return (key, value, present, used_by) for each known database key.

    For a path-valued key, ``present`` means the path exists. For a setting-valued key (see
    :data:`NON_PATH_DB_KEYS`) it just means a value was supplied.
    """
    results = []
    for key, used_by in DB_KEYS.items():
        value = settings.db(key)
        if key in NON_PATH_DB_KEYS:
            present = bool(value)
        else:
            present = bool(value) and os.path.exists(value)
        results.append((key, value, present, used_by))
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
