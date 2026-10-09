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

from metawrap2.limits import LIMITS_KEYS, LimitsConfig, LimitsError, config_from_table

if sys.version_info >= (3, 11):
    import tomllib as _toml
else:  # pragma: no cover - exercised on <3.11
    import tomli as _toml

# One conda env per *workflow*, not per module or per tool. Several modules share an env whenever
# they run the same tools or run in sequence, which is most of them: the alternative - one env per
# module - meant solving the same aligner/samtools graph five times and keeping five copies of it
# current, for no benefit a user ever sees.
#
# The grouping was chosen by trying it: every combination below was solved against conda-forge +
# bioconda before being written down, and the two places the pipeline is split in an awkward-looking
# way (CheckM2 away from CheckM1, GTDB-Tk away from taxator-tk) are splits the solver forced, not
# preferences. See envs/<name>.yaml, where each file says what it holds and why.
MODULE_ENVS: Dict[str, str] = {
    # Everything that touches reads: QC, assembly, re-assembly, and the blob plot. They share the
    # aligner stack and SPAdes, and in a normal study they run one after another.
    "read_qc": "metawrap2-reads",
    "assembly": "metawrap2-reads",
    "reassemble_bins": "metawrap2-reads",
    "blobology": "metawrap2-reads",
    # Everything that turns contigs into a MAG catalogue: binning, cross-sample coverage,
    # dereplication, abundance. coverage shares binning's aligner and samtools; dereplicate and
    # quant_bins add three single static binaries (skani, sylph, salmon) that cost nothing to carry.
    "binning": "metawrap2-bins",
    "coverage": "metawrap2-bins",
    "dereplicate": "metawrap2-bins",
    "quant_bins": "metawrap2-bins",
    # Scoring and consolidating bins: CheckM1, GUNC, DAS Tool. GUNC runs by default here, so it was
    # never genuinely optional.
    "bin_refinement": "metawrap2-binqc",
    # Read- and bin-level taxonomy.
    "classify_bins": "metawrap2-taxonomy",
    "kraken2": "metawrap2-taxonomy",
    # Gene calling: Prokka and Bakta solve together, so --bakta needs no second install.
    "annotate_bins": "metawrap2-annotate",
}

#: Every env a full install builds, each once. MODULE_ENVS has more keys than values now, so anything
#: that iterates it to build or check envs has to go through this or it does the same solve five times.
ALL_MODULE_ENVS = sorted(set(MODULE_ENVS.values()))

# Envs that are still separate, and why each one is. Two of these are forced by the solver and the
# third by install size; none of them is a second copy of a tool that could have shared a home.
OPTIONAL_ENVS: Dict[str, str] = {
    # Benchmark-only: AMBER is reached by `metawrap2 benchmark score`, never by a module. It also
    # pins matplotlib below 3.9, which nothing else here wants.
    "amber": "metawrap2-benchmark-amber",
    # Forced. CheckM2 has no common solution with the CheckM1 + GUNC + DAS Tool set.
    "checkm2": "metawrap2-checkm2",
    # The three learned binners, together, because they are all pytorch and solve as one. Separate
    # from metawrap2-bins on size: pytorch is ~1.5 GB and nobody who only wants metaBAT2 should wait
    # for it. The env's pytorch is CUDA-capable, so `--device cuda` needs no second install.
    "semibin2": "metawrap2-learned",
    "vamb": "metawrap2-learned",
    "comebin": "metawrap2-learned",
    # Forced. GTDB-Tk has no common solution with taxator-tk, and it is the one tool here with a
    # 100+ GB reference database, so a user who wants it is opting in already.
    "gtdbtk": "metawrap2-gtdbtk",
}

# Tools that live in an OPTIONAL_ENVS environment rather than a module's own, keyed the same way
# OPTIONAL_ENVS is. Kept separate because `doctor` probes a tool inside the env that is supposed to
# hold it: listing gunc under bin_refinement would report every healthy install as broken, since
# bin_refinement's env holds CheckM and nothing else.
OPTIONAL_ENV_TOOLS: Dict[str, List[str]] = {
    "amber": ["amber.py"],
    "checkm2": ["checkm2"],
    # All three learned binners share one env, so `doctor` probes all three whichever key it is
    # asked about: a half-built metawrap2-learned should report as broken, not as fine.
    # "SemiBin2" is capitalised and is the only entry point the 2.5.0 package installs; COMEBin's
    # entry point is a bash script, name and extension included - there is no "comebin" on PATH.
    "semibin2": ["SemiBin2", "vamb", "run_comebin.sh"],
    "vamb": ["SemiBin2", "vamb", "run_comebin.sh"],
    "comebin": ["SemiBin2", "vamb", "run_comebin.sh"],
    "gtdbtk": ["gtdbtk"],
}

# The same two tables keyed by the yaml name an opt-in env is built from (`bin_refinement-gunc`)
# rather than by the module flag. That name is what a user types at `install-env` and `doctor`, and
# it is not a key of MODULE_ENVS - so without these, both commands fall back to "a module I have
# never heard of", which for doctor meant reporting a freshly built env as "0/0 tools OK" without
# probing anything in it.
#: Which module's flag selects each opt-in env. Stated rather than parsed out of the env name: the
#: names used to be ``metawrap2-<module>-<tool>``, so `doctor` split on the dash to find the owner,
#: and consolidating the envs left ``metawrap2-learned`` with no second dash - which did not fail,
#: it silently reported no opt-in envs at all. An env with no module (AMBER, reached only by
#: `benchmark score`) maps to None and is skipped.
OPTIONAL_ENV_OWNERS: Dict[str, Optional[str]] = {
    "amber": None,
    "checkm2": "bin_refinement",
    "semibin2": "binning",
    "vamb": "binning",
    "comebin": "binning",
    "gtdbtk": "classify_bins",
}

OPTIONAL_ENVS_BY_MODULE: Dict[str, str] = {
    name[len("metawrap2-") :]: name for name in OPTIONAL_ENVS.values()
}
OPTIONAL_TOOLS_BY_MODULE: Dict[str, List[str]] = {
    OPTIONAL_ENVS[flag][len("metawrap2-") :]: tools
    for flag, tools in OPTIONAL_ENV_TOOLS.items()
    if flag in OPTIONAL_ENVS
}

# External tools each module invokes (used by the preflight check).
MODULE_TOOLS: Dict[str, List[str]] = {
    "read_qc": ["trim_galore", "fastp", "fastqc", "bmtagger.sh", "bowtie2"],
    "assembly": ["metaspades.py", "megahit", "quast", "bwa", "bwa-mem2", "minimap2"],
    # bowtie2 was listed here but the module never invoked it, and the aligner benchmark gave no
    # reason to start: it costs what bwa costs and agrees with it less well than minimap2.
    "binning": ["metabat2", "run_MaxBin.pl", "concoct", "minimap2", "bwa", "bwa-mem2", "samtools"],
    # The binning tools it shares, minus the binners: it measures coverage and nothing else.
    # fairy moved in with the aligners: it is one static binary and the route it serves
    # (--coverage fairy) should not need a second install command to try.
    "coverage": ["minimap2", "bwa", "bwa-mem2", "samtools", "fairy"],
    # gunc and DAS_Tool moved in with CheckM1: GUNC runs by default unless --skip-gunc, and
    # DAS Tool is an alternative consolidator for this same step.
    "bin_refinement": ["checkm", "gunc", "DAS_Tool"],
    "reassemble_bins": ["spades.py", "bwa", "samtools", "minimap2"],
    "quant_bins": ["sylph", "salmon", "coverm"],
    "dereplicate": ["skani"],
    "classify_bins": ["blastn", "taxator", "binner", "taxknife"],
    "annotate_bins": ["prokka", "bakta"],
    "kraken2": ["kraken2", "ktImportText"],
    # Same aligner set as binning, for the same reason: blobology maps reads back only to get
    # per-contig abundance for the plot's y-axis, so it shares the registry in metawrap2.aligners.
    "blobology": ["blastn", "minimap2", "bwa", "bwa-mem2", "samtools"],
}

# Database keys and which module(s) use them.
DB_KEYS = {
    "KRAKEN2_DB": "kraken2",
    "BMTAGGER_DB": "read_qc (host removal)",
    "BLASTDB": "blobology / classify_bins",
    "BLASTDB_NAME": "blobology / classify_bins (db name inside BLASTDB; default nt)",
    "TAXDUMP": "blobology / classify_bins",
    "CHECKM_DB": "bin_refinement (CheckM1 reference data)",
    "CHECKM2DB": "bin_refinement (--checkm2 DIAMOND reference database)",
    "GUNC_DB": "bin_refinement (GUNC reference database, unless --skip-gunc)",
    "BAKTA_DB": "annotate_bins (--bakta)",
    "GTDBTK_DATA_PATH": "classify_bins (--gtdbtk)",
}

# Databases only an opt-in code path needs; absent ones are not counted as a failure by
# `metawrap2 check` unless that path was actually requested.
OPTIONAL_DB_KEYS = {
    "BAKTA_DB",
    "GTDBTK_DATA_PATH",
    "CHECKM_DB",
    "CHECKM2DB",
    "BLASTDB_NAME",
    "GUNC_DB",
}

# Keys whose value is a setting, not a filesystem path - checking them for existence is
# meaningless and reporting them as MISSING is just confusing.
NON_PATH_DB_KEYS = {"BLASTDB_NAME"}

#: Database keys whose value names a *file*, not a directory holding files.
#:
#: GUNC takes ``--db_file`` and nothing else, so install-db correctly points GUNC_DB at the
#: ``.dmnd`` itself. Everything that checked a configured database assumed a directory, so a perfectly
#: good GUNC install reported as "MISSING ... not a directory" in `doctor` for as long as it existed -
#: a false alarm on the one database whose correct value could not be a directory.
FILE_DB_KEYS = {"GUNC_DB"}


class ConfigError(Exception):
    """The config file cannot be used as written."""


@dataclass
class Settings:
    databases: Dict[str, str] = field(default_factory=dict)
    #: The per-module default thread count, or None when the config does not set one - in which
    #: case the answer is the ``[limits]`` budget. This has to be None rather than 1: with 1 as
    #: the default, "the user did not configure threads" and "the user asked for one thread" are
    #: the same value, so the budget could never be used as the fallback.
    threads: Optional[int] = None
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
    #: The whole run's resource budget - threads, memory, disk ceiling, GPUs. ``threads`` above is
    #: still the per-module default a user asks for; this is the ceiling it is capped to, and the
    #: pool a cohort driver divides among concurrent samples. See metawrap2.limits.
    limits: LimitsConfig = field(default_factory=LimitsConfig)
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
    "use_conda_envs": "run tools inside the metawrap2-* conda envs (true) or off PATH (false)",
    "scratch_dir": "where modules put temporary space",
    "skip_space_check": "skip the pre-flight disk-space estimate",
    "strict_space_check": "let a failed disk-space estimate stop the run, not just warn",
}

#: The top-level tables a config may contain.
CONFIG_TABLES = {
    "settings": "global defaults - see SETTINGS_KEYS",
    "databases": "database paths - see DB_KEYS",
    "limits": "the run's resource budget - see metawrap2.limits.LIMITS_KEYS",
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

    * An unknown **table**, or an unknown ``[settings]`` or ``[limits]`` key, is a problem. That
      namespace is small, entirely MetaWrap2's own, and a typo there means the setting silently
      does not apply - indistinguishable from never having written it. A did-you-mean makes it a
      five-second fix. ``[limits]`` values are checked too, not just their names, because a
      malformed budget is worth catching in the first second rather than the first hour.
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
    limits = data.get("limits")
    if isinstance(limits, dict):
        for key in limits:
            if key not in LIMITS_KEYS:
                problems.append(
                    "unknown key 'limits.%s'%s%s. Known keys: %s"
                    % (key, where, _did_you_mean(key, LIMITS_KEYS), ", ".join(LIMITS_KEYS))
                )
        try:
            config_from_table(limits)
        except LimitsError as exc:
            # The message already names the key, so the file goes in front rather than trailing
            # after the offending value, where it reads as part of it.
            problems.append("%s%s" % (("%s: " % path) if path else "", exc))
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
    try:
        limits = config_from_table(dict(data.get("limits", {})))
    except LimitsError:
        # Only reachable with strict=False, where the problem was already reported above. Falling
        # back to the default budget is safer than propagating: the caller asked to run anyway,
        # and the defaults are conservative.
        limits = LimitsConfig()
    configured_threads = data.get("settings", {}).get("threads")
    return Settings(
        databases=dict(data.get("databases", {})),
        threads=int(configured_threads) if configured_threads is not None else None,
        use_conda_envs=bool(data.get("settings", {}).get("use_conda_envs", True)),
        scratch_dir=str(data.get("settings", {}).get("scratch_dir", "") or ""),
        skip_space_check=bool(data.get("settings", {}).get("skip_space_check", False)),
        strict_space_check=bool(data.get("settings", {}).get("strict_space_check", False)),
        limits=limits,
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
