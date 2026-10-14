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
from typing import Any, Dict, Iterable, List, Mapping, Optional, Set, Tuple

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
    # Gene calling: Prokka and Bakta solve together, so --bakta needs no second install. `annotate`
    # calls genes in-process with Pyrodigal and only needs this env when --bakta is given.
    "annotate_bins": "metawrap2-annotate",
    "annotate": "metawrap2-annotate",
    # Only for Krona's ktImportText, which draws the hierarchy plots. Everything pathways computes
    # is in-process over vendored KEGG data, so this env is for the picture, not the numbers.
    "pathways": "metawrap2-taxonomy",
    # AMRFinderPlus, hAMRonization and argNorm. RGI is not in here; see envs/amr.yaml.
    "amr": "metawrap2-amr",
    # geNomad and CheckV. One env because they are two halves of one answer and are always run
    # together; see envs/viruses.yaml for the numpy pin they need.
    "viruses": "metawrap2-viruses",
    # SingleM gets its own env: it pins a narrow dependency set and is the only tool here that
    # needs the GTDB marker metapackage, so a user who never asks "what did my bins miss" should
    # not carry it.
    "appraise": "metawrap2-singlem",
    # Nonpareil lives with the reads: it is the only other thing in the pipeline that reads a raw
    # read file and nothing else, and in a study it runs beside read_qc. It brings r-base along
    # (bioconda builds it against R, which carries its model-fitting package), which costs install
    # size and nothing else - the estimates here are computed in Python, see metawrap2.nonpareil.
    "nonpareil": "metawrap2-reads",
    # Pangenomes and phylogenies over a cohort's MAGs: they share a pipeline (one writes the
    # alignment the other reads), so they share an env.
    "pangenome": "metawrap2-comparative",
    "phylogeny": "metawrap2-comparative",
    # inStrain, and the aligners that feed it. Its own env on a pin: inStrain resolves to Python
    # 3.8, which nothing else here can be built against.
    "strains": "metawrap2-strains",
    # MMseqs2 comes with the annotate env already - Bakta pulls it in - so the bin-free path needs
    # no extra install.
    "catalogue": "metawrap2-annotate",
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
    # TRACS, the second strain backend. Its own env on a pin, like RGI's: it resolves an older
    # scientific-python stack than inStrain's env holds.
    "tracs": "metawrap2-tracs",
    # Forced. RGI pins older dependencies than AMRFinderPlus and hAMRonization resolve to, so the
    # three cannot share an env without downgrading the default backend to satisfy the opt-in one.
    "rgi": "metawrap2-rgi",
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
    "rgi": ["rgi"],
    "tracs": ["tracs"],
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
    "rgi": "amr",
    "tracs": "strains",
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
    # flye is listed: --long-reads is a documented path of this module, and a check that passes
    # without it would send a long-read user into a run that dies at the assembler.
    "assembly": ["metaspades.py", "megahit", "quast", "bwa", "bwa-mem2", "minimap2", "flye"],
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
    # `annotate` needs nothing external for gene calling - Pyrodigal is a Python dependency - so
    # only the optional --bakta path has a tool to check for.
    # run_dbcan (--cazy) and gecco (--bgc) ship in this env too. GECCO needs no database, which is
    # why it has no DB_KEYS entry: an install of this env can run it immediately.
    "annotate": ["bakta", "run_dbcan", "gecco"],
    "nonpareil": ["nonpareil"],
    "strains": ["inStrain", "minimap2", "samtools"],
    "pangenome": ["ppanggolin"],
    "phylogeny": ["iqtree", "fasttree"],
    "appraise": ["singlem"],
    "pathways": ["ktImportText"],
    "amr": ["amrfinder"],
    # cctyper is in this env too but is not listed: it only runs for --crispr-hosts, and a check
    # that fails for a user who never asked for host prediction is a check that gets ignored.
    "viruses": ["genomad", "checkv"],
    "catalogue": ["mmseqs"],
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
    "BAKTA_DB": "annotate / annotate_bins (--bakta)",
    # The KOfam HMM profiles (profiles.tar.gz, 1.5 GB) and their per-profile thresholds (ko_list),
    # both from `metawrap2 install-db kofam`. The thresholds decide what counts as a hit.
    "KOFAM_DB": "annotate (--kofam), the KOfam directory: profiles/ and ko_list",
    # BRITE ko00001, projected to one row per (KO, pathway), from `metawrap2 install-db kegg`.
    "KEGG_BRITE": "pathways (the functional hierarchy and Krona plots)",
    # dbCAN's HMMs, DIAMOND databases and substrate mapping. Fetched from dbCAN's pinned S3 release
    # rather than its `db_current` HTTP path, which is a moving snapshot: a CAZyme family call is
    # only reproducible against the release that made it.
    "DBCAN_DB": "annotate (--cazy), the dbCAN database directory",
    # Pfam-A's HMMs, for `annotate --pfam`. The release matters: Pfam renumbers and merges families
    # between releases, so a domain call is only reproducible against the release that made it -
    # which is why the installer reads the version out of relnotes.txt and the table records it.
    "PFAM_DB": "annotate (--pfam), the directory holding Pfam-A.hmm",
    # Point this at a dated release rather than the `latest` symlink amrfinder_update writes: an ARG
    # call is only as datable as its database, and a rerun against a moved symlink is a different
    # experiment wearing the same name.
    "AMRFINDER_DB": "amr (AMRFinderPlus reference database)",
    # geNomad's database also carries its plasmid markers, so `amr`'s mobility work and `viruses`
    # share this one download rather than needing two.
    "GENOMAD_DB": "viruses (geNomad marker database; also plasmid calls for amr)",
    "CHECKV_DB": "viruses (CheckV reference database)",
    # CARD's data, for `amr --rgi`. Versioned by CARD itself (4.0.2 as of writing) and released
    # every few months, so the directory name should carry the version rather than say `latest`.
    "CARD_DB": "amr (--rgi; CARD data directory holding card.json)",
    "GTDBTK_DATA_PATH": "classify_bins (--gtdbtk)",
    # The GTDB marker metapackage SingleM profiles against: 2.2 GB compressed, 9.7 GB unpacked.
    # Its release is part of the answer - a per-cell fraction is only as current as the reference
    # it was measured against - so point this at the versioned directory, not at a `latest` link.
    "SINGLEM_METAPACKAGE": "appraise (SingleM GTDB marker metapackage)",
}

# Databases only an opt-in code path needs; absent ones are not counted as a failure by
# `metawrap2 check` unless that path was actually requested.
OPTIONAL_DB_KEYS = {
    "BAKTA_DB",
    "KOFAM_DB",
    "KEGG_BRITE",
    "DBCAN_DB",
    "AMRFINDER_DB",
    "GENOMAD_DB",
    "CHECKV_DB",
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


#: Environment names from MetaWrap2 2.0.x, which had one env per tool rather than one per workflow.
#: Mapped rather than rejected because they are what the 2.0.x documentation says to type, what a
#: 2.0.x error message quotes back, and therefore what an upgrading user reaches for. Two issue
#: reports from 2.0.1 users named the first two of these exactly.
LEGACY_ENV_ALIASES: Dict[str, str] = {
    "annotate_bins-bakta": "metawrap2-annotate",
    "annotate_bins-prokka": "metawrap2-annotate",
    "bin_refinement-checkm2": "metawrap2-checkm2",
    "bin_refinement-gunc": "metawrap2-binqc",
    "bin_refinement-dastool": "metawrap2-binqc",
    "binning-semibin2": "metawrap2-learned",
    "binning-vamb": "metawrap2-learned",
    "binning-comebin": "metawrap2-learned",
    "classify_bins-gtdbtk": "metawrap2-gtdbtk",
}


def legacy_module_for(name: str) -> Optional[str]:
    """The module a 2.0.x env name belonged to, when that module still exists.

    The 2.0.x names are ``<module>-<tool>`` - ``bin_refinement-checkm2``,
    ``annotate_bins-bakta`` - so the module is the part before the first dash. Worth recovering
    because `check` and `doctor` are asked about *modules*: given only the environment, they have no
    tool list and no databases to look at, and report "OK" having checked nothing. Checking the
    module the name came from is the answer the user wanted.
    """
    if name not in LEGACY_ENV_ALIASES:
        return None
    module = name.split("-", 1)[0]
    return module if module in MODULE_ENVS else None


#: Every registry that maps something a user or a module can name onto a conda env. Listed in one
#: place because the bug this exists to prevent was caused by there being several, keyed differently,
#: with each consumer reading a different subset:
#:
#: * ``MODULE_ENVS`` is keyed by module (``binning``).
#: * ``OPTIONAL_ENVS`` is keyed by **tool** (``semibin2``, ``vamb``, ``comebin``, ``amber``) and is
#:   what ``optional_env_for`` reads at run time.
#: * ``OPTIONAL_ENVS_BY_MODULE`` is keyed by **env base name** (``learned``, ``benchmark-amber``).
#: * ``LEGACY_ENV_ALIASES`` is keyed by the 2.0.x per-tool env name.
#:
#: ``OPTIONAL_ENVS`` and ``OPTIONAL_ENVS_BY_MODULE`` overlap but do not agree: ``semibin2`` is a key
#: of the first and not the second, so a resolver reading only the second refused a name the running
#: pipeline itself uses. That is the same defect as the one 2.0.1 shipped - a name the runtime needs
#: and the installer cannot reach - and it is why resolution reads all of them rather than whichever
#: one the author of a given call site happened to remember.
def tools_for_name(name: str) -> List[str]:
    """The executables to verify when a user asks about *name*.

    A module's tools come from MODULE_TOOLS; a tool or opt-in env's come from OPTIONAL_ENV_TOOLS; a
    2.0.x alias defers to the module it belonged to. Empty only when nothing is known to check,
    which the caller should report as such rather than as success.
    """
    if name in MODULE_TOOLS:
        return list(MODULE_TOOLS[name])
    if name in OPTIONAL_ENV_TOOLS:
        return list(OPTIONAL_ENV_TOOLS[name])
    if name in OPTIONAL_TOOLS_BY_MODULE:
        return list(OPTIONAL_TOOLS_BY_MODULE[name])
    module = legacy_module_for(name)
    if module:
        return list(MODULE_TOOLS.get(module, []))

    # An env's own base name - "bins", "reads", "binqc", "taxonomy". Several modules share each of
    # those, so the tools to verify are the union over every module that uses it. Without this,
    # `check bins` resolved the env, found no tool list and reported OK having verified nothing -
    # the same vacuous answer a 2.0.x alias used to give.
    env = resolve_env_name(name)
    if env:
        shared: List[str] = []
        for module, module_env in MODULE_ENVS.items():
            if module_env != env:
                continue
            for tool in MODULE_TOOLS.get(module, []):
                if tool not in shared:
                    shared.append(tool)
        if shared:
            return shared
    return []


#: Order is precedence: a module name wins over a tool name, which wins over an env base name,
#: which wins over a 2.0.x alias. Declared here, after all four dicts above exist.
ENV_REGISTRIES: Tuple[Mapping[str, str], ...] = (
    MODULE_ENVS,
    OPTIONAL_ENVS,
    OPTIONAL_ENVS_BY_MODULE,
    LEGACY_ENV_ALIASES,
)


def resolve_env_name(name: str) -> Optional[str]:
    """The conda env *name* refers to, or None if it refers to nothing.

    Five things a user or a module can name, all of which appear in the docs, in a module's own
    code, or in an older release's output:

    * a module - ``binning`` -> ``metawrap2-bins``
    * a tool with its own env - ``semibin2`` -> ``metawrap2-learned``
    * an env's own base name - ``bins`` -> ``metawrap2-bins``
    * an opt-in env - ``learned`` -> ``metawrap2-learned``
    * a 2.0.x per-tool env name - ``bin_refinement-checkm2`` -> ``metawrap2-checkm2``

    Every command that turns a name into an environment goes through this. Three used to do it
    themselves and disagreed in three different ways: `install-env` raised KeyError, `doctor`
    invented ``metawrap2-<whatever was typed>`` and reported it missing, and `check` declared
    "Preflight OK" having matched nothing. Returning None is what lets all of them refuse a name
    rather than invent an answer for it.
    """
    for registry in ENV_REGISTRIES:
        found = registry.get(name)
        if found:
            return found
    prefixed = "metawrap2-%s" % name
    if prefixed in known_envs():
        return prefixed
    return None


def known_envs() -> Set[str]:
    """Every conda env name this package can create or use."""
    return {env for registry in ENV_REGISTRIES for env in registry.values()}


def installable_names() -> List[str]:
    """Every name the env commands accept, for an error message that can be acted on."""
    names = {key for registry in ENV_REGISTRIES for key in registry}
    names |= {env.removeprefix("metawrap2-") for env in known_envs()}
    return sorted(names)


def explain_unknown_name(name: str) -> List[str]:
    """Lines to print for a name that resolves to nothing, including a near-miss suggestion."""
    options = installable_names()
    lines = ["'%s' is not a module or an environment." % name]
    close = difflib.get_close_matches(name, options, n=1)
    if close:
        lines.append("Did you mean '%s'?" % close[0])
    lines.append("Known names: %s" % ", ".join(options))
    return lines


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


def env_prefix(env: str) -> Optional[str]:
    """Where a conda env lives, or None if it does not exist."""
    for directory in env_dirs():
        candidate = os.path.join(directory, env)
        if os.path.isdir(candidate):
            return candidate
    return None


def tool_in_env(tool: str, env: str) -> Optional[bool]:
    """Whether *tool* is executable inside *env*, or None if the env does not exist.

    Looks for ``<prefix>/bin/<tool>`` rather than shelling out to ``mamba run -n env which`` - the
    answer is the same and a `check` over every module would otherwise pay a conda activation per
    tool.
    """
    prefix = env_prefix(env)
    if prefix is None:
        return None
    return os.access(os.path.join(prefix, "bin", tool), os.X_OK)


def check_tools(modules: List[str]) -> List[Tuple[str, str, bool]]:
    """Return (tool, module, reachable) for the tools used by the given modules.

    *reachable* means reachable **the way the run will reach it**: inside the module's conda env
    when envs are in use, and on PATH when they are not. That distinction is the whole point of
    this function rather than a refinement of it.

    Checking PATH unconditionally is what it used to do, and it produced the worst possible answer
    for a user with a tool installed system-wide and the env missing: `check` reported the tool
    present, and the run then died inside ``mamba run`` with "The given prefix does not exist". The
    check passed on a question nobody had asked.

    An env that does not exist yields False here rather than falling back to PATH, because a run
    with conda envs enabled will not fall back either.
    """
    settings = load_settings()
    results = []
    seen = set()
    for module in modules:
        env = MODULE_ENVS.get(module) if settings.use_conda_envs else None
        for tool in MODULE_TOOLS.get(module, []):
            if (tool, module) in seen:
                continue
            seen.add((tool, module))
            if env:
                in_env = tool_in_env(tool, env)
                reachable = bool(in_env)
            else:
                reachable = shutil.which(tool) is not None
            results.append((tool, module, reachable))
    return results
