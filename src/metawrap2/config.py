"""MetaWrap2 configuration and environment preflight.

Configuration lives in a small TOML file (``metawrap2.toml``): database paths, a few
settings, and (optionally) overrides for any surfaced module constant. Resolution order
for the config file: an explicit ``--config`` path, then ``$METAWRAP2_CONFIG``, then
``~/.metawrap2/config.toml``, then the packaged default. Values a user does not set fall
back to the in-code module defaults.

Also provides the preflight check, which reports missing tools, databases and conda envs before
a long run starts.
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

# One conda env per workflow. Modules that run the same tools or run in sequence share an env, so
# the aligner/samtools stack is solved once. Every grouping below was solved against conda-forge +
# bioconda; CheckM2 is split from CheckM1 and GTDB-Tk from taxator-tk because the solver found no
# common solution. Each envs/<name>.yaml describes its contents.
MODULE_ENVS: Dict[str, str] = {
    # Read-level modules: QC, assembly, re-assembly and the blob plot. They share the aligners and
    # SPAdes and usually run in sequence.
    "read_qc": "metawrap2-reads",
    "assembly": "metawrap2-reads",
    "reassemble_bins": "metawrap2-reads",
    "blobology": "metawrap2-reads",
    # Contigs to MAG catalogue: binning, cross-sample coverage, dereplication, abundance. coverage
    # shares binning's aligner and samtools; dereplicate and quant_bins add three small static
    # binaries (skani, sylph, salmon).
    "binning": "metawrap2-bins",
    "coverage": "metawrap2-bins",
    "dereplicate": "metawrap2-bins",
    "quant_bins": "metawrap2-bins",
    # Bin consolidation and the CheckM1 scorer (--checkm1): CheckM1, GUNC, DAS Tool. GUNC runs by
    # default. The default scorer, CheckM2, is in metawrap2-checkm2 (OPTIONAL_ENVS).
    "bin_refinement": "metawrap2-binqc",
    # Read- and bin-level taxonomy.
    "classify_bins": "metawrap2-taxonomy",
    "kraken2": "metawrap2-taxonomy",
    # Gene calling: Prokka and Bakta solve together, so --bakta needs no second install. `annotate`
    # calls genes in-process with Pyrodigal and needs this env only for --bakta.
    "annotate_bins": "metawrap2-annotate",
    "annotate": "metawrap2-annotate",
    # Only for Krona's ktImportText (hierarchy plots). The pathway computation runs in-process over
    # vendored KEGG data.
    "pathways": "metawrap2-taxonomy",
    # AMRFinderPlus, hAMRonization and argNorm. RGI is not in here; see envs/amr.yaml.
    "amr": "metawrap2-amr",
    # geNomad and CheckV, which always run together. See envs/viruses.yaml for their numpy pin.
    "viruses": "metawrap2-viruses",
    # SingleM has its own env: it pins a narrow dependency set and is the only user of the GTDB
    # marker metapackage.
    "appraise": "metawrap2-singlem",
    # Nonpareil reads only raw reads and runs beside read_qc. Bioconda builds it against R, so it
    # brings r-base; the estimates themselves are computed in Python (metawrap2.nonpareil).
    "nonpareil": "metawrap2-reads",
    # Pangenomes and phylogenies over a cohort's MAGs. One writes the alignment the other reads.
    "pangenome": "metawrap2-comparative",
    "phylogeny": "metawrap2-comparative",
    # inStrain and its aligners. Separate because inStrain resolves to Python 3.8.
    "strains": "metawrap2-strains",
    # MMseqs2 is already in the annotate env (a Bakta dependency).
    "catalogue": "metawrap2-annotate",
}

#: Every env a full install builds, each listed once. Iterate this, not MODULE_ENVS, when building
#: or checking envs.
ALL_MODULE_ENVS = sorted(set(MODULE_ENVS.values()))

# Envs kept separate from the module envs, with the reason for each. Most are opt-in; the ones in
# DEFAULT_OPTIONAL_ENVS are used by a default run.
OPTIONAL_ENVS: Dict[str, str] = {
    # Benchmark-only: AMBER is used by `metawrap2 benchmark score`, not by a module. It pins
    # matplotlib below 3.9.
    "amber": "metawrap2-benchmark-amber",
    # CheckM2 has no common solution with the CheckM1 + GUNC + DAS Tool set. It is the default bin
    # scorer, so it is in DEFAULT_OPTIONAL_ENVS.
    "checkm2": "metawrap2-checkm2",
    # The three pytorch-based binners, which solve together. Separate from metawrap2-bins because
    # pytorch is ~1.5 GB. The pytorch build is CUDA-capable, so `--device cuda` needs no extra
    # install.
    "semibin2": "metawrap2-learned",
    "vamb": "metawrap2-learned",
    "comebin": "metawrap2-learned",
    # TRACS, the second strain backend. It resolves an older scientific-python stack than
    # inStrain's env.
    "tracs": "metawrap2-tracs",
    # eggNOG-mapper 2.1.x imports distutils (removed in Python 3.12). The Python ceiling it needs
    # is kept out of the annotate env used by bakta, dbCAN and GECCO.
    "eggnog": "metawrap2-eggnog",
    # RGI pins older dependencies than AMRFinderPlus and hAMRonization; sharing an env would
    # downgrade the default backend.
    "rgi": "metawrap2-rgi",
    # GTDB-Tk has no common solution with taxator-tk. Its reference database is 100+ GB.
    "gtdbtk": "metawrap2-gtdbtk",
}

#: OPTIONAL_ENVS keys that a default run uses, and the modules whose default run uses each:
#: bin_refinement scores bins with CheckM2 unless --checkm1, reassemble_bins always. `quickstart` builds
#: these with those modules, and `doctor` counts a missing one as a failure. (The learned binners of
#: the default binning tier are added by quickstart from binning.DEFAULT_PRESET.)
DEFAULT_OPTIONAL_ENVS: Dict[str, Tuple[str, ...]] = {
    "checkm2": ("bin_refinement", "reassemble_bins"),
}

# Tools that live in an OPTIONAL_ENVS env, keyed like OPTIONAL_ENVS. `doctor` probes each tool
# inside the env that holds it, so these must not be listed under the module's own env.
OPTIONAL_ENV_TOOLS: Dict[str, List[str]] = {
    "amber": ["amber.py"],
    "checkm2": ["checkm2"],
    # The learned binners share one env, so every key probes all three and a partly built env
    # reports as broken. "SemiBin2" is the only entry point the 2.5.0 package installs. COMEBin's
    # entry point is the bash script run_comebin.sh; there is no "comebin" executable.
    "semibin2": ["SemiBin2", "vamb", "run_comebin.sh"],
    "vamb": ["SemiBin2", "vamb", "run_comebin.sh"],
    "comebin": ["SemiBin2", "vamb", "run_comebin.sh"],
    "gtdbtk": ["gtdbtk"],
    "rgi": ["rgi"],
    "tracs": ["tracs"],
    # install-db uses the downloader and annotate uses emapper, so both are required.
    "eggnog": ["emapper.py", "download_eggnog_data.py"],
}

# OPTIONAL_ENVS_BY_MODULE and OPTIONAL_TOOLS_BY_MODULE below re-key the two tables above by env
# base name, which is what a user types at `install-env` and `doctor`. Without them `doctor`
# reported an opt-in env as "0/0 tools OK" without probing it.
#: The module whose flag selects each opt-in env. Listed explicitly because the env names no
#: longer encode the module. An env with no module (AMBER, used only by `benchmark score`) maps to
#: None and is skipped.
OPTIONAL_ENV_OWNERS: Dict[str, Optional[str]] = {
    "amber": None,
    "checkm2": "bin_refinement",
    "semibin2": "binning",
    "vamb": "binning",
    "comebin": "binning",
    "gtdbtk": "classify_bins",
    "rgi": "amr",
    "tracs": "strains",
    "eggnog": "annotate",
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
    # flye is needed for --long-reads.
    "assembly": ["metaspades.py", "megahit", "quast", "bwa", "bwa-mem2", "minimap2", "flye"],
    # bowtie2 is not used: in the aligner benchmark it cost as much as bwa and agreed with it less
    # well than minimap2.
    "binning": ["metabat2", "run_MaxBin.pl", "concoct", "minimap2", "bwa", "bwa-mem2", "samtools"],
    # binning's aligners without the binners. fairy (--coverage fairy) is a single static binary
    # installed with them.
    "coverage": ["minimap2", "bwa", "bwa-mem2", "samtools", "fairy"],
    # GUNC runs by default unless --skip-gunc; DAS Tool is an alternative consolidator.
    "bin_refinement": ["checkm", "gunc", "DAS_Tool"],
    "reassemble_bins": ["spades.py", "bwa", "samtools", "minimap2"],
    "quant_bins": ["sylph", "salmon", "coverm"],
    "dereplicate": ["skani"],
    "classify_bins": ["blastn", "taxator", "binner", "taxknife"],
    "annotate_bins": ["prokka", "bakta"],
    # Gene calling uses Pyrodigal in-process, so only the optional paths need tools: bakta
    # (--bakta), run_dbcan (--cazy) and gecco (--bgc). GECCO needs no database, so it has no
    # DB_KEYS entry.
    "annotate": ["bakta", "run_dbcan", "gecco"],
    "nonpareil": ["nonpareil"],
    "strains": ["inStrain", "minimap2", "samtools"],
    "pangenome": ["ppanggolin"],
    "phylogeny": ["iqtree", "fasttree"],
    "appraise": ["singlem"],
    "pathways": ["ktImportText"],
    "amr": ["amrfinder"],
    # cctyper is also in this env but is not listed, because it runs only for --crispr-hosts.
    "viruses": ["genomad", "checkv"],
    "catalogue": ["mmseqs"],
    "kraken2": ["kraken2", "ktImportText"],
    # Same aligners as binning: blobology maps reads back for per-contig abundance (the plot's
    # y-axis), using the registry in metawrap2.aligners.
    "blobology": ["blastn", "minimap2", "bwa", "bwa-mem2", "samtools"],
}

# Database keys and which module(s) use them.
DB_KEYS = {
    "KRAKEN2_DB": "kraken2",
    # The databaseNmers.kmer_distrib files, one per read length. `bracken-build` writes them into
    # the KRAKEN2_DB directory and they are valid only for that database, so this key usually
    # repeats it. The pre-built k2_standard download ships 50 to 300 bp.
    "BRACKEN_DB": "kraken2 (--bracken; Bracken k-mer distributions)",
    "BMTAGGER_DB": "read_qc (host removal)",
    "BLASTDB": "blobology / classify_bins",
    "BLASTDB_NAME": "blobology / classify_bins (db name inside BLASTDB; default nt)",
    "TAXDUMP": "blobology / classify_bins",
    "CHECKM_DB": "bin_refinement / binning --run-checkm with --checkm1 (CheckM1 reference data)",
    "CHECKM2DB": "bin_refinement / reassemble_bins (CheckM2 DIAMOND reference database)",
    "GUNC_DB": "bin_refinement (GUNC reference database, unless --skip-gunc)",
    "BAKTA_DB": "annotate / annotate_bins (--bakta)",
    # The KOfam HMM profiles (profiles.tar.gz, 1.5 GB) and per-profile hit thresholds (ko_list),
    # both from `metawrap2 install-db kofam`.
    "KOFAM_DB": "annotate (--kofam), the KOfam directory: profiles/ and ko_list",
    # BRITE ko00001, projected to one row per (KO, pathway), from `metawrap2 install-db kegg`.
    "KEGG_BRITE": "pathways (the functional hierarchy and Krona plots)",
    # dbCAN's HMMs, DIAMOND databases and substrate mapping. Fetched from a pinned dbCAN release
    # because the `db_current` HTTP path is a moving snapshot and CAZyme calls depend on the
    # release.
    "DBCAN_DB": "annotate (--cazy), the dbCAN database directory",
    # Pfam-A's HMMs, for `annotate --pfam`. Pfam renumbers and merges families between releases,
    # so the installer reads the version from relnotes.txt and the output table records it.
    "PFAM_DB": "annotate (--pfam), the directory holding Pfam-A.hmm",
    # eggNOG's annotation SQLite and DIAMOND protein database, for `annotate --eggnog`. About 50 GB,
    # fetched by `install-db eggnog` from eggnog5.embl.de; eggNOG-mapper's download_eggnog_data.py
    # uses a hardcoded host that no longer responds. eggNOG renumbers orthologous groups between
    # releases, so every row records the database version.
    "EGGNOG_DB": "annotate (--eggnog), the eggNOG data directory: eggnog.db and "
    "eggnog_proteins.dmnd",
    # Point this at a dated release, not the `latest` symlink amrfinder_update writes, so reruns
    # use the same database.
    "AMRFINDER_DB": "amr (AMRFinderPlus reference database)",
    # geNomad's database includes plasmid markers, so `amr` mobility calls and `viruses` share it.
    "GENOMAD_DB": "viruses (geNomad marker database; also plasmid calls for amr)",
    "CHECKV_DB": "viruses (CheckV reference database)",
    # CARD's data, for `amr --rgi`. CARD releases every few months (4.0.2 at the time of writing);
    # name the directory by version, not `latest`.
    "CARD_DB": "amr (--rgi; CARD data directory holding card.json)",
    "GTDBTK_DATA_PATH": "classify_bins (--gtdbtk)",
    # The GTDB marker metapackage SingleM profiles against: 2.2 GB compressed, 9.7 GB unpacked.
    # Results depend on its release, so point this at the versioned directory, not a `latest` link.
    "SINGLEM_METAPACKAGE": "appraise (SingleM GTDB marker metapackage)",
}

# Databases only an opt-in code path needs. `metawrap2 check` counts an absent one as a failure
# only when that path was requested.
OPTIONAL_DB_KEYS = {
    "BAKTA_DB",
    "KOFAM_DB",
    "KEGG_BRITE",
    "DBCAN_DB",
    "EGGNOG_DB",
    "AMRFINDER_DB",
    "GENOMAD_DB",
    "CHECKV_DB",
    "GTDBTK_DATA_PATH",
    "CHECKM_DB",
    "CHECKM2DB",
    "BLASTDB_NAME",
    "GUNC_DB",
    "BRACKEN_DB",
}

# Keys whose value is a setting, not a filesystem path, so they are not checked for existence.
NON_PATH_DB_KEYS = {"BLASTDB_NAME"}

#: Database keys whose value is a file, not a directory. GUNC takes only ``--db_file``, so
#: install-db points GUNC_DB at the ``.dmnd`` itself.
FILE_DB_KEYS = {"GUNC_DB"}


class ConfigError(Exception):
    """The config file cannot be used as written."""


@dataclass
class Settings:
    databases: Dict[str, str] = field(default_factory=dict)
    #: The per-module default thread count, or None when unset, in which case the ``[limits]``
    #: budget applies. None (not 1) keeps "unset" distinct from "one thread".
    threads: Optional[int] = None
    use_conda_envs: bool = True
    #: Where modules put temporary space. Empty means inside the run's output directory. On a
    #: cluster /tmp is often small and node-local, and metaSPAdes intermediates can reach hundreds
    #: of GB. See metawrap2.scratch.
    scratch_dir: str = ""
    #: Set true to skip the pre-flight disk-space estimate.
    skip_space_check: bool = False
    #: Set true to make a failed disk-space estimate stop the run instead of warning. Useful on
    #: shared filesystems.
    strict_space_check: bool = False
    #: Seconds of silence before a step prints a progress heartbeat; 0 turns heartbeats off. None
    #: means unset, and :data:`metawrap2.logfilter.HEARTBEAT_SECONDS` applies.
    heartbeat_seconds: Optional[float] = None
    #: How long a step must run before a low CPU-utilisation warning, and the fraction of its thread
    #: grant below which it warns. 0 seconds turns the warning off. None means the default in
    #: :mod:`metawrap2.usage`.
    low_utilisation_after_seconds: Optional[float] = None
    low_utilisation_fraction: Optional[float] = None
    #: The run's resource budget: threads, memory, disk ceiling, GPUs. ``threads`` is capped to
    #: it, and a cohort driver divides it among concurrent samples. See metawrap2.limits.
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


#: The keys ``[settings]`` understands. Any other key is reported as an error.
SETTINGS_KEYS = {
    "threads": "default thread count for every module",
    "use_conda_envs": "run tools inside the metawrap2-* conda envs (true) or off PATH (false)",
    "scratch_dir": "where modules put temporary space",
    "skip_space_check": "skip the pre-flight disk-space estimate",
    "strict_space_check": "let a failed disk-space estimate stop the run, not just warn",
    "heartbeat_seconds": "seconds of silence before a step says it is still running (0 for never)",
    "low_utilisation_after_seconds": "how long a step must run before its core usage is worth "
    "warning about (0 for never)",
    "low_utilisation_fraction": "the share of its thread grant a step has to keep busy to avoid "
    "that warning",
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


#: Database keys that are recognised but no longer used, with what to use instead. They come from
#: older configs, so they produce a note and do not stop a run.
RETIRED_DB_KEYS = {
    "KRAKEN_DB": "the KRAKEN1 module was removed; use KRAKEN2_DB with the kraken2 module",
}


def check_config(data: Dict[str, Any], path: str = "") -> Tuple[List[str], List[str]]:
    """Check a parsed config. Returns ``(problems, notes)``.

    *problems* stop a run; *notes* are printed and ignored.

    * An unknown table, or an unknown ``[settings]`` or ``[limits]`` key, is a problem, reported
      with a did-you-mean suggestion. ``[limits]`` values are validated as well.
    * An unknown or retired ``[databases]`` key is a note. That table is open so that locally
      patched modules can add keys, and retired keys are harmless.

    All findings are returned together.
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
            # The message names the key; put the file name in front of it.
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

    With *strict* (the default), unknown tables or ``[settings]`` keys raise ConfigError. Pass
    ``strict=False`` to load the config anyway.
    """
    path = find_config(explicit)
    if not path:
        return Settings()
    with open(path, "rb") as fh:
        data = _toml.load(fh)
    problems, notes = check_config(data, path)
    for note in notes:
        # Notes are printed once and do not fail the load.
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
        # Only reachable with strict=False, after the problem was reported. Use the conservative
        # default budget.
        limits = LimitsConfig()
    configured_threads = data.get("settings", {}).get("threads")
    settings_table = data.get("settings", {})

    def optional_float(key: str) -> Optional[float]:
        """A numeric setting, or None when it is absent.

        0 turns these settings off, so an absent key stays None and the in-code default applies.
        """
        value = settings_table.get(key)
        return None if value is None else float(value)

    return Settings(
        databases=dict(data.get("databases", {})),
        threads=int(configured_threads) if configured_threads is not None else None,
        use_conda_envs=bool(data.get("settings", {}).get("use_conda_envs", True)),
        scratch_dir=str(data.get("settings", {}).get("scratch_dir", "") or ""),
        skip_space_check=bool(data.get("settings", {}).get("skip_space_check", False)),
        strict_space_check=bool(data.get("settings", {}).get("strict_space_check", False)),
        heartbeat_seconds=optional_float("heartbeat_seconds"),
        low_utilisation_after_seconds=optional_float("low_utilisation_after_seconds"),
        low_utilisation_fraction=optional_float("low_utilisation_fraction"),
        limits=limits,
        overrides={k: dict(v) for k, v in data.get("modules", {}).items()},
    )


def env_dirs() -> List[str]:
    """Directories conda/mamba keep environments in, most authoritative first.

    ``conda env list`` is the authoritative source, because envs can live anywhere ``envs_dirs``
    points. The default install locations are added only when no env manager is on PATH; adding
    them otherwise can pick up envs from a second conda installation that ``mamba run -n`` cannot
    enter.
    """
    dirs: List[str] = []
    #: True once an env manager has listed its environments, even if it listed none. The path
    #: guesses are used only when this stays False.
    listed = False
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
        listed = True
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
    if not listed:
        # No env manager on PATH: fall back to the default install locations.
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


#: Environment names from MetaWrap2 2.0.x, which had one env per tool. They are mapped to the
#: current envs because the 2.0.x documentation and error messages use them.
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

    The 2.0.x names are ``<module>-<tool>`` (``bin_refinement-checkm2``), so the module is the
    part before the first dash. `check` and `doctor` use it to find the module's tools and
    databases.
    """
    if name not in LEGACY_ENV_ALIASES:
        return None
    module = name.split("-", 1)[0]
    return module if module in MODULE_ENVS else None


#: ENV_REGISTRIES (below) lists every registry that maps a name onto a conda env:
#:
#: * ``MODULE_ENVS`` is keyed by module (``binning``).
#: * ``OPTIONAL_ENVS`` is keyed by tool (``semibin2``, ``vamb``, ``comebin``, ``amber``) and is
#:   what ``optional_env_for`` reads at run time.
#: * ``OPTIONAL_ENVS_BY_MODULE`` is keyed by env base name (``learned``, ``benchmark-amber``).
#: * ``LEGACY_ENV_ALIASES`` is keyed by the 2.0.x per-tool env name.
#:
#: The keys of ``OPTIONAL_ENVS`` and ``OPTIONAL_ENVS_BY_MODULE`` differ (``semibin2`` is only in
#: the first), so name resolution must read all of them.
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

    # An env base name ("bins", "reads", "binqc", "taxonomy"): return the union of tools over
    # every module that uses that env.
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

    Accepted names:

    * a module - ``binning`` -> ``metawrap2-bins``
    * a tool with its own env - ``semibin2`` -> ``metawrap2-learned``
    * an env's own base name - ``bins`` -> ``metawrap2-bins``
    * an opt-in env - ``learned`` -> ``metawrap2-learned``
    * a 2.0.x per-tool env name - ``bin_refinement-checkm2`` -> ``metawrap2-checkm2``

    Every command that maps a name to an env uses this, so all of them reject unknown names the
    same way.
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
    """True if a conda env named *env* exists in any of conda's environment directories.

    This includes other conda installations: ``conda env list`` reads the machine-wide
    ``~/.conda/environments.txt``. Such an env is usable only by prefix, so callers that run in it
    must use :func:`env_prefix` and ``-p <prefix>``; ``-n <name>`` resolves only against the
    active root.
    """
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


def active_conda_root() -> Optional[str]:
    """The conda/mamba installation that ``-n <name>`` resolves names against, or None.

    ``MAMBA_ROOT_PREFIX`` and ``CONDA_ROOT`` name it outright; otherwise it is derived from
    ``CONDA_PREFIX``, the active env: ``<root>/envs/<name>`` for a named env, the root itself for
    base. Used to flag envs that belong to a different conda installation.
    """
    for var in ("MAMBA_ROOT_PREFIX", "CONDA_ROOT"):
        value = os.environ.get(var, "")
        if value:
            return os.path.normpath(os.path.expanduser(value))
    prefix = os.environ.get("CONDA_PREFIX", "")
    if not prefix:
        return None
    prefix = os.path.normpath(prefix)
    parent = os.path.dirname(prefix)
    if os.path.basename(parent) == "envs":
        return os.path.dirname(parent)
    return prefix


def prefix_in_active_root(prefix: str) -> bool:
    """Whether *prefix* lives under the active conda root, so that ``-n`` could have found it.

    True when there is no active root to compare against.
    """
    root = active_conda_root()
    if not root:
        return True
    prefix = os.path.normpath(prefix)
    return prefix == root or prefix.startswith(root + os.sep)


def tool_in_env(tool: str, env: str) -> Optional[bool]:
    """Whether *tool* is executable inside *env*, or None if the env does not exist.

    Looks for ``<prefix>/bin/<tool>`` instead of running ``mamba run -n env which``, which would
    cost one conda activation per tool.
    """
    prefix = env_prefix(env)
    if prefix is None:
        return None
    return os.access(os.path.join(prefix, "bin", tool), os.X_OK)


def check_tools(modules: List[str]) -> List[Tuple[str, str, bool]]:
    """Return (tool, module, reachable) for the tools used by the given modules.

    *reachable* means reachable the way the run will look for it: inside the module's conda env
    when envs are in use, on PATH otherwise. A missing env yields False with no PATH fallback,
    matching what the run does.
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
