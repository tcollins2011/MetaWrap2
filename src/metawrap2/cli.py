"""MetaWrap2 command-line entry point: ``metawrap2 <module> [options]``.

Routes ``metawrap2 <module> ...`` to that module's ``main()``. Each module parses its own
arguments and can also be run directly (``python -m metawrap2.modules.binning ...``).
"""

from __future__ import annotations

import importlib
import os
import sys
from typing import List, Optional

from . import __version__
from .command import EnvManagerMissing
from .config import ConfigError

# module name -> (import path, one-line description shown in help)
MODULES = {
    "read_qc": ("metawrap2.modules.read_qc", "Read QC (trimming and host/contaminant removal)"),
    "assembly": ("metawrap2.modules.assembly", "Metagenomic assembly (metaSPAdes or MEGAHIT)"),
    "kraken2": ("metawrap2.modules.kraken2", "Taxonomy profiling of reads/contigs with KRAKEN2"),
    "blobology": ("metawrap2.modules.blobology", "GC-vs-abundance blobplots of contigs and bins"),
    "binning": (
        "metawrap2.modules.binning",
        "Initial binning (metaBAT2, MaxBin2, CONCOCT, SemiBin2, VAMB, COMEBin)",
    ),
    "coverage": (
        "metawrap2.modules.coverage",
        "One cross-sample depth table for a cohort, split per assembly",
    ),
    "bin_refinement": (
        "metawrap2.modules.bin_refinement",
        "Consolidate bin sets into a refined set",
    ),
    "reassemble_bins": (
        "metawrap2.modules.reassemble_bins",
        "Reassemble bins from recruited reads",
    ),
    "quant_bins": ("metawrap2.modules.quant_bins", "Estimate bin abundance across samples"),
    "dereplicate": (
        "metawrap2.modules.dereplicate",
        "Collapse several samples' bins into one non-redundant catalogue",
    ),
    "classify_bins": ("metawrap2.modules.classify_bins", "Assign taxonomy to genomic bins"),
    "annotate": (
        "metawrap2.modules.annotate",
        "Call genes once over bins or an assembly, and annotate them",
    ),
    "diff_abundance": (
        "metawrap2.modules.diff_abundance",
        "Differential abundance between sample groups, by ANCOM-BC",
    ),
    "catalogue": (
        "metawrap2.modules.catalogue",
        "Collapse genes from many samples into one non-redundant catalogue",
    ),
    "quantify_function": (
        "metawrap2.modules.quantify_function",
        "Sample-level gene and pathway abundance from per-genome annotation",
    ),
    "viruses": (
        "metawrap2.modules.viruses",
        "Recover viral contigs, estimate completeness, group into vOTUs",
    ),
    "amr": (
        "metawrap2.modules.amr",
        "Resistance and virulence genes with their genomic context",
    ),
    "appraise": (
        "metawrap2.modules.appraise",
        "What share of each sample's community the bins represent, per cell",
    ),
    "pangenome": (
        "metawrap2.modules.pangenome",
        "A pangenome per species in a cohort, from dereplication's clusters",
    ),
    "phylogeny": (
        "metawrap2.modules.phylogeny",
        "A tree from a marker or core-gene alignment, with support values",
    ),
    "strains": (
        "metawrap2.modules.strains",
        "Within-population variation and strain sharing between samples",
    ),
    "nonpareil": (
        "metawrap2.modules.nonpareil",
        "Diversity coverage per sample, and whether samples are comparable",
    ),
    "pathways": (
        "metawrap2.modules.pathways",
        "KEGG module completeness and the functional hierarchy, per genome and per sample",
    ),
    # Kept because published scripts use it. The Prokka/Bakta wrapper, replaced by `annotate`.
    "annotate_bins": ("metawrap2.modules.annotate_bins", "Deprecated: use `annotate`"),
}

# Non-module utility subcommands.
UTILITIES = {
    # Listed first because it is the starting point after cloning.
    "quickstart": (
        "metawrap2.commands.quickstart",
        "Build the conda environments and verify them, in one command (start here)",
    ),
    "run": ("metawrap2.commands.run", "Run the whole pipeline for a study from a sample sheet"),
    "cohort": (
        "metawrap2.commands.cohort",
        "Run a study per sample and dereplicate into one genome catalogue",
    ),
    "check": ("metawrap2.commands.check", "Verify tools/databases (and conda envs) for a module"),
    "install-env": (
        "metawrap2.commands.install_env",
        "Create (and verify) the conda env(s) for module(s)",
    ),
    "install-db": ("metawrap2.commands.install_db", "Download, index, and configure the databases"),
    "history": (
        "metawrap2.commands.history_cmd",
        "Show every MetaWrap2 run on this machine, and when",
    ),
    "status": (
        "metawrap2.commands.status_cmd",
        "Show what has run in a study, and what --resume would do next",
    ),
    "report": (
        "metawrap2.commands.report_cmd",
        "One PDF report, and one per-genome table, for a finished run",
    ),
    "score": (
        "metawrap2.commands.score_cmd",
        "Score any directories of genomes with one scorer, and compare them",
    ),
    "calibrate": (
        "metawrap2.commands.calibrate",
        "Check the built-in disk/memory estimates against real runs",
    ),
    "benchmark": (
        "metawrap2.commands.benchmark",
        "Measure tools and settings against known-truth datasets",
    ),
    # `doctor` is the documented name (with `--fix`); `test` is an older alias. `check` is for
    # scripts (pass/fail exit code); `doctor` gives a readable report.
    "doctor": (
        "metawrap2.commands.doctor",
        "Check every module's env, tools, databases and Python extras; `--fix` repairs envs",
    ),
    # The older spelling of `doctor`. Both work; neither is deprecated.
    "test": ("metawrap2.commands.doctor", "Older name for `doctor`"),
    "config": (
        "metawrap2.commands.config_cmd",
        "Scaffold (init) or inspect (show) the metawrap2.toml config",
    ),
    "completion": ("metawrap2.commands.completion", "Print a shell completion script (bash|zsh)"),
}


#: Module and flag names removed since earlier versions, with what to use instead. Reported
#: in place of argparse's "unrecognized argument".
RETIRED = {
    "kraken": "the KRAKEN1 module was dropped. Use `kraken2`, which takes the same arguments.",
    "--metabat1": "the original metaBAT binner was dropped. Use --metabat2.",
    "--metabat": "the original metaBAT binner was dropped. Use --metabat2.",
    "phylosift": "the PhyloSift module was never enabled and no longer exists. "
    "For bin taxonomy use `classify_bins`.",
}


def _retirement_notice(argv: List[str]) -> Optional[str]:
    """The explanation for the first retired name in *argv*, if there is one."""
    for token in argv:
        if token in RETIRED:
            return "%s: %s" % (token, RETIRED[token])
    return None


#: Timeout in seconds for `git rev-parse`. This runs on `-h` and `-v`, so a stalled network
#: filesystem must not hang them.
GIT_REVISION_TIMEOUT = 2.0


def git_revision() -> Optional[str]:
    """The short commit this package is running from, or None if that cannot be established.

    A checkout is usually past the last tag, so `-v` output in bug reports needs the commit.

    Returns None instead of raising when there is no repository (wheel install), no git binary,
    or git fails, so `-v` always prints. Output that does not look like an abbreviated sha is
    rejected.
    """
    import subprocess

    here = os.path.dirname(os.path.abspath(__file__))
    try:
        proc = subprocess.run(
            ["git", "-C", here, "rev-parse", "--short", "HEAD"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=GIT_REVISION_TIMEOUT,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode != 0:
        return None
    rev = (proc.stdout or "").strip()
    if not (7 <= len(rev) <= 40) or any(c not in "0123456789abcdef" for c in rev):
        return None
    return rev


def version_string() -> str:
    """What `metawrap2 -v` prints: ``MetaWrap2 2.1.2``, plus the commit in a checkout.

    Formatted to be pasted into a bug report as is.
    """
    rev = git_revision()
    if rev:
        return "MetaWrap2 %s (git %s)" % (__version__, rev)
    return "MetaWrap2 %s" % __version__


def help_message() -> str:
    lines = [
        "",
        version_string(),
        "Usage: metawrap2 [module] [options]",
        "",
        "  Modules:",
    ]
    # Name column sized to the longest name.
    width = max(len(name) for name in list(MODULES) + list(UTILITIES)) + 2
    for name, (_path, desc) in MODULES.items():
        lines.append("\t%-*s%s" % (width, name, desc))
    lines.append("")
    lines.append("  Utilities:")
    for name, (_path, desc) in UTILITIES.items():
        lines.append("\t%-*s%s" % (width, name, desc))
    lines += [
        "",
        "  Global options (before the module name):",
        "\t%-20s print/record the commands without running them" % "--dry-run",
        "\t%-20s overwrite an existing MetaWrap2 output directory" % "--force",
        "\t%-20s reuse an existing output directory, skipping finished steps" % "--resume",
        "\t%-20s don't collapse repeated tool messages in the captured logs" % "--verbose-logs",
        "\t%-20s seconds of tool silence before a progress heartbeat (0 off)"
        % "--heartbeat-seconds N",
        "\t%-20s start even if the disk-space estimate says there is no room"
        % "--skip-space-check",
        "\t%-20s skip input/intermediate file validation (use if it flags a valid file)"
        % "--skip-validation",
        "",
        # Wider column than above: these flag names do not fit in 20 characters.
        "  Resource budget (override [limits] in metawrap2.toml for this run):",
        "\t%-24s thread ceiling for the whole run" % "--max-threads N",
        "\t%-24s or, as a share of this machine's cores (default 0.50)" % "--max-thread-fraction F",
        "\t%-24s memory ceiling for the whole run" % "--max-memory-gb N",
        "\t%-24s or, as a share of this machine's RAM (default 0.80)" % "--max-memory-fraction F",
        "\t%-24s stop when a filesystem is this full (default 0.70)" % "--max-disk-fraction F",
        "\t%-24s how many GPUs the run may use (0 for none)" % "--max-gpus N",
        "\t%-24s refuse to start a GPU step on a card another process holds" % "--gpu-strict",
        # Listed because it is a documented [limits] key like the others.
        "\t%-24s seconds between disk checks during a long step (default 60)"
        % "--disk-check-interval S",
        "",
        "\t%-20s show this help message" % "--help | -h",
        "\t%-20s show MetaWrap2 version" % "--version | -v",
        "",
    ]
    return "\n".join(lines)


#: Global flags that take a value, mapped to the ``[limits]`` field they override. Both
#: ``--flag value`` and ``--flag=value`` are accepted.
_LIMIT_FLAGS = {
    "--max-threads": "max_threads",
    "--max-thread-fraction": "max_thread_fraction",
    "--max-memory-gb": "max_memory_gb",
    "--max-memory-fraction": "max_memory_fraction",
    "--max-disk-fraction": "max_disk_fraction",
    "--max-gpus": "max_gpus",
    "--disk-check-interval": "disk_check_interval",
}


#: Value-carrying global flags that are not ``[limits]`` keys. Separate from
#: :data:`_LIMIT_FLAGS` because these are ``[settings]`` and go through their own setter rather
#: than the budget's override table and validator.
_SETTING_FLAGS = {"--heartbeat-seconds"}


def _apply_setting_flag(flag: str, value: str) -> None:
    """Apply one ``[settings]`` global flag, rejecting a value the setting cannot take."""
    from . import command, limits

    try:
        seconds = float(value)
    except ValueError:
        raise limits.LimitsError("%s needs a number of seconds, got %r" % (flag, value)) from None
    if seconds < 0:
        raise limits.LimitsError("%s cannot be negative (got %s); 0 turns it off" % (flag, value))
    command.set_heartbeat_seconds(seconds)


def _apply_global_flags(argv: List[str]) -> List[str]:
    """Strip global flags (--dry-run, --force, the --max-* budget) and apply them.

    Parsed by hand because they come before the module name, and each module parses the rest.
    """
    from . import command, limits

    kept: List[str] = []
    pending: Optional[str] = None
    for tok in argv:
        if pending is not None:
            if pending in _SETTING_FLAGS:
                _apply_setting_flag(pending, tok)
            else:
                limits.set_override(_LIMIT_FLAGS[pending], tok)
            pending = None
        elif tok == "--dry-run":
            command.set_dry_run(True)
        elif tok == "--force":
            command.set_force(True)
        elif tok == "--resume":
            command.set_resume(True)
        elif tok == "--verbose-logs":
            command.set_verbose_logs(True)
        elif tok == "--skip-space-check":
            command.set_skip_space_check(True)
        elif tok == "--skip-validation":
            command.set_skip_validation(True)
        elif tok == "--gpu-strict":
            # A boolean [limits] key, so it takes no value and is not in _LIMIT_FLAGS. The string
            # goes through the same validator as the config file.
            limits.set_override("gpu_strict", "true")
        elif tok in _LIMIT_FLAGS or tok in _SETTING_FLAGS:
            pending = tok
        elif "=" in tok and tok.split("=", 1)[0] in _LIMIT_FLAGS:
            flag, value = tok.split("=", 1)
            limits.set_override(_LIMIT_FLAGS[flag], value)
        elif "=" in tok and tok.split("=", 1)[0] in _SETTING_FLAGS:
            flag, value = tok.split("=", 1)
            _apply_setting_flag(flag, value)
        else:
            kept.append(tok)
    if pending is not None:
        raise limits.LimitsError("%s needs a value" % pending)
    return kept


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(help_message())
        return 0
    if argv[0] in ("-v", "--version"):
        print(version_string())
        return 0

    from .limits import LimitsError

    try:
        argv = _apply_global_flags(argv)
    except LimitsError as exc:
        # Report a bad budget flag before any module starts.
        print("\n%s\n" % exc, file=sys.stderr)
        return 2
    if not argv:
        print(help_message())
        return 0
    name, rest = argv[0], argv[1:]
    entry = MODULES.get(name) or UTILITIES.get(name)
    if entry is None:
        notice = _retirement_notice([name])
        if notice:
            print("\n%s\n" % notice, file=sys.stderr)
            return 2
        print("Unknown module '%s'.\n%s" % (name, help_message()))
        return 1

    # Explain retired flags here, before the module's parser rejects them as unrecognised.
    notice = _retirement_notice(rest)
    if notice:
        print("\n%s\n" % notice, file=sys.stderr)
        return 2

    mod = importlib.import_module(entry[0])
    try:
        return int(mod.main(rest) or 0)
    except (ConfigError, EnvManagerMissing) as exc:
        # Config errors and a missing mamba on PATH are user errors: print the message, not a
        # traceback.
        print("\n%s\n" % exc, file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        # Also covers interrupt.Stopped, so a scheduler's SIGTERM exits 143 without a traceback.
        from . import interrupt

        print("\nStopped before finishing. Re-run with --resume to continue.", file=sys.stderr)
        return interrupt.exit_code()


if __name__ == "__main__":
    sys.exit(main())
