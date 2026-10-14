"""MetaWrap2 command-line entry point: ``metawrap2 <module> [options]``.

Thin dispatcher. Each module owns its own argument parsing and is fully usable on its own
(``python -m metawrap2.modules.binning ...``); this just routes to the right module's
``main()`` so ``metawrap2 binning ...`` works too. Kept deliberately simple.
"""

from __future__ import annotations

import importlib
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
    # Kept as a name because it is in published scripts. It is the Prokka/Bakta wrapper, which
    # `annotate` replaces; see that module's docstring for what it could not answer.
    "annotate_bins": ("metawrap2.modules.annotate_bins", "Deprecated: use `annotate`"),
}

# Non-module utility subcommands.
UTILITIES = {
    # First in the list on purpose: it is the answer to "I just cloned this, now what?"
    "quickstart": (
        "metawrap2.commands.quickstart",
        "Install everything and verify it, in one command (start here)",
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
    "calibrate": (
        "metawrap2.commands.calibrate",
        "Check the built-in disk/memory estimates against real runs",
    ),
    "benchmark": (
        "metawrap2.commands.benchmark",
        "Measure tools and settings against known-truth datasets",
    ),
    "test": (
        "metawrap2.commands.doctor",
        "Test module conda envs + software, run unit tests, print a status map",
    ),
    # Same command under the name people reach for when something is broken - and the name its
    # own repair flag uses (`--fix`). Both spellings work; neither is deprecated.
    "doctor": ("metawrap2.commands.doctor", "Alias for `test`; `doctor --fix` repairs what it can"),
    "config": (
        "metawrap2.commands.config_cmd",
        "Scaffold (init) or inspect (show) the metawrap2.toml config",
    ),
    "completion": ("metawrap2.commands.completion", "Print a shell completion script (bash|zsh)"),
}


#: Module and flag names that existed in earlier versions and no longer do. An "unrecognized
#: argument" from argparse is technically correct and tells the user nothing; each of these says
#: what happened and what to use instead.
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


def help_message() -> str:
    lines = [
        "",
        "MetaWrap2 v=%s" % __version__,
        "Usage: metawrap2 [module] [options]",
        "",
        "  Modules:",
    ]
    # The name column fits the longest name: a fixed width of 16 ran `quantify_function` straight
    # into its own description.
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
        "\t%-20s start even if the disk-space estimate says there is no room"
        % "--skip-space-check",
        "\t%-20s skip input/intermediate file validation (use if it flags a valid file)"
        % "--skip-validation",
        "",
        # Wider column than the block above: these flag names do not fit in 20 characters, and a
        # description that starts mid-word is worse than a column that does not line up with it.
        "  Resource budget (override [limits] in metawrap2.toml for this run):",
        "\t%-24s thread ceiling for the whole run" % "--max-threads N",
        "\t%-24s or, as a share of this machine's cores (default 0.50)" % "--max-thread-fraction F",
        "\t%-24s memory ceiling for the whole run" % "--max-memory-gb N",
        "\t%-24s or, as a share of this machine's RAM (default 0.80)" % "--max-memory-fraction F",
        "\t%-24s stop when a filesystem is this full (default 0.70)" % "--max-disk-fraction F",
        "\t%-24s how many GPUs the run may use (0 for none)" % "--max-gpus N",
        "",
        "\t%-20s show this help message" % "--help | -h",
        "\t%-20s show MetaWrap2 version" % "--version | -v",
        "",
    ]
    return "\n".join(lines)


#: Global flags that carry a value, mapped to the ``[limits]`` field they override. Kept as a
#: table rather than more branches below because every one of them behaves identically, and
#: because the same table is what makes ``--flag=value`` work without repeating the parsing.
_LIMIT_FLAGS = {
    "--max-threads": "max_threads",
    "--max-thread-fraction": "max_thread_fraction",
    "--max-memory-gb": "max_memory_gb",
    "--max-memory-fraction": "max_memory_fraction",
    "--max-disk-fraction": "max_disk_fraction",
    "--max-gpus": "max_gpus",
    "--disk-check-interval": "disk_check_interval",
}


def _apply_global_flags(argv: List[str]) -> List[str]:
    """Strip global flags (--dry-run, --force, the --max-* budget) and apply them.

    Done by hand rather than with argparse because these come *before* the module name, and each
    module owns its own parser for everything after it.
    """
    from . import command, limits

    kept: List[str] = []
    pending: Optional[str] = None
    for tok in argv:
        if pending is not None:
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
        elif tok in _LIMIT_FLAGS:
            pending = tok
        elif "=" in tok and tok.split("=", 1)[0] in _LIMIT_FLAGS:
            flag, value = tok.split("=", 1)
            limits.set_override(_LIMIT_FLAGS[flag], value)
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
        print("MetaWrap2 v=%s" % __version__)
        return 0

    from .limits import LimitsError

    try:
        argv = _apply_global_flags(argv)
    except LimitsError as exc:
        # A bad budget flag is a user error, and one worth catching before any module starts.
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

    # A retired flag reaches the module's parser as an unrecognised argument, which says nothing
    # about why it is gone. Intercept it here so every module gets the explanation for free.
    notice = _retirement_notice(rest)
    if notice:
        print("\n%s\n" % notice, file=sys.stderr)
        return 2

    mod = importlib.import_module(entry[0])
    try:
        return int(mod.main(rest) or 0)
    except (ConfigError, EnvManagerMissing) as exc:
        # A mistyped config key is a user error, not a bug, so it gets a message rather than a
        # traceback - and it is caught here so every module and utility benefits.
        #
        # EnvManagerMissing is the same kind of thing: mamba not being on PATH is something the user
        # fixes, and the exception exists purely to say so in two sentences with both ways out. It
        # was still escaping as a traceback, which buried its own message under a stack and made a
        # PATH problem look like a crash - the exact confusion it was written to prevent.
        print("\n%s\n" % exc, file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        # Covers interrupt.Stopped, so a scheduler's SIGTERM exits 143 rather than tracebacking.
        from . import interrupt

        print("\nStopped before finishing. Re-run with --resume to continue.", file=sys.stderr)
        return interrupt.exit_code()


if __name__ == "__main__":
    sys.exit(main())
