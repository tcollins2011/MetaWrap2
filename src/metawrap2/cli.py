"""MetaWrap2 command-line entry point: ``metawrap2 <module> [options]``.

Thin dispatcher. Each module owns its own argument parsing and is fully usable on its own
(``python -m metawrap2.modules.binning ...``); this just routes to the right module's
``main()`` so ``metawrap2 binning ...`` works too. Kept deliberately simple.
"""

from __future__ import annotations

import importlib
import sys
from typing import List

from . import __version__

# module name -> (import path, one-line description shown in help)
MODULES = {
    "read_qc": ("metawrap2.modules.read_qc", "Read QC (trimming and host/contaminant removal)"),
    "assembly": ("metawrap2.modules.assembly", "Metagenomic assembly (metaSPAdes or MEGAHIT)"),
    "kraken2": ("metawrap2.modules.kraken2", "Taxonomy profiling of reads/contigs with KRAKEN2"),
    "blobology": ("metawrap2.modules.blobology", "GC-vs-abundance blobplots of contigs and bins"),
    "binning": ("metawrap2.modules.binning", "Initial binning (metaBAT2, MaxBin2, CONCOCT)"),
    "bin_refinement": ("metawrap2.modules.bin_refinement", "Consolidate bin sets into a refined set"),
    "reassemble_bins": ("metawrap2.modules.reassemble_bins", "Reassemble bins from recruited reads"),
    "quant_bins": ("metawrap2.modules.quant_bins", "Estimate bin abundance across samples"),
    "classify_bins": ("metawrap2.modules.classify_bins", "Assign taxonomy to genomic bins"),
    "annotate_bins": ("metawrap2.modules.annotate_bins", "Functional annotation of draft genomes"),
}

# Non-module utility subcommands.
UTILITIES = {
    "check": ("metawrap2.commands.check", "Verify tools/databases (and conda envs) for a module"),
    "install-env": ("metawrap2.commands.install_env", "Create (and verify) the conda env(s) for module(s)"),
    "test": ("metawrap2.commands.doctor", "Test module conda envs + software, run unit tests, print a status map"),
    "config": ("metawrap2.commands.config_cmd", "Scaffold (init) or inspect (show) the metawrap2.toml config"),
    "completion": ("metawrap2.commands.completion", "Print a shell completion script (bash|zsh)"),
}


def help_message() -> str:
    lines = ["", "MetaWrap2 v=%s" % __version__, "Usage: metawrap2 [module] [options]", "", "  Modules:"]
    for name, (_path, desc) in MODULES.items():
        lines.append("\t%-16s%s" % (name, desc))
    lines.append("")
    lines.append("  Utilities:")
    for name, (_path, desc) in UTILITIES.items():
        lines.append("\t%-16s%s" % (name, desc))
    lines += [
        "",
        "  Global options (before the module name):",
        "\t--dry-run        print/record the commands without running them",
        "\t--force          overwrite an existing MetaWrap2 output directory",
        "\t--resume         reuse an existing output directory, skipping finished steps",
        "",
        "\t--help | -h      show this help message",
        "\t--version | -v   show MetaWrap2 version",
        "",
    ]
    return "\n".join(lines)


def _apply_global_flags(argv: List[str]) -> List[str]:
    """Strip global flags (--dry-run, --force) from argv and apply them to the runner."""
    from . import command
    kept = []
    for tok in argv:
        if tok == "--dry-run":
            command.set_dry_run(True)
        elif tok == "--force":
            command.set_force(True)
        elif tok == "--resume":
            command.set_resume(True)
        else:
            kept.append(tok)
    return kept


def main(argv: List[str] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(help_message())
        return 0
    if argv[0] in ("-v", "--version"):
        print("MetaWrap2 v=%s" % __version__)
        return 0

    argv = _apply_global_flags(argv)
    if not argv:
        print(help_message())
        return 0
    name, rest = argv[0], argv[1:]
    entry = MODULES.get(name) or UTILITIES.get(name)
    if entry is None:
        print("Unknown module '%s'.\n%s" % (name, help_message()))
        return 1

    mod = importlib.import_module(entry[0])
    return int(mod.main(rest) or 0)


if __name__ == "__main__":
    sys.exit(main())
