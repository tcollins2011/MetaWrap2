"""`metawrap2 config init|show|migrate` - scaffold, inspect, and convert the MetaWrap2 config.

- ``config init``     writes a starter ``~/.metawrap2/config.toml`` you can edit.
- ``config show``     prints the resolved settings (which file is in effect, database paths,
  thread default, whether conda envs are used).
- ``config migrate``  converts an older shell-style configuration - a file of
  ``KEY=/path`` lines meant to be ``source``d - into the TOML MetaWrap2 reads, keeping the paths
  and saying what happened to anything that no longer applies.
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, List, Optional, Tuple

from ..config import RETIRED_DB_KEYS, find_config, load_settings

#: The starter config. This is the single source of truth: ``metawrap2.toml.example`` at the
#: repository root is generated from it, and a test asserts the two have not drifted apart -
#: they were previously two hand-maintained copies, and the example had fallen behind.
DEFAULT_CONFIG = """\
# MetaWrap2 configuration. Copy to ~/.metawrap2/config.toml (or pass with --config).
# `metawrap2 install-db` writes the [databases] paths here for you.
# Anything you don't set falls back to the in-code defaults documented in each module.

[settings]
# Default thread count when a module isn't given -t. Modules use this instead of silently
# running single-threaded, so it is worth setting to something sensible for your machine.
threads = 8
# Set false to run tools directly off your PATH instead of per-module conda envs.
use_conda_envs = true

[databases]
# Point these at your own databases, or let `metawrap2 install-db` fill them in.
# kraken2 module:
KRAKEN2_DB = "/path/to/my/kraken2_database"
# read_qc host removal (the directory holding <host>.bitmask and <host>.srprism*):
BMTAGGER_DB = "/path/to/my/bmtagger_database"
# blobology and classify_bins:
BLASTDB = "/path/to/my/NCBI_nt"
# The database's name inside BLASTDB. Only needed if it is not called "nt" (e.g. "core_nt").
# BLASTDB_NAME = "nt"
TAXDUMP = "/path/to/my/NCBI_tax"
# CheckM1 reference data, used by bin_refinement / reassemble_bins / binning --run-checkm.
# CheckM stores this path itself; `metawrap2 install-db checkm` runs `checkm data setRoot`.
# CHECKM_DB = "/path/to/my/checkm_data"

# Only needed for the opt-in modern tools:
# annotate_bins --bakta
# BAKTA_DB = "/path/to/my/bakta_db"
# classify_bins --gtdbtk
# GTDBTK_DATA_PATH = "/path/to/my/gtdb_release"

# Optional: override any surfaced module constant without editing the source.
# [modules.binning]
# MIN_CONTIG_LEN = 2000
"""


#: Keys an older shell-style config could set, mapped to the key MetaWrap2 reads. Most are
#: unchanged; the point of the table is that it is explicit about which ones are not.
SHELL_CONFIG_KEYS = {
    "KRAKEN2_DB": "KRAKEN2_DB",
    "BMTAGGER_DB": "BMTAGGER_DB",
    "BLASTDB": "BLASTDB",
    "TAXDUMP": "TAXDUMP",
    "CHECKM_DB": "CHECKM_DB",
    "BAKTA_DB": "BAKTA_DB",
    "GTDBTK_DATA_PATH": "GTDBTK_DATA_PATH",
}

#: Variables such a file may contain that are not database paths at all - they located the
#: installation's own scripts, which MetaWrap2 resolves from the installed package instead.
SHELL_CONFIG_IGNORED = {"SOFT", "PIPES", "mw_path", "bin_path"}


def parse_shell_config(text: str) -> Tuple[Dict[str, str], List[str], List[str]]:
    """Read ``KEY=value`` lines. Returns ``(translated, retired, unrecognised)``.

    Deliberately a line scanner rather than anything that executes the file: the input is a shell
    script, and running a file to read its settings is how a config becomes an attack surface.
    Only plain assignments are understood, which is all such a file ever held; anything with
    command substitution in it is reported as unrecognised rather than guessed at.
    """
    translated: Dict[str, str] = {}
    retired: List[str] = []
    unrecognised: List[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        name, _, value = line.partition("=")
        name = name.strip()
        value = value.split("#", 1)[0].strip().strip("\"'")
        if not name or name in SHELL_CONFIG_IGNORED:
            continue
        if "$" in value or "`" in value:
            unrecognised.append("%s (value is computed, not a literal path)" % name)
            continue
        if name in SHELL_CONFIG_KEYS:
            translated[SHELL_CONFIG_KEYS[name]] = value
        elif name in RETIRED_DB_KEYS:
            retired.append(name)
        else:
            unrecognised.append(name)
    return translated, retired, unrecognised


def render_config(databases: Dict[str, str], threads: Optional[int] = None) -> str:
    """A minimal config holding exactly *databases*, with comments saying where it came from."""
    lines = [
        "# MetaWrap2 configuration, converted from a shell-style config.",
        "# Check the paths below, then move this to ~/.metawrap2/config.toml (or pass --config).",
        "",
        "[settings]",
        "# Default thread count when a module is not given -t.",
        "threads = %d" % (threads or 8),
        "use_conda_envs = true",
        "",
        "[databases]",
    ]
    for key in sorted(databases):
        lines.append('%s = "%s"' % (key, databases[key]))
    lines.append("")
    return "\n".join(lines)


def _migrate(args) -> int:
    if not os.path.isfile(args.source):
        print("No such file: %s" % args.source)
        return 1
    with open(args.source, errors="replace") as fh:
        translated, retired, unrecognised = parse_shell_config(fh.read())

    if not translated and not retired:
        print(
            "\nFound nothing to convert in %s.\n"
            "Expected lines like KRAKEN2_DB=/path/to/db. If this is not that kind of file,\n"
            "`metawrap2 config init` writes a fresh config instead.\n" % args.source
        )
        return 1

    text = render_config(translated, threads=args.threads)
    if args.output:
        if os.path.exists(args.output) and not args.force:
            print("%s already exists. Pass --force to overwrite it." % args.output)
            return 1
        os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True)
        with open(args.output, "w") as fh:
            fh.write(text)
        print("\nWrote %s" % args.output)
    else:
        print(text)

    print("\n  Converted %d database path(s):" % len(translated))
    for key in sorted(translated):
        exists = "" if os.path.exists(translated[key]) else "   <- this path does not exist here"
        print("    %-18s %s%s" % (key, translated[key], exists))
    for key in retired:
        print("\n  Dropped %s: %s" % (key, RETIRED_DB_KEYS[key]))
    if unrecognised:
        print("\n  Left out (not database paths MetaWrap2 reads): %s" % ", ".join(unrecognised))
    if not args.output:
        print(
            "\n  Save it with:  metawrap2 config migrate %s -o ~/.metawrap2/config.toml\n"
            % args.source
        )
    else:
        print("\n  Check it with: metawrap2 config show --config %s\n" % args.output)
    return 0


def _init(args) -> int:
    path = args.path or os.path.expanduser("~/.metawrap2/config.toml")
    if os.path.exists(path) and not args.force:
        print("%s already exists (use --force to overwrite)." % path)
        return 1
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(DEFAULT_CONFIG)
    print("Wrote starter config to %s\nEdit it, then run: metawrap2 check <module>" % path)
    return 0


def _show(args) -> int:
    path = find_config(args.config)
    settings = load_settings(args.config)
    print("\nMetaWrap2 configuration")
    print("  config file:    %s" % (path or "<none found; using built-in defaults>"))
    print("  threads:        %d" % settings.threads)
    print("  use_conda_envs: %s" % settings.use_conda_envs)
    print("  databases:")
    if settings.databases:
        for key, val in settings.databases.items():
            exists = "ok" if val and os.path.exists(val) else "missing"
            print("    %-12s %s  [%s]" % (key, val, exists))
    else:
        print("    <none set>")
    print()
    return 0


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(prog="metawrap2 config")
    sub = ap.add_subparsers(dest="action", required=True)

    p_init = sub.add_parser("init", help="write a starter config file")
    p_init.add_argument("--path", help="where to write (default ~/.metawrap2/config.toml)")
    p_init.add_argument("--force", action="store_true", help="overwrite an existing file")

    p_show = sub.add_parser("show", help="print the resolved configuration")
    p_show.add_argument("--config", help="path to a metawrap2.toml")

    p_mig = sub.add_parser(
        "migrate", help="convert an older shell-style config (KEY=/path lines) to TOML"
    )
    p_mig.add_argument("source", help="the file to convert")
    p_mig.add_argument("-o", "--output", help="write here instead of printing to stdout")
    p_mig.add_argument("--force", action="store_true", help="overwrite the output file")
    p_mig.add_argument(
        "-t", "--threads", type=int, help="thread default to write into [settings] (default 8)"
    )

    args = ap.parse_args(argv)
    return {"init": _init, "show": _show, "migrate": _migrate}[args.action](args)
