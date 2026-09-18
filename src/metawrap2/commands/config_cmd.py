"""`metawrap2 config init|show` - scaffold and inspect the MetaWrap2 config.

- ``config init``  writes a starter ``~/.metawrap2/config.toml`` you can edit.
- ``config show``  prints the resolved settings (which file is in effect, database paths,
  thread default, whether conda envs are used).
"""

from __future__ import annotations

import argparse
import os
from typing import List

from ..config import find_config, load_settings

DEFAULT_CONFIG = """\
# MetaWrap2 configuration. Edit the values below; anything you don't set falls back to the
# in-code defaults documented in each module.

[settings]
# Default thread count when a module isn't given -t.
threads = 1
# Set false to run tools directly off your PATH instead of per-module conda envs.
use_conda_envs = true

[databases]
# Point these at your own downloaded databases.
KRAKEN2_DB = "/path/to/my/kraken2_database"
BMTAGGER_DB = "/path/to/my/bmtagger_database"
BLASTDB = "/path/to/my/NCBI_nt"
TAXDUMP = "/path/to/my/NCBI_tax"

# Optional: override any surfaced module constant without editing the source.
# [modules.binning]
# MIN_CONTIG_LEN = 2000
"""


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
    print("")
    return 0


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(prog="metawrap2 config")
    sub = ap.add_subparsers(dest="action", required=True)

    p_init = sub.add_parser("init", help="write a starter config file")
    p_init.add_argument("--path", help="where to write (default ~/.metawrap2/config.toml)")
    p_init.add_argument("--force", action="store_true", help="overwrite an existing file")

    p_show = sub.add_parser("show", help="print the resolved configuration")
    p_show.add_argument("--config", help="path to a metawrap2.toml")

    args = ap.parse_args(argv)
    return _init(args) if args.action == "init" else _show(args)
