"""`metawrap2 status <study_dir>` - what has run in this study, and what a resume would do.

The manifest already records every step's command, inputs, outputs, timing and outcome. Until
now nothing read it back, so the information existed but was only visible to ``--resume``. This
command prints it:

    metawrap2 status study/            # per-step table, and what --resume would skip or redo
    metawrap2 status study/ --audit    # check every recorded output is still present, non-empty
    metawrap2 status study/ --step binning --verbose    # one step's inputs, outputs, command
    metawrap2 status study/ --json     # the manifest, for scripting

The most useful column is the last one. A step is shown as ``would skip`` only under the same
verification ``--resume`` applies - completed, outputs intact, inputs unchanged - so this answers
"what will happen if I re-run this?" before you re-run it.
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
from typing import Any, Dict, List, Optional

from .. import usage as _usage
from ..manifest import COMPLETED, MANIFEST_NAME, Manifest

#: Must match metawrap2.commands.run._STATE_DIR. Imported lazily in manifest_path() rather than
#: at module scope, because importing the driver to print a table is a lot of machinery.
_STATE_DIR = os.path.join(".metawrap2", "run_steps")

_MARK = {
    COMPLETED: "ok",
    "failed": "FAILED",
    "interrupted": "interrupted",
    "running": "running?",
}


def manifest_path(study_dir: str) -> str:
    """The manifest inside *study_dir* - or *study_dir* itself if it already is one."""
    if os.path.isfile(study_dir):
        return study_dir
    return os.path.join(study_dir, _STATE_DIR, MANIFEST_NAME)


def _duration(seconds: Optional[float]) -> str:
    if seconds is None:
        return "-"
    seconds = float(seconds)
    if seconds < 90:
        return "%.0fs" % seconds
    if seconds < 5400:
        return "%.0fm" % (seconds / 60)
    return "%.1fh" % (seconds / 3600)


def _resume_note(manifest: Manifest, name: str, strict: bool) -> str:
    skip, reason = manifest.can_skip(name, whole_file=strict)
    return ("would skip: %s" if skip else "would re-run: %s") % reason


def _print_table(manifest: Manifest, steps: List[Dict[str, Any]], strict: bool) -> None:
    print("\n  %-28s %-13s %9s  %s" % ("step", "status", "duration", "on --resume"))
    print("  " + "-" * 110)
    for step in steps:
        print(
            "  %-28s %-13s %9s  %s"
            % (
                step["name"],
                _MARK.get(step.get("status", ""), step.get("status", "?")),
                _duration(step.get("duration_seconds")),
                _resume_note(manifest, step["name"], strict),
            )
        )
    print()


def _print_step(manifest: Manifest, steps: List[Dict[str, Any]]) -> None:
    for step in steps:
        print("\n  %s" % step["name"])
        print("    status     %s" % step.get("status", "?"))
        print("    started    %s" % step.get("started", "?"))
        print("    duration   %s" % _duration(step.get("duration_seconds")))
        if step.get("log"):
            print("    log        %s" % step["log"])
        cost = step.get("usage") or {}
        if cost:
            print("    cost       %s" % _usage.summarise(cost))
        command = step.get("command") or []
        if command:
            print("    command    %s" % " ".join(shlex.quote(a) for a in command))
        for label in ("inputs", "outputs"):
            entries = step.get(label) or []
            print("    %-10s %d" % (label, len(entries)))
            for entry in entries:
                detail = entry.get("kind", "?")
                if entry.get("kind") == "file":
                    detail = "%s, %d bytes" % (detail, entry.get("size", 0))
                elif entry.get("kind") == "directory":
                    count = entry.get("entries", 0)
                    detail = "%s, %d %s" % (detail, count, "entry" if count == 1 else "entries")
                print("      - %s (%s)" % (entry.get("path", "?"), detail))
    print()


def _print_audit(manifest: Manifest, root: str) -> int:
    from .run import Layout

    problems = manifest.audit(root, ignore=Layout.DRIVER_OWNED)
    bad = problems["missing"] + problems["empty"] + problems["no_outputs"]
    print()
    for path in problems["missing"]:
        print("  MISSING    %s" % path)
    for path in problems["empty"]:
        print("  EMPTY      %s" % path)
    for name in problems["no_outputs"]:
        print("  UNVERIFIED %s declared no outputs, so --resume can never verify it" % name)
    for path in problems["undeclared"]:
        print("  note       %s is not declared as any step's output" % path)
    if not bad:
        print("  Every recorded output is present and non-empty.")
    print()
    return 1 if bad else 0


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(
        prog="metawrap2 status",
        description="Show what has run in a study directory, and what --resume would do next.",
    )
    ap.add_argument("study", help="a study output directory (or a manifest.json directly)")
    ap.add_argument("--step", help="show only this step, in full")
    ap.add_argument(
        "-v", "--verbose", action="store_true", help="show each step's inputs, outputs and command"
    )
    ap.add_argument(
        "--audit",
        action="store_true",
        help="verify every recorded output still exists and is non-empty (exit 1 if not)",
    )
    ap.add_argument(
        "--strict-fingerprints",
        action="store_true",
        help="hash whole files when deciding whether an input changed, as --resume would",
    )
    ap.add_argument("--json", action="store_true", help="print the raw manifest")
    args = ap.parse_args(argv)

    path = manifest_path(args.study)
    if not os.path.isfile(path):
        print(
            "\nNo MetaWrap2 run manifest at %s.\nEither %s is not a `metawrap2 run` study "
            "directory, or it was made before manifests existed.\n" % (path, args.study)
        )
        return 1

    manifest = Manifest(path)
    steps = manifest.summary()
    if args.step:
        steps = [s for s in steps if s["name"] == args.step]
        if not steps:
            print(
                "\nNo step named %r in %s.\nKnown steps: %s\n"
                % (args.step, path, ", ".join(s["name"] for s in manifest.summary()) or "none")
            )
            return 1

    if args.json:
        print(json.dumps(steps if args.step else manifest.data, indent=2, sort_keys=True))
        return 0
    if args.audit:
        return _print_audit(manifest, os.path.abspath(args.study))
    if not steps:
        print("\nThe manifest at %s records no steps yet.\n" % path)
        return 0
    if args.verbose or args.step:
        _print_step(manifest, steps)
        return 0

    _print_table(manifest, steps, args.strict_fingerprints)
    done = sum(1 for s in steps if s.get("status") == COMPLETED)
    print("  %d of %d step(s) completed. Manifest: %s\n" % (done, len(steps), path))
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
