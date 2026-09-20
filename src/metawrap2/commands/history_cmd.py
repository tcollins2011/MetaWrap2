"""`metawrap2 history` - what has been run on this machine, and when.

Every module invocation appends to ``~/.metawrap2/history.jsonl`` (see
:mod:`metawrap2.history`), outside any run's output directory, and stays there for as long as
MetaWrap2 is installed. This command reads it back.

    metawrap2 history                      # the 20 most recent runs
    metawrap2 history -n 100               # more of them
    metawrap2 history --module binning     # just one module
    metawrap2 history --failed             # only runs that did not finish
    metawrap2 history --output STUDY/      # what produced a particular directory
    metawrap2 history --commands           # the exact command lines, to re-run or paste
    metawrap2 history --json               # the raw records, for scripting
    metawrap2 history --stats              # totals per module, and time spent
    metawrap2 history --diff A B            # what differed between two runs
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
from typing import Any, Dict, List, Optional

from .. import history as _history

_STATUS_MARK = {
    _history.COMPLETED: "ok",
    _history.FAILED: "FAILED",
    _history.INTERRUPTED: "interrupted",
    _history.STARTED: "started?",
}


def _duration(entry: Dict[str, Any]) -> str:
    seconds = entry.get("duration_seconds")
    if seconds is None:
        return "-"
    seconds = float(seconds)
    if seconds < 90:
        return "%.0fs" % seconds
    if seconds < 5400:
        return "%.0fm" % (seconds / 60)
    return "%.1fh" % (seconds / 3600)


def _print_table(rows: List[Dict[str, Any]]) -> None:
    print("\n  %-19s %-16s %-12s %8s  %s" % ("when", "module", "status", "duration", "output"))
    print("  " + "-" * 96)
    for entry in rows:
        print(
            "  %-19s %-16s %-12s %8s  %s"
            % (
                entry.get("time", "?"),
                entry.get("module", "?"),
                _STATUS_MARK.get(entry.get("status", ""), entry.get("status", "?")),
                _duration(entry),
                entry.get("output", "-"),
            )
        )
    print()


def _print_stats(rows: List[Dict[str, Any]]) -> None:
    per_module: Dict[str, Dict[str, float]] = {}
    for entry in rows:
        stats = per_module.setdefault(
            entry.get("module", "?"), {"runs": 0, "completed": 0, "failed": 0, "seconds": 0.0}
        )
        stats["runs"] += 1
        if entry.get("status") == _history.COMPLETED:
            stats["completed"] += 1
        elif entry.get("status") in (_history.FAILED, _history.INTERRUPTED):
            stats["failed"] += 1
        stats["seconds"] += float(entry.get("duration_seconds") or 0)

    print("\n  %-16s %6s %10s %8s %12s" % ("module", "runs", "completed", "failed", "total time"))
    print("  " + "-" * 60)
    for module in sorted(per_module, key=lambda m: -per_module[m]["runs"]):
        stats = per_module[module]
        print(
            "  %-16s %6d %10d %8d %12s"
            % (
                module,
                stats["runs"],
                stats["completed"],
                stats["failed"],
                _duration({"duration_seconds": stats["seconds"]}),
            )
        )
    total = sum(s["seconds"] for s in per_module.values())
    print("  " + "-" * 60)
    print(
        "  %-16s %6d %10s %8s %12s"
        % ("all", len(rows), "", "", _duration({"duration_seconds": total}))
    )
    print()


#: Fields worth comparing between two runs, in the order a reader wants them.
_DIFF_FIELDS = (
    ("module", "module"),
    ("version", "MetaWrap2 version"),
    ("conda_env", "conda env"),
    ("python", "Python"),
    ("host", "host"),
    ("user", "user"),
    ("status", "outcome"),
    ("duration_seconds", "duration (s)"),
)


def _find_run(rows: List[Dict[str, Any]], token: str) -> Optional[Dict[str, Any]]:
    """Locate a run by run_id (or a unique prefix of one), or by exact timestamp."""
    exact = [r for r in rows if r.get("run_id") == token or r.get("time") == token]
    if exact:
        return exact[-1]
    prefixed = [r for r in rows if str(r.get("run_id", "")).startswith(token)]
    return prefixed[-1] if len(prefixed) == 1 else None


def _diff_values(a: Any, b: Any) -> bool:
    return a != b


def _print_diff(rows: List[Dict[str, Any]], left: str, right: str) -> int:
    """Compare two runs field by field. Returns an exit code."""
    one = _find_run(rows, left)
    two = _find_run(rows, right)
    for token, found in ((left, one), (right, two)):
        if found is None:
            print(
                "\nNo run matching %r (or the prefix matches more than one).\n"
                "List run ids with: metawrap2 history --json\n" % token
            )
            return 1
    assert one is not None and two is not None  # for the type checker; guarded above

    print("\n  %-20s %-34s %s" % ("", one.get("run_id", "?")[:32], two.get("run_id", "?")[:32]))
    print("  " + "-" * 92)
    differences = 0
    for key, label in _DIFF_FIELDS:
        a, b = one.get(key), two.get(key)
        if a is None and b is None:
            continue
        mark = " "
        if _diff_values(a, b):
            mark = "*"
            differences += 1
        print("  %s %-18s %-34s %s" % (mark, label, a, b))

    # Databases and the command line are the two things that most often explain a difference in
    # results, and both are structured, so they are compared per key rather than as blobs.
    db_one = one.get("databases") or {}
    db_two = two.get("databases") or {}
    for key in sorted(set(db_one) | set(db_two)):
        a, b = db_one.get(key, "-"), db_two.get(key, "-")
        if _diff_values(a, b):
            differences += 1
            print("  * %-18s %-34s %s" % ("db " + key, a, b))

    cmd_one = " ".join(shlex.quote(x) for x in one.get("command", []))
    cmd_two = " ".join(shlex.quote(x) for x in two.get("command", []))
    if cmd_one != cmd_two:
        differences += 1
        print("\n  * command differs:")
        print("      %s" % cmd_one)
        print("      %s" % cmd_two)

    print("  " + "-" * 92)
    if differences:
        print("  %d difference(s), marked *." % differences)
    else:
        print("  No differences in anything the history records.")
    # The history does not hold per-tool versions; those live with the data that used them.
    print(
        "  Tool versions and the full package list are per-run, in each output directory's\n"
        "  run_config.json and run_environment.json.\n"
    )
    return 0


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(
        prog="metawrap2 history", description="Show the history of MetaWrap2 runs on this machine."
    )
    ap.add_argument(
        "-n",
        "--limit",
        type=int,
        default=20,
        help="how many of the most recent runs to show (0 for all; default 20)",
    )
    ap.add_argument("--module", help="only runs of this module")
    ap.add_argument("--output", help="only runs that wrote to this output directory")
    ap.add_argument("--since", metavar="YYYY-MM-DD", help="only runs at or after this date")
    ap.add_argument(
        "--failed", action="store_true", help="only runs that failed or were interrupted"
    )
    ap.add_argument(
        "--commands",
        action="store_true",
        help="print the command line of each run instead of a table",
    )
    ap.add_argument("--json", action="store_true", help="print the raw records as JSON")
    ap.add_argument("--stats", action="store_true", help="summarise by module")
    ap.add_argument(
        "--diff",
        nargs=2,
        metavar=("RUN", "RUN"),
        help="compare two runs field by field (run ids, or unique prefixes of them)",
    )
    ap.add_argument("--path", help="read a specific history file")
    args = ap.parse_args(argv)

    target = _history.history_path(args.path)
    rows = _history.runs(args.path, module=args.module, since=args.since, output=args.output)
    if args.failed:
        rows = [r for r in rows if r.get("status") in (_history.FAILED, _history.INTERRUPTED)]

    if args.diff:
        # Deliberately before the "no rows" message and the limit: a diff names its runs, so
        # filters and -n are irrelevant to it.
        return _print_diff(_history.runs(args.path), args.diff[0], args.diff[1])

    if not rows:
        if not os.path.isfile(target):
            print(
                "\nNo history yet: %s does not exist.\nIt is created the first time you run "
                "a module.\n" % target
            )
        else:
            print("\nNo runs in %s match those filters.\n" % target)
        return 0

    rows = rows[-args.limit :] if args.limit else rows

    if args.json:
        print(json.dumps(rows, indent=2, sort_keys=True))
        return 0
    if args.stats:
        _print_stats(rows)
        return 0
    if args.commands:
        print()
        for entry in rows:
            print(
                "# %s  %s  (%s)"
                % (entry.get("time", "?"), entry.get("module", "?"), entry.get("status", "?"))
            )
            print(" ".join(shlex.quote(a) for a in entry.get("command", [])))
            print()
        return 0

    _print_table(rows)
    print("  %d run(s) shown from %s\n" % (len(rows), target))
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
