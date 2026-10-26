"""`metawrap2 calibrate <study_dir>...` - compare the built-in estimates with recorded runs.

Two sets of numbers in MetaWrap2 are configured, not measured:

* :data:`metawrap2.scratch.SPACE_MULTIPLIERS` - peak disk need as a multiple of input size;
* the memory and thread budget from ``[limits]``, divided between concurrent steps.

Every run records, per step, input and output sizes, peak RSS, CPU time, and the share of the
budget the step was granted (see :mod:`metawrap2.limits`). This command compares those
observations with the constants and with each step's grant.

    metawrap2 calibrate study1/ study2/          # table per module
    metawrap2 calibrate study*/ --suggest        # constants that fit the observations
    metawrap2 calibrate study/ --json

The output-to-input ratio is a lower bound on the disk multiplier, because scratch is deleted when
a step finishes and its peak size is not recorded. A multiplier below the observed ratio is wrong;
one above it may still be too low. Peak memory is measured directly, and the report labels which
figures are bounds.

The budget table compares two measured values (grant and peak). It cannot show whether a tool
would use more memory if granted more.
"""

from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Optional

from ..cliformat import make_parser
from ..manifest import COMPLETED, Manifest
from ..scratch import SPACE_MULTIPLIERS, human
from ..usage import human as human_bytes
from .status_cmd import manifest_path

#: Steps with smaller total input are skipped: fixed overheads dominate and the ratios are noise.
MIN_INPUT_BYTES = 10 * 1024 * 1024  # 10 MB


def _module_of(step_name: str) -> str:
    """Steps are named ``module`` or ``module:sample``."""
    return step_name.split(":", 1)[0]


def _total(entries: List[Dict[str, Any]]) -> int:
    total = 0
    for entry in entries:
        if entry.get("kind") == "file":
            total += int(entry.get("size") or 0)
        elif entry.get("kind") == "directory":
            # Only the entry count is recorded, not the bytes, so directories are not counted.
            continue
    return total


def _cores_busy(usage: Dict[str, Any], wall: Optional[float]) -> Optional[float]:
    """Return the mean number of busy cores over the step's wall time, or None if not measured.

    Compared with the thread grant, this shows whether a step used its threads.
    """
    user, system = usage.get("user_seconds"), usage.get("system_seconds")
    if user is None or system is None or not wall:
        return None
    return round((user + system) / wall, 2)


def observations(study_dirs: List[str]) -> Dict[str, List[Dict[str, Any]]]:
    """Return per-module observations gathered from each study's manifest."""
    found: Dict[str, List[Dict[str, Any]]] = {}
    for study in study_dirs:
        path = manifest_path(study)
        if not os.path.isfile(path):
            continue
        for step in Manifest(path).summary():
            if step.get("status") != COMPLETED:
                continue
            inputs = _total(step.get("inputs") or [])
            if inputs < MIN_INPUT_BYTES:
                continue
            usage = step.get("usage") or {}
            grant = step.get("grant") or {}
            wall = usage.get("wall_seconds") or step.get("duration_seconds")
            cores_busy = _cores_busy(usage, wall)
            found.setdefault(_module_of(step["name"]), []).append(
                {
                    "study": study,
                    "step": step["name"],
                    "input_bytes": inputs,
                    "output_bytes": _total(step.get("outputs") or []),
                    "peak_rss_bytes": usage.get("peak_rss_bytes"),
                    "wall_seconds": step.get("duration_seconds"),
                    # The step's grant, to compare the measurement against.
                    "granted_memory_gb": grant.get("memory_gb"),
                    "granted_threads": grant.get("threads"),
                    "cores_busy": cores_busy,
                }
            )
    return found


def _ratio(entry: Dict[str, Any]) -> float:
    return entry["output_bytes"] / entry["input_bytes"] if entry["input_bytes"] else 0.0


def summarise(found: Dict[str, List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
    """Return one row per module: observed ratios and memory, and the configured multiplier."""
    rows = []
    for module in sorted(found):
        entries = found[module]
        ratios = [_ratio(e) for e in entries]
        memories = [e["peak_rss_bytes"] for e in entries if e["peak_rss_bytes"]]
        # Only steps with both a grant and a measured peak are compared. Older runs have no grant.
        against_grant = [
            e for e in entries if e.get("granted_memory_gb") and e.get("peak_rss_bytes")
        ]
        used = [e["peak_rss_bytes"] / (e["granted_memory_gb"] * 1024**3) for e in against_grant]
        threaded = [
            e["cores_busy"] / e["granted_threads"]
            for e in entries
            if e.get("cores_busy") and e.get("granted_threads")
        ]
        rows.append(
            {
                "module": module,
                "runs": len(entries),
                "configured_multiplier": SPACE_MULTIPLIERS.get(module),
                "max_output_ratio": max(ratios) if ratios else None,
                "mean_output_ratio": sum(ratios) / len(ratios) if ratios else None,
                "max_peak_rss_bytes": max(memories) if memories else None,
                "max_rss_per_input_gb": (
                    max(m / (e["input_bytes"] / 1024**3) for m, e in zip(memories, entries))
                    if memories
                    else None
                ),
                "runs_with_grant": len(against_grant),
                "max_granted_memory_gb": (
                    max(e["granted_memory_gb"] for e in against_grant) if against_grant else None
                ),
                #: Peak RSS as a fraction of the memory granted. Above 1 means the tool did not
                #: honour -m; well below 1 means memory was not the constraint.
                "max_memory_used_fraction": max(used) if used else None,
                #: Busy cores as a fraction of the threads granted. A low value means more threads
                #: will not speed the step up.
                "max_thread_used_fraction": max(threaded) if threaded else None,
            }
        )
    return rows


def _verdict(row: Dict[str, Any]) -> str:
    """Return whether the observations contradict the configured multiplier."""
    configured = row["configured_multiplier"]
    observed = row["max_output_ratio"]
    if configured is None:
        return "no multiplier configured"
    if observed is None:
        return "no ratio observed"
    if configured < observed:
        return "TOO LOW: outputs alone exceeded it"
    if configured > observed * 8:
        return "possibly generous"
    return "consistent"


def _print_table(rows: List[Dict[str, Any]]) -> None:
    print(
        "\n  %-16s %5s %7s %9s %9s %11s  %s"
        % ("module", "runs", "config", "out/in max", "mean", "peak memory", "verdict")
    )
    print("  " + "-" * 104)
    for row in rows:
        print(
            "  %-16s %5d %7s %9s %9s %11s  %s"
            % (
                row["module"],
                row["runs"],
                ("%.1fx" % row["configured_multiplier"]) if row["configured_multiplier"] else "-",
                ("%.2fx" % row["max_output_ratio"]) if row["max_output_ratio"] else "-",
                ("%.2fx" % row["mean_output_ratio"]) if row["mean_output_ratio"] else "-",
                human_bytes(row["max_peak_rss_bytes"]) if row["max_peak_rss_bytes"] else "-",
                _verdict(row),
            )
        )
    print("  " + "-" * 104)
    print(
        "  out/in is a lower bound on the disk multiplier: scratch is deleted when a step\n"
        "  finishes, so its peak size is not recorded. A multiplier below the observed ratio\n"
        "  is too low; one above it may still be too low. Peak memory is measured.\n"
    )


def _budget_verdict(row: Dict[str, Any]) -> str:
    """Return how the steps' peak memory compares with their grant."""
    used = row["max_memory_used_fraction"]
    if used is None:
        return "no grant recorded"
    if used > 1.0:
        return "EXCEEDED the grant by %.0f%%" % ((used - 1) * 100)
    if used > 0.85:
        return "close to the grant"
    if used < 0.25:
        return "used under a quarter of it"
    return "comfortable"


def _print_budget_table(rows: List[Dict[str, Any]]) -> None:
    """Print granted versus used memory and threads per module."""
    measured = [r for r in rows if r["runs_with_grant"]]
    if not measured:
        print(
            "  No step recorded its granted budget, so peak memory cannot be compared.\n"
            "  Runs made by `metawrap2 run` record the budget.\n"
        )
        return
    print(
        "  %-16s %5s %9s %9s %7s %7s  %s"
        % ("module", "runs", "granted", "peak", "used", "threads", "verdict")
    )
    print("  " + "-" * 104)
    for row in measured:
        print(
            "  %-16s %5d %9s %9s %7s %7s  %s"
            % (
                row["module"],
                row["runs_with_grant"],
                "%d GB" % row["max_granted_memory_gb"],
                human_bytes(row["max_peak_rss_bytes"]) if row["max_peak_rss_bytes"] else "-",
                "%.0f%%" % (row["max_memory_used_fraction"] * 100),
                (
                    "%.0f%%" % (row["max_thread_used_fraction"] * 100)
                    if row["max_thread_used_fraction"]
                    else "-"
                ),
                _budget_verdict(row),
            )
        )
    print("  " + "-" * 104)
    print(
        "  'granted' is the largest -m given to any of these steps, 'peak' the largest peak RSS,\n"
        "  'used' peak as a fraction of granted, and 'threads' the fraction of the thread grant\n"
        "  that was busy. For a module well under its grant, lower [limits] max_memory_fraction\n"
        "  (or raise -j to run more steps at once). A module over its grant is not honouring\n"
        "  -m; give it a smaller grant or run fewer steps at once.\n"
    )


def _print_suggestions(rows: List[Dict[str, Any]]) -> None:
    """Print a SPACE_MULTIPLIERS block that no observation contradicts."""
    print("\n  SPACE_MULTIPLIERS block consistent with these observations:\n")
    print("  SPACE_MULTIPLIERS = {")
    for row in rows:
        configured = row["configured_multiplier"] or 2.0
        observed = row["max_output_ratio"] or 0.0
        # Keep the configured value unless it is below the observed ratio; then double the
        # observed ratio, since the real peak also includes scratch.
        suggested = configured if configured >= observed else round(observed * 2, 1)
        note = "" if suggested == configured else "  # raised: observed %.2fx output" % observed
        print('      "%s": %.1f,%s' % (row["module"], suggested, note))
    print("  }")
    print(
        "\n  Only the modules observed in these studies appear. Paste over the block in\n"
        "  src/metawrap2/scratch.py, keeping any module not listed here.\n"
    )


def main(argv: List[str]) -> int:
    ap = make_parser(
        prog="metawrap2 calibrate",
        description="Compare MetaWrap2's built-in disk and memory estimates against real runs.",
    )
    ap.add_argument("studies", nargs="+", help="study output directories to read manifests from")
    ap.add_argument(
        "--suggest",
        action="store_true",
        help="print a SPACE_MULTIPLIERS block consistent with the observations",
    )
    ap.add_argument("--json", action="store_true", help="print the raw observations")
    ap.add_argument(
        "--per-step", action="store_true", help="list every observed step, not just the summary"
    )
    args = ap.parse_args(argv)

    found = observations(args.studies)
    if not found:
        print(
            "\nNothing to calibrate from. Looked for run manifests in: %s\n"
            "Only completed steps whose inputs total at least %s are used; smaller steps\n"
            "do not show how the estimates scale.\n"
            % (", ".join(args.studies), human(MIN_INPUT_BYTES))
        )
        return 1

    rows = summarise(found)
    if args.json:
        print(json.dumps({"summary": rows, "observations": found}, indent=2, sort_keys=True))
        return 0

    _print_table(rows)
    _print_budget_table(rows)
    if args.per_step:
        print(
            "  %-28s %12s %12s %11s %9s %7s"
            % ("step", "inputs", "outputs", "peak memory", "granted", "cores")
        )
        print("  " + "-" * 88)
        for module in sorted(found):
            for entry in found[module]:
                print(
                    "  %-28s %12s %12s %11s %9s %7s"
                    % (
                        entry["step"],
                        human(entry["input_bytes"]),
                        human(entry["output_bytes"]),
                        human_bytes(entry["peak_rss_bytes"]) if entry["peak_rss_bytes"] else "-",
                        (
                            "%d GB" % entry["granted_memory_gb"]
                            if entry["granted_memory_gb"]
                            else "-"
                        ),
                        (
                            "%.1f/%d" % (entry["cores_busy"], entry["granted_threads"])
                            if entry["cores_busy"] and entry["granted_threads"]
                            else "-"
                        ),
                    )
                )
        print()
    if args.suggest:
        _print_suggestions(rows)

    # Exit nonzero when a configured multiplier is contradicted, so this is usable in CI.
    wrong = [r["module"] for r in rows if _verdict(r).startswith("TOO LOW")]
    if wrong:
        print(
            "  %d multiplier(s) contradicted by observation: %s\n" % (len(wrong), ", ".join(wrong))
        )
        return 1
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
