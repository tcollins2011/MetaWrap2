"""`metawrap2 calibrate <study_dir>...` - check the built-in estimates against real runs.

Two sets of numbers in MetaWrap2 are asserted rather than measured:

* :data:`metawrap2.scratch.SPACE_MULTIPLIERS` - peak disk need as a multiple of input size;
* the memory budget a user passes with ``-m``, divided between workers.

Every run already records, per step, the size of each input, the size of each output, and (since
per-step measurement was added) peak RSS and CPU time. So the observed values can be compared
against the constants instead of the constants being trusted.

    metawrap2 calibrate study1/ study2/          # table per module
    metawrap2 calibrate study*/ --suggest        # constants that fit the observations
    metawrap2 calibrate study/ --json

**What this can and cannot measure.** Output-to-input ratio is measurable and is a *lower bound*
on the disk multiplier - scratch is deleted when a step finishes, so no manifest can say how large
it got at its peak. A multiplier below the observed output ratio is therefore definitely wrong; one
above it may still be too low. Peak memory, by contrast, is measured exactly. The report says which
is which rather than presenting both as equally solid, because a confidently wrong disk estimate is
how a run fills a filesystem four hours in.
"""

from __future__ import annotations

import argparse
import json
import os
from typing import Any, Dict, List

from ..manifest import COMPLETED, Manifest
from ..scratch import SPACE_MULTIPLIERS, human
from ..usage import human as human_bytes
from .status_cmd import manifest_path

#: A step whose inputs total less than this tells us nothing about scaling: fixed overheads
#: dominate, and dividing by a tiny denominator produces meaningless ratios.
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
            # Only the entry count was recorded, not the bytes. Skipping is honest; guessing a
            # per-file size would put invented numbers into a calibration report.
            continue
    return total


def observations(study_dirs: List[str]) -> Dict[str, List[Dict[str, Any]]]:
    """Per-module observations gathered from each study's manifest."""
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
            found.setdefault(_module_of(step["name"]), []).append(
                {
                    "study": study,
                    "step": step["name"],
                    "input_bytes": inputs,
                    "output_bytes": _total(step.get("outputs") or []),
                    "peak_rss_bytes": usage.get("peak_rss_bytes"),
                    "wall_seconds": step.get("duration_seconds"),
                }
            )
    return found


def _ratio(entry: Dict[str, Any]) -> float:
    return entry["output_bytes"] / entry["input_bytes"] if entry["input_bytes"] else 0.0


def summarise(found: Dict[str, List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
    """One row per module: observed ratios and memory, against the configured multiplier."""
    rows = []
    for module in sorted(found):
        entries = found[module]
        ratios = [_ratio(e) for e in entries]
        memories = [e["peak_rss_bytes"] for e in entries if e["peak_rss_bytes"]]
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
            }
        )
    return rows


def _verdict(row: Dict[str, Any]) -> str:
    """Whether the configured multiplier is contradicted by what was observed."""
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
        "  out/in is a LOWER BOUND on the disk multiplier: scratch is deleted when a step\n"
        "  finishes, so no manifest records how large it got at its peak. A multiplier below\n"
        "  the observed ratio is definitely wrong; one above it may still be too low.\n"
        "  Peak memory is measured exactly.\n"
    )


def _print_suggestions(rows: List[Dict[str, Any]]) -> None:
    """Print a SPACE_MULTIPLIERS block that is at least consistent with what was seen."""
    print("\n  A SPACE_MULTIPLIERS block no observation here contradicts:\n")
    print("  SPACE_MULTIPLIERS = {")
    for row in rows:
        configured = row["configured_multiplier"] or 2.0
        observed = row["max_output_ratio"] or 0.0
        # Keep the configured value unless it is below what was actually produced; then raise it
        # with headroom, because the real peak includes scratch this cannot see.
        suggested = configured if configured >= observed else round(observed * 2, 1)
        note = "" if suggested == configured else "  # raised: observed %.2fx output" % observed
        print('      "%s": %.1f,%s' % (row["module"], suggested, note))
    print("  }")
    print(
        "\n  Only the modules observed in these studies appear. Paste over the block in\n"
        "  src/metawrap2/scratch.py, keeping any module not listed here.\n"
    )


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(
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
            "A study contributes only completed steps whose inputs total at least %s - a step\n"
            "smaller than that says nothing about how the estimates scale.\n"
            % (", ".join(args.studies), human(MIN_INPUT_BYTES))
        )
        return 1

    rows = summarise(found)
    if args.json:
        print(json.dumps({"summary": rows, "observations": found}, indent=2, sort_keys=True))
        return 0

    _print_table(rows)
    if args.per_step:
        print("  %-28s %12s %12s %11s" % ("step", "inputs", "outputs", "peak memory"))
        print("  " + "-" * 70)
        for module in sorted(found):
            for entry in found[module]:
                print(
                    "  %-28s %12s %12s %11s"
                    % (
                        entry["step"],
                        human(entry["input_bytes"]),
                        human(entry["output_bytes"]),
                        human_bytes(entry["peak_rss_bytes"]) if entry["peak_rss_bytes"] else "-",
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
