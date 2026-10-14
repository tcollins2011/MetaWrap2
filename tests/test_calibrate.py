"""`metawrap2 calibrate` - the estimates and the budget, read against what runs really did.

The disk multipliers and the `[limits]` budget are both numbers someone chose. These tests cover
the reading of them back: that a step's grant reaches the manifest, that calibrate compares the
grant with the measured peak rather than with itself, and that it says so plainly enough to act on.
"""

from __future__ import annotations

import os

import pytest

from metawrap2.commands import calibrate
from metawrap2.manifest import COMPLETED, MANIFEST_NAME, Manifest

_STATE = os.path.join(".metawrap2", "run_steps")

BIG = 200 * 1024 * 1024  # comfortably over MIN_INPUT_BYTES


def study(tmp_path, steps):
    """Build a study directory whose manifest holds *steps* (name -> entry)."""
    path = tmp_path / _STATE / MANIFEST_NAME
    manifest = Manifest(str(path))
    for name, entry in steps.items():
        manifest.data["steps"][name] = entry
    manifest.save()
    return str(tmp_path)


def step(
    *,
    inputs=BIG,
    outputs=BIG,
    peak_gb=None,
    granted_gb=None,
    granted_threads=None,
    user_seconds=None,
    wall=100.0,
    status=COMPLETED,
):
    entry = {
        "status": status,
        "inputs": [{"kind": "file", "path": "r1.fastq", "size": inputs}],
        "outputs": [{"kind": "file", "path": "out.fa", "size": outputs}],
        "duration_seconds": wall,
    }
    usage = {"wall_seconds": wall}
    if peak_gb is not None:
        usage["peak_rss_bytes"] = int(peak_gb * 1024**3)
    if user_seconds is not None:
        usage["user_seconds"] = user_seconds
        usage["system_seconds"] = 0.0
    entry["usage"] = usage
    grant = {}
    if granted_gb is not None:
        grant["memory_gb"] = granted_gb
    if granted_threads is not None:
        grant["threads"] = granted_threads
    if grant:
        entry["grant"] = grant
    return entry


# --- the comparison and the report ----------------------------------------------------------


def test_the_report_prints_the_budget_table(tmp_path, capsys):
    """calibrate compares the grant with the measured peak, of the worst run, and never with zero.

    Several runs of a module: the budget has to cover the worst of them, not the average. A study
    from before the budget existed must produce no comparison at all, rather than a division by zero
    or a fabricated 0 GB grant.
    """
    found = calibrate.observations(
        [
            study(
                tmp_path,
                {
                    "assembly:a": step(peak_gb=10, granted_gb=64, granted_threads=16),
                    "assembly:b": step(peak_gb=58, granted_gb=64, granted_threads=16),
                    "binning": step(peak_gb=40),
                },
            )
        ]
    )
    by_module = {row["module"]: row for row in calibrate.summarise(found)}
    assembly = by_module["assembly"]
    assert assembly["runs_with_grant"] == 2
    assert assembly["max_granted_memory_gb"] == 64
    assert assembly["max_memory_used_fraction"] == pytest.approx(58 / 64)
    assert "close to" in calibrate._budget_verdict(assembly)

    ungranted = by_module["binning"]
    assert ungranted["runs_with_grant"] == 0
    assert ungranted["max_memory_used_fraction"] is None
    assert calibrate._budget_verdict(ungranted) == "no grant recorded"

    rc = calibrate.main([str(tmp_path)])
    out = capsys.readouterr().out
    assert rc == 0
    assert "granted" in out and "64 GB" in out
    assert "max_memory_fraction" in out  # tells the reader what to change
