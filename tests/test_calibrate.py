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


# --- the manifest side ---------------------------------------------------------------------


def test_the_manifest_records_what_a_step_was_granted(tmp_path):
    manifest = Manifest(str(tmp_path / "manifest.json"))
    started = manifest.start("assembly", ["metawrap2", "assembly"], [], grant={"memory_gb": 64})
    manifest.finish("assembly", started, COMPLETED)
    assert Manifest(str(tmp_path / "manifest.json")).step("assembly").data["grant"] == {
        "memory_gb": 64
    }


def test_a_step_with_no_grant_records_no_grant_key(tmp_path):
    # Absent rather than zero: calibrate must be able to tell "not recorded" from "given none",
    # or every pre-budget run shows up as a step that was granted nothing and exceeded it.
    manifest = Manifest(str(tmp_path / "manifest.json"))
    manifest.start("assembly", ["metawrap2", "assembly"], [])
    assert "grant" not in manifest.step("assembly").data


def test_the_grant_survives_an_interrupted_step(tmp_path):
    # The OOM case: the step that died is the one whose grant you want to look at.
    manifest = Manifest(str(tmp_path / "manifest.json"))
    started = manifest.start("assembly", ["x"], [], grant={"memory_gb": 8, "threads": 4})
    manifest.finish("assembly", started, "interrupted")
    assert manifest.step("assembly").data["grant"]["memory_gb"] == 8


# --- reading it back -----------------------------------------------------------------------


def test_observations_carry_the_grant_through(tmp_path):
    found = calibrate.observations(
        [study(tmp_path, {"assembly": step(peak_gb=40, granted_gb=64, granted_threads=16)})]
    )
    entry = found["assembly"][0]
    assert entry["granted_memory_gb"] == 64
    assert entry["granted_threads"] == 16


def test_cores_busy_is_cpu_time_over_wall_time(tmp_path):
    found = calibrate.observations(
        [study(tmp_path, {"assembly": step(peak_gb=1, user_seconds=800.0, wall=100.0)})]
    )
    assert found["assembly"][0]["cores_busy"] == 8.0


def test_cores_busy_is_none_when_cpu_time_was_not_measured(tmp_path):
    found = calibrate.observations([study(tmp_path, {"assembly": step(peak_gb=1)})])
    assert found["assembly"][0]["cores_busy"] is None


# --- the comparison ------------------------------------------------------------------------


def test_used_fraction_is_the_peak_over_the_grant(tmp_path):
    found = calibrate.observations([study(tmp_path, {"assembly": step(peak_gb=32, granted_gb=64)})])
    row = calibrate.summarise(found)[0]
    assert row["max_memory_used_fraction"] == pytest.approx(0.5)
    assert row["max_granted_memory_gb"] == 64
    assert row["runs_with_grant"] == 1


def test_a_step_over_its_grant_is_called_out(tmp_path):
    # The kernel does not enforce -m; a tool that ignores it goes over. That is worth naming,
    # because no amount of raising the budget fixes a tool that does not read it.
    found = calibrate.observations([study(tmp_path, {"assembly": step(peak_gb=96, granted_gb=64)})])
    assert "EXCEEDED" in calibrate._budget_verdict(calibrate.summarise(found)[0])


def test_a_step_well_under_its_grant_says_so(tmp_path):
    found = calibrate.observations([study(tmp_path, {"binning": step(peak_gb=4, granted_gb=64)})])
    assert "under a quarter" in calibrate._budget_verdict(calibrate.summarise(found)[0])


def test_a_step_near_its_grant_says_so(tmp_path):
    found = calibrate.observations([study(tmp_path, {"assembly": step(peak_gb=60, granted_gb=64)})])
    assert "close to" in calibrate._budget_verdict(calibrate.summarise(found)[0])


def test_runs_without_a_grant_are_not_counted_as_zero(tmp_path):
    # A study from before the budget existed must produce no comparison at all, rather than a
    # division by zero or a fabricated 0 GB grant.
    found = calibrate.observations([study(tmp_path, {"assembly": step(peak_gb=40)})])
    row = calibrate.summarise(found)[0]
    assert row["runs_with_grant"] == 0
    assert row["max_memory_used_fraction"] is None
    assert calibrate._budget_verdict(row) == "no grant recorded"


def test_thread_fraction_shows_a_step_that_did_not_thread(tmp_path):
    found = calibrate.observations(
        [
            study(
                tmp_path,
                {"binning": step(peak_gb=4, granted_gb=8, granted_threads=24, user_seconds=120.0)},
            )
        ]
    )
    row = calibrate.summarise(found)[0]
    assert row["max_thread_used_fraction"] == pytest.approx(1.2 / 24)


def test_the_worst_run_is_the_one_reported(tmp_path):
    # Several runs of a module: the budget has to cover the worst of them, not the average.
    found = calibrate.observations(
        [
            study(
                tmp_path,
                {
                    "assembly:a": step(peak_gb=10, granted_gb=64),
                    "assembly:b": step(peak_gb=58, granted_gb=64),
                },
            )
        ]
    )
    row = calibrate.summarise(found)[0]
    assert row["runs_with_grant"] == 2
    assert row["max_memory_used_fraction"] == pytest.approx(58 / 64)


# --- the report ----------------------------------------------------------------------------


def test_the_report_prints_the_budget_table(tmp_path, capsys):
    rc = calibrate.main(
        [study(tmp_path, {"assembly": step(peak_gb=40, granted_gb=64, granted_threads=16)})]
    )
    out = capsys.readouterr().out
    assert rc == 0
    assert "granted" in out and "64 GB" in out
    assert "max_memory_fraction" in out  # tells the reader what to change


def test_the_report_says_when_there_is_nothing_to_compare(tmp_path, capsys):
    calibrate.main([study(tmp_path, {"assembly": step(peak_gb=40)})])
    assert "No step recorded the budget it was granted" in capsys.readouterr().out


def test_per_step_lists_the_grant(tmp_path, capsys):
    calibrate.main(
        [
            study(
                tmp_path,
                {
                    "assembly": step(
                        peak_gb=40, granted_gb=64, granted_threads=16, user_seconds=800.0
                    )
                },
            ),
            "--per-step",
        ]
    )
    out = capsys.readouterr().out
    assert "8.0/16" in out  # cores busy against threads granted


def test_json_output_includes_the_grant(tmp_path, capsys):
    import json

    calibrate.main([study(tmp_path, {"assembly": step(peak_gb=40, granted_gb=64)}), "--json"])
    parsed = json.loads(capsys.readouterr().out)
    assert parsed["summary"][0]["max_granted_memory_gb"] == 64
