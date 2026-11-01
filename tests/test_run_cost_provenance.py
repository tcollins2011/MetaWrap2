"""Run cost (CPU-hours, wall-clock, peak memory) in provenance.txt.

A module's provenance.txt gives the cost of its own run and of each MetaWrap2 step it inherits.
The study's provenance.txt, written by `run` and `cohort`, gives the catalogue strategy, each
step's cost and the total over the completed steps in the manifest.
"""

from __future__ import annotations

import json

import pytest

from metawrap2 import command as _cmd
from metawrap2.commands import cohort as cohort_cmd
from metawrap2.commands import run as driver
from metawrap2.manifest import COMPLETED, FAILED, MANIFEST_NAME, Manifest
from metawrap2.modules import _common
from metawrap2.provenance import RunRecorder, describe_cost, total_usage

GB = 1024**3


def _recorder(module, inputs=()):
    return RunRecorder(
        module=module,
        version="2.1.5",
        parameters={"threads": 4},
        config_file=None,
        conda_env=None,
        databases={},
        inputs=list(inputs),
        command_line=["metawrap2", module],
    )


def test_describe_cost_gives_cpu_hours_wall_clock_and_peak_memory():
    usage = {
        "wall_seconds": 723.0,
        "user_seconds": 4000.0,
        "system_seconds": 500.0,
        "peak_rss_bytes": 3 * GB,
    }
    assert describe_cost(usage) == "1.25 CPU-hours, 12m03s wall-clock, 3.0 GB peak memory"
    # A platform that measures wall time only.
    assert describe_cost({"wall_seconds": 5}) == "5s wall-clock"
    assert describe_cost({}) == ""


def test_total_sums_time_and_takes_the_largest_peak():
    total = total_usage(
        [
            {"wall_seconds": 60, "user_seconds": 3000, "system_seconds": 600, "peak_rss_bytes": GB},
            {"wall_seconds": 30, "user_seconds": 1800, "peak_rss_bytes": 4 * GB},
            {},
        ]
    )
    assert total == {"wall_seconds": 90.0, "cpu_seconds": 5400.0, "peak_rss_bytes": 4 * GB}
    assert total_usage([{}, {}]) == {}


def test_a_module_provenance_gives_its_own_cost(tmp_path):
    rec = _recorder("binning")
    out = tmp_path / "binning_out"
    usage = {"wall_seconds": 90.0, "user_seconds": 7200.0, "system_seconds": 0.0}
    usage["peak_rss_bytes"] = 2 * GB
    rec.finalize(str(out), capture_versions=False, usage=usage)

    text = (out / "provenance.txt").read_text()
    assert "    cost:     2.00 CPU-hours, 1m30s wall-clock, 2.0 GB peak memory" in text
    # One step, so no total line.
    assert "Total over" not in text
    config = json.loads((out / "run_config.json").read_text())
    assert config["usage"] == usage
    assert config["lineage"][-1]["usage"] == usage


def test_chained_runs_add_up_in_the_last_provenance(tmp_path):
    """A run whose input is a MetaWrap2 output totals the cost of every costed step."""
    qc = tmp_path / "qc"
    _recorder("read_qc").finalize(
        str(qc),
        capture_versions=False,
        usage={"wall_seconds": 60, "user_seconds": 3600, "peak_rss_bytes": GB},
    )
    asm = tmp_path / "asm"
    _recorder("assembly", [str(qc / "final_pure_reads_1.fastq")]).finalize(
        str(asm),
        capture_versions=False,
        usage={"wall_seconds": 120, "user_seconds": 7200, "peak_rss_bytes": 8 * GB},
    )
    text = (asm / "provenance.txt").read_text()
    assert text.count("    cost:") == 2
    assert (
        "Total over the 2 of 2 steps with a recorded cost: "
        "3.00 CPU-hours, 3m00s wall-clock, 8.0 GB peak memory"
    ) in text


def test_finish_run_records_the_measured_usage(tmp_path, monkeypatch):
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))
    monkeypatch.setattr(_cmd.runner, "dry_run", False)
    measured = {"wall_seconds": 10.0, "user_seconds": 36.0, "system_seconds": 0.0}
    monkeypatch.setattr(_common, "self_usage", lambda rec: dict(measured))
    real_finalize = RunRecorder.finalize
    monkeypatch.setattr(
        RunRecorder,
        "finalize",
        lambda self, out, **kw: real_finalize(self, out, **{**kw, "capture_versions": False}),
    )
    out = tmp_path / "out"
    _common.finish_run(_recorder("read_qc"), str(out), status="completed")
    assert "    cost:     0.01 CPU-hours, 10s wall-clock" in (out / "provenance.txt").read_text()


# -- the study's provenance.txt -------------------------------------------------------------

SHEET = """
[settings]
coassemble = %s

[[samples]]
name = "A"
r1 = "raw/A_1.fastq.gz"
r2 = "raw/A_2.fastq.gz"

[[samples]]
name = "B"
r1 = "raw/B_1.fastq.gz"
r2 = "raw/B_2.fastq.gz"
"""


def _step(manifest, name, started, finished, status, usage):
    manifest.data["steps"][name] = {
        "status": status,
        "started": started,
        "finished": finished,
        "duration_seconds": usage.get("wall_seconds", 1),
        "usage": usage,
    }


def _fake_execute(state_box, fail=False):
    """Stand in for `execute`: record three steps (two completed, one failed) in the manifest."""

    def execute(steps, layout, study, threads, jobs, config, flags, state, dry_run, **kwargs):
        state_box.append(state)
        m = state.manifest
        # read_qc:A and read_qc:B overlap, as with -j 2.
        _step(
            m,
            "read_qc:A",
            "2026-01-01T10:00:00",
            "2026-01-01T11:00:00",
            COMPLETED,
            {
                "wall_seconds": 3600,
                "user_seconds": 7000,
                "system_seconds": 200,
                "peak_rss_bytes": GB,
            },
        )
        _step(
            m,
            "read_qc:B",
            "2026-01-01T10:30:00",
            "2026-01-01T12:00:00",
            COMPLETED,
            {"wall_seconds": 5400, "user_seconds": 10800, "peak_rss_bytes": 6 * GB},
        )
        _step(
            m,
            "assembly:coassembly",
            "2026-01-01T12:00:00",
            "2026-01-01T15:00:00",
            FAILED,
            {"wall_seconds": 10800, "user_seconds": 99999, "peak_rss_bytes": 99 * GB},
        )
        m.save()
        if fail:
            driver.error("1 of 1 step(s) failed; stopping.")
        return 0

    return execute


def _run_study(tmp_path, monkeypatch, command=driver.main, coassemble=True, fail=False, extra=()):
    states: list = []
    monkeypatch.setattr(driver, "execute", _fake_execute(states, fail=fail))
    monkeypatch.setattr(driver, "preflight_databases", lambda *args: None)
    sheet = tmp_path / "samples.toml"
    sheet.write_text(SHEET % ("true" if coassemble else "false"))
    out = tmp_path / "study"
    rc = command([str(sheet), "-o", str(out), *extra])
    return rc, out


def test_run_writes_the_study_total_and_strategy(tmp_path, monkeypatch):
    rc, out = _run_study(tmp_path, monkeypatch)
    assert rc == 0
    text = (out / "provenance.txt").read_text()
    assert "catalogue strategy: co-assembly, 2 sample(s), dereplication off" in text
    assert "Total over 2 completed step(s):" in text
    # (7000 + 200 + 10800) / 3600; the failed assembly is not counted.
    assert "  CPU-hours:    5.00" in text
    # 10:00 to 12:00, not the 2h30m the two overlapping steps sum to.
    assert "  wall-clock:   2h00m (first step start to last step finish)" in text
    assert "  step time:    2h30m (summed over steps)" in text
    assert "  peak memory:  6.0 GB (read_qc:B)" in text
    # Each step is listed, the failed one without a cost.
    assert "read_qc:A" in text and "1.0 GB peak memory" in text
    failed_line = next(line for line in text.splitlines() if "assembly:coassembly" in line)
    assert "failed" in failed_line and "CPU-hours" not in failed_line


def test_cohort_records_the_per_sample_strategy(tmp_path, monkeypatch):
    rc, out = _run_study(tmp_path, monkeypatch, command=cohort_cmd.main, coassemble=False)
    assert rc == 0
    text = (out / "provenance.txt").read_text()
    assert "catalogue strategy: per-sample, 2 sample(s), dereplication on" in text


def test_a_failed_study_still_writes_its_provenance(tmp_path, monkeypatch):
    with pytest.raises(SystemExit):
        _run_study(tmp_path, monkeypatch, fail=True)
    text = (tmp_path / "study" / "provenance.txt").read_text()
    assert "Total over 2 completed step(s):" in text


def test_a_dry_run_writes_no_study_provenance(tmp_path, monkeypatch):
    rc, out = _run_study(tmp_path, monkeypatch, extra=("--dry-run",))
    assert rc == 0
    assert not (out / "provenance.txt").exists()


def test_a_study_with_no_completed_step_says_so(tmp_path):
    layout = driver.Layout(str(tmp_path / "study"))
    state = driver.StepState(layout, resume=False)
    study = driver.load_sheet(str(_sheet(tmp_path)))
    path = driver.write_study_provenance(state, layout, study)
    with open(path) as fh:
        text = fh.read()
    assert "No step has completed yet." in text
    assert "catalogue strategy: co-assembly" in text


def _sheet(tmp_path):
    sheet = tmp_path / "samples.toml"
    sheet.write_text(SHEET % "true")
    (tmp_path / "study").mkdir(exist_ok=True)
    return sheet


def test_a_step_recorded_without_usage_counts_its_duration(tmp_path):
    manifest = Manifest(str(tmp_path / MANIFEST_NAME))
    manifest.data["steps"]["old"] = {"status": COMPLETED, "duration_seconds": 42.0}
    assert manifest.completed_usage() == [("old", {"wall_seconds": 42.0})]
    assert manifest.completed_span_seconds() is None
