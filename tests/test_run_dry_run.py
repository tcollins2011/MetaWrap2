"""`run --dry-run` and `cohort --dry-run` print the plan and write nothing under -o."""

from __future__ import annotations

import json
import os

import pytest

from metawrap2 import command as _cmd
from metawrap2.commands import cohort as cohort_cmd
from metawrap2.commands import run as driver
from metawrap2.manifest import COMPLETED, MANIFEST_NAME, SCHEMA_VERSION

SHEET = """
[settings]
coassemble = %s

[[samples]]
name = "A"
r1 = "raw/A_1.fastq"
r2 = "raw/A_2.fastq"

[[samples]]
name = "B"
r1 = "raw/B_1.fastq"
r2 = "raw/B_2.fastq"
"""


@pytest.fixture(autouse=True)
def _runner(tmp_path, monkeypatch):
    monkeypatch.setenv("METAWRAP2_HISTORY", str(tmp_path / "history.jsonl"))
    monkeypatch.setattr(driver, "preflight_databases", lambda *args: None)
    _cmd.runner.configure(resume=False, force=False, dry_run=False)
    yield
    _cmd.runner.configure(resume=False, force=False, dry_run=False)


def _sheet(tmp_path, coassemble=True):
    raw = tmp_path / "raw"
    raw.mkdir(exist_ok=True)
    for name in ("A_1", "A_2", "B_1", "B_2"):
        (raw / ("%s.fastq" % name)).write_text("@r\nACGT\n+\nIIII\n")
    sheet = tmp_path / "samples.toml"
    sheet.write_text(SHEET % ("true" if coassemble else "false"))
    return str(sheet)


def _snapshot(root):
    """Every path under *root* with its mtime and content."""
    seen = {}
    for parent, dirs, files in os.walk(root):
        for name in dirs + files:
            path = os.path.join(parent, name)
            stat = os.lstat(path)
            content = None
            if os.path.isfile(path) and not os.path.islink(path):
                with open(path, "rb") as fh:
                    content = fh.read()
            seen[os.path.relpath(path, root)] = (stat.st_mtime_ns, content)
    return seen


@pytest.mark.parametrize(
    "command, coassemble",
    [(driver.main, True), (driver.main, False), (cohort_cmd.main, False)],
)
def test_a_dry_run_into_a_new_directory_creates_nothing(tmp_path, capsys, command, coassemble):
    out = tmp_path / "study"
    rc = command([_sheet(tmp_path, coassemble), "-o", str(out), "--dry-run"])
    assert rc == 0
    assert not out.exists()
    printed = capsys.readouterr().out
    # The full plan is still printed: every step's module command.
    assert printed.count("[plan]") >= 4
    assert "read_qc" in printed and "assembly" in printed and "binning" in printed


def test_a_global_dry_run_creates_nothing_either(tmp_path):
    """`metawrap2 --dry-run run ...` sets the shared runner instead of the flag after `run`."""
    _cmd.runner.configure(resume=False, force=False, dry_run=True)
    out = tmp_path / "study"
    assert driver.main([_sheet(tmp_path), "-o", str(out)]) == 0
    assert not out.exists()


def _existing_study(tmp_path):
    out = tmp_path / "study"
    (out / "LOGS").mkdir(parents=True)
    (out / "LOGS" / "read_qc_A.log").write_text("$ metawrap2 read_qc ...\n")
    (out / "CLEAN_READS").mkdir()
    (out / "CLEAN_READS" / "ALL_READS_1.fastq").write_text("@r\nACGT\n+\nIIII\n")
    manifest = {
        "schema": SCHEMA_VERSION,
        "steps": {"read_qc:A": {"status": COMPLETED, "started": 1.0, "finished": 2.0}},
    }
    state_dir = driver.Layout(str(out)).state_dir
    os.makedirs(state_dir)
    with open(os.path.join(state_dir, MANIFEST_NAME), "w") as fh:
        json.dump(manifest, fh)
    (out / "provenance.txt").write_text("old provenance\n")
    # Old timestamps, so a rewrite with the same content still shows up.
    for path in _snapshot(out):
        os.utime(os.path.join(out, path), (1_000_000, 1_000_000), follow_symlinks=False)
    os.utime(out, (1_000_000, 1_000_000))
    return out


@pytest.mark.parametrize("extra", [(), ("--resume",), ("--force",)])
@pytest.mark.parametrize("command", [driver.main, cohort_cmd.main])
def test_a_dry_run_leaves_an_existing_study_untouched(tmp_path, command, extra):
    out = _existing_study(tmp_path)
    before = _snapshot(out)
    root_mtime = os.stat(out).st_mtime_ns
    rc = command([_sheet(tmp_path), "-o", str(out), "--dry-run", *extra])
    assert rc == 0
    assert _snapshot(out) == before
    assert os.stat(out).st_mtime_ns == root_mtime
