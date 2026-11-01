"""Database readiness: tool environments, `check`, `doctor` and the sample-sheet preflight.

Covers:

* the configured database paths reach each tool's environment (GTDB-Tk reads only
  GTDBTK_DATA_PATH, which a rebuilt conda env loses), with the config winning over the shell;
* tool output that is not UTF-8 does not stop a run;
* `check <module>` lists only that module's databases, calls one optional only when a flag skips
  it, and fails with the install command when a required one is unset;
* `doctor` states which modules its database table covers, and its counts match its rows;
* `viruses` names --genomad-db and the install command when geNomad's database is unset;
* `run`/`cohort` stop before the first step when the plan needs a database that is not ready.

No conda env, tool or database is used: env and tool lookups are stubbed and databases are
temporary directories.
"""

from __future__ import annotations

import sys

import pytest

from metawrap2 import command, dbcheck
from metawrap2.config import Settings

# -- database paths in the tool's environment ------------------------------------------------------


def test_configured_database_paths_reach_the_tool_environment(monkeypatch):
    monkeypatch.setenv("GTDBTK_DATA_PATH", "/from/the/shell")
    monkeypatch.setenv("CHECKVDB", "/shell/checkv")
    runner = command.CommandRunner(
        databases={"GTDBTK_DATA_PATH": "/db/gtdbtk", "CHECKV_DB": "/db/checkv", "KRAKEN2_DB": "/k"}
    )
    child = runner._child_env(isolate_user_site=True)
    assert child is not None
    # The config wins over the shell, because it is what `check` reports.
    assert child["GTDBTK_DATA_PATH"] == "/db/gtdbtk"
    assert child["CHECKVDB"] == "/db/checkv"
    assert child["PYTHONNOUSERSITE"] == "1"
    # Kraken2 takes its database as an argument; no variable is invented for it.
    assert "KRAKEN2_DB" not in child


def test_an_unset_database_leaves_the_users_variable_alone(monkeypatch):
    monkeypatch.setenv("BAKTA_DB", "/shell/bakta")
    child = command.CommandRunner(databases={})._child_env(isolate_user_site=True)
    assert child is not None and child["BAKTA_DB"] == "/shell/bakta"


def test_a_per_call_variable_wins_over_the_database_mapping():
    runner = command.CommandRunner(databases={"CHECKM_DB": "/db/checkm"})
    child = runner._child_env(False, {"CHECKM_DATA_PATH": "/explicit"})
    assert child is not None and child["CHECKM_DATA_PATH"] == "/explicit"


def test_the_running_modules_config_is_used_when_the_runner_has_none():
    class Recorder:
        def __init__(self):
            self.databases = {"GTDBTK_DATA_PATH": "/recorded/gtdbtk"}

        def record_command(self, cmd):
            pass

    runner = command.CommandRunner(recorder=Recorder())
    child = runner._child_env(isolate_user_site=False)
    assert child is not None and child["GTDBTK_DATA_PATH"] == "/recorded/gtdbtk"


def test_a_tool_sees_the_variable_when_it_runs(tmp_path):
    out = tmp_path / "seen.txt"
    runner = command.CommandRunner(databases={"GTDBTK_DATA_PATH": "/db/gtdbtk"})
    runner.configure(heartbeat_seconds=0)
    runner.run(
        [sys.executable, "-c", "import os; print(os.environ.get('GTDBTK_DATA_PATH', ''))"],
        tool="probe",
        stdout_path=str(out),
    )
    assert out.read_text().strip() == "/db/gtdbtk"


# -- tool output that is not UTF-8 --------------------------------------------------------------


def test_a_byte_that_is_not_utf8_does_not_stop_the_run(capsys):
    runner = command.CommandRunner(databases={})
    runner.configure(heartbeat_seconds=0)
    script = "import sys; sys.stdout.buffer.write(b'before \\xff after\\n'); sys.stdout.flush()"
    assert runner.run([sys.executable, "-c", script], tool="latin1-tool") == 0
    out = capsys.readouterr().out
    assert "before � after" in out


# -- `check <module>` ---------------------------------------------------------------------------


@pytest.fixture
def check_cmd(monkeypatch):
    """`check` with every env and tool present, so only the database verdict varies."""
    from metawrap2.commands import check

    monkeypatch.setattr(check, "conda_env_exists", lambda env: True)
    monkeypatch.setattr(check, "tool_in_env", lambda tool, env: "/env/bin/" + tool)
    monkeypatch.setattr(check.shutil, "which", lambda tool: "/bin/" + tool)
    return check


def _config(tmp_path, **databases):
    path = tmp_path / "metawrap2.toml"
    lines = ["[databases]"] + ['%s = "%s"' % (k, v) for k, v in databases.items()]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_check_viruses_fails_when_genomad_db_is_unset(check_cmd, tmp_path, capsys):
    rc = check_cmd.main(["viruses", "--config", _config(tmp_path)])
    out = capsys.readouterr().out
    assert rc == 1
    assert out.rstrip().splitlines()[-1] == (
        "Preflight FAILED for viruses: GENOMAD_DB is not set "
        "(install with 'metawrap2 install-db genomad')"
    )
    # CheckV has a skip flag, so it is optional and the label names the flag.
    checkv = [line for line in out.splitlines() if "CHECKV_DB" in line]
    assert len(checkv) == 1 and "(optional - skip with --skip-checkv)" in checkv[0]
    genomad = [line for line in out.splitlines() if "GENOMAD_DB" in line and "[" in line]
    assert "optional" not in genomad[0]


def test_check_lists_only_the_modules_databases(check_cmd, tmp_path, capsys):
    genomad = tmp_path / "genomad_db"
    genomad.mkdir()
    rc = check_cmd.main(["viruses", "--config", _config(tmp_path, GENOMAD_DB=genomad)])
    out = capsys.readouterr().out
    assert rc == 0
    assert "Preflight OK for: viruses" in out
    rows = [line for line in out.split("Databases:")[1].splitlines() if line.strip()[:1] == "["]
    assert [row.split("]")[1].split()[0] for row in rows] == ["GENOMAD_DB", "CHECKV_DB"]
    assert "not needed by viruses; 'metawrap2 doctor --all' lists them." in out
    assert "KRAKEN2_DB" not in out


def test_a_module_with_no_database_says_so(check_cmd, capsys):
    assert check_cmd.main(["dereplicate"]) == 0
    assert "reads no database" in capsys.readouterr().out


def test_opt_in_databases_name_the_flag_that_needs_them(check_cmd, tmp_path, capsys):
    assert check_cmd.main(["annotate", "--config", _config(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "(optional - only with --kofam)" in out
    assert "(optional - only with --cazy)" in out


# -- `doctor` -----------------------------------------------------------------------------------


def _db_rows(out):
    table = out.split("Database status")[1]
    body = table.split("-" * 70)[1]
    return [line for line in body.splitlines() if line.strip() and not line.startswith(" " * 22)]


def test_doctor_states_its_database_scope_and_counts_its_rows(tmp_path, capsys):
    from metawrap2.commands import doctor
    from metawrap2.commands.run import DEFAULT_STEPS

    modules, scope = doctor._database_scope([], named=False, every=False)
    assert modules == list(DEFAULT_STEPS)
    doctor._report_databases(modules, _config(tmp_path), scope=scope)
    out = capsys.readouterr().out
    assert "Database status (default pipeline; --all for every module)" in out
    rows = _db_rows(out)
    assert "%d of %d usable" % (0, len(rows)) in out
    # Every default-pipeline database has a skip flag, so none fails doctor.
    assert "0 required and unusable, %d optional and unusable" % len(rows) in out


def test_doctor_all_covers_every_modules_databases(tmp_path, capsys):
    from metawrap2.commands import doctor

    modules, scope = doctor._database_scope(["read_qc"], named=True, every=True)
    bad = doctor._report_databases(modules, _config(tmp_path), scope=scope)
    out = capsys.readouterr().out
    assert "(every module" in out
    keys = {row.split()[0] for row in _db_rows(out)}
    assert {"GENOMAD_DB", "KRAKEN2_DB", "SINGLEM_METAPACKAGE", "GTDBTK_DATA_PATH"} <= keys
    # Required and missing: KRAKEN2_DB, BLASTDB, TAXDUMP, GENOMAD_DB, SINGLEM_METAPACKAGE.
    assert bad == 5
    assert "%d required and unusable" % bad in out


# -- `viruses` ----------------------------------------------------------------------------------


def test_viruses_names_the_flag_and_install_command_without_genomad(capsys):
    from metawrap2.modules import viruses

    with pytest.raises(SystemExit):
        viruses._database(None, Settings(), "GENOMAD_DB", "geNomad", "--genomad-db DIR", "genomad")
    err = capsys.readouterr()
    text = err.out + err.err
    assert "--genomad-db" in text and "metawrap2 install-db genomad" in text


# -- preflight of a planned run -----------------------------------------------------------------


def test_sheet_flags_decide_which_databases_are_needed():
    settings = Settings()
    errors, _ = dbcheck.preflight(settings, [("bin_refinement", ["--skip-gunc"])])
    assert errors == []
    errors, _ = dbcheck.preflight(settings, [("bin_refinement", [])])
    assert errors == [
        "GUNC_DB is not set (install with 'metawrap2 install-db gunc') (bin_refinement)"
    ]
    errors, _ = dbcheck.preflight(settings, [("classify_bins", ["--gtdbtk"])])
    assert [e.split()[0] for e in errors] == ["GTDBTK_DATA_PATH"]
    errors, _ = dbcheck.preflight(settings, [("classify_bins", [])])
    assert [e.split()[0] for e in errors] == ["BLASTDB", "TAXDUMP"]


def test_a_database_needed_by_many_steps_is_reported_once():
    steps = [("read_qc", ["-o", "a"]), ("read_qc", ["-o", "b"])]
    errors, _ = dbcheck.preflight(Settings(), steps)
    assert len(errors) == 1


def test_a_path_flag_is_checked_instead_of_the_config(tmp_path):
    errors, _ = dbcheck.preflight(
        Settings(), [("viruses", ["--genomad-db", str(tmp_path), "--skip-checkv"])]
    )
    assert errors == []
    errors, _ = dbcheck.preflight(
        Settings(), [("viruses", ["--genomad-db=/no/such/dir", "--skip-checkv"])]
    )
    assert errors[0].startswith("--genomad-db points at /no/such/dir")


def test_an_unset_checkm2db_only_warns():
    errors, warnings = dbcheck.preflight(Settings(), [("reassemble_bins", [])])
    assert errors == [] and [w.split()[0] for w in warnings] == ["CHECKM2DB"]


def test_run_stops_before_the_first_step_and_dry_run_warns(capsys):
    from metawrap2.commands import run

    steps = [run.Step(name="read_qc:s1", module="read_qc", argv=["-1", "r.fq", "-o", "x"])]
    with pytest.raises(SystemExit):
        run.preflight_databases(steps, Settings(), dry_run=False)
    out = capsys.readouterr()
    assert "metawrap2 install-db host" in out.out + out.err
    run.preflight_databases(steps, Settings(), dry_run=True)
    out = capsys.readouterr()
    assert "A run without --dry-run stops here." in out.out + out.err
    skipped = [run.Step(name="read_qc:s1", module="read_qc", argv=["--skip-bmtagger"])]
    run.preflight_databases(skipped, Settings(), dry_run=False)
