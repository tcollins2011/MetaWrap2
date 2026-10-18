"""What a first install reads: step numbering, environment counts, and messages that say nothing.

These are the output contracts of `quickstart` and `install-env` - the two commands a new user
runs before anything else works. Each one here was a real report from a first install: a step
number that skipped from 1 to 3, an environment count that counted modules, a skip message that
did not say which conda installation it had found, and half an hour of silence in the slowest
step.

Nothing here touches conda, the network, or a database. Every test fakes the step it is not
about, so a failure names one behaviour.
"""

from __future__ import annotations

import pytest

from metawrap2 import cli
from metawrap2.commands import config_cmd, install_env, quickstart

try:  # Python 3.11+; the package still supports 3.10, where tomllib does not exist
    import tomllib as _toml
except ModuleNotFoundError:  # pragma: no cover
    import tomli as _toml


def flat(text: str) -> str:
    """Captured output with the banner borders and the line wrapping taken back out.

    Messages are wrapped to the banner's text column, so asserting on a sentence means either
    knowing where it wraps or removing the wrapping. The second is the stable one.
    """
    kept = []
    for raw in text.split("\n"):
        line = raw.strip()
        if line and set(line) <= set("-*#="):
            continue  # a border row
        kept.append(line.strip("-*#").strip())
    return " ".join(" ".join(kept).split())


@pytest.fixture
def quickstart_calls(monkeypatch):
    """Let `quickstart.main` run end to end while recording what it would have invoked."""
    seen = {"env": None, "db": None, "doctor": None, "doctor_modules": None}

    def fake_install_env(modules, threads, from_lock, log_dir=None):
        seen["env"] = {
            "modules": list(modules),
            "threads": threads,
            "from_lock": from_lock,
            "log_dir": log_dir,
        }
        return 0

    monkeypatch.setattr(quickstart, "_package_manager", lambda: "mamba")
    monkeypatch.setattr(quickstart, "_run_install_env", fake_install_env)
    monkeypatch.setattr(
        quickstart,
        "_run_install_db",
        lambda names, directory, small, jobs, config: seen.__setitem__("db", list(names)) or 0,
    )

    def fake_doctor(modules, config, databases):
        seen["doctor"] = databases
        seen["doctor_modules"] = list(modules)
        return 0

    monkeypatch.setattr(quickstart, "_run_doctor", fake_doctor)
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    return seen


# --- step numbering -------------------------------------------------------------------------


def test_a_skipped_step_is_numbered_and_says_why(quickstart_calls, capsys):
    """`--databases none` printed STEP 1/3 then STEP 3/3, so step 2 looked like a silent failure."""
    assert quickstart.main(["--binning-only", "--databases", "none"]) == 0
    out = flat(capsys.readouterr().out)
    assert "STEP 1/3" in out
    assert "STEP 2/3: databases (skipped: --databases none)" in out
    assert "STEP 3/3" in out


def test_skipping_environments_still_starts_at_step_one(quickstart_calls, capsys):
    """`--skip-envs` used to print STEP 3/3 as the first thing the command said."""
    assert quickstart.main(["--binning-only", "--skip-envs", "--databases", "none"]) == 0
    out = flat(capsys.readouterr().out)
    assert "STEP 1/3: conda environments (skipped: --skip-envs)" in out
    assert quickstart_calls["env"] is None, "--skip-envs must not create environments"


# --- the environment count ------------------------------------------------------------------


def test_the_step_banner_counts_environments_not_modules(quickstart_calls, capsys):
    """ "creating 23 conda environment(s)" counted module names; 11 environments were built.

    The genome-recovery path names seven things that resolve to four environments, so the count
    in the banner and the number of names handed to install-env have to differ.
    """
    assert quickstart.main(["--binning-only", "--databases", "none"]) == 0
    out = flat(capsys.readouterr().out)
    envs = quickstart._unique_envs(quickstart_calls["env"]["modules"])
    assert len(envs) < len(quickstart_calls["env"]["modules"]), "this path must share envs"
    assert "creating %d conda environment(s)" % len(envs) in out


def test_the_slow_step_states_how_long_it_took_when_measured(quickstart_calls, capsys):
    """Silence with no expected duration is an unbounded wait; the measurement bounds it."""
    assert quickstart.main(["--binning-only", "--databases", "none"]) == 0
    out = flat(capsys.readouterr().out)
    assert "29 minutes" in out


# --- install-env: one line per environment, naming where it is -------------------------------


@pytest.fixture
def existing_envs(monkeypatch):
    """Make every environment look already built, at a predictable prefix."""
    monkeypatch.setattr(install_env, "ensure_mamba", lambda: True)
    monkeypatch.setattr(install_env, "conda_env_exists", lambda env: True)
    monkeypatch.setattr(install_env, "env_prefix", lambda env: "/opt/conda/envs/%s" % env)
    monkeypatch.setattr(install_env, "prefix_in_active_root", lambda prefix: True)
    monkeypatch.setattr(
        install_env, "_create_one", lambda *a, **k: pytest.fail("nothing should be built")
    )


def test_an_existing_environment_is_reported_once_with_its_location(existing_envs, capsys):
    """Three modules share metawrap2-reads: the skip path printed the same line three times.

    The build path already deduplicated; the "already exists" branch returned before reaching the
    bookkeeping that does it. The location is part of the message because `conda_env_exists`
    answers across every conda installation on the machine, so the name alone cannot tell a user
    whether the environment found is the one they will run against.
    """
    assert install_env.main(["read_qc", "assembly", "reassemble_bins", "--no-test"]) == 0
    out = capsys.readouterr().out
    expected = (
        "metawrap2-reads already exists at /opt/conda/envs/metawrap2-reads; "
        "skipping (use --force to recreate)."
    )
    exists = [line for line in out.split("\n") if "already exists" in line]
    assert exists == [expected]
    assert out.count("already covered by this run") == 2


def test_an_environment_from_another_conda_installation_says_so(existing_envs, capsys, monkeypatch):
    """A prefix outside the active root is usable but is not the user's; the name hides that."""
    monkeypatch.setattr(install_env, "prefix_in_active_root", lambda prefix: False)
    assert install_env.main(["read_qc", "--no-test"]) == 0
    out = capsys.readouterr().out
    assert "/opt/conda/envs/metawrap2-reads" in out
    assert "built by a different conda" in out


def test_each_environment_announces_its_start_and_its_elapsed_time(monkeypatch, capsys, tmp_path):
    """The slowest step in the install printed nothing at all between its banner and its result.

    A conda solve's output goes to a per-env log that is only shown on failure, so success was 29
    minutes of silence. The counter is over the real number of environments, and the env name is
    on both lines because parallel workers interleave them.
    """
    monkeypatch.setattr(install_env, "ensure_mamba", lambda: True)
    monkeypatch.setattr(install_env, "conda_env_exists", lambda env: False)
    monkeypatch.setattr(
        install_env,
        "_create_one",
        lambda module, env, yaml, force, log_dir, lock=None: (module, True, "%s/x.log" % log_dir),
    )
    argv = ["read_qc", "binning", "--no-test", "-t", "1", "--log-dir", str(tmp_path)]
    assert install_env.main(argv) == 0
    out = capsys.readouterr().out
    assert "[1/2] metawrap2-reads" in out
    assert "[2/2] metawrap2-bins" in out
    assert out.count("solving...") == 2
    assert out.count("done in") == 2


def test_the_elapsed_time_reads_as_minutes_and_seconds():
    """ "done in 192.0 seconds" is arithmetic; "done in 3m12s" is an answer."""
    assert install_env.format_duration(45) == "45s"
    assert install_env.format_duration(192) == "3m12s"
    assert install_env.format_duration(3840) == "1h04m"


# --- quickstart forwards --log-dir -----------------------------------------------------------


def test_log_dir_is_forwarded_to_install_env(quickstart_calls, tmp_path):
    """`install-env --log-dir` existed and quickstart could not ask for it.

    Without it, an install that succeeded with warnings left its only record in a temp directory
    whose path was printed once and never again.
    """
    argv = ["--binning-only", "--databases", "none", "--log-dir", str(tmp_path)]
    assert quickstart.main(argv) == 0
    assert quickstart_calls["env"]["log_dir"] == str(tmp_path)


# --- the non-interactive database default ----------------------------------------------------


def test_a_non_interactive_run_says_it_is_downloading_nothing(quickstart_calls, capsys):
    """It downloaded nothing, correctly, and never said so - which reads as a silent failure.

    The behaviour must not change: under sbatch or in CI there is nobody to answer a prompt, and
    a 300 GB download must never start unattended. Only the sentence is new.
    """
    assert quickstart.main(["--binning-only"]) == 0
    out = flat(capsys.readouterr().out)
    assert "No --databases given and stdin is not a terminal; downloading nothing." in out
    assert "--databases small" in out
    assert quickstart_calls["db"] is None


def test_skip_envs_reports_database_status_not_environment_status(
    quickstart_calls, capsys, tmp_path
):
    """`--skip-envs` is documented as "databases only" and reported on environments.

    It printed the module table - every row missing, because creating the environments had just
    been skipped on purpose - and exited 1 for a run that did exactly what was asked.
    """
    config = tmp_path / "metawrap2.toml"
    config.write_text("[databases]\n")
    argv = ["--binning-only", "--skip-envs", "--databases", "none", "--config", str(config)]
    assert quickstart.main(argv) == 0
    out = flat(capsys.readouterr().out)
    assert "MetaWrap2 database status" in out
    assert quickstart_calls["doctor"] is None, "the module table is not what --skip-envs is for"


# --- the two modules that cannot have a row --------------------------------------------------


def test_a_full_quickstart_asks_about_every_module(quickstart_calls, capsys):
    """The readiness table has to cover what `metawrap2 -h` lists, or the two cannot be compared.

    It used to pass the conda environment names, which held the table at 22 rows of 24 - and the
    two it dropped were the modules with no environment, one of which is gated on a Python extra
    that can actually be missing. So the module missing from the table was precisely the one with
    an unmet dependency. A full quickstart now names no modules at all and lets `doctor` report
    everything it knows about.
    """
    assert quickstart.main(["--databases", "none"]) == 0
    assert (
        quickstart_calls["doctor_modules"] == []
    ), "a full quickstart must not narrow the table to a module list; " "got %r" % (
        quickstart_calls["doctor_modules"],
    )
    absent = [name for name in cli.MODULES if name not in quickstart.MODULE_ENVS]
    assert absent, "this test is about modules with no env; there have to be some"


# --- the version string ----------------------------------------------------------------------


def test_the_version_names_the_commit_in_a_checkout(monkeypatch, capsys):
    """`MetaWrap2 v=2.1.2` from a clone one commit past the tag made every bug report ambiguous."""
    monkeypatch.setattr(cli, "git_revision", lambda: "16f3d14")
    assert cli.main(["-v"]) == 0
    assert capsys.readouterr().out.strip() == "MetaWrap2 %s (git 16f3d14)" % cli.__version__


def test_the_version_falls_back_to_the_plain_number_without_git(monkeypatch, capsys):
    """A wheel install has no repository and no git binary, and `-v` must still answer."""
    monkeypatch.setattr(cli, "git_revision", lambda: None)
    assert cli.main(["-v"]) == 0
    printed = capsys.readouterr().out
    assert printed.strip() == "MetaWrap2 %s" % cli.__version__
    assert "v=" not in printed, "`v=` is not the convention and is awkward to paste"


@pytest.mark.parametrize(
    "outcome",
    [
        OSError("no git on PATH"),
        (128, ""),  # git in a directory that is not a repository
        (0, "not a sha at all"),  # some future git saying something unexpected
    ],
)
def test_a_version_string_never_fails_however_git_behaves(monkeypatch, outcome, capsys):
    """`-v` has to answer on any machine: no repository, no git binary, no surprises quoted back."""
    import subprocess

    def fake_run(_cmd, **_kwargs):
        if isinstance(outcome, Exception):
            raise outcome
        code, text = outcome
        return subprocess.CompletedProcess(_cmd, code, stdout=text)

    monkeypatch.setattr(subprocess, "run", fake_run)
    assert cli.git_revision() is None
    assert cli.version_string() == "MetaWrap2 %s" % cli.__version__


def test_every_global_flag_in_the_limit_table_is_in_the_help_text():
    """`--disk-check-interval` was accepted, documented in the config, and absent from -h."""
    help_text = cli.help_message()
    for flag in cli._LIMIT_FLAGS:
        assert flag in help_text, flag


# --- the starter config has no live placeholder paths ----------------------------------------


def test_a_fresh_config_has_no_active_database_key(tmp_path):
    """Four placeholders were uncommented, so a brand-new config had four "missing" databases.

    "Configured but missing" and "not configured" need different actions from the user, and an
    active `/path/to/my/kraken2_database` collapses them. It also left `install-db` to decide
    whether to replace the key it finds or write a second one.
    """
    path = tmp_path / "config.toml"
    assert config_cmd.main(["init", "--path", str(path)]) == 0
    with open(path, "rb") as fh:
        assert _toml.load(fh).get("databases", {}) == {}


def test_config_show_reports_a_fresh_config_as_not_set_rather_than_missing(tmp_path, capsys):
    """The symptom: `config show` on an untouched config listed four databases as [missing]."""
    path = tmp_path / "config.toml"
    assert config_cmd.main(["init", "--path", str(path)]) == 0
    capsys.readouterr()
    assert config_cmd.main(["show", "--config", str(path)]) == 0
    out = capsys.readouterr().out
    assert "[missing]" not in out
    assert "<none set>" in out
