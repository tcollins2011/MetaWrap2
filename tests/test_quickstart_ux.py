"""Tests for the output of `quickstart` and `install-env`.

Covers step numbering when steps are skipped, environment counts, skip messages that name the conda
installation, progress lines during environment creation, the version string and the starter
config.

No test touches conda, the network or a database; each fakes the steps it does not test.
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

    Messages are wrapped to the banner width, so the wrapping is removed before matching.
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
    """A skipped step 2 is printed as STEP 2/3 with the reason."""
    assert quickstart.main(["--binning-only", "--databases", "none"]) == 0
    out = flat(capsys.readouterr().out)
    assert "STEP 1/3" in out
    assert "STEP 2/3: databases (skipped: --databases none)" in out
    assert "STEP 3/3" in out


def test_skipping_environments_still_starts_at_step_one(quickstart_calls, capsys):
    """With `--skip-envs` the output still starts at STEP 1/3, marked skipped."""
    assert quickstart.main(["--binning-only", "--skip-envs", "--databases", "none"]) == 0
    out = flat(capsys.readouterr().out)
    assert "STEP 1/3: conda environments (skipped: --skip-envs)" in out
    assert quickstart_calls["env"] is None, "--skip-envs must not create environments"


# --- the environment count ------------------------------------------------------------------


def test_the_step_banner_counts_environments_not_modules(quickstart_calls, capsys):
    """The banner counts distinct environments, not module names.

    The genome-recovery path has seven names that resolve to four environments.
    """
    assert quickstart.main(["--binning-only", "--databases", "none"]) == 0
    out = flat(capsys.readouterr().out)
    envs = quickstart._unique_envs(quickstart_calls["env"]["modules"])
    assert len(envs) < len(quickstart_calls["env"]["modules"]), "this path must share envs"
    assert "creating %d conda environment(s)" % len(envs) in out


def test_the_slow_step_states_how_long_it_took_when_measured(quickstart_calls, capsys):
    """The slow environment step prints its measured duration (29 minutes)."""
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
    """An environment shared by three modules is reported once as existing, with its prefix.

    The prefix is included because `conda_env_exists` checks every conda installation on the
    machine.
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
    """An environment outside the active conda root is reported as built by a different conda."""
    monkeypatch.setattr(install_env, "prefix_in_active_root", lambda prefix: False)
    assert install_env.main(["read_qc", "--no-test"]) == 0
    out = capsys.readouterr().out
    assert "/opt/conda/envs/metawrap2-reads" in out
    assert "built by a different conda" in out


def test_each_environment_announces_its_start_and_its_elapsed_time(monkeypatch, capsys, tmp_path):
    """Each environment prints a start line and an elapsed-time line.

    The solver's output goes to a per-env log shown only on failure. The counter is over distinct
    environments, and both lines name the env because parallel workers interleave.
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
    """Durations are formatted as seconds, minutes and seconds, or hours and minutes."""
    assert install_env.format_duration(45) == "45s"
    assert install_env.format_duration(192) == "3m12s"
    assert install_env.format_duration(3840) == "1h04m"


# --- quickstart forwards --log-dir -----------------------------------------------------------


def test_log_dir_is_forwarded_to_install_env(quickstart_calls, tmp_path):
    """quickstart passes --log-dir to install-env, so the solver logs can be kept."""
    argv = ["--binning-only", "--databases", "none", "--log-dir", str(tmp_path)]
    assert quickstart.main(argv) == 0
    assert quickstart_calls["env"]["log_dir"] == str(tmp_path)


# --- the non-interactive database default ----------------------------------------------------


def test_a_non_interactive_run_says_it_is_downloading_nothing(quickstart_calls, capsys):
    """Without a terminal or --databases, nothing is downloaded and the output says so.

    Under sbatch or CI there is no one to answer a prompt, and a 300 GB download must not start
    unattended.
    """
    assert quickstart.main(["--binning-only"]) == 0
    out = flat(capsys.readouterr().out)
    assert "No --databases given and stdin is not a terminal; downloading nothing." in out
    assert "--databases small" in out
    assert quickstart_calls["db"] is None


def test_skip_envs_reports_database_status_not_environment_status(
    quickstart_calls, capsys, tmp_path
):
    """`--skip-envs` reports database status, not the module table, and exits 0."""
    config = tmp_path / "metawrap2.toml"
    config.write_text("[databases]\n")
    argv = ["--binning-only", "--skip-envs", "--databases", "none", "--config", str(config)]
    assert quickstart.main(argv) == 0
    out = flat(capsys.readouterr().out)
    assert "MetaWrap2 database status" in out
    assert quickstart_calls["doctor"] is None, "the module table is not what --skip-envs is for"


# --- the two modules that cannot have a row --------------------------------------------------


def test_a_full_quickstart_asks_about_every_module(quickstart_calls, capsys):
    """A full quickstart passes no module list to `doctor`, so every module is reported.

    Passing environment names left out the modules with no environment, one of which depends on an
    optional Python extra.
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
    """In a git checkout, `-v` includes the commit."""
    monkeypatch.setattr(cli, "git_revision", lambda: "16f3d14")
    assert cli.main(["-v"]) == 0
    assert capsys.readouterr().out.strip() == "MetaWrap2 %s (git 16f3d14)" % cli.__version__


def test_the_version_falls_back_to_the_plain_number_without_git(monkeypatch, capsys):
    """Without git, `-v` prints the plain version."""
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
    """`-v` works with no git binary, outside a repository, or with unexpected git output."""
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
    """Every flag in the limits table appears in -h."""
    help_text = cli.help_message()
    for flag in cli._LIMIT_FLAGS:
        assert flag in help_text, flag


# --- the starter config has no live placeholder paths ----------------------------------------


def test_a_fresh_config_has_no_active_database_key(tmp_path):
    """A new config has no active database keys, only commented placeholders.

    An active placeholder path would be reported as a missing database instead of an unset one.
    """
    path = tmp_path / "config.toml"
    assert config_cmd.main(["init", "--path", str(path)]) == 0
    with open(path, "rb") as fh:
        assert _toml.load(fh).get("databases", {}) == {}


def test_config_show_reports_a_fresh_config_as_not_set_rather_than_missing(tmp_path, capsys):
    """`config show` on a new config lists no database as [missing]."""
    path = tmp_path / "config.toml"
    assert config_cmd.main(["init", "--path", str(path)]) == 0
    capsys.readouterr()
    assert config_cmd.main(["show", "--config", str(path)]) == 0
    out = capsys.readouterr().out
    assert "[missing]" not in out
    assert "<none set>" in out


# --- the default install covers the default bin scorer, CheckM2 ------------------------------


def test_the_binning_install_builds_checkm2_and_fetches_its_database(quickstart_calls):
    """bin_refinement scores with CheckM2 by default, so its env and database are in the set."""
    assert quickstart.main(["--binning-only", "--databases", "small"]) == 0
    assert "checkm2" in quickstart_calls["env"]["modules"]
    assert "checkm2" in quickstart_calls["db"]
    # CheckM1's data is only for --checkm1.
    assert "checkm" not in quickstart_calls["db"]
    assert (
        "checkm2" in quickstart.BINNING_DATABASES and "checkm" not in quickstart.BINNING_DATABASES
    )


def test_a_full_install_builds_every_module_env_and_checkm2(quickstart_calls):
    """A full quickstart names every module env, plus CheckM2 and the learned binners."""
    from metawrap2.config import ALL_MODULE_ENVS

    assert quickstart.main(["--databases", "none"]) == 0
    built = quickstart._unique_envs(quickstart_calls["env"]["modules"])
    assert set(ALL_MODULE_ENVS) <= set(built)
    assert {"metawrap2-checkm2", "metawrap2-learned"} <= set(built)


def test_installing_only_reassemble_bins_also_builds_checkm2(quickstart_calls):
    assert quickstart.main(["--modules", "reassemble_bins", "--databases", "none"]) == 0
    assert "checkm2" in quickstart_calls["env"]["modules"]
