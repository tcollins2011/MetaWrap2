"""Tests for the things a first-time or migrating user hits: quickstart, config migration,
retired names, and the disk-space check's severity.

The property that matters most here is the one about databases. They run from 550 MB to ~300 GB
and take hours to days, so `quickstart` must never start one without being told to - and in
particular must never be able to do so on a non-interactive run, where nobody is there to say no.
"""

from __future__ import annotations

import os

import pytest

from metawrap2.commands import config_cmd, quickstart

# --- quickstart: databases are opt-in ------------------------------------------------------


@pytest.fixture
def calls(monkeypatch):
    """Record what quickstart would invoke, without invoking anything."""
    seen = {"env": None, "db": None, "doctor": None}

    monkeypatch.setattr(quickstart, "_package_manager", lambda: "mamba")
    monkeypatch.setattr(
        quickstart,
        "_run_install_env",
        lambda modules, threads, from_lock: seen.__setitem__("env", (list(modules), threads)) or 0,
    )
    monkeypatch.setattr(
        quickstart,
        "_run_install_db",
        lambda names, directory, small, jobs, config: seen.__setitem__(
            "db", (list(names), small, directory)
        )
        or 0,
    )
    monkeypatch.setattr(
        quickstart,
        "_run_doctor",
        lambda modules, config, databases: seen.__setitem__("doctor", databases) or 0,
    )
    return seen


def test_no_database_is_downloaded_by_default(calls, monkeypatch):
    """The whole point: an unattended quickstart installs environments and nothing else."""
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    assert quickstart.main(["--binning-only"]) == 0
    assert calls["env"] is not None, "environments should still be created"
    assert calls["db"] is None, "no database download may be started without being asked"


def test_a_non_interactive_run_never_prompts(monkeypatch, capsys):
    """A batch script or CI job must not be able to trigger a 300 GB download."""
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    monkeypatch.setattr(
        "builtins.input", lambda *a: pytest.fail("must not prompt without a terminal")
    )
    assert quickstart._ask_about_databases(["checkm"], 3) == "none"


def test_an_interactive_run_defaults_to_no(monkeypatch):
    """Pressing Enter must mean 'no', not 'yes, download hundreds of gigabytes'."""
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    monkeypatch.setattr("sys.stdout.isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda *a: "")
    assert quickstart._ask_about_databases(["checkm"], 3) == "none"


@pytest.mark.parametrize("answer,expected", [("1", "none"), ("2", "small"), ("3", "full")])
def test_the_prompt_maps_answers(monkeypatch, answer, expected):
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    monkeypatch.setattr("sys.stdout.isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda *a: answer)
    assert quickstart._ask_about_databases(["checkm"], 3) == expected


def test_an_unparseable_answer_means_no(monkeypatch):
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    monkeypatch.setattr("sys.stdout.isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda *a: "yes please")
    assert quickstart._ask_about_databases(["checkm"], 3) == "none"


def test_ctrl_c_at_the_prompt_means_no(monkeypatch):
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    monkeypatch.setattr("sys.stdout.isatty", lambda: True)

    def interrupt(*_a):
        raise KeyboardInterrupt

    monkeypatch.setattr("builtins.input", interrupt)
    assert quickstart._ask_about_databases(["checkm"], 3) == "none"


def test_databases_small_installs_the_capped_variants(calls, monkeypatch):
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    assert quickstart.main(["--binning-only", "--databases", "small"]) == 0
    names, small, _dir = calls["db"]
    assert small is True
    assert names == list(quickstart.BINNING_DATABASES)


def test_databases_full_installs_the_production_ones(calls, monkeypatch):
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    assert quickstart.main(["--databases", "full"]) == 0
    names, small, _dir = calls["db"]
    assert small is False
    assert names == ["--all"]  # a full install needs every database, not just the binning ones


def test_binning_only_installs_the_genome_recovery_modules(calls, monkeypatch):
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    quickstart.main(["--binning-only", "--databases", "none"])
    modules, _threads = calls["env"]
    assert modules == list(quickstart.BINNING_PATH)
    assert "kraken2" not in modules  # its database is one of the enormous ones


def test_doctor_is_not_asked_about_databases_that_were_skipped(calls, monkeypatch):
    """Otherwise a successful env-only install would report itself as failed."""
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    quickstart.main(["--binning-only", "--databases", "none"])
    assert calls["doctor"] is False
    quickstart.main(["--binning-only", "--databases", "small"])
    assert calls["doctor"] is True


def test_an_unknown_module_is_rejected(monkeypatch, capsys):
    monkeypatch.setattr(quickstart, "_package_manager", lambda: "mamba")
    assert quickstart.main(["--modules", "not_a_module"]) == 1
    assert "Unknown module" in capsys.readouterr().out


def test_missing_conda_is_explained_not_tracebacked(monkeypatch, capsys):
    monkeypatch.setattr(quickstart, "_package_manager", lambda: None)
    assert quickstart.main([]) == 1
    out = capsys.readouterr().out
    assert "conda" in out and "use_conda_envs = false" in out


def test_a_failed_step_names_the_command_to_resume_with(calls, monkeypatch, capsys):
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    monkeypatch.setattr(quickstart, "_run_install_env", lambda *a, **k: 1)
    assert quickstart.main(["--binning-only"]) == 1
    assert "metawrap2 install-env" in capsys.readouterr().out


# --- config migration ----------------------------------------------------------------------

SHELL_CONFIG = """\
# Paths to custom pipelines and scripts
mw_path=$(which metawrap)
SOFT=${bin_path}/metawrap-scripts

KRAKEN_DB=/scratch/old/MY_KRAKEN_DB
KRAKEN2_DB=/scratch/db/kraken2
BMTAGGER_DB=/scratch/db/bmtagger
BLASTDB=/scratch/db/NCBI_nt
TAXDUMP=/scratch/db/NCBI_tax
SOMETHING_ELSE=/scratch/db/whatever
"""


def test_a_shell_config_is_translated():
    translated, retired, unrecognised = config_cmd.parse_shell_config(SHELL_CONFIG)
    assert translated == {
        "KRAKEN2_DB": "/scratch/db/kraken2",
        "BMTAGGER_DB": "/scratch/db/bmtagger",
        "BLASTDB": "/scratch/db/NCBI_nt",
        "TAXDUMP": "/scratch/db/NCBI_tax",
    }
    assert retired == ["KRAKEN_DB"]
    assert "SOMETHING_ELSE" in unrecognised


def test_installation_paths_are_not_treated_as_databases():
    """Those variables located the old installation's own scripts, which is resolved from the
    installed package now - carrying them over would write nonsense into the config."""
    translated, _retired, unrecognised = config_cmd.parse_shell_config(SHELL_CONFIG)
    assert "SOFT" not in translated and not any("SOFT" in u for u in unrecognised)
    assert "mw_path" not in translated


def test_a_computed_value_is_reported_rather_than_guessed():
    """Nothing is executed to read this file, so a value built by the shell cannot be resolved."""
    _t, _r, unrecognised = config_cmd.parse_shell_config("MY_DB=$HOME/db\n")
    assert any("MY_DB" in u and "computed" in u for u in unrecognised)


def test_comments_and_blank_lines_are_ignored():
    translated, _r, _u = config_cmd.parse_shell_config(
        "# BLASTDB=/commented/out\n\nTAXDUMP=/real/path   # trailing comment\n"
    )
    assert translated == {"TAXDUMP": "/real/path"}
    assert "BLASTDB" not in translated


def test_quoted_values_are_unquoted():
    translated, _r, _u = config_cmd.parse_shell_config('TAXDUMP="/quoted/path"\n')
    assert translated == {"TAXDUMP": "/quoted/path"}


def test_the_rendered_config_is_valid_and_loads(tmp_path):
    """A migration that produces a config MetaWrap2 then rejects is worse than no migration."""
    from metawrap2.config import check_config, load_settings

    try:
        import tomllib as toml
    except ModuleNotFoundError:  # pragma: no cover
        import tomli as toml

    translated, _r, _u = config_cmd.parse_shell_config(SHELL_CONFIG)
    path = tmp_path / "migrated.toml"
    path.write_text(config_cmd.render_config(translated, threads=16))

    with open(path, "rb") as fh:
        data = toml.load(fh)
    problems, notes = check_config(data, str(path))
    assert problems == [] and notes == []

    settings = load_settings(str(path))
    assert settings.threads == 16
    assert settings.db("TAXDUMP") == "/scratch/db/NCBI_tax"


def test_migrate_writes_a_file_and_refuses_to_clobber(tmp_path, capsys):
    source = tmp_path / "old-config"
    source.write_text(SHELL_CONFIG)
    target = tmp_path / "config.toml"

    argv = ["migrate", str(source), "-o", str(target)]
    assert config_cmd.main(argv) == 0
    assert "KRAKEN2_DB" in target.read_text()
    assert "KRAKEN_DB =" not in target.read_text()  # the retired key is dropped, not carried

    assert config_cmd.main(argv) == 1
    assert "--force" in capsys.readouterr().out
    assert config_cmd.main(argv + ["--force"]) == 0


def test_migrate_says_when_a_path_does_not_exist(tmp_path, capsys):
    source = tmp_path / "old-config"
    source.write_text("TAXDUMP=/definitely/not/here\n")
    assert config_cmd.main(["migrate", str(source)]) == 0
    assert "does not exist here" in capsys.readouterr().out


def test_migrate_on_something_that_is_not_a_config(tmp_path, capsys):
    source = tmp_path / "random.txt"
    source.write_text("this file has nothing to do with anything\n")
    assert config_cmd.main(["migrate", str(source)]) == 1
    assert "config init" in capsys.readouterr().out


def test_migrate_on_a_missing_file(tmp_path, capsys):
    assert config_cmd.main(["migrate", str(tmp_path / "nope")]) == 1
    assert "No such file" in capsys.readouterr().out


# --- a retired key must not block an unrelated run ----------------------------------------


def test_a_retired_database_key_is_a_note_not_a_failure():
    """It is harmless and only a removed module ever read it; refusing to run binning over it
    would be absurd."""
    from metawrap2.config import check_config

    problems, notes = check_config({"databases": {"KRAKEN_DB": "/x", "TAXDUMP": "/y"}}, "cfg.toml")
    assert problems == []
    assert len(notes) == 1 and "KRAKEN_DB" in notes[0] and "delete the line" in notes[0]


def test_a_mistyped_database_key_is_also_only_a_note():
    from metawrap2.config import check_config

    problems, notes = check_config({"databases": {"CHECKM_DBB": "/x"}}, "cfg.toml")
    assert problems == []
    assert any("CHECKM_DB" in n for n in notes)  # the suggestion is still made


def test_a_mistyped_setting_is_still_a_failure():
    """That namespace is small and entirely ours, and a typo there silently does nothing."""
    from metawrap2.config import check_config

    problems, _notes = check_config({"settings": {"use_conda_env": False}}, "cfg.toml")
    assert len(problems) == 1 and "use_conda_envs" in problems[0]


def test_a_config_with_only_a_retired_key_still_loads(tmp_path):
    from metawrap2.config import load_settings

    path = tmp_path / "cfg.toml"
    path.write_text('[databases]\nKRAKEN_DB = "/old/kraken"\n')
    settings = load_settings(str(path))  # must not raise
    assert settings.db("KRAKEN2_DB") == ""


# --- retired CLI names explain themselves --------------------------------------------------


@pytest.mark.parametrize("name", ["kraken", "phylosift"])
def test_a_retired_module_name_explains_itself(name, capsys):
    from metawrap2 import cli

    assert cli.main([name, "-o", "out"]) == 2
    err = capsys.readouterr().err
    assert name in err and ("kraken2" in err or "classify_bins" in err)


@pytest.mark.parametrize("flag", ["--metabat1", "--metabat"])
def test_a_retired_flag_explains_itself(flag, capsys):
    from metawrap2 import cli

    assert cli.main(["binning", flag, "-a", "x.fa", "-o", "o"]) == 2
    assert "--metabat2" in capsys.readouterr().err


def test_a_genuinely_unknown_module_still_prints_the_help(capsys):
    from metawrap2 import cli

    assert cli.main(["not_a_thing"]) == 1
    assert "Unknown module" in capsys.readouterr().out


# --- the disk-space check warns rather than refuses -----------------------------------------


def test_a_tight_disk_warns_but_does_not_stop(tmp_path, monkeypatch, capsys):
    """The estimate is a chosen multiplier, not a measurement; refusing on it is the worse error."""
    from metawrap2 import scratch
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    monkeypatch.setattr(scratch, "free_bytes", lambda path: 1024)

    _common.check_disk_space("assembly", Settings(), [str(reads)], str(tmp_path))
    out = capsys.readouterr().out
    assert "may not have room" in out
    assert "strict_space_check" in out  # and says how to make it a hard stop


def test_strict_space_check_makes_it_stop(tmp_path, monkeypatch):
    from metawrap2 import scratch
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    monkeypatch.setattr(scratch, "free_bytes", lambda path: 1024)

    with pytest.raises(SystemExit):
        _common.check_disk_space(
            "assembly", Settings(strict_space_check=True), [str(reads)], str(tmp_path)
        )


def test_skip_space_check_bypasses_it_entirely(tmp_path, monkeypatch, capsys):
    from metawrap2 import scratch
    from metawrap2.config import Settings
    from metawrap2.modules import _common

    monkeypatch.setattr(
        scratch, "check_space", lambda *a: pytest.fail("the check should not have run")
    )
    _common.check_disk_space("assembly", Settings(skip_space_check=True), [], str(tmp_path))
    assert os.path.isdir(tmp_path)
