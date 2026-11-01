"""Flags and messages shared between commands, and the starter config's comments."""

from __future__ import annotations

import re
import sys

import pytest

from metawrap2 import gunc
from metawrap2.commands import completion, config_cmd
from metawrap2.commands import run as driver
from metawrap2.modules import bin_refinement, diff_abundance

if sys.version_info >= (3, 11):
    import tomllib as _toml
else:  # pragma: no cover - exercised on <3.11
    import tomli as _toml


# --- GUNC flag names -------------------------------------------------------------------------


def test_bin_refinement_accepts_scores_no_gunc():
    base = ["-o", "out", "-A", "bins"]
    assert bin_refinement._parse_args(base + ["--no-gunc"]).run_gunc is False
    assert bin_refinement._parse_args(base + ["--skip-gunc"]).run_gunc is False
    assert bin_refinement._parse_args(base).run_gunc is True


def test_the_alias_is_not_advertised(capsys):
    with pytest.raises(SystemExit):
        bin_refinement._parse_args(["-h"])
    assert "--no-gunc" not in capsys.readouterr().out


def test_the_gunc_hint_names_the_callers_flag():
    assert gunc.db_hint("--no-gunc").endswith("or pass --no-gunc.")
    assert "--skip-gunc" not in gunc.db_hint("--no-gunc")
    assert gunc.DB_HINT == gunc.db_hint("--skip-gunc")


def test_run_gunc_passes_the_flag_into_its_failure_hint(monkeypatch, tmp_path):
    seen = {}

    def fake_run(_command, **kwargs):
        seen.update(kwargs)
        return 0

    monkeypatch.setattr(gunc, "run", fake_run)
    monkeypatch.setattr(gunc.runner, "dry_run", True)
    gunc.run_gunc(str(tmp_path), threads=1, db="db", env=None, skip_flag="--no-gunc")
    assert seen["hint"] == gunc.db_hint("--no-gunc")


# --- completion ------------------------------------------------------------------------------


def test_completion_prints_bash_and_zsh(capsys):
    assert completion.main(["bash"]) == 0
    assert capsys.readouterr().out.startswith("_metawrap2()")
    assert completion.main(["zsh"]) == 0
    assert capsys.readouterr().out.startswith("#compdef metawrap2")
    assert completion.main([]) == 0
    assert capsys.readouterr().out.startswith("_metawrap2()")


@pytest.mark.parametrize("shell", ["fish", "notashell", "powershell"])
def test_completion_refuses_other_shells(capsys, shell):
    with pytest.raises(SystemExit) as stopped:
        completion.main([shell])
    assert stopped.value.code == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "invalid choice" in captured.err


def test_completion_help_is_help(capsys):
    with pytest.raises(SystemExit) as stopped:
        completion.main(["-h"])
    assert stopped.value.code == 0
    out = capsys.readouterr().out
    assert "usage: metawrap2 completion" in out and "_metawrap2()" not in out


# --- diff_abundance -t -----------------------------------------------------------------------


@pytest.mark.parametrize("flag", ["-t", "--threads"])
def test_diff_abundance_takes_and_drops_threads(flag):
    args = diff_abundance._parse_args(["-o", "out", "-i", "t.tsv", "-g", "g.tsv", flag, "8"])
    assert not hasattr(args, "threads")


def test_diff_abundance_does_not_advertise_threads(capsys):
    with pytest.raises(SystemExit):
        diff_abundance._parse_args(["-h"])
    assert "--threads" not in capsys.readouterr().out


# --- the example sheet -----------------------------------------------------------------------


@pytest.mark.parametrize("cohort", [False, True])
def test_the_example_sheet_says_when_it_is_read(cohort):
    header = driver.example_sheet(cohort=cohort).split("[settings]")[0]
    flat = " ".join(line.lstrip("# ") for line in header.splitlines())
    assert "read once, when the run starts" in flat
    assert "next --resume" in flat
    assert "[limits] in metawrap2.toml is re-read by every step" in flat


# --- config init ------------------------------------------------------------------------------

#: A commented-out setting, with an optional comment after its value.
_KEY_LINE = re.compile(r"^# ?([a-z_]+) = (\S+)\s*(?:#\s*(.*))?$")


def test_config_keys_carry_only_a_short_unit_after_the_value():
    """Prose goes above a key; uncommenting and editing the value leaves nothing contradictory."""
    keys = [m for m in map(_KEY_LINE.match, config_cmd.DEFAULT_CONFIG.splitlines()) if m]
    assert {"max_disk_fraction", "max_gpus", "max_thread_fraction"} <= {m.group(1) for m in keys}
    for match in keys:
        trailing = match.group(3) or ""
        assert len(trailing.split()) <= 3, match.group(0)
        assert not re.search(r"\d", trailing), match.group(0)  # no value restated as a number


def test_an_uncommented_and_edited_limit_still_parses():
    text = config_cmd.DEFAULT_CONFIG.replace("# max_disk_fraction = 0.7", "max_disk_fraction = 0.9")
    assert _toml.loads(text)["limits"]["max_disk_fraction"] == 0.9
