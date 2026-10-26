"""Tests for the sheet `metawrap2 run --example` prints, the documented first run.

The sheet once listed ``--metabat2 --maxbin2`` as the defaults; that is the ``fast`` tier, which
gives about a quarter fewer clean MAGs, and it used ``-c 50`` instead of the reported thresholds.
Assertions use the binning module's ``PRESETS`` and ``DEFAULT_PRESET``, not literals.
"""

from __future__ import annotations

from typing import List

import pytest

from metawrap2.commands import cohort as cohort_cmd
from metawrap2.commands import run as driver
from metawrap2.modules import binning as _binning


@pytest.fixture
def example(tmp_path, capsys):
    """The example sheet as printed by --example, loaded back."""
    assert driver.main(["--example"]) == 0
    path = tmp_path / "samples.toml"
    path.write_text(capsys.readouterr().out)
    return driver.load_sheet(str(path))


def test_the_example_sheet_runs_the_binning_tier_it_says_it_does(example):
    """The example sheet selects the binners of the default preset."""
    assert sorted(driver.selected_binners(example.extra("binning"))) == sorted(
        _binning.PRESETS[_binning.DEFAULT_PRESET]
    )
    # The sheet names the preset, not a list of binners, so the default is defined in one place.
    assert example.extra("binning") == ["--preset", _binning.DEFAULT_PRESET]


def test_the_example_sheet_uses_the_documented_quality_thresholds(example):
    """The example sheet uses the documented -c 70 -x 10 thresholds."""
    assert driver._quality_thresholds(example) == (70, 10)


def test_each_command_s_own_example_sheet_runs_without_a_warning(tmp_path, monkeypatch, capsys):
    """`run` and `cohort` accept their own example sheets without a warning."""
    said: List[str] = []
    monkeypatch.setattr(driver, "execute", lambda *args, **kwargs: 0)
    monkeypatch.setattr(driver, "warning", said.append)
    # `<command> --example > samples.toml` then `<command> samples.toml`, as in the README and
    # installation.md.
    for name, command in (("run", driver.main), ("cohort", cohort_cmd.main)):
        capsys.readouterr()  # discard the previous command's output
        assert command(["--example"]) == 0
        sheet = tmp_path / ("samples_%s.toml" % name)
        sheet.write_text(capsys.readouterr().out)
        assert command([str(sheet), "-o", str(tmp_path / ("study_%s" % name))]) == 0
        assert said == [], (name, said)


def test_the_cohort_example_sheet_asks_for_per_sample_assembly(tmp_path, capsys):
    """`cohort --example` prints a sheet with per-sample assembly."""
    assert cohort_cmd.main(["--example"]) == 0
    path = tmp_path / "samples.toml"
    path.write_text(capsys.readouterr().out)
    study = driver.load_sheet(str(path))
    assert study.coassemble is False
    assert study.coassemble_chosen is True
    # Otherwise the same as the `run` sheet.
    assert study.extra("binning") == ["--preset", _binning.DEFAULT_PRESET]


def test_the_example_sheet_does_not_pin_steps(example):
    """The example sheet does not pin steps, so run and cohort each use their own default."""
    assert example.steps_chosen is False
    assert example.steps == driver.DEFAULT_STEPS
    # In cohort mode the same sheet gets the cohort steps, including the catalogue.
    example.coassemble = False
    driver.resolve_steps(example)
    assert example.steps == driver.COHORT_STEPS
