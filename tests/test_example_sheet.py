"""Tests for the sheet `metawrap2 run --example` prints, which is the documented first run.

That sheet named ``--metabat2 --maxbin2`` and called them the defaults. They are the ``fast``
tier, so a user starting where the README tells them to start got about a quarter fewer clean MAGs
than the table they had just read promised, at thresholds (``-c 50``) that were not the ones any
reported count uses. The assertions here are made against the binning module's own ``PRESETS`` and
``DEFAULT_PRESET`` rather than against literals, so re-pointing a tier there fails here instead of
silently making the shipped sheet wrong for a third time.
"""

from __future__ import annotations

from typing import List

import pytest

from metawrap2.commands import cohort as cohort_cmd
from metawrap2.commands import run as driver
from metawrap2.modules import binning as _binning


@pytest.fixture
def example(tmp_path, capsys):
    """The example sheet exactly as a user gets it: printed by --example, then loaded back."""
    assert driver.main(["--example"]) == 0
    path = tmp_path / "samples.toml"
    path.write_text(capsys.readouterr().out)
    return driver.load_sheet(str(path))


def test_the_example_sheet_runs_the_binning_tier_it_says_it_does(example):
    """What regresses: the sheet promising one tier and buying another's MAG count."""
    assert sorted(driver.selected_binners(example.extra("binning"))) == sorted(
        _binning.PRESETS[_binning.DEFAULT_PRESET]
    )
    # Named as a tier, not as a list of binners: a list is a third place for "the default" to be
    # written down, and the other two already disagreed with it.
    assert example.extra("binning") == ["--preset", _binning.DEFAULT_PRESET]


def test_the_example_sheet_uses_the_documented_quality_thresholds(example):
    """What regresses: a sheet at -c 50 under docs that report everything at -c 70."""
    assert driver._quality_thresholds(example) == (70, 10)


def test_each_command_s_own_example_sheet_runs_without_a_warning(tmp_path, monkeypatch, capsys):
    """What regresses: the documented first command opening with a warning about its own sheet."""
    said: List[str] = []
    monkeypatch.setattr(driver, "execute", lambda *args, **kwargs: 0)
    monkeypatch.setattr(driver, "warning", said.append)
    # `<command> --example > samples.toml` then `<command> samples.toml`, which is the on-ramp
    # both the README and installation.md give, for each of the two commands that offer it.
    for name, command in (("run", driver.main), ("cohort", cohort_cmd.main)):
        capsys.readouterr()  # the previous run's console output is not part of the next sheet
        assert command(["--example"]) == 0
        sheet = tmp_path / ("samples_%s.toml" % name)
        sheet.write_text(capsys.readouterr().out)
        assert command([str(sheet), "-o", str(tmp_path / ("study_%s" % name))]) == 0
        assert said == [], (name, said)


def test_the_cohort_example_sheet_asks_for_per_sample_assembly(tmp_path, capsys):
    """What regresses: `cohort --example` printing a sheet `cohort` then has to override."""
    assert cohort_cmd.main(["--example"]) == 0
    path = tmp_path / "samples.toml"
    path.write_text(capsys.readouterr().out)
    study = driver.load_sheet(str(path))
    assert study.coassemble is False
    assert study.coassemble_chosen is True
    # Still the same sheet otherwise: the mode is the only thing that differs between the two.
    assert study.extra("binning") == ["--preset", _binning.DEFAULT_PRESET]


def test_the_example_sheet_does_not_pin_steps(example):
    """What regresses: one pinned list silently deciding the plan for both run and cohort."""
    assert example.steps_chosen is False
    assert example.steps == driver.DEFAULT_STEPS
    # The point of not pinning: the same file handed to `metawrap2 cohort` gets the stages that
    # only exist in that mode, instead of the co-assembly's step list with no catalogue at the end.
    example.coassemble = False
    driver.resolve_steps(example)
    assert example.steps == driver.COHORT_STEPS
