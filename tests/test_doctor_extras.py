"""The readiness table's row count, and the dependency it never looked at.

From a fresh install. `quickstart` finished with ``22 OK, 0 missing, 0 broken (of 22 modules)`` and
the next command said:

    ***** differential abundance needs scikit-bio >= 0.7.1: pip install 'metawrap2[stats]'

`metawrap2 -h` lists 24 modules. The table covered the 22 that own a conda env, because `doctor`
probed conda envs and databases only - so the two modules that run in MetaWrap2's own Python were
*absent* from it, and the one of those with an unmet dependency was absent from a report that said
everything was ready.

So the row count is pinned against the CLI's own module table (a hardcoded list is how the two came
apart), and a module gated on a declared Python extra is probed like any other dependency: OK when
it is satisfied, MISSING with the install command when it is not, BROKEN when it is installed and
does not import.

Nothing here probes a conda env or imports an extra: ``_check_module`` and the three functions that
answer "what is installed" are replaced, so the result does not depend on the machine.
"""

from __future__ import annotations

from typing import List, Optional

import pytest

from metawrap2 import cli
from metawrap2.commands import doctor, install_env

REQUIREMENT = "scikit-bio>=0.7.1"


def _healthy_env_row(module: str) -> dict:
    return {
        "module": module,
        "env": "metawrap2-" + module,
        "status": doctor.OK,
        "detail": "3/3 tools OK",
        "tools": [],
    }


@pytest.fixture
def table(monkeypatch):
    """Answer every conda probe with a healthy env, so only the extras rows are under test."""
    monkeypatch.setattr(doctor, "_check_module", _healthy_env_row)
    monkeypatch.setattr(doctor, "_check_optional_envs", lambda modules: [])
    monkeypatch.setattr(doctor, "_extra_requirements", lambda extra: [REQUIREMENT])
    monkeypatch.setattr(doctor, "_installed_version", lambda dist: "0.7.4")
    monkeypatch.setattr(doctor, "_import_probe", lambda dist: (0, ""))


def _row(out: str, module: str) -> str:
    """The table line for *module*, so a test can assert about its cells rather than the page."""
    lines = [line for line in out.splitlines() if line.strip().startswith(module + " ")]
    assert len(lines) == 1, "expected exactly one row for %s, got %r" % (module, lines)
    return lines[0]


def _extras_module() -> str:
    """The module the suite exercises the extras path with."""
    modules = sorted(doctor.MODULE_EXTRAS)
    assert modules, "this file is about modules gated on a Python extra; there has to be one"
    return modules[0]


def test_the_table_lists_every_module_the_cli_lists(table, capsys):
    """One row per module `metawrap2 -h` offers - the count the user compares it against.

    If this regresses, a module drops out of the readiness table silently, which is how
    `diff_abundance` came to be missing from the one report that was supposed to say it could not
    run.
    """
    assert doctor.main(["--tools-only"]) == 0
    out = capsys.readouterr().out
    for module in cli.MODULES:
        assert _row(out, module), module
    assert "(of %d modules)" % len(cli.MODULES) in out


def test_an_unmet_extra_is_missing_with_the_command_that_installs_it(table, monkeypatch, capsys):
    """An extra that is not installed is MISSING, named, and carries its own pip command.

    This is the reported defect: the dependency was real, the message at run time was exactly
    right, and the readiness table said nothing at all.
    """
    monkeypatch.setattr(doctor, "_installed_version", lambda dist: None)
    assert doctor.main(["--tools-only"]) == 1
    row = _row(capsys.readouterr().out, _extras_module())
    assert doctor.PYTHON_EXTRA in row
    assert doctor.MISSING in row
    assert "needs %s" % REQUIREMENT in row
    assert "pip install 'metawrap2[stats]'" in row


def test_a_satisfied_extra_is_ok(table, capsys):
    """A module whose extra is installed is ready, and must not be nagged about.

    The counterpart to the test above: a check that reports every extra as missing would be
    ignored, which is the same outcome as not checking at all.
    """
    assert doctor.main(["--tools-only"]) == 0
    row = _row(capsys.readouterr().out, _extras_module())
    assert doctor.OK in row
    assert "pip install" not in row


def test_an_extra_below_its_floor_is_missing(table, monkeypatch, capsys):
    """An old scikit-bio imports perfectly and still cannot run ANCOM-BC: the floor is the check.

    scikit-bio ships ancombc as of 0.7.1. Probing by import alone would report 0.6 as ready and
    leave the user to find out from the module.
    """
    monkeypatch.setattr(doctor, "_installed_version", lambda dist: "0.6.0")
    assert doctor.main(["--tools-only"]) == 1
    out = capsys.readouterr().out
    assert doctor.MISSING in _row(out, _extras_module())
    # The version that is there, under the row, because "install it" is the wrong instruction for
    # a user who has it installed already.
    assert "scikit-bio 0.6.0 is installed" in out


def test_an_extra_that_does_not_import_is_broken_not_missing(table, monkeypatch, capsys):
    """Installed-but-unimportable is a different state from absent, and needs a different fix.

    An extra built against another numpy's ABI is the usual way this happens. The import runs in a
    subprocess so that it cannot take `doctor` itself down - a command that reports nothing about
    anything is the worst of the three answers.
    """
    monkeypatch.setattr(
        doctor, "_import_probe", lambda dist: (1, "ImportError: numpy.core.multiarray failed")
    )
    assert doctor.main(["--tools-only"]) == 1
    row = _row(capsys.readouterr().out, _extras_module())
    assert doctor.BROKEN in row


def test_fix_does_not_ask_install_env_to_build_an_environment_for_an_extra(
    table, monkeypatch, capsys
):
    """`--fix` must leave extras to pip, and must not report "nothing to fix" over a MISSING row.

    `install-env diff_abundance` would be asked to build an environment that has never existed,
    and an exit code of 0 beside a table that says MISSING is the contradiction this whole change
    is about.
    """
    monkeypatch.setattr(doctor, "_installed_version", lambda dist: None)
    rebuilt: List[Optional[List[str]]] = []
    monkeypatch.setattr(install_env, "main", lambda argv: rebuilt.append(argv) or 0)

    assert doctor.main(["--tools-only", "--fix"]) == 1
    assert rebuilt == [], "nothing here is an environment that can be rebuilt"
    assert "cannot be fixed by rebuilding an environment" in capsys.readouterr().out


# ── `check` had the same blind spot ──────────────────────────────────────────────────────────


def test_check_accepts_a_module_that_has_no_conda_environment():
    """`check diff_abundance` refused a real module as "not a module or an environment".

    Two of the 24 modules run in this package's own interpreter, so `resolve_env_name` returns None
    for them - and the unknown-name guard read that as "no such name". The module whose dependency
    can actually be missing was the one the preflight would not look at.
    """
    from metawrap2.commands import check as check_cmd

    assert check_cmd.main(["quantify_function"]) == 0


def test_check_reports_a_python_extra_rather_than_passing_a_module_it_did_not_probe(capsys):
    """An extra is a dependency like any other; a preflight that ignores it is not a preflight."""
    import importlib.util

    from metawrap2.commands import check as check_cmd

    installed = importlib.util.find_spec("skbio") is not None
    assert check_cmd.main(["diff_abundance"]) == (0 if installed else 1)
    out = capsys.readouterr().out
    assert "Python extras:" in out
    assert "diff_abundance" in out
    if not installed:
        assert "pip install 'metawrap2[stats]'" in out
