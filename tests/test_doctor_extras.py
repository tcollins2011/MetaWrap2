"""`doctor` readiness table: one row per CLI module, and Python-extra dependencies.

Guards against a case where `quickstart` reported ``22 OK, 0 missing, 0 broken (of 22 modules)``
and the next command failed with:

    ***** differential abundance needs scikit-bio >= 0.7.1: pip install 'metawrap2[stats]'

`metawrap2 -h` lists 24 modules; `doctor` only probed conda envs and databases, so the two modules
that run in MetaWrap2's own Python were missing from the table.

Checked: the row count matches the CLI's module table, and a module gated on a Python extra is OK
when the extra is satisfied, MISSING (with the install command) when not, and BROKEN when it is
installed but does not import.

``_check_module`` and the three installed-package probes are stubbed, so no conda env is probed and
no extra is imported.
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
    """Stub every conda probe as healthy, so only the extras rows are tested."""
    monkeypatch.setattr(doctor, "_check_module", _healthy_env_row)
    monkeypatch.setattr(doctor, "_check_optional_envs", lambda modules: [])
    monkeypatch.setattr(doctor, "_extra_requirements", lambda extra: [REQUIREMENT])
    monkeypatch.setattr(doctor, "_installed_version", lambda dist: "0.7.4")
    monkeypatch.setattr(doctor, "_import_probe", lambda dist: (0, ""))


def _row(out: str, module: str) -> str:
    """The table line for *module*."""
    lines = [line for line in out.splitlines() if line.strip().startswith(module + " ")]
    assert len(lines) == 1, "expected exactly one row for %s, got %r" % (module, lines)
    return lines[0]


def _extras_module() -> str:
    """The module the suite exercises the extras path with."""
    modules = sorted(doctor.MODULE_EXTRAS)
    assert modules, "this file is about modules gated on a Python extra; there has to be one"
    return modules[0]


def test_the_table_lists_every_module_the_cli_lists(table, capsys):
    """The table has one row per module that `metawrap2 -h` lists."""
    assert doctor.main(["--tools-only"]) == 0
    out = capsys.readouterr().out
    for module in cli.MODULES:
        assert _row(out, module), module
    assert "(of %d modules)" % len(cli.MODULES) in out


def test_an_unmet_extra_is_missing_with_the_command_that_installs_it(table, monkeypatch, capsys):
    """An extra that is not installed is MISSING, named, and shows its pip command."""
    monkeypatch.setattr(doctor, "_installed_version", lambda dist: None)
    assert doctor.main(["--tools-only"]) == 1
    row = _row(capsys.readouterr().out, _extras_module())
    assert doctor.PYTHON_EXTRA in row
    assert doctor.MISSING in row
    assert "needs %s" % REQUIREMENT in row
    assert "pip install 'metawrap2[stats]'" in row


def test_a_satisfied_extra_is_ok(table, capsys):
    """A module whose extra is installed is OK, with no install hint."""
    assert doctor.main(["--tools-only"]) == 0
    row = _row(capsys.readouterr().out, _extras_module())
    assert doctor.OK in row
    assert "pip install" not in row


def test_an_extra_below_its_floor_is_missing(table, monkeypatch, capsys):
    """An extra below its version floor is MISSING even though it imports.

    scikit-bio has ancombc from 0.7.1, so 0.6 must not be reported as ready.
    """
    monkeypatch.setattr(doctor, "_installed_version", lambda dist: "0.6.0")
    assert doctor.main(["--tools-only"]) == 1
    out = capsys.readouterr().out
    assert doctor.MISSING in _row(out, _extras_module())
    # The installed version is shown, since the fix is an upgrade, not an install.
    assert "scikit-bio 0.6.0 is installed" in out


def test_an_extra_that_does_not_import_is_broken_not_missing(table, monkeypatch, capsys):
    """An installed extra that fails to import is BROKEN, not MISSING.

    Usually caused by a build against another numpy ABI. The import runs in a subprocess so it
    cannot crash `doctor`.
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
    """`--fix` leaves extras to pip and exits non-zero while a row is MISSING.

    No `install-env` call is made, since an extra has no environment to build.
    """
    monkeypatch.setattr(doctor, "_installed_version", lambda dist: None)
    rebuilt: List[Optional[List[str]]] = []
    monkeypatch.setattr(install_env, "main", lambda argv: rebuilt.append(argv) or 0)

    assert doctor.main(["--tools-only", "--fix"]) == 1
    assert rebuilt == [], "nothing here is an environment that can be rebuilt"
    assert "cannot be fixed by rebuilding an environment" in capsys.readouterr().out


# ── `check` ──────────────────────────────────────────────────────────────────────────────────


def test_check_accepts_a_module_that_has_no_conda_environment():
    """`check` accepts a module without a conda env.

    Two of the 24 modules run in this package's interpreter, so `resolve_env_name` returns None for
    them; the unknown-name guard used to reject them as "not a module or an environment".
    """
    from metawrap2.commands import check as check_cmd

    assert check_cmd.main(["quantify_function"]) == 0


def test_check_reports_a_python_extra_rather_than_passing_a_module_it_did_not_probe(capsys):
    """`check` probes a module's Python extra and reports the pip command when it is missing."""
    import importlib.util

    from metawrap2.commands import check as check_cmd

    installed = importlib.util.find_spec("skbio") is not None
    assert check_cmd.main(["diff_abundance"]) == (0 if installed else 1)
    out = capsys.readouterr().out
    assert "Python extras:" in out
    assert "diff_abundance" in out
    if not installed:
        assert "pip install 'metawrap2[stats]'" in out
