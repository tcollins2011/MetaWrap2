"""Tests for the shared help formatter and the argument types that state their own range.

Defaults are written in one place so help text cannot drift from the code.
"""

from __future__ import annotations

import pytest


def test_every_shipped_command_renders_its_help():
    """Every module and command renders --help and exits 0."""
    import importlib
    import pkgutil

    import metawrap2.commands
    import metawrap2.modules

    checked = 0
    for package in (metawrap2.modules, metawrap2.commands):
        for info in pkgutil.iter_modules(package.__path__):
            if info.name.startswith("_"):
                continue
            module = importlib.import_module("%s.%s" % (package.__name__, info.name))
            parse = getattr(module, "_parse_args", None)
            if parse is None:
                continue
            with pytest.raises(SystemExit) as exit_info:
                parse(["--help"])
            assert exit_info.value.code == 0
            checked += 1
    assert checked >= 10, "expected to have exercised most of the entry points, got %d" % checked
