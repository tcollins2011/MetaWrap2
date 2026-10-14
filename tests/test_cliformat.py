"""The one help format every command wears, and the argument types that state their own range.

The point of all of this is that a default appears in exactly one place. Two dozen entry points used to
describe their defaults by hand, with the failure mode you would expect: the code changes and the help
string keeps quoting the old number.
"""

from __future__ import annotations

import pytest


def test_every_shipped_command_renders_its_help():
    """A formatter that raises on one parser would only be found when someone typed --help."""
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
