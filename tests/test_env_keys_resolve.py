"""Every environment key the source looks up must exist.

The env tables were consolidated from 22 entries to 9, and `assess.gunc_settings` kept looking up
"gunc" in OPTIONAL_ENVS, which raised partway through a benchmark. This test scans the source for
literal lookups.
"""

from __future__ import annotations

import os
import re

from metawrap2 import config

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src", "metawrap2")

#: ``TABLE["literal"]`` with a plain string literal, which fails at runtime rather than at import.
LOOKUP = re.compile(
    r"\b(OPTIONAL_ENVS|MODULE_ENVS|OPTIONAL_ENV_TOOLS|MODULE_TOOLS)\[\"([^\"]+)\"\]"
)


def _sources():
    for root, _dirs, files in os.walk(SRC):
        for name in files:
            if name.endswith(".py"):
                yield os.path.join(root, name)


def test_every_literal_env_lookup_in_the_source_resolves():
    tables = {
        "OPTIONAL_ENVS": config.OPTIONAL_ENVS,
        "MODULE_ENVS": config.MODULE_ENVS,
        "OPTIONAL_ENV_TOOLS": config.OPTIONAL_ENV_TOOLS,
        "MODULE_TOOLS": config.MODULE_TOOLS,
    }
    missing = []
    checked = 0
    for path in _sources():
        with open(path) as fh:
            text = fh.read()
        for table, key in LOOKUP.findall(text):
            checked += 1
            if key not in tables[table]:
                missing.append("%s: %s[%r]" % (os.path.relpath(path, SRC), table, key))
    assert not missing, "lookups that would raise KeyError at runtime:\n  " + "\n  ".join(missing)
    # There are currently no such lookups (the one that existed now resolves through MODULE_ENVS),
    # so the count is printed, not asserted. The test asserts instead that the pattern still
    # matches the shape it looks for.
    print("literal env lookups checked: %d" % checked)
    assert LOOKUP.findall('x = _config.OPTIONAL_ENVS["gunc"]') == [("OPTIONAL_ENVS", "gunc")]
    assert LOOKUP.findall('MODULE_TOOLS["binning"]') == [("MODULE_TOOLS", "binning")]
