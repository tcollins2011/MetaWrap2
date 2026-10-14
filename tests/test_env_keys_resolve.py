"""Every environment key the source looks up must exist.

The env tables were consolidated from 22 entries to 9, and one caller kept indexing a key that had
gone - `assess.gunc_settings` asking OPTIONAL_ENVS for "gunc". It raised only on a code path that runs
partway through a benchmark, so a chain failed after the expensive part had already succeeded. This
sweeps the source for literal lookups instead of waiting for one to be exercised.
"""

from __future__ import annotations

import os
import re

from metawrap2 import config

SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src", "metawrap2")

#: ``TABLE["literal"]`` where the literal is a plain string - the form that fails at runtime rather
#: than at import, and the form a rename silently breaks.
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
    # Zero is the right answer today - the one such lookup that existed has been replaced by a
    # resolution through MODULE_ENVS - so the count is reported rather than asserted. What is asserted
    # instead is that the pattern still recognises the shape it is looking for, since a regex that
    # matches nothing would pass this test forever without reading anything.
    print("literal env lookups checked: %d" % checked)
    assert LOOKUP.findall('x = _config.OPTIONAL_ENVS["gunc"]') == [("OPTIONAL_ENVS", "gunc")]
    assert LOOKUP.findall('MODULE_TOOLS["binning"]') == [("MODULE_TOOLS", "binning")]
