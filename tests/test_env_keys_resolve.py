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


def test_every_env_named_in_the_tables_has_a_yaml_to_build_it_from():
    """An env nothing can build is an env `install-env` reports as missing and cannot fix."""
    envs_dir = os.path.join(SRC, "envs")
    missing = []
    for name in sorted(set(config.MODULE_ENVS.values()) | set(config.OPTIONAL_ENVS.values())):
        yaml = os.path.join(envs_dir, "%s.yaml" % name[len("metawrap2-") :])
        if not os.path.isfile(yaml):
            missing.append(os.path.basename(yaml))
    assert not missing, "env(s) with no yaml: %s" % ", ".join(missing)


def test_every_yaml_declares_the_env_name_it_is_resolved_by():
    """`install-env` finds a yaml by name, then mamba uses the name *inside* it. They must agree."""
    envs_dir = os.path.join(SRC, "envs")
    wrong = []
    for name in sorted(os.listdir(envs_dir)):
        if not name.endswith(".yaml"):
            continue
        with open(os.path.join(envs_dir, name)) as fh:
            declared = next(
                (ln.split(":", 1)[1].strip() for ln in fh if ln.startswith("name:")), None
            )
        expected = "metawrap2-%s" % name[: -len(".yaml")]
        if declared != expected:
            wrong.append("%s declares %r, expected %r" % (name, declared, expected))
    assert not wrong, "\n  ".join(wrong)


def test_every_optional_env_has_tools_for_doctor_to_probe():
    """Without this, doctor reports a freshly built env as "0/0 tools OK" - a pass that probed nothing."""
    for flag in config.OPTIONAL_ENVS:
        assert config.OPTIONAL_ENV_TOOLS.get(flag), "no tools listed for the %r env" % flag


def test_all_module_envs_is_the_deduplicated_set_of_module_envs():
    """Several modules share an env now, so anything iterating MODULE_ENVS to build would repeat work."""
    assert config.ALL_MODULE_ENVS == sorted(set(config.MODULE_ENVS.values()))
    assert len(config.ALL_MODULE_ENVS) < len(config.MODULE_ENVS)


def test_no_source_file_spells_out_a_conda_env_name_that_no_longer_exists():
    """A hardcoded env name is a line that breaks silently when the envs are reorganised.

    `install-db checkm` pointed CheckM1 at its data by running `checkm data setRoot` inside
    "metawrap2-bin_refinement". Consolidating the environments replaced that env, and because the
    failure path there is only a printed NOTE, the step stopped working without saying so - CheckM then
    died in its tree stage with an EOFError from a worker process, which names nothing useful.
    """
    import re

    # Strings that share the prefix but are not conda environments. Listed rather than pattern-matched
    # so that adding one is a deliberate act: the whole value of this test is that an unfamiliar
    # "metawrap2-something" has to be accounted for.
    not_envs = {
        "metawrap2-disk-monitor",  # the disk-watcher thread's name
        "metawrap2-install-env-",  # prefix of install-env's temporary log directory
    }
    known = set(config.MODULE_ENVS.values()) | set(config.OPTIONAL_ENVS.values()) | not_envs
    pattern = re.compile(r"[\"']metawrap2-[a-z0-9_.-]+[\"']")
    offenders = []
    for path in _sources():
        with open(path) as fh:
            for number, line in enumerate(fh, 1):
                if line.lstrip().startswith("#"):
                    continue
                for hit in pattern.findall(line):
                    name = hit.strip("\"'")
                    if name not in known:
                        offenders.append("%s:%d %s" % (os.path.relpath(path, SRC), number, name))
    assert not offenders, "env name(s) that do not exist:\n  " + "\n  ".join(offenders)
