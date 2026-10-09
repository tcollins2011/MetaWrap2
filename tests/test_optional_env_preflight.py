"""A missing opt-in conda env is reported before the work, not at the point of use.

This is the difference between two failure modes rather than between failing and succeeding. Module
envs were always checked up front by env_for, so "metawrap2-bins does not exist" arrives in the
first second. Opt-in envs were looked up straight out of OPTIONAL_ENVS with no check, so the first
sign of trouble was mamba complaining about a missing prefix at whatever point the tool was
reached - after every read had been aligned for a learned binner, and after CheckM1 and
consolidation had both finished for GUNC.
"""

import argparse
import ast
import importlib
import os
import re

import pytest

from metawrap2.config import OPTIONAL_ENVS, Settings
from metawrap2.modules import _common, binning


@pytest.fixture
def env_manager_present(monkeypatch):
    """Satisfy the other precondition for running in an env, so only the env itself is under test.

    Running a tool in a conda env needs two things: the env, and the manager that enters it. Both are
    checked here, and a test about one of them should not depend on the developer's PATH for the
    other - mamba is typically not on PATH in a non-login shell, which is how pytest runs.
    """
    monkeypatch.setattr(_common, "_require_env_manager", lambda env: None)


def test_a_missing_env_names_the_install_command_rather_than_the_prefix(
    monkeypatch, env_manager_present
):
    monkeypatch.setattr(_common, "conda_env_exists", lambda env: False)
    said = []

    def fake_error(message, exit_code=1):
        said.append(message)
        raise SystemExit(exit_code)

    # The metawrap2 logger does not propagate to the root logger, so caplog cannot see it.
    monkeypatch.setattr(_common, "error", fake_error)
    with pytest.raises(SystemExit):
        _common.optional_env_for("comebin", Settings(), "--comebin")
    # The yaml name, not the env name and not the module name: that is what install-env takes.
    assert "install-env learned" in said[0]
    assert "--comebin needs it" in said[0]


def test_an_existing_env_is_returned_unchanged(monkeypatch, env_manager_present):
    monkeypatch.setattr(_common, "conda_env_exists", lambda env: True)
    assert _common.optional_env_for("vamb", Settings(), "--vamb") == OPTIONAL_ENVS["vamb"]


def test_conda_envs_off_skips_the_check_entirely(monkeypatch):
    """With use_conda_envs = false the tool is expected on PATH, so there is no env to miss."""
    monkeypatch.setattr(
        _common, "conda_env_exists", lambda env: pytest.fail("should not have been consulted")
    )
    assert (
        _common.optional_env_for("semibin2", Settings(use_conda_envs=False), "--semibin2") is None
    )


def _option_strings(module):
    """Every flag *module*'s parser accepts, without running the module.

    The parsers are built inside _parse_args, so the only way to the parser object is to interrupt
    the call just before it consumes argv.
    """

    class _Caught(Exception):
        def __init__(self, parser):
            self.parser = parser

    def capture(self, argv=None, namespace=None):
        raise _Caught(self)

    original = argparse.ArgumentParser.parse_args
    argparse.ArgumentParser.parse_args = capture
    try:
        module._parse_args([])
    except _Caught as caught:
        parser = caught.parser
    finally:
        argparse.ArgumentParser.parse_args = original
    return {opt for action in parser._actions for opt in action.option_strings}


MODULE_DIR = os.path.dirname(importlib.import_module("metawrap2.modules").__file__)


def _needed_by_phrases():
    """Every (module, needed_by) pair in the tree, read out of the source rather than listed here.

    Discovered rather than hard-coded so a call site added later is covered without anyone having to
    remember to add it. Only literal strings can be checked: binning builds its three as "--%s" %
    from the flag name itself, which cannot go wrong in the way this test is looking for, and
    quant_bins interpolates the chosen quantifier into a literal prefix - the prefix is what matters
    and is recovered from the format string.
    """
    found: list = []
    for name in sorted(os.listdir(MODULE_DIR)):
        if not name.endswith(".py"):
            continue
        with open(os.path.join(MODULE_DIR, name)) as fh:
            tree = ast.parse(fh.read(), filename=name)
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if getattr(func, "id", None) != "optional_env_for" or len(node.args) < 3:
                continue
            phrase = node.args[2]
            if isinstance(phrase, ast.Constant) and isinstance(phrase.value, str):
                text = phrase.value
            elif isinstance(phrase, ast.BinOp) and isinstance(phrase.left, ast.Constant):
                # "--quantifier %s" % args.quantifier - the literal half is the part to check.
                text = phrase.left.value
            else:
                continue
            if re.search(r"--[a-z0-9]", text):
                found.append((name[: -len(".py")], text))
    return found


NEEDED_BY = _needed_by_phrases()


@pytest.mark.parametrize("module_name,phrase", NEEDED_BY)
def test_no_message_tells_the_user_to_pass_a_flag_that_does_not_exist(module_name, phrase):
    """Why needed_by is the caller's argument instead of being derived from the OPTIONAL_ENVS key.

    Deriving "--" + key reads correctly for --semibin2 and --checkm2 and is wrong elsewhere: GUNC
    runs by default and is turned *off* with --skip-gunc, and CoverM is selected as
    --quantifier coverm, not --coverm. So every flag a message mentions is checked against the
    parser that would have to accept it - which is also how the first draft of this test was caught
    asserting classify_bins had no --gtdbtk flag when it has one.
    """
    module = importlib.import_module("metawrap2.modules.%s" % module_name)
    accepted = _option_strings(module)
    mentioned = re.findall(r"--[a-z0-9][a-z0-9-]*", phrase)
    assert mentioned, "%r names no flag, so it cannot tell the user what to pass" % phrase
    for flag in mentioned:
        assert flag in accepted, "%s does not accept %s" % (module_name, flag)


def test_binnings_three_derived_flags_are_all_real():
    """binning builds its needed_by as "--%s" % flag, so there is no literal for the test above to
    pick up - but the keys it derives from still have to match the flags it defines."""
    accepted = _option_strings(binning)
    for flag in ("semibin2", "vamb", "comebin"):
        assert "--%s" % flag in accepted
        assert flag in OPTIONAL_ENVS


@pytest.mark.conda_envs
def test_binning_checks_the_learned_binners_env_before_aligning_anything(
    tmp_path, monkeypatch, env_manager_present
):
    """The reason this is resolved in main() and not where the binner runs.

    A dry run would otherwise reach the SemiBin2 command and return 0; the assertion is that it does
    not get as far as the alignment banner.
    """
    from metawrap2.command import runner

    monkeypatch.setattr(runner, "dry_run", True)
    monkeypatch.setattr(_common, "conda_env_exists", lambda env: env != "metawrap2-learned")
    assembly = tmp_path / "asm.fa"
    assembly.write_text(">c1\n%s\n" % ("ACGT" * 500))
    reads = []
    for mate in ("1", "2"):
        read = tmp_path / ("r_%s.fastq" % mate)
        read.write_text("@r1\nACGT\n+\n!!!!\n")
        reads.append(str(read))
    with pytest.raises(SystemExit):
        binning.main(
            ["--semibin2", "-a", str(assembly), "-o", str(tmp_path / "out"), "-t", "2"] + reads
        )
    # Nothing was aligned: the work directory the alignment stage creates does not exist yet.
    assert not (tmp_path / "out" / "work_files").exists()
