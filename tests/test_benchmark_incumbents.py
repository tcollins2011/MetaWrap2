"""The variant marked `default = true` is the one the module actually ships.

`default = true` exists so a proposed default is always read next to the one it would replace. It can
only do that if it names the right one, and nothing was checking: refinement_completeness carried the
marker on c50_x10 for as long as it existed, while bin_refinement's DEFAULT_MIN_COMPLETION has always
been 70. `benchmark list` printed "-c 50 -x 10 <- current default" and a results table would have
shown the incumbent on a row no user's run produces.

The check is to parse the incumbent's own flags with the module's own parser and compare the result
against parsing nothing: the variant that ships is exactly the one whose flags are a no-op. That
makes the assertion independent of how a default is spelled, which matters because a default can live
in a constant, in `add_argument(default=...)`, or in a settings lookup.
"""

import argparse
import importlib

import pytest

from metawrap2.benchmark.registry import load_configurations

# Comparisons whose incumbent cannot be checked this way, with the reason. Kept explicit so a new
# comparison is checked by default and skipping one is a decision someone wrote down.
UNCHECKABLE = {
    # The binner set has no shipped default expressible as a no-op: `metawrap2 binning` requires the
    # user to name at least one binner, so every variant's flags are load-bearing. What ships is a
    # recommendation in the docs, not a parser default.
    "binner_set": "binning has no default binner - at least one flag is mandatory",
}


#: Modules whose default is resolved in code rather than by the parser, with the resolution written
#: out so the comparison is against what the module will actually do. assembly's
#: ``do_megahit = args.megahit or not args.metaspades`` means --megahit parses to a non-default
#: namespace and is still exactly the shipped behaviour; comparing raw flags would call that a bug.
#: Each entry has to mirror its module, so it is kept next to the line it mirrors by naming it.
def _assembly_effective(ns):
    """assembly.py: do_metaspades = args.metaspades; do_megahit = args.megahit or not metaspades."""
    out = {k: v for k, v in ns.items() if k not in ("metaspades", "megahit")}
    out["do_metaspades"] = ns["metaspades"]
    out["do_megahit"] = ns["megahit"] or not ns["metaspades"]
    return out


RESOLVERS = {"assembly": _assembly_effective}


def _parsed(module, argv):
    """What *module*'s parser makes of *argv*, without running the module.

    The parsers are built inside _parse_args, so the only way to the namespace is to let it parse and
    stop the module before it acts on the result.
    """

    class _Caught(Exception):
        def __init__(self, namespace):
            self.namespace = namespace

    original = argparse.ArgumentParser.parse_args

    def capture(self, args=None, namespace=None):
        raise _Caught(original(self, args, namespace))

    argparse.ArgumentParser.parse_args = capture
    try:
        module._parse_args(argv)
    except _Caught as caught:
        return vars(caught.namespace)
    finally:
        argparse.ArgumentParser.parse_args = original
    raise AssertionError("%s._parse_args did not reach parse_args" % module.__name__)


def _parser(module):
    """*module*'s parser, fully built, before it has consumed any argv."""

    class _Caught(Exception):
        def __init__(self, parser):
            self.parser = parser

    original = argparse.ArgumentParser.parse_args
    argparse.ArgumentParser.parse_args = lambda self, a=None, n=None: (_ for _ in ()).throw(
        _Caught(self)
    )
    try:
        module._parse_args([])
    except _Caught as caught:
        return caught.parser
    finally:
        argparse.ArgumentParser.parse_args = original
    raise AssertionError("%s._parse_args did not reach parse_args" % module.__name__)


def _required_args(module):
    """The minimum argv *module*'s parser accepts, derived from the parser rather than listed.

    Required flags say nothing about defaults - they are what the module cannot run without - but
    every comparison's flags have to be parsed alongside them or argparse refuses before any default
    is visible. Reading them off the parser means a module that gains a required flag does not
    silently turn this test into a SystemExit.
    """
    argv = []
    for action in _parser(module)._actions:
        positional = not action.option_strings
        if not (action.required or (positional and action.nargs not in ("?", "*"))):
            continue
        if positional and action.nargs in (None, 1):
            argv.append("placeholder")
        elif positional:
            argv.append("placeholder")  # nargs="+": one value is enough to satisfy it
        else:
            argv.append(action.option_strings[0])
            if action.nargs != 0:
                argv.append("placeholder")
    return argv


def _effective(module_name, namespace):
    """What the module will do, which is the namespace unless its default is resolved in code."""
    resolver = RESOLVERS.get(module_name)
    return namespace if resolver is None else resolver(namespace)


COMPARISONS = sorted(load_configurations().items())


@pytest.mark.parametrize("name,comparison", COMPARISONS, ids=[n for n, _ in COMPARISONS])
def test_the_incumbent_variants_flags_are_a_no_op(name, comparison):
    if name in UNCHECKABLE:
        pytest.skip(UNCHECKABLE[name])
    incumbent = comparison.incumbent
    assert incumbent is not None, "%s declares no incumbent" % name
    module = importlib.import_module("metawrap2.modules.%s" % comparison.module)
    base = _required_args(module)
    shipped = _effective(comparison.module, _parsed(module, base))
    # common_args deliberately are not defaults - they switch off work that is identical across
    # variants so it cannot bury the difference - so only the variant's own flags are checked.
    claimed = _effective(comparison.module, _parsed(module, base + list(incumbent.args)))
    differing = {k: (shipped[k], claimed[k]) for k in shipped if shipped[k] != claimed[k]}
    assert not differing, "%s: variant %r is marked default but changes %s" % (
        name,
        incumbent.name,
        ", ".join("%s %r->%r" % (k, a, b) for k, (a, b) in sorted(differing.items())),
    )


@pytest.mark.parametrize("name,comparison", COMPARISONS, ids=[n for n, _ in COMPARISONS])
def test_every_other_variant_actually_changes_something(name, comparison):
    """The converse, and the cheaper half of the same guarantee.

    A variant whose flags parse to the shipped defaults is not a challenger, it is a duplicate of the
    incumbent under another name - and it would be recorded as its own row, so the same configuration
    would appear twice in a results table with two different costs.
    """
    if name in UNCHECKABLE:
        pytest.skip(UNCHECKABLE[name])
    module = importlib.import_module("metawrap2.modules.%s" % comparison.module)
    base = _required_args(module)
    shipped = _effective(comparison.module, _parsed(module, base))
    for variant in comparison.variants:
        if variant.default:
            continue
        got = _effective(comparison.module, _parsed(module, base + list(variant.args)))
        assert got != shipped, "%s: variant %r parses identically to the incumbent" % (
            name,
            variant.name,
        )
