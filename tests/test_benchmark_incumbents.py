"""The variant marked `default = true` matches what the module ships.

`default = true` marks the incumbent, which tables show beside the proposed alternatives. This was
not checked: refinement_completeness marked c50_x10 while bin_refinement's DEFAULT_MIN_COMPLETION
is 70, so `benchmark list` printed "-c 50 -x 10 <- current default".

The test parses the incumbent's flags with the module's parser and compares the result with
parsing no flags: the shipped variant is the one whose flags are a no-op. This works whether a
default is set in a constant, in `add_argument(default=...)`, or in a settings lookup.
"""

import argparse
import importlib

from metawrap2.benchmark.registry import load_configurations

# Comparisons whose incumbent cannot be checked this way, with the reason. New comparisons are
# checked unless listed here.
UNCHECKABLE = {
    # `metawrap2 binning` requires at least one binner, so no variant's flags are a no-op. The
    # recommended set is in the docs, not a parser default.
    "binner_set": "binning has no default binner - at least one flag is mandatory",
}


#: Modules whose default is resolved in code after parsing, with that resolution reproduced here.
#: assembly's ``do_megahit = args.megahit or not args.metaspades`` means --megahit parses to a
#: non-default namespace but gives the shipped behaviour. Each entry names the line it mirrors.
def _assembly_effective(ns):
    """assembly.py: do_metaspades = args.metaspades; do_megahit = args.megahit or not metaspades."""
    out = {k: v for k, v in ns.items() if k not in ("metaspades", "megahit")}
    out["do_metaspades"] = ns["metaspades"]
    out["do_megahit"] = ns["megahit"] or not ns["metaspades"]
    return out


RESOLVERS = {"assembly": _assembly_effective}


def _parsed(module, argv):
    """What *module*'s parser makes of *argv*, without running the module.

    The parsers are built inside _parse_args, so the module is allowed to parse and is stopped
    before it acts on the result.
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

    Every comparison's flags must be parsed together with the required ones, or argparse exits
    before any default is visible. Reading them from the parser keeps the test working when a
    module gains a required flag.
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


def test_the_incumbent_variants_flags_are_a_no_op():
    """The incumbent parses to the shipped defaults, and every other variant does not.

    A non-incumbent variant whose flags parse to the shipped defaults would duplicate the
    incumbent under another name and appear twice in a results table.
    """
    for name, comparison in sorted(load_configurations().items()):
        if name in UNCHECKABLE:
            continue
        incumbent = comparison.incumbent
        assert incumbent is not None, "%s declares no incumbent" % name
        module = importlib.import_module("metawrap2.modules.%s" % comparison.module)
        base = _required_args(module)
        shipped = _effective(comparison.module, _parsed(module, base))
        # common_args are not defaults (they disable work that is identical across variants), so
        # only the variant's own flags are checked.
        claimed = _effective(comparison.module, _parsed(module, base + list(incumbent.args)))
        differing = {k: (shipped[k], claimed[k]) for k in shipped if shipped[k] != claimed[k]}
        assert not differing, "%s: variant %r is marked default but changes %s" % (
            name,
            incumbent.name,
            ", ".join("%s %r->%r" % (k, a, b) for k, (a, b) in sorted(differing.items())),
        )
        for variant in comparison.variants:
            if variant.default:
                continue
            if variant.is_external:
                # An external arm has no MetaWrap2 flags to parse. The registry already rejects an
                # external arm that declares args or is the default.
                continue
            got = _effective(comparison.module, _parsed(module, base + list(variant.args)))
            assert got != shipped, "%s: variant %r parses identically to the incumbent" % (
                name,
                variant.name,
            )
