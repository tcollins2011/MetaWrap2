"""The one help format every command wears, and the argument types that state their own range.

The point of all of this is that a default appears in exactly one place. Two dozen entry points used to
describe their defaults by hand, with the failure mode you would expect: the code changes and the help
string keeps quoting the old number.
"""

from __future__ import annotations

import argparse

import pytest

from metawrap2.cliformat import (
    MetaWrap2HelpFormatter,
    bounded_float,
    bounded_int,
    make_parser,
    percent,
)


def _help(parser: argparse.ArgumentParser) -> str:
    return parser.format_help()


def test_a_default_is_read_off_the_parser_rather_than_written_into_the_help():
    p = make_parser(prog="x")
    p.add_argument("--threads", type=int, default=8, help="how many threads")
    assert "[default: 8]" in _help(p)


def test_a_help_string_that_explains_its_own_default_is_left_alone():
    """ "default minimap2; see ALIGNERS for why" says more than "[default: minimap2]"."""
    p = make_parser(prog="x")
    p.add_argument("--aligner", default="minimap2", help="aligner (default minimap2; see ALIGNERS)")
    text = _help(p)
    assert "[default: minimap2]" not in text
    assert "default minimap2" in text


def test_a_store_true_flag_gets_no_default_annotation():
    """ "[default: False]" on every switch is noise: not passing the flag already means off."""
    p = make_parser(prog="x")
    p.add_argument("--keep-bams", action="store_true", help="keep the BAMs")
    assert "default" not in _help(p).split("--keep-bams")[1].split("\n")[0].lower()


def test_a_required_option_gets_no_default_annotation():
    p = make_parser(prog="x")
    p.add_argument("-o", "--output", required=True, help="output directory")
    assert "[default:" not in _help(p)


def test_an_option_with_no_useful_default_gets_no_annotation():
    p = make_parser(prog="x")
    p.add_argument("--config", help="path to metawrap2.toml")
    p.add_argument("--extra", nargs="*", default=[], help="extras")
    assert "[default:" not in _help(p)


def test_a_positional_gets_no_default_annotation():
    """nargs="*" defaults to None or [], neither of which means anything to a reader."""
    p = make_parser(prog="x")
    p.add_argument("reads", nargs="*", help="read files")
    assert "[default:" not in _help(p)


def test_a_bounded_option_prints_its_range_and_enforces_it():
    p = make_parser(prog="x")
    p.add_argument("-c", "--completeness", type=percent(), default=70, help="minimum completeness")
    text = _help(p)
    assert "[0-100]" in text and "[default: 70]" in text

    assert p.parse_args(["-c", "0"]).completeness == 0.0
    assert p.parse_args(["-c", "100"]).completeness == 100.0
    for bad in ("-5", "101"):
        with pytest.raises(SystemExit):
            p.parse_args(["-c", bad])


def test_a_one_sided_bound_prints_the_side_it_has():
    p = make_parser(prog="x")
    p.add_argument("-l", "--min-len", type=bounded_int(1), default=1000, help="minimum length")
    p.add_argument("--ceiling", type=bounded_float(None, 10), default=5, help="at most ten")
    text = _help(p)
    assert "[>=1]" in text and "[<=10]" in text


def test_the_bound_is_inclusive_at_both_ends():
    check = bounded_int(1, 5)
    assert check("1") == 1 and check("5") == 5
    for bad in ("0", "6"):
        with pytest.raises(argparse.ArgumentTypeError):
            check(bad)


def test_a_value_that_is_not_a_number_is_rejected_by_type_not_by_range():
    check = bounded_float(0, 1)
    with pytest.raises(argparse.ArgumentTypeError, match="not a number"):
        check("banana")


def test_an_int_bound_keeps_the_int_type():
    """The cast is part of the type, so a bounded int does not hand the module a float."""
    assert isinstance(bounded_int(0)("3"), int)
    assert isinstance(bounded_float(0)("3"), float)


def test_choices_are_not_printed_twice():
    """argparse already shows them in the option itself unless a metavar hides them."""
    shown = make_parser(prog="x")
    shown.add_argument("--route", choices=("a", "b"), default="a", help="which route")
    assert _help(shown).count("'a'") + _help(shown).count("{a,b}") >= 1
    assert "[choices:" not in _help(shown)

    hidden = make_parser(prog="x")
    hidden.add_argument(
        "--route", choices=("a", "b"), default="a", metavar="{a,b,more}", help="which route"
    )
    assert "[choices: a, b]" in _help(hidden)


def test_make_parser_uses_the_formatter_but_lets_a_caller_override_it():
    assert make_parser(prog="x").formatter_class is MetaWrap2HelpFormatter
    other = make_parser(prog="x", formatter_class=argparse.RawDescriptionHelpFormatter)
    assert other.formatter_class is argparse.RawDescriptionHelpFormatter


def test_a_percent_sign_in_a_help_string_still_renders():
    """argparse %-formats help text, so the module's own "%%" habit must survive the formatter."""
    p = make_parser(prog="x")
    p.add_argument(
        "-x", "--contamination", type=percent(), default=10, help="maximum %% contamination"
    )
    assert "maximum % contamination" in _help(p)


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
