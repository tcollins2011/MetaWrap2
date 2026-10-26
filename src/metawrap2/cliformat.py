"""Shared help formatter for ``metawrap2`` commands, and argument types with declared ranges.

`MetaWrap2HelpFormatter` reads each option's default from the parser, so help text cannot go out of
date. `bounded_int` and `bounded_float` declare a range once; it is both checked and printed.

A help string that already mentions its default is left unchanged, so authors can explain a default
in their own words when there is a reason to.
"""

from __future__ import annotations

import argparse
from typing import Callable, Optional, Sequence


class _Bounded:
    """An argparse ``type`` that knows its own limits, so the help can print them.

    Out-of-range values are rejected at parse time instead of failing later in a downstream tool.
    """

    def __init__(
        self,
        cast: Callable[[str], object],
        low: Optional[float],
        high: Optional[float],
        name: str,
    ) -> None:
        self._cast = cast
        self.low = low
        self.high = high
        self.__name__ = name

    def __call__(self, text: str):
        try:
            value = self._cast(text)
        except (TypeError, ValueError):
            raise argparse.ArgumentTypeError("%r is not a %s" % (text, self.__name__))
        number = float(value)  # type: ignore[arg-type]
        if self.low is not None and number < self.low:
            raise argparse.ArgumentTypeError(
                "%s is below the minimum of %s" % (text, _num(self.low))
            )
        if self.high is not None and number > self.high:
            raise argparse.ArgumentTypeError(
                "%s is above the maximum of %s" % (text, _num(self.high))
            )
        return value

    def describe(self) -> str:
        if self.low is not None and self.high is not None:
            return "%s-%s" % (_num(self.low), _num(self.high))
        if self.low is not None:
            return ">=%s" % _num(self.low)
        assert self.high is not None  # one bound at least, or the type would not be bounded
        return "<=%s" % _num(self.high)


def _num(value: float) -> str:
    return "%g" % value


def bounded_int(low: Optional[float] = None, high: Optional[float] = None) -> _Bounded:
    """An ``int`` option restricted to *low*..*high* inclusive, either end optional."""
    return _Bounded(int, low, high, "integer")


def bounded_float(low: Optional[float] = None, high: Optional[float] = None) -> _Bounded:
    """A ``float`` option restricted to *low*..*high* inclusive, either end optional."""
    return _Bounded(float, low, high, "number")


#: A percentage, as CheckM and every completeness/contamination threshold in this tool report one.
def percent() -> _Bounded:
    return _Bounded(float, 0.0, 100.0, "percentage")


class MetaWrap2HelpFormatter(argparse.HelpFormatter):
    """Appends the accepted range and the default to each option's help, when neither is stated.

    Up to three annotations, in this order:

    * ``[0-100]`` - the range, from a `bounded_int` / `bounded_float` type.
    * ``[choices: a, b, c]`` - only when argparse is not already showing them in the option itself,
      which it does whenever the option has no explicit ``metavar``.
    * ``[default: x]`` - skipped when the option is required, when there is no useful default
      (``None``, ``False``, an empty list), and when the help already contains the word "default".

    Nothing is appended to ``store_true`` / ``store_false`` flags.
    """

    def _get_help_string(self, action: argparse.Action) -> str:
        help_text = action.help or ""
        if "%(default)" in help_text:
            return help_text
        parts = []

        if isinstance(action.type, _Bounded):
            parts.append("[%s]" % action.type.describe())

        # argparse shows choices in the option string unless a metavar replaces them.
        if action.choices and action.metavar:
            parts.append("[choices: %s]" % ", ".join(str(c) for c in action.choices))

        if self._should_show_default(action, help_text):
            parts.append("[default: %s]" % self._format_default(action.default))

        if not parts:
            return help_text
        return (help_text + " " if help_text else "") + " ".join(parts)

    @staticmethod
    def _should_show_default(action: argparse.Action, help_text: str) -> bool:
        if action.required or "default" in help_text.lower():
            return False
        if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
            return False
        if isinstance(action, (argparse._HelpAction, argparse._VersionAction)):
            return False
        # Positionals have no `required` flag set and no meaningful default, so skip them here.
        if not action.option_strings:
            return False
        return action.default not in (None, False, [], (), "")

    @staticmethod
    def _format_default(default: object) -> str:
        if isinstance(default, (list, tuple)):
            return " ".join(str(item) for item in default)
        return str(default)


def make_parser(*args, **kwargs) -> argparse.ArgumentParser:
    """`argparse.ArgumentParser` with this tool's help format. Same signature otherwise."""
    kwargs.setdefault("formatter_class", MetaWrap2HelpFormatter)
    return argparse.ArgumentParser(*args, **kwargs)


def add_choice_argument(
    parser: argparse.ArgumentParser,
    *flags: str,
    choices: Sequence[str],
    default: str,
    help: str,
    **kwargs,
) -> None:
    """Add a ``--flag {a,b,c}`` option, so all such options are declared the same way."""
    parser.add_argument(*flags, choices=list(choices), default=default, help=help, **kwargs)
