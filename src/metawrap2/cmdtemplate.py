"""Filling command templates so each value stays one shell word.

:func:`metawrap2.command.run` splits a command string with ``shlex``, so a value inserted with plain
``str.format`` breaks at a space. :func:`fill_command` quotes each value instead. This module imports
nothing from MetaWrap2, so any module can use it; :mod:`metawrap2.modules._common` re-exports it.
"""

from __future__ import annotations

import shlex
import string
from typing import Any, Sequence

__all__ = ["Raw", "fill_command", "quote_value", "quoted_words"]


class Raw(str):
    """A command fragment that :func:`fill_command` inserts unquoted.

    For a value that is already quoted or is meant to be several words, such as a list of files
    from :func:`quoted_words` or a whole command run under ``bash -c``.
    """


def quote_value(value: Any, format_spec: str = "") -> str:
    """*value* as one shell word: shell-quoted, unless it is :class:`Raw` or empty.

    An empty value adds nothing, as with ``str.format``, so an unset optional flag such as
    ``{cuda}`` does not become an empty argument.
    """
    text = format(value, format_spec)
    if not text or isinstance(value, Raw):
        return text
    return shlex.quote(text)


class _QuotingFormatter(string.Formatter):
    def format_field(self, value: Any, format_spec: str) -> str:
        return quote_value(value, format_spec)


_QUOTING_FORMATTER = _QuotingFormatter()


def fill_command(template: str, **values: Any) -> str:
    """``template.format(**values)``, with each value shell-quoted unless it is :class:`Raw`.

    :func:`metawrap2.command.run` splits a command string with ``shlex``, so an unquoted path with
    a space becomes two arguments (FastQC was given ``--outdir`` truncated at the space). Quoting
    leaves ordinary paths and numbers unchanged, and an empty value inserts nothing. The result is
    also safe inside ``bash -c``.
    """
    return _QUOTING_FORMATTER.format(template, **values)


def quoted_words(items: Sequence[Any]) -> Raw:
    """*items* shell-quoted and joined with spaces, for one ``{files}``-style template field."""
    return Raw(" ".join(shlex.quote(str(item)) for item in items))
