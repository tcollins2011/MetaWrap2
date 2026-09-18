"""Structured logging that preserves metaWRAP's on-screen banner style.

The original ``print_comment.py`` drew a 120-char banner box around each message.
We reproduce that output byte-for-byte (so ``comm``/``announcement``/``warning``/
``error`` look exactly as users expect) while optionally teeing every message to a
machine-readable run log for provenance.

The banners map to the same delimiters the shell modules already use:
    comm         -> "-"
    warning      -> "*"
    error        -> "*"   (then exit 1, handled by the caller)
    announcement -> "#"
"""

from __future__ import annotations

import logging
import os
from typing import Optional

from .constants import BANNER_MAX_LINE as MAX_LINE
from .constants import BANNER_WIDTH as WIDTH

_LOGGER = logging.getLogger("metawrap")
_run_log_configured = False


def format_comment(message: str, delim: str) -> str:
    """Return the banner-boxed rendering of *message*, identical to the legacy script.

    Historically this lived in ``print_comment.py`` (Python 2), where ``/2`` was integer
    division. We use ``//`` here so the spacing is identical under Python 3 instead of
    crashing on ``" " * float``.
    """
    out = ["\n" + delim * WIDTH]

    # Legacy behavior: `line` is seeded empty, so the first word keeps a leading space.
    rendered = []
    line = ""
    for word in message.split(" "):
        if (len(line) + 1 + len(word)) > MAX_LINE:
            rendered.append(line)
            line = word
        else:
            line = line + " " + word
    rendered.append(line)

    for line in rendered:
        edge1 = (WIDTH - len(line)) // 2 - 5
        edge2 = WIDTH - edge1 - len(line) - 10
        # guard against pathological long lines so we never emit negative padding
        edge1 = max(edge1, 0)
        edge2 = max(edge2, 0)
        out.append(delim * 5 + " " * edge1 + line + " " * edge2 + delim * 5)

    out.append(delim * WIDTH + "\n")
    return "\n".join(out)


def configure_run_log(path: Optional[str]) -> None:
    """Tee all banner messages to *path* (appended), in addition to stdout."""
    global _run_log_configured
    if not path:
        return
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    handler = logging.FileHandler(path)
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    _LOGGER.addHandler(handler)
    _LOGGER.setLevel(logging.INFO)
    _run_log_configured = True


def _emit(message: str, delim: str, level: int) -> None:
    print(format_comment(message, delim))
    if _run_log_configured:
        _LOGGER.log(level, message)


def comm(message: str) -> None:
    _emit(message, "-", logging.INFO)


def announcement(message: str) -> None:
    _emit(message, "#", logging.INFO)


def warning(message: str) -> None:
    _emit(message, "*", logging.WARNING)


def error(message: str, exit_code: int = 1) -> None:
    """Print an error banner and exit, matching the shell ``error`` helper."""
    _emit(message, "*", logging.ERROR)
    raise SystemExit(exit_code)
