"""MetaWrap2's logging: one :class:`logging.Logger`, two destinations, banner output on screen.

Everything MetaWrap2 says goes through the ``metawrap2`` logger. It carries two handlers:

**the console** (stdout)
    renders banner-boxed messages in the MetaWrap2 style;

**the run log** (``<output>/metawrap2.log``)
    the same messages as plain timestamped lines, plus every line of every external tool's
    stdout and stderr, in order.

``run.stdout`` / ``run.stderr`` are still written separately, so the raw tool streams are
available unmixed.

A Logger is used instead of ``print`` so that levels (tool output is DEBUG, progress INFO,
problems WARNING/ERROR) keep the console readable while the file keeps everything, and other
destinations can be added as handlers.

The banner helpers are the module's public interface and are used everywhere:

    comm("...")          a step-level message      -> INFO,    "-" banner
    announcement("...")  a stage heading           -> INFO,    "#" banner
    warning("...")       something to look at      -> WARNING, "*" banner
    error("...")         fatal; raises SystemExit  -> ERROR,   "*" banner
"""

from __future__ import annotations

import logging
import logging.handlers
import os
import sys
import textwrap
from typing import List, NoReturn, Optional, TextIO

from .constants import BANNER_MAX_LINE as MAX_LINE
from .constants import BANNER_WIDTH as WIDTH

LOGGER_NAME = "metawrap2"

#: Name of the combined log written into a run's output directory.
RUN_LOG_NAME = "metawrap2.log"

#: Size at which the run log rolls over to metawrap2.log.1, .2, ... A normal run writes a few
#: hundred KB, so this only applies to a tool with runaway output, which could otherwise fill the
#: filesystem.
MAX_LOG_BYTES = 256 * 1024 * 1024  # 256 MB

#: How many rolled-over logs to keep.
LOG_BACKUPS = 3

#: A message's banner delimiter travels as a logging "extra" so the console handler can draw
#: the box and the file handler can ignore it.
_DELIM_KEY = "metawrap2_delim"

#: Level used for external tools' output: kept out of the console (the runner already streams
#: it live) but recorded in the run log.
TOOL_LEVEL = logging.DEBUG

_logger = logging.getLogger(LOGGER_NAME)
_logger.setLevel(logging.DEBUG)
_logger.propagate = False  # our handlers only; never the root logger's

_console_handler: Optional[logging.Handler] = None
_file_handler: Optional[logging.Handler] = None


#: How many delimiters frame each side of a content line, and therefore how much of the box
#: width is border rather than text.
_EDGE = 5

#: The space between the two borders, and so the upper limit on the wrap column.
_INNER = WIDTH - 2 * _EDGE

#: The column the text is wrapped to: the style's line length, clamped to what fits.
_WRAP = min(MAX_LINE, _INNER)


def wrap_banner_text(message: str) -> List[str]:
    """Break *message* into the lines that go inside the banner, each at most :data:`_WRAP` wide.

    Explicit newlines are kept as line breaks, and each segment is then wrapped to the text
    column.

    Wrapping happens before padding. Padding first, measured on a segment that still contained
    newlines, misplaced the right border (it broke three of four lines in `quickstart`'s closing
    banner).
    """
    lines: List[str] = []
    for segment in message.split("\n"):
        text = segment.strip()
        if not text:
            # Keep intentional blank lines between paragraphs.
            lines.append("")
            continue
        # break_long_words and break_on_hyphens are off: long tokens are paths, URLs and commands
        # the user copies, which must not be split.
        lines.extend(
            textwrap.wrap(text, width=_WRAP, break_long_words=False, break_on_hyphens=False)
        )
    return lines or [""]


def format_comment(message: str, delim: str) -> str:
    """Return the banner-boxed rendering of *message*.

    Ported from ``print_comment.py`` (Python 2), where ``/2`` was integer division;
    ``str.center`` gives the same spacing (the odd space goes to the right).

    A single word longer than the text column overruns the right border instead of being split,
    because such tokens are paths or commands the user may copy.
    """
    border = delim * _EDGE
    out = ["\n" + delim * WIDTH]
    for line in wrap_banner_text(message):
        out.append(border + line.center(_INNER) + border)
    out.append(delim * WIDTH + "\n")
    return "\n".join(out)


class ConsoleHandler(logging.StreamHandler):
    """A stream handler that resolves ``sys.stdout`` at emit time, not at construction.

    ``logging.StreamHandler(sys.stdout)`` keeps the stream object, so a later replacement of
    ``sys.stdout`` (test capture, redirection) leaves it writing to an old or closed stream.
    Looking it up per record avoids that.
    """

    @property
    def stream(self) -> TextIO:
        return sys.stdout

    @stream.setter
    def stream(self, value: object) -> None:
        # StreamHandler.__init__ assigns to self.stream; ignore it and stay dynamic.
        pass


class BannerFormatter(logging.Formatter):
    """Console formatter: draws the banner box when a record carries a delimiter."""

    def format(self, record: logging.LogRecord) -> str:
        message = record.getMessage()
        delim = getattr(record, _DELIM_KEY, None)
        if delim:
            return format_comment(message, delim)
        return message


class PlainFormatter(logging.Formatter):
    """File formatter: one timestamped line per message, banners flattened.

    Tool output arrives already tagged with its stream (see :func:`log_tool_output`), so the
    level and the tag together tell a reader whether MetaWrap2 or a tool said something.
    """

    def __init__(self) -> None:
        super().__init__(fmt="%(asctime)s %(levelname)-7s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")


def get_logger() -> logging.Logger:
    """The MetaWrap2 logger. Prefer the helpers below for user-facing messages."""
    return _logger


def configure_console() -> None:
    """Attach (once) the console handler that renders banners to stdout."""
    global _console_handler
    if _console_handler is not None:
        return
    handler = ConsoleHandler()
    handler.setFormatter(BannerFormatter())
    handler.setLevel(logging.INFO)  # tool chatter stays out of the console
    _logger.addHandler(handler)
    _console_handler = handler


def configure_run_log(path: Optional[str]) -> None:
    """Send everything - messages and tool output - to *path* as well as the screen.

    Called by each module's ``start_run`` with ``<output>/metawrap2.log``. Replaces any
    previously configured run log, so one process running several modules writes each one's
    log into its own output directory.
    """
    global _file_handler
    if _file_handler is not None:
        _logger.removeHandler(_file_handler)
        _file_handler.close()
        _file_handler = None
    if not path:
        return
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    # Rotating: one tool with runaway output (a progress bar the noise filter misses, a solver
    # looping on a warning) can write tens of GB. A normal run never reaches the cap.
    handler = logging.handlers.RotatingFileHandler(
        path, maxBytes=MAX_LOG_BYTES, backupCount=LOG_BACKUPS
    )
    handler.setFormatter(PlainFormatter())
    handler.setLevel(logging.DEBUG)
    _logger.addHandler(handler)
    _file_handler = handler


def log_tool_output(line: str, stream: str = "stdout") -> None:
    """Record one line of an external tool's output in the run log (not on the console).

    The runner already streams tool output to the screen; sending it through the console handler
    would print every line twice.
    """
    text = line.rstrip("\n")
    if text:
        _logger.log(TOOL_LEVEL, "[%s] %s", stream, text)


def _emit(message: str, delim: str, level: int) -> None:
    configure_console()
    _logger.log(level, "%s", message, extra={_DELIM_KEY: delim})


def comm(message: str) -> None:
    _emit(message, "-", logging.INFO)


def announcement(message: str) -> None:
    _emit(message, "#", logging.INFO)


def warning(message: str) -> None:
    _emit(message, "*", logging.WARNING)


def error(message: str, exit_code: int = 1) -> NoReturn:
    """Print an error banner and exit, matching the shell ``error`` helper."""
    _emit(message, "*", logging.ERROR)
    raise SystemExit(exit_code)
