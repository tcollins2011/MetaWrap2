"""MetaWrap2's logging: one :class:`logging.Logger`, two destinations, banner output on screen.

Everything MetaWrap2 says goes through the ``metawrap2`` logger. It carries two handlers:

**the console** (stdout)
    renders banner-boxed messages in the established MetaWrap2 style, so the on-screen output a
    user recognises is unchanged;

**the run log** (``<output>/metawrap2.log``)
    the same messages as plain timestamped lines, plus every line of every external tool's
    stdout and stderr - so one file tells the whole story of a run in order, which neither
    ``run.stdout`` nor ``run.stderr`` does on its own.

``run.stdout`` / ``run.stderr`` are still written separately, because a tool's stdout is
sometimes voluminous and it is useful to have the raw streams unmixed.

Why a Logger rather than ``print``: levels (a tool's chatter is DEBUG, a module's progress is
INFO, a problem is WARNING/ERROR) mean the console can stay readable while the file keeps
everything; and anything that wants to add a destination - a syslog, a GUI, a test capture -
attaches a handler instead of patching print.

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
from typing import NoReturn, Optional, TextIO

from .constants import BANNER_MAX_LINE as MAX_LINE
from .constants import BANNER_WIDTH as WIDTH

LOGGER_NAME = "metawrap2"

#: Name of the combined log written into a run's output directory.
RUN_LOG_NAME = "metawrap2.log"

#: Cap on the run log before it rolls over to metawrap2.log.1, .2, ... A normal run produces a
#: few hundred KB, so this is reached only by a tool whose output has gone wrong - which is
#: exactly when an unbounded log turns a failed step into a full filesystem.
MAX_LOG_BYTES = 256 * 1024 * 1024  # 256 MB

#: How many rolled-over logs to keep. The oldest output of a runaway tool is the interesting
#: part (where it started going wrong), but keeping the whole thing defeats the cap.
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


class ConsoleHandler(logging.StreamHandler):
    """A stream handler that resolves ``sys.stdout`` at emit time, not at construction.

    ``logging.StreamHandler(sys.stdout)`` captures the stream object. Anything that later
    replaces ``sys.stdout`` - a test harness capturing output, a caller redirecting it - leaves
    the handler writing to the old object, which at best goes unseen and at worst raises
    because the stream has been closed. Looking it up per record costs nothing and means the
    logger always writes wherever stdout currently points.
    """

    @property  # type: ignore[override]
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
    # Rotating, not plain: the log holds every line of every tool, and one pathological tool
    # (a progress bar that defeats the noise filter, a solver looping on a warning) can write
    # tens of GB - filling the very filesystem the run needs to finish. The cap is generous
    # enough that a normal run never rotates, so the common case is still one readable file.
    handler = logging.handlers.RotatingFileHandler(
        path, maxBytes=MAX_LOG_BYTES, backupCount=LOG_BACKUPS
    )
    handler.setFormatter(PlainFormatter())
    handler.setLevel(logging.DEBUG)
    _logger.addHandler(handler)
    _file_handler = handler


def log_tool_output(line: str, stream: str = "stdout") -> None:
    """Record one line of an external tool's output in the run log (not on the console).

    The runner streams tool output to the screen itself, live and unbuffered; duplicating it
    through the console handler would double every line.
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
