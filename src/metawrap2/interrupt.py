"""Handle SIGTERM and SIGHUP the same way as Ctrl-C.

On Ctrl-C the running step is recorded as ``interrupted`` so ``--resume`` redoes it. On a cluster a
run is more often stopped by SIGTERM (a scheduler enforcing a wall-clock or memory limit) or SIGHUP
(the terminal closed), which killed the process and left the step recorded as ``running``.

:func:`install_handlers` maps these signals onto :class:`KeyboardInterrupt`, so the existing
interrupt handling applies: ``finally`` blocks run, the manifest records ``interrupted``, the
history entry is closed, and the exit code is 128 + signal number.

SIGINT keeps Python's own handler, except when it was inherited as ignored; see
:func:`_inherited_ignored_sigint`.
"""

from __future__ import annotations

import signal
from typing import Dict, List, Optional

#: Signals that stop the run. SIGHUP covers a dropped SSH session for a run started without
#: nohup or a scheduler.
STOP_SIGNALS: List[str] = ["SIGTERM", "SIGHUP"]

#: The signal number, set when a handler fires. Read it with :func:`received`.
_received: Optional[int] = None

_installed = False


class Stopped(KeyboardInterrupt):
    """Raised in place of a fatal signal.

    Subclasses KeyboardInterrupt so ``except KeyboardInterrupt`` handles it, while code can still
    tell a signal from Ctrl-C.
    """

    def __init__(self, signum: int):
        self.signum = signum
        super().__init__("stopped by %s" % signal_name(signum))


def signal_name(signum: int) -> str:
    try:
        return signal.Signals(signum).name
    except ValueError:
        return "signal %d" % signum


def received() -> Optional[int]:
    """The signal that stopped this process, if one did."""
    return _received


def exit_code(default: int = 130) -> int:
    """Exit code for how this process was stopped: 128 + signal number, else *default*."""
    return 128 + _received if _received is not None else default


def _inherited_ignored_sigint() -> List[str]:
    """``["SIGINT"]`` if this process inherited SIGINT as ignored, otherwise ``[]``.

    A shell sets SIGINT and SIGQUIT to SIG_IGN for background jobs, and CPython installs its
    SIGINT handler only when it finds SIG_DFL. A run started with ``&`` or under nohup therefore
    ignores ``kill -INT``. In that case SIGINT is handled here; otherwise Python's handler is kept.
    """
    if getattr(signal, "SIGINT", None) is None:
        return []
    try:
        return ["SIGINT"] if signal.getsignal(signal.SIGINT) is signal.SIG_IGN else []
    except (OSError, ValueError):
        return []


def install_handlers(signals: Optional[List[str]] = None) -> Dict[str, bool]:
    """Route each of *signals* to :class:`Stopped`. Returns which ones were installed.

    SIGINT is added only when it was inherited as ignored (see :func:`_inherited_ignored_sigint`).
    Safe to call more than once on any platform: a signal the OS lacks, or that cannot be handled
    from the current thread, is skipped without raising.
    """
    global _installed
    installed: Dict[str, bool] = {}

    def handler(signum, _frame):
        global _received
        _received = signum
        raise Stopped(signum)

    for name in list(signals or STOP_SIGNALS) + _inherited_ignored_sigint():
        signum = getattr(signal, name, None)
        if signum is None:
            installed[name] = False
            continue
        try:
            signal.signal(signum, handler)
            installed[name] = True
        except (OSError, ValueError, RuntimeError):
            # ValueError: not the main thread. OSError: signal cannot be caught.
            installed[name] = False
    _installed = any(installed.values())
    return installed


def handlers_installed() -> bool:
    return _installed
