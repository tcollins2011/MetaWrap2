"""Treat being killed the same as being interrupted.

Ctrl-C was already handled: the running step is recorded as ``interrupted`` so ``--resume``
redoes it rather than trusting half-written outputs. But Ctrl-C is the *rare* case. The common
one on a cluster is a scheduler sending SIGTERM because the job hit its wall clock or memory
limit, or SIGHUP because the terminal went away - and those killed the process outright, leaving
the step recorded as ``running`` forever and nothing to explain why.

:func:`install_handlers` maps those signals onto :class:`KeyboardInterrupt`, which means every
piece of interrupt handling already written applies to them unchanged: the ``finally`` blocks
run, the manifest records ``interrupted``, the history entry is closed out, and the exit code is
the conventional 128 + signal number.

This is deliberately *not* a general signal framework. Only the signals that mean "stop now,
your outputs are not finished" are touched; SIGINT keeps Python's own behaviour, which already
raises KeyboardInterrupt - except in the one case where there is no such behaviour to keep, which
:func:`install_handlers` explains.
"""

from __future__ import annotations

import signal
from typing import Dict, List, Optional

#: Signals that mean "you are being stopped". SIGHUP is included because a dropped SSH session
#: is how a long run most often dies when it was not started under nohup or a scheduler.
STOP_SIGNALS: List[str] = ["SIGTERM", "SIGHUP"]

#: Set when a handler fires, so an exit code and a message can name the cause. Read it with
#: :func:`received`.
_received: Optional[int] = None

_installed = False


class Stopped(KeyboardInterrupt):
    """Raised in place of a fatal signal.

    Subclasses KeyboardInterrupt on purpose: every existing ``except KeyboardInterrupt`` keeps
    working, while code that wants to distinguish a scheduler kill from a user's Ctrl-C can.
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
    """The conventional exit code for how this process was stopped: 128 + signal number."""
    return 128 + _received if _received is not None else default


def _inherited_ignored_sigint() -> List[str]:
    """``["SIGINT"]`` if this process inherited SIGINT as ignored, otherwise ``[]``.

    A shell puts a background job's SIGINT and SIGQUIT to SIG_IGN, and the disposition is
    inherited; CPython only installs its own SIGINT handler when it finds SIG_DFL, so a run
    started with ``&`` or under nohup - exactly what a long run is started with - cannot be
    interrupted by SIGINT at all. "kill -INT" on it does nothing, silently, and the only way left
    to stop it is SIGTERM.

    Taking SIGINT over in that case restores the behaviour the rest of this module assumes rather
    than overriding it: where Python's handler exists it is left alone, and where the shell threw
    it away it is put back.
    """
    if getattr(signal, "SIGINT", None) is None:
        return []
    try:
        return ["SIGINT"] if signal.getsignal(signal.SIGINT) is signal.SIG_IGN else []
    except (OSError, ValueError):
        return []


def install_handlers(signals: Optional[List[str]] = None) -> Dict[str, bool]:
    """Route each of *signals* to :class:`Stopped`. Returns which ones were installed.

    SIGINT is added to them only when it was inherited as ignored - see
    :func:`_inherited_ignored_sigint`.

    Safe to call more than once and on any platform: a signal this OS does not have, or that
    cannot be handled from the current thread, is skipped rather than raising. Losing the
    handler must never be the reason a run fails to start.
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
