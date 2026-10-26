"""Tests for treating a scheduler kill the same as Ctrl-C.

:class:`Stopped` must be catchable as ``KeyboardInterrupt``, so the ``finally`` and ``except``
blocks written for Ctrl-C also record the outcome of a killed run.
"""

from __future__ import annotations

import os
import signal

import pytest

from metawrap2 import interrupt


@pytest.fixture(autouse=True)
def _restore_handlers():
    """Put the process's signal handlers back, so one test cannot affect another."""
    saved = {}
    for name in interrupt.STOP_SIGNALS + ["SIGINT"]:
        signum = getattr(signal, name, None)
        if signum is not None:
            saved[signum] = signal.getsignal(signum)
    yield
    for signum, handler in saved.items():
        signal.signal(signum, handler)


def test_a_background_runs_ignored_sigint_is_taken_back(monkeypatch):
    """A shell sets a background job's SIGINT to SIG_IGN, and Python then leaves it ignored.

    Runs started with "&", nohup or a batch scheduler start this way, and "kill -INT" then had no
    effect. Checks that the handler is reinstalled.
    """
    monkeypatch.setattr(interrupt, "_received", None, raising=False)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    assert interrupt.install_handlers(["SIGTERM"])["SIGINT"] is True

    with pytest.raises(interrupt.Stopped):
        os.kill(os.getpid(), signal.SIGINT)
    assert interrupt.received() == signal.SIGINT

    # A kill runs the cleanup (the manifest depends on this).
    cleaned = []
    with pytest.raises(KeyboardInterrupt):
        try:
            os.kill(os.getpid(), signal.SIGTERM)
        finally:
            cleaned.append("recorded the outcome")
    assert cleaned == ["recorded the outcome"]
    assert interrupt.received() == signal.SIGTERM
    # 128 + signal number is the shell convention, so `echo $?` is readable.
    assert interrupt.exit_code() == 128 + int(signal.SIGTERM)
