"""Tests for treating a scheduler kill the same as Ctrl-C.

The whole design rests on one property: :class:`Stopped` must be catchable by code that only
knows about ``KeyboardInterrupt``, because every ``finally`` and ``except`` already written for
Ctrl-C is what makes a killed run record its outcome.
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
    for name in interrupt.STOP_SIGNALS:
        signum = getattr(signal, name, None)
        if signum is not None:
            saved[signum] = signal.getsignal(signum)
    yield
    for signum, handler in saved.items():
        signal.signal(signum, handler)


def test_stopped_is_a_keyboard_interrupt():
    """The point of the whole module: existing Ctrl-C handling applies unchanged."""
    assert issubclass(interrupt.Stopped, KeyboardInterrupt)
    with pytest.raises(KeyboardInterrupt):
        raise interrupt.Stopped(signal.SIGTERM)


def test_sigterm_unwinds_a_finally_block(monkeypatch):
    """The behaviour the manifest depends on: a kill runs the cleanup, it does not skip it."""
    monkeypatch.setattr(interrupt, "_received", None, raising=False)
    interrupt.install_handlers(["SIGTERM"])
    cleaned = []
    with pytest.raises(interrupt.Stopped):
        try:
            os.kill(os.getpid(), signal.SIGTERM)
        finally:
            cleaned.append("recorded the outcome")
    assert cleaned == ["recorded the outcome"]


def test_the_signal_that_stopped_us_is_recorded(monkeypatch):
    monkeypatch.setattr(interrupt, "_received", None, raising=False)
    interrupt.install_handlers(["SIGTERM"])
    with pytest.raises(interrupt.Stopped):
        os.kill(os.getpid(), signal.SIGTERM)
    assert interrupt.received() == signal.SIGTERM
    # 128 + signal number is the shell convention, so `echo $?` is readable.
    assert interrupt.exit_code() == 128 + int(signal.SIGTERM)
