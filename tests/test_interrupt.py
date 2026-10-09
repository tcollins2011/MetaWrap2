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
    for name in interrupt.STOP_SIGNALS + ["SIGINT"]:
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


def test_a_background_runs_ignored_sigint_is_taken_back(monkeypatch):
    """A shell sets a background job's SIGINT to SIG_IGN, and Python then leaves it ignored.

    That is the state every long run starts in - "&", nohup, a batch scheduler - so "kill -INT" on
    one did nothing at all, silently: no KeyboardInterrupt, no message, the tool kept running.
    """
    monkeypatch.setattr(interrupt, "_received", None, raising=False)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    assert interrupt.install_handlers(["SIGTERM"])["SIGINT"] is True

    with pytest.raises(interrupt.Stopped):
        os.kill(os.getpid(), signal.SIGINT)
    assert interrupt.received() == signal.SIGINT


def test_pythons_own_sigint_handling_is_left_alone(monkeypatch):
    """Where Ctrl-C already works, this module has nothing to add - and KeyboardInterrupt from the
    interpreter is not a Stopped, so nothing reports an interactive Ctrl-C as a scheduler kill."""
    monkeypatch.setattr(interrupt, "_received", None, raising=False)
    signal.signal(signal.SIGINT, signal.default_int_handler)
    assert "SIGINT" not in interrupt.install_handlers(["SIGTERM"])
    assert signal.getsignal(signal.SIGINT) is signal.default_int_handler


def test_the_signal_that_stopped_us_is_recorded(monkeypatch):
    monkeypatch.setattr(interrupt, "_received", None, raising=False)
    interrupt.install_handlers(["SIGTERM"])
    with pytest.raises(interrupt.Stopped):
        os.kill(os.getpid(), signal.SIGTERM)
    assert interrupt.received() == signal.SIGTERM
    # 128 + signal number is the shell convention, so `echo $?` is readable.
    assert interrupt.exit_code() == 128 + int(signal.SIGTERM)
