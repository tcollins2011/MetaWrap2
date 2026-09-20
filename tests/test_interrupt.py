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


def test_stopped_names_the_signal():
    exc = interrupt.Stopped(signal.SIGTERM)
    assert exc.signum == signal.SIGTERM
    assert "SIGTERM" in str(exc)


def test_install_handlers_reports_what_it_installed():
    installed = interrupt.install_handlers()
    assert installed["SIGTERM"] is True
    assert interrupt.handlers_installed() is True


def test_install_handlers_is_safe_to_call_twice():
    interrupt.install_handlers()
    assert interrupt.install_handlers()["SIGTERM"] is True


def test_a_signal_this_platform_lacks_is_skipped_not_raised():
    assert interrupt.install_handlers(["SIGNAL_THAT_DOES_NOT_EXIST"]) == {
        "SIGNAL_THAT_DOES_NOT_EXIST": False
    }


def test_an_uncatchable_signal_is_skipped(monkeypatch):
    """Losing the handler must never be the reason a run fails to start."""

    def refuse(signum, handler):
        raise OSError("cannot catch this one")

    monkeypatch.setattr(signal, "signal", refuse)
    assert interrupt.install_handlers(["SIGTERM"]) == {"SIGTERM": False}


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


def test_exit_code_falls_back_when_no_signal_arrived(monkeypatch):
    monkeypatch.setattr(interrupt, "_received", None, raising=False)
    assert interrupt.exit_code() == 130  # the Ctrl-C convention
    assert interrupt.exit_code(default=1) == 1


def test_signal_name_handles_an_unknown_number():
    assert interrupt.signal_name(int(signal.SIGTERM)) == "SIGTERM"
    assert "9999" in interrupt.signal_name(9999)
