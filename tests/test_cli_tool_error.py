"""A failed tool ends the run with its name, exit code and last output, not a traceback."""

from __future__ import annotations

import types

import pytest

from metawrap2 import cli, command
from metawrap2.command import ToolError


@pytest.fixture
def failing_module(monkeypatch):
    tail = "\n".join("line %d" % i for i in range(1, 26))

    def main(argv):
        raise ToolError(
            tool="fastqc",
            returncode=2,
            cmd=["fastqc", "--outdir", "/o"],
            log_tail=tail,
            hint="check the input",
        )

    fake = types.SimpleNamespace(main=main)
    monkeypatch.setitem(cli.MODULES, "fake_qc", ("fake.module", "test"))
    monkeypatch.setattr(cli.importlib, "import_module", lambda name: fake)
    yield
    command.runner.configure(verbose_logs=False)


def test_tool_error_prints_the_tool_and_its_last_lines(failing_module, capsys):
    assert cli.main(["fake_qc"]) == 1
    err = capsys.readouterr().err
    assert "fastqc failed (exit code 2)." in err
    assert "Command: fastqc --outdir /o" in err
    assert "Hint: check the input" in err
    assert "line 25" in err and "line 11" in err
    assert "line 10\n" not in err  # only the last 15 lines
    assert "Traceback" not in err
    assert "--verbose-logs" in err


def test_verbose_logs_adds_the_traceback(failing_module, capsys):
    assert cli.main(["--verbose-logs", "fake_qc"]) == 1
    err = capsys.readouterr().err
    assert "Traceback" in err
    assert "fastqc failed (exit code 2)." in err
