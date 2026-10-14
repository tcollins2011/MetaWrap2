"""Confining a tool to the parallelism it was asked for.

A numerical library picks its thread count from the core count unless told otherwise, so a tool that
forks N workers and does array maths in each gets N x cores threads contending. CONCOCT does exactly
that, and it cost 30 CPU-hours per CAMI2 sample where 5.6 was enough - found by running the same
CONCOCT version on the same assembly under a pipeline whose environment carried MKL rather than a
pthreads OpenBLAS.
"""

from __future__ import annotations

import os

from metawrap2 import command
from metawrap2.command import THREAD_LIMITS
from metawrap2.modules import binning


def test_concoct_is_run_with_the_limits_and_the_other_binners_are_not(tmp_path, monkeypatch):
    """CONCOCT forks and does array maths; metaBAT2 and MaxBin2 do neither, so they are left alone.

    The learned binners are deliberately excluded too: torch is internally threaded and does not fork
    per worker, so capping it to one thread would make them slower rather than cheaper.
    """
    calls = []

    def fake_run(cmd, **kw):
        calls.append((cmd if isinstance(cmd, str) else " ".join(cmd), kw.get("extra_env")))
        return 0

    monkeypatch.setattr(binning, "run", fake_run)
    monkeypatch.setattr(binning, "announcement", lambda *a, **k: None)
    monkeypatch.setattr(binning, "comm", lambda *a, **k: None)
    monkeypatch.setattr(binning, "ensure_dir", lambda *a, **k: None)
    monkeypatch.setattr(binning, "_adopt_bins", lambda *a, **k: None)
    monkeypatch.setattr(binning, "_split_concoct_bins", lambda *a, **k: 1)
    monkeypatch.setattr(binning, "dry_run", lambda: False)

    class _Args:
        min_len = 1000
        threads = 96
        memory = 100
        output = str(tmp_path / "out")
        universal = False

    binning._run_concoct(_Args(), "e", str(tmp_path), "asm.fa")
    limited = [env for cmd, env in calls if cmd.startswith("concoct ")]
    assert limited and limited[0] == THREAD_LIMITS, "concoct must be capped"
    # The cut-up and merge helpers around it are I/O, not maths; they do not need capping and must not
    # silently acquire it, or a future reader will think the cap is cargo cult.
    others = [env for cmd, env in calls if not cmd.startswith("concoct ")]
    assert all(env is None for env in others), others

    # The limits reach the child without dropping the rest of the environment: PYTHONNOUSERSITE
    # still applied, and a child handed a five-variable environment would lose PATH and fail to
    # find its own tools.
    child = command.runner._child_env(True, THREAD_LIMITS)
    assert child["OMP_NUM_THREADS"] == "1"
    assert child["PYTHONNOUSERSITE"] == "1"
    assert child.get("PATH") == os.environ.get("PATH")
