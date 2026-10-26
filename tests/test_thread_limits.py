"""Tests for limiting the threads of numerical libraries in tools that fork workers.

Numerical libraries default to one thread per core, so a tool that forks N workers doing array maths
runs N x cores threads. CONCOCT does this; it used 30 CPU-hours per CAMI2 sample where 5.6 were
enough (same version and assembly, MKL in the environment instead of a pthreads OpenBLAS).
"""

from __future__ import annotations

import os

from metawrap2 import command
from metawrap2.command import THREAD_LIMITS
from metawrap2.modules import binning


def test_concoct_is_run_with_the_limits_and_the_other_binners_are_not(tmp_path, monkeypatch):
    """Only the concoct call gets THREAD_LIMITS.

    metaBAT2 and MaxBin2 do not fork array-maths workers. The learned binners are also excluded:
    torch is threaded internally without forking, so a one-thread cap would only slow them down.
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
    # The cut-up and merge helpers are I/O-bound and run without the limits.
    others = [env for cmd, env in calls if not cmd.startswith("concoct ")]
    assert all(env is None for env in others), others

    # The limits are added to the child environment, which keeps PYTHONNOUSERSITE and PATH.
    child = command.runner._child_env(True, THREAD_LIMITS)
    assert child["OMP_NUM_THREADS"] == "1"
    assert child["PYTHONNOUSERSITE"] == "1"
    assert child.get("PATH") == os.environ.get("PATH")
