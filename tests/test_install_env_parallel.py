"""Parallel `install-env`: the first environment is built alone, and failures are retried serially."""

import os
import threading

import pytest

from metawrap2.commands import install_env as ie

ENVS = ("reads", "bins", "taxonomy")
MODULES = ("read_qc", "binning", "kraken2")


@pytest.fixture
def queued(tmp_path, monkeypatch):
    """Three envs that reach the build step, with _create_one replaced by a recorder."""
    target = tmp_path / "envs"
    (target / "locks").mkdir(parents=True)
    monkeypatch.setenv("METAWRAP2_ENVS_DIR", str(target))
    for name in ENVS:
        (target / ("%s.yaml" % name)).write_text(
            "name: metawrap2-%s\ndependencies: [samtools]\n" % name
        )
    monkeypatch.setattr(ie, "ensure_mamba", lambda: True)
    monkeypatch.setattr(ie, "conda_env_exists", lambda env: False)
    monkeypatch.setattr(ie, "env_prefix", lambda env: None)
    state = {
        "calls": [],
        "running": 0,
        "first_alone": None,
        "fail_once": set(),
        "fail_always": set(),
    }
    lock = threading.Lock()

    def fake_create(module, env, yaml, force, log_dir, lock_file=None):
        with lock:
            state["running"] += 1
            if not state["calls"]:
                state["first_alone"] = state["running"] == 1
            state["calls"].append((env, force))
            fail = env in state["fail_always"] or (
                env in state["fail_once"] and [c[0] for c in state["calls"]].count(env) == 1
            )
        with lock:
            state["running"] -= 1
        return module, not fail, os.path.join(log_dir, "x.log")

    monkeypatch.setattr(ie, "_create_one", fake_create)
    return state


def test_the_first_environment_is_built_before_the_others_start(queued):
    assert ie.main(list(MODULES) + ["--no-test", "-t", "4"]) == 0
    assert queued["first_alone"] is True
    assert sorted(c[0] for c in queued["calls"]) == sorted("metawrap2-%s" % e for e in ENVS)


def test_a_failed_parallel_build_is_retried_on_its_own(queued, capsys):
    queued["fail_once"].add("metawrap2-bins")
    assert ie.main(list(MODULES) + ["--no-test", "-t", "4"]) == 0
    assert [c[0] for c in queued["calls"]].count("metawrap2-bins") == 2
    assert "Retrying 1 environment(s) one at a time" in capsys.readouterr().out


def test_an_environment_that_fails_twice_fails_the_command(queued):
    queued["fail_always"].add("metawrap2-taxonomy")
    assert ie.main(list(MODULES) + ["--no-test", "-t", "4"]) == 1


def test_a_serial_build_is_not_retried(queued):
    queued["fail_once"].add("metawrap2-bins")
    assert ie.main(list(MODULES) + ["--no-test", "-t", "1"]) == 1
    assert [c[0] for c in queued["calls"]].count("metawrap2-bins") == 1


def test_the_envs_still_building_are_listed_after_each_finishes(queued, monkeypatch, capsys):
    """After bins finishes, taxonomy is still building and is named; the last finish prints none."""
    started = threading.Event()
    listed = threading.Event()
    real_line = ie.still_building_line

    def line(building, total):
        text = real_line(building, total)
        if text:
            listed.set()
        return text

    def fake_create(module, env, yaml, force, log_dir, lock_file=None):
        if env == "metawrap2-taxonomy":
            started.set()
            listed.wait(5)
        elif env == "metawrap2-bins":
            started.wait(5)
        return module, True, os.path.join(log_dir, "x.log")

    monkeypatch.setattr(ie, "still_building_line", line)
    monkeypatch.setattr(ie, "_create_one", fake_create)
    assert ie.main(list(MODULES) + ["--no-test", "-t", "4"]) == 0
    lines = [ln for ln in capsys.readouterr().out.splitlines() if "still building" in ln]
    assert lines == ["  still building: taxonomy (1 of 3)"]


def test_the_still_building_line():
    assert ie.still_building_line([], 3) == ""
    assert (
        ie.still_building_line(["metawrap2-learned", "metawrap2-annotate"], 12)
        == "  still building: annotate, learned (2 of 12)"
    )
