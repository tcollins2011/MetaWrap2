"""`run`/`cohort -j N` gives concurrent GPU binning steps one card each.

Without it every learned binner on bare ``--device cuda`` trains on the framework's default card.
"""

from __future__ import annotations

import threading
from typing import Dict, List

import pytest

from metawrap2.commands import run as driver


def _binning(name: str, *options: str) -> driver.Step:
    return driver.Step(
        name="binning:%s" % name, module="binning", argv=list(options), parallel_group="binning"
    )


@pytest.mark.parametrize(
    "options, device",
    [
        (("--preset", "deep"), "cuda"),
        (("--preset", "gpu"), "cuda"),  # the older name of deep
        (("--preset", "deep", "--device", "cpu"), None),
        (("--preset", "deep", "--device=cuda:1"), "cuda:1"),
        (("--semibin2", "--device", "auto"), "auto"),
        (("--metabat2", "--device", "cuda"), None),  # no learned binner, so no GPU
        (("--preset", "default"), None),  # learned binners train on the CPU by default
        ((), None),
    ],
)
def test_step_device_follows_the_binning_module(options, device):
    assert driver.step_device(_binning("A", *options)) == device


def test_only_steps_that_leave_the_card_open_are_assigned():
    assert driver.assignable(_binning("A", "--preset", "deep"))
    assert driver.assignable(_binning("A", "--semibin2", "--device", "auto"))
    assert not driver.assignable(_binning("A", "--preset", "deep", "--device", "cuda:1"))
    assert not driver.assignable(driver.Step(name="assembly:A", module="assembly", argv=[]))


def test_cards_come_from_an_inherited_cuda_visible_devices(monkeypatch):
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    assert driver.gpu_cards(2) == ["0", "1"]
    assert driver.gpu_cards(0) == []
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "3,1,5")
    assert driver.gpu_cards(2) == ["3", "1"]
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "-1")
    assert driver.gpu_cards(2) == []


def test_allocator_hands_out_the_least_used_card():
    cards = driver.CardAllocator(["0", "1"])
    assert [cards.acquire(), cards.acquire(), cards.acquire()] == ["0", "1", "0"]
    cards.release("1")
    assert cards.acquire() == "1"
    assert driver.CardAllocator([]).acquire() is None


def test_pin_rewrites_bare_cuda_and_sets_the_env_for_auto():
    env: Dict[str, str] = {}
    argv = ["--preset", "deep"]
    driver.pin_to_card(argv, env, "cuda", "1")
    assert argv == ["--preset", "deep", "--device", "cuda:1"] and env == {}

    argv = ["--device", "cuda", "-t", "4"]
    driver.pin_to_card(argv, env, "cuda", "1")
    assert argv == ["--device", "cuda:1", "-t", "4"]

    argv = ["--device=cuda"]
    driver.pin_to_card(argv, env, "cuda", "0")
    assert argv == ["--device=cuda:0"]

    # auto keeps its CPU fallback, so the card is chosen through the environment.
    argv = ["--semibin2", "--device", "auto"]
    driver.pin_to_card(argv, env, "auto", "1")
    assert argv == ["--semibin2", "--device", "auto"] and env["CUDA_VISIBLE_DEVICES"] == "1"


def test_banner_names_the_cards():
    assert driver.group_banner(2, "binning", 2, 16, 60, ["0", "1"]) == (
        "RUNNING 2 x binning (2 at a time, 16 threads, 60GB and 1 GPU each: cuda:0, cuda:1)"
    )
    # More steps at a time than cards: every card is listed and shared.
    assert driver.group_banner(3, "binning", 3, 16, 60, ["0"]) == (
        "RUNNING 3 x binning (3 at a time, 16 threads, 60GB and 1 of 1 GPUs each: cuda:0)"
    )
    # Without GPUs the banner is unchanged.
    assert driver.group_banner(2, "read_qc", 2, 16, 60, []) == (
        "RUNNING 2 x read_qc (2 at a time, 16 threads and 60GB each)"
    )
    assert driver.group_banner(2, "kraken2", 2, 16, None, []) == (
        "RUNNING 2 x kraken2 (2 at a time, 16 threads each)"
    )


def test_gpu_share_needs_concurrency_and_a_budget(monkeypatch):
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    deep = [_binning(s, "--preset", "deep") for s in "AB"]
    assert driver.gpu_share(deep, 2, 2) == ["0", "1"]
    assert driver.gpu_share(deep, 1, 2) == []  # one at a time: nothing to divide
    assert driver.gpu_share(deep, 2, 0) == []  # --max-gpus 0
    named = [_binning(s, "--preset", "deep", "--device", "cuda:0") for s in "AB"]
    assert driver.gpu_share(named, 2, 2) == []


def _study() -> driver.Study:
    return driver.Study(
        samples=[driver.Sample(name=s, r1="raw/%s_1.fq" % s, r2="raw/%s_2.fq" % s) for s in "AB"],
        coassemble=False,
    )


def test_concurrent_steps_run_on_different_cards(tmp_path, monkeypatch):
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    layout = driver.Layout(str(tmp_path / "study"))
    steps = [_binning(s, "--preset", "deep") for s in "AB"]
    state = driver.StepState(layout, resume=False)
    both_started = threading.Barrier(2, timeout=10)
    seen: List[str] = []
    lock = threading.Lock()

    def fake_run_measured(argv, env=None, **_kwargs):
        both_started.wait()  # both steps hold a card at the same time
        with lock:
            seen.append(argv[argv.index("--device") + 1])
        assert "CUDA_VISIBLE_DEVICES" not in env
        return 0, {}

    banners: List[str] = []
    monkeypatch.setattr(driver._usage, "run_measured", fake_run_measured)
    monkeypatch.setattr(driver, "announcement", banners.append)
    monkeypatch.setattr(driver, "comm", lambda *_a: None)
    assert driver.execute(steps, layout, _study(), 32, 2, None, [], state, False, 60, gpus=2) == 0
    assert sorted(seen) == ["cuda:0", "cuda:1"]
    assert banners[0].endswith("and 1 GPU each: cuda:0, cuda:1)")
    assert state.manifest.step("binning:A").data["grant"]["gpu"] in ("cuda:0", "cuda:1")


def test_dry_run_plan_shows_the_round_robin(tmp_path, monkeypatch, capsys):
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    layout = driver.Layout(str(tmp_path / "study"))
    steps = [_binning(s, "--semibin2", "--device", "auto") for s in "AB"]
    state = driver.StepState(layout, resume=False)
    monkeypatch.setattr(driver, "announcement", lambda *_a: None)
    monkeypatch.setattr(driver, "comm", lambda *_a: None)
    assert driver.execute(steps, layout, _study(), 32, 2, None, [], state, True, 60, gpus=2) == 0
    plans = sorted(line for line in capsys.readouterr().out.splitlines() if "[plan]" in line)
    assert [p.split()[1] for p in plans] == ["CUDA_VISIBLE_DEVICES=0", "CUDA_VISIBLE_DEVICES=1"]


def test_a_named_card_is_left_alone(tmp_path, monkeypatch):
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    layout = driver.Layout(str(tmp_path / "study"))
    steps = [_binning(s, "--preset", "deep", "--device", "cuda:1") for s in "AB"]
    state = driver.StepState(layout, resume=False)
    seen: List[str] = []

    def fake_run_measured(argv, env=None, **_kwargs):
        seen.append(argv[argv.index("--device") + 1])
        return 0, {}

    banners: List[str] = []
    monkeypatch.setattr(driver._usage, "run_measured", fake_run_measured)
    monkeypatch.setattr(driver, "announcement", banners.append)
    monkeypatch.setattr(driver, "comm", lambda *_a: None)
    assert driver.execute(steps, layout, _study(), 32, 2, None, [], state, False, 60, gpus=2) == 0
    assert seen == ["cuda:1", "cuda:1"]
    assert "GPU" not in banners[0]
