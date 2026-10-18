"""Peak GPU memory, which rusage cannot see and which decides whether a step runs at all.

Three tools here can use a GPU - SemiBin2, VAMB and COMEBin - and a card either holds the model or
the step dies in the allocator, so the VRAM figures in the README's hardware table are measured per
step rather than copied from an upstream README. These tests pin the parts that have to be right
for those figures to mean anything, all without a driver: the parse of what nvidia-smi writes, the
per-device peak, and the rule that a step which touched no GPU reports nothing rather than zero.
"""

import os

import pytest

from metawrap2 import usage
from metawrap2.modules import _common


class _Rec:
    """Only what self_usage reads off a RunRecorder."""

    def elapsed(self):
        return 12.0


class _Watcher:
    def __init__(self, peak_mib, devices):
        self.peak_mib = peak_mib
        self.devices = devices


def test_parses_a_driver_row():
    rows = usage.parse_compute_apps("12345, GPU-aaa, 2048\n67890, GPU-bbb, 512\n")
    assert rows == [(12345, "GPU-aaa", 2048), (67890, "GPU-bbb", 512)]


def test_drops_rows_the_driver_could_not_fill():
    """A row with no readable memory is dropped, never counted as zero.

    The regression: a zero pulls a *peak* down when it is averaged or summed with real samples,
    while a dropped row can only fail to raise it. nvidia-smi writes `[Not Found]` for a process
    that exited between the enumeration and the query, and `[N/A]` on some driver versions in
    containers.
    """
    text = "[Not Found], GPU-aaa, [N/A]\n999, GPU-aaa, [Not Found]\n1, GPU-aaa, 64\n"
    assert usage.parse_compute_apps(text) == [(1, "GPU-aaa", 64)]


def test_peak_is_per_device_not_summed_across_them():
    """Two cards holding 900 MiB each is a 900 MiB requirement, not 1800.

    The figure answers "will this fit on one card", so devices are never added together. Within a
    device the processes *are* added, because several workers on one card have to fit together.
    """
    samples = [
        [(1, "GPU-a", 500), (2, "GPU-a", 400), (3, "GPU-b", 900)],
        [(1, "GPU-a", 100)],
    ]
    assert usage.peak_device_memory(samples) == 900


def test_peak_counts_only_the_processes_asked_for():
    """Another tenant's training job on the same card must not become this step's requirement."""
    samples = [[(1, "GPU-a", 300), (4242, "GPU-a", 20000)]]
    assert usage.peak_device_memory(samples, belongs=lambda pid: pid == 1) == 300


def test_no_samples_is_zero_rather_than_an_error():
    assert usage.peak_device_memory([]) == 0
    assert usage.peak_device_memory([[]]) == 0


@pytest.mark.skipif(not os.path.isdir("/proc"), reason="reads /proc, which only Linux has")
def test_process_tree_walk_finds_a_parent_and_stops_at_a_stranger():
    assert usage.in_process_tree(os.getpid(), os.getpid())
    assert usage.in_process_tree(os.getpid(), os.getppid())
    # pid 1 is nobody's descendant here, so a GPU process belonging to another user's run is not
    # attributed to this step.
    assert not usage.in_process_tree(1, os.getpid())


def test_a_step_that_touched_no_gpu_reports_no_vram_key():
    """Absence rather than zero, for the same reason a blank cell is not a zero elsewhere.

    `peak_gpu_memory_bytes: 0` on every CPU step of a machine that happens to have a card would
    read as "measured, used nothing", which is a claim about the tool rather than about the step.
    """
    rc, measured = usage.run_measured(["true"])
    assert rc == 0
    assert "peak_gpu_memory_bytes" not in measured
    assert "wall_seconds" in measured


def test_a_standalone_module_reports_vram_too(monkeypatch):
    """The regression: VRAM appeared under `run` and nowhere else.

    A module invoked directly measures itself from the inside, and rusage has no GPU figure to
    reconstruct one from - so the three binners that can train on a GPU printed peak VRAM as a
    study step and nothing at all when run the way the docs run them.
    """
    monkeypatch.setattr(_common, "_active_gpu_watcher", _Watcher(2048, {"GPU-a"}))
    measured = _common.self_usage(_Rec())
    assert measured["peak_gpu_memory_bytes"] == 2048 * usage.MIB
    assert measured["gpu_devices"] == 1
    assert "2.0 GB peak VRAM on 1 GPU(s)" in usage.summarise(measured)


def test_a_standalone_cpu_run_on_a_gpu_machine_reports_none(monkeypatch):
    """A watcher that saw nothing contributes nothing, even though it ran."""
    monkeypatch.setattr(_common, "_active_gpu_watcher", _Watcher(0, set()))
    measured = _common.self_usage(_Rec())
    assert "peak_gpu_memory_bytes" not in measured
    assert "VRAM" not in usage.summarise(measured)


def test_summary_states_vram_and_how_many_devices_held_it():
    line = usage.summarise(
        {
            "wall_seconds": 100,
            "peak_rss_bytes": 2 * 1024**3,
            "user_seconds": 200,
            "system_seconds": 0,
            "peak_gpu_memory_bytes": 6 * 1024**3,
            "gpu_devices": 2,
        }
    )
    assert "6.0 GB peak VRAM on 2 GPU(s)" in line
    # The CPU figures are still there: VRAM is an addition to the cost line, not a replacement.
    assert "2.0 cores busy" in line and "2.0 GB peak memory" in line
