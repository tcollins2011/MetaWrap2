"""Peak GPU memory measurement, which rusage does not provide.

SemiBin2, VAMB and COMEBin can use a GPU, and a step fails if the model does not fit, so the VRAM
figures in the README's hardware table are measured per step. Tested without a driver: parsing
nvidia-smi output, the per-device peak, and that a step that used no GPU reports nothing, not zero.
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

    A zero would lower a peak when combined with real samples. nvidia-smi writes `[Not Found]` for
    a process that exited between enumeration and query, and `[N/A]` on some drivers in
    containers.
    """
    text = "[Not Found], GPU-aaa, [N/A]\n999, GPU-aaa, [Not Found]\n1, GPU-aaa, 64\n"
    assert usage.parse_compute_apps(text) == [(1, "GPU-aaa", 64)]


def test_peak_is_per_device_not_summed_across_them():
    """Two cards holding 900 MiB each is a 900 MiB requirement, not 1800.

    The figure is per card, so devices are not summed. Processes on one device are summed, since
    they must fit together.
    """
    samples = [
        [(1, "GPU-a", 500), (2, "GPU-a", 400), (3, "GPU-b", 900)],
        [(1, "GPU-a", 100)],
    ]
    assert usage.peak_device_memory(samples) == 900


def test_peak_counts_only_the_processes_asked_for():
    """GPU processes outside this step's process tree are not counted."""
    samples = [[(1, "GPU-a", 300), (4242, "GPU-a", 20000)]]
    assert usage.peak_device_memory(samples, belongs=lambda pid: pid == 1) == 300


def test_no_samples_is_zero_rather_than_an_error():
    assert usage.peak_device_memory([]) == 0
    assert usage.peak_device_memory([[]]) == 0


@pytest.mark.skipif(not os.path.isdir("/proc"), reason="reads /proc, which only Linux has")
def test_process_tree_walk_finds_a_parent_and_stops_at_a_stranger():
    assert usage.in_process_tree(os.getpid(), os.getpid())
    assert usage.in_process_tree(os.getpid(), os.getppid())
    # pid 1 is not a descendant of this step, so its GPU use is not attributed to the step.
    assert not usage.in_process_tree(1, os.getpid())


def test_a_step_that_touched_no_gpu_reports_no_vram_key():
    """A step that used no GPU has no VRAM fields, instead of a zero.

    `peak_gpu_memory_bytes: 0` on CPU steps of a GPU machine would read as a measurement.
    """
    rc, measured = usage.run_measured(["true"])
    assert rc == 0
    assert "peak_gpu_memory_bytes" not in measured
    assert "wall_seconds" in measured


def test_a_standalone_module_reports_vram_too(monkeypatch):
    """A module invoked directly reports peak VRAM, as it does under `run`.

    A standalone module measures itself from the inside, and rusage has no GPU figure.
    """
    monkeypatch.setattr(_common, "_active_gpu_watcher", _Watcher(2048, {"GPU-a"}))
    measured = _common.self_usage(_Rec())
    assert measured["peak_gpu_memory_bytes"] == 2048 * usage.MIB
    assert measured["gpu_devices"] == 1
    assert "2.0 GB peak VRAM on 1 GPU(s)" in usage.summarise(measured)


def test_a_standalone_cpu_run_on_a_gpu_machine_reports_none(monkeypatch):
    """A watcher that recorded no GPU use adds no VRAM fields."""
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
    # The CPU figures are still present alongside VRAM.
    assert "2.0 cores busy" in line and "2.0 GB peak memory" in line
