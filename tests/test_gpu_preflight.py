"""Tests for the GPU preflight that checks a card before a GPU step starts.

`nvidia-smi` is faked, so no GPU or driver is needed.
"""

from __future__ import annotations

import types

import pytest

from metawrap2 import limits, usage
from metawrap2.modules import binning


def _card(index: int, utilisation: int = 0, free_mb: int = 40000, processes=()) -> limits.GpuState:
    return limits.GpuState(
        index=index,
        uuid="GPU-%d" % index,
        name="NVIDIA GPU",
        utilisation_percent=utilisation,
        memory_total_mb=40960,
        memory_free_mb=free_mb,
        processes=tuple(processes),
    )


def _probe(monkeypatch, *states: limits.GpuState) -> None:
    monkeypatch.setattr(limits, "query_gpus", lambda: limits.GpuProbe(states=tuple(states)))


def _busy_card(index: int) -> limits.GpuState:
    """A card held by another process at high utilisation."""
    return _card(
        index,
        utilisation=94,
        free_mb=2100,
        processes=[limits.GpuProcess(pid=31274, memory_mb=38000)],
    )


# --- what the driver is asked, and what is made of the answer ------------------------------


def test_the_query_is_parsed_field_by_field_so_one_missing_number_is_not_a_zero():
    """`[N/A]` is parsed as unknown, not as 0."""
    text = "0, GPU-aaa, NVIDIA GPU, 94, 40960, 2100\n1, GPU-bbb, NVIDIA GPU, [N/A], 40960, [N/A]\n"
    first, second = limits.parse_gpu_query(text)
    assert (first.index, first.utilisation_percent, first.memory_free_mb) == (0, 94, 2100)
    assert second.utilisation_percent is None and second.memory_free_mb is None
    assert "utilisation unknown" in second.describe()
    assert "free memory unknown" in second.describe()


def test_processes_are_joined_onto_cards_by_uuid(monkeypatch):
    """Processes from --query-compute-apps are matched to cards by uuid."""
    csv = "0, GPU-aaa, GPU, 3, 40960, 40000\n1, GPU-bbb, GPU, 97, 40960, 1200\n"
    monkeypatch.setattr(limits.shutil, "which", lambda name: "/usr/bin/" + name)
    monkeypatch.setattr(
        limits.subprocess, "run", lambda *a, **kw: types.SimpleNamespace(returncode=0, stdout=csv)
    )
    monkeypatch.setattr(usage, "gpu_compute_apps", lambda: [(31274, "GPU-bbb", 38000)])
    probe = limits.query_gpus()
    assert probe.problem is None
    assert probe.states[0].processes == ()
    assert probe.states[1].processes == (limits.GpuProcess(pid=31274, memory_mb=38000),)
    assert probe.states[0].busy is False and probe.states[1].busy is True


def test_a_card_held_at_low_utilisation_still_counts_as_in_use():
    """A card with a process on it is busy even at near 0% utilisation."""
    assert _card(0, utilisation=1, processes=[limits.GpuProcess(pid=9, memory_mb=30000)]).busy


def test_a_busy_card_with_no_visible_process_still_counts_as_in_use():
    """High utilisation marks a card busy even with no listed process, as in a container."""
    assert _card(0, utilisation=90, processes=[]).busy
    assert not _card(0, utilisation=3, processes=[]).busy


# --- the decisions ------------------------------------------------------------------------


def test_an_index_above_the_device_count_is_refused_and_says_how_many_there_are(monkeypatch):
    _probe(monkeypatch, _card(0), _card(1))
    report = limits.gpu_preflight("cuda:3")
    assert report.ok is False
    assert "this host does not have" in report.message
    assert "2 cards, cuda:0 to cuda:1" in report.message
    assert "no GPU in it" in report.message


def test_a_free_card_is_reported_in_one_line_and_starts(monkeypatch):
    _probe(monkeypatch, _card(1, utilisation=0, free_mb=40000))
    report = limits.gpu_preflight("cuda:1")
    assert report.ok and not report.contended
    assert report.message.count("\n") == 0
    assert "GPU preflight: cuda:1 (NVIDIA GPU)" in report.message
    assert "utilisation 0%" in report.message
    assert "39.1 GB of 40.0 GB free" in report.message
    assert "no other process on it" in report.message


def test_a_busy_card_warns_with_the_cost_of_sharing_it_spelled_out(monkeypatch):
    _probe(monkeypatch, _card(0), _busy_card(1))
    report = limits.gpu_preflight("cuda:1")
    assert report.ok is True and report.contended is True
    assert "pid 31274 holding 37.1 GB" in report.message
    assert "wall time inflates" in report.message
    assert "cannot be separated out" in report.message
    assert "--gpu-strict" in report.message
    assert "gpu_strict" in report.message


def test_the_strict_flag_turns_the_same_contention_into_a_refusal(monkeypatch):
    _probe(monkeypatch, _busy_card(1))
    report = limits.gpu_preflight("cuda:1", strict=True)
    assert report.ok is False and report.contended is True
    assert "Refusing to start" in report.message
    assert "wall time inflates" in report.message


def test_strict_does_not_refuse_a_free_card(monkeypatch):
    """Strict mode still starts on a free card."""
    _probe(monkeypatch, _card(1))
    assert limits.gpu_preflight("cuda:1", strict=True).ok is True


def test_bare_cuda_with_every_card_busy_says_what_will_happen_instead(monkeypatch):
    """With every card busy, bare `--device cuda` warns that it will share one (no CPU fallback)."""
    _probe(monkeypatch, _busy_card(0), _busy_card(1))
    report = limits.gpu_preflight("cuda")
    assert report.ok is True and report.contended is True
    assert "all 2 cards in use (cuda:0 at 94%, cuda:1 at 94%)" in report.message
    assert "does not fall back to the CPU" in report.message
    assert "shares it" in report.message
    assert limits.gpu_preflight("cuda", strict=True).ok is False


def test_bare_cuda_names_the_free_cards_and_that_nothing_pins_one(monkeypatch):
    _probe(monkeypatch, _busy_card(0), _card(1), _card(2))
    report = limits.gpu_preflight("cuda")
    assert report.ok and not report.contended
    assert "2 of 3 cards free (cuda:1, cuda:2)" in report.message
    assert "pins nothing" in report.message


def test_no_driver_is_said_plainly_and_does_not_stop_the_step(monkeypatch):
    monkeypatch.setattr(limits.shutil, "which", lambda name: None)
    report = limits.gpu_preflight("cuda:1")
    assert report.ok is True and report.contended is False
    assert "no nvidia-smi on PATH" in report.message
    assert "Starting anyway" in report.message


def test_a_driver_that_does_not_answer_is_said_plainly_too(monkeypatch):
    monkeypatch.setattr(limits.shutil, "which", lambda name: "/usr/bin/" + name)

    def explode(*args, **kwargs):
        raise OSError("driver not loaded")

    monkeypatch.setattr(limits.subprocess, "run", explode)
    report = limits.gpu_preflight("cuda:1")
    assert report.ok is True
    assert "nvidia-smi did not answer" in report.message


def test_a_cpu_run_probes_nothing_at_all(monkeypatch):
    """CPU devices and a malformed index do not query the driver."""

    def forbidden():
        raise AssertionError("the driver was queried for a CPU run")

    monkeypatch.setattr(limits, "query_gpus", forbidden)
    for device in ("cpu", "", "auto"):
        report = limits.gpu_preflight(device)
        assert report.ok and report.message == ""
    assert limits.gpu_preflight("cuda:x").message == ""


# --- the strict setting comes from either the flag or the config ---------------------------


def test_gpu_strict_reads_the_config_and_is_beaten_by_an_override():
    try:
        assert limits.gpu_strict(limits.LimitsConfig()) is False
        assert limits.gpu_strict(limits.LimitsConfig(gpu_strict=True)) is True
        # A flag passes the value as a string.
        limits.set_override("gpu_strict", "true")
        assert limits.gpu_strict(limits.LimitsConfig()) is True
    finally:
        limits.clear_overrides()


def test_gpu_strict_is_a_validated_limits_key():
    assert "gpu_strict" in limits.LIMITS_KEYS
    assert limits.config_from_table({"gpu_strict": True}).gpu_strict is True
    assert limits.config_from_table({"gpu_strict": "no"}).gpu_strict is False
    with pytest.raises(limits.LimitsError):
        limits.config_from_table({"gpu_strict": "maybe"})


# --- wired into the steps that can run for hours -------------------------------------------


class _Args:
    def __init__(self, output, device, strict_gpu=False, threads=4, min_len=1000):
        self.output = output
        self.device = device
        self.strict_gpu = strict_gpu
        self.threads = threads
        self.min_len = min_len


def _capture(monkeypatch):
    calls = []
    monkeypatch.setattr(binning, "run", lambda cmd, **kw: calls.append((cmd, kw)))
    monkeypatch.setattr(binning, "resolve_device", lambda device, env: device)
    monkeypatch.setattr(binning, "_adopt_bins", lambda *a, **kw: None)
    return calls


@pytest.mark.parametrize("binner", ("_run_semibin2", "_run_vamb", "_run_comebin"))
def test_every_gpu_binner_checks_the_card_before_it_starts(binner, tmp_path, monkeypatch):
    """Each GPU binner runs the preflight; a missing card index stops it before any command."""
    calls = _capture(monkeypatch)
    monkeypatch.setattr(binning, "_comebin_assembly", lambda args, work, assembly: assembly)
    _probe(monkeypatch, _card(0), _card(1))
    work = tmp_path / "work"
    work.mkdir()
    (work / "sample_0.bam").write_text("")
    with pytest.raises(SystemExit):
        getattr(binning, binner)(_Args(str(tmp_path / "out"), "cuda:3"), "e", str(work), "asm.fa")
    assert calls == []


def test_a_shared_card_lets_the_step_run_unless_strict(tmp_path, monkeypatch):
    calls = _capture(monkeypatch)
    _probe(monkeypatch, _busy_card(3))
    work = tmp_path / "work"
    work.mkdir()
    (work / "sample_0.bam").write_text("")
    binning._run_semibin2(_Args(str(tmp_path / "out"), "cuda:3"), "e", str(work), "asm.fa")
    assert len(calls) == 1
    assert calls[0][1]["extra_env"] == {"CUDA_VISIBLE_DEVICES": "3"}

    calls.clear()
    with pytest.raises(SystemExit):
        binning._run_semibin2(
            _Args(str(tmp_path / "out"), "cuda:3", strict_gpu=True), "e", str(work), "asm.fa"
        )
    assert calls == []


def test_a_cpu_step_does_not_touch_the_driver(tmp_path, monkeypatch):
    calls = _capture(monkeypatch)

    def forbidden():
        raise AssertionError("the driver was queried for a --device cpu step")

    monkeypatch.setattr(limits, "query_gpus", forbidden)
    work = tmp_path / "work"
    work.mkdir()
    (work / "sample_0.bam").write_text("")
    binning._run_semibin2(_Args(str(tmp_path / "out"), "cpu"), "e", str(work), "asm.fa")
    assert len(calls) == 1


def test_a_dry_run_makes_no_claim_about_this_host(tmp_path, monkeypatch):
    """`--dry-run` does not query the driver."""
    calls = _capture(monkeypatch)
    monkeypatch.setattr(binning, "dry_run", lambda: True)

    def forbidden():
        raise AssertionError("the driver was queried under --dry-run")

    monkeypatch.setattr(limits, "query_gpus", forbidden)
    binning.gpu_preflight(_Args(str(tmp_path / "out"), "cuda:3"), "cuda:3")
    assert calls == []
