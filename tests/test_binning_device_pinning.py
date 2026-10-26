"""Tests that `--device cuda:N` runs each GPU binner on device N.

SemiBin2 takes only an accelerator kind (`--engine gpu`) and VAMB a bare `--cuda`, so without
``CUDA_VISIBLE_DEVICES`` both run on device 0 whatever index was asked for. The results are still
correct, so the error only shows as contention on a shared multi-GPU host. The tests check the
environment passed to each command and need no GPU.
"""

from metawrap2.modules import binning


class _Args:
    def __init__(self, output, device, threads=4, min_len=1000):
        self.output = output
        self.device = device
        self.threads = threads
        self.min_len = min_len


def _capture(monkeypatch):
    calls = []
    monkeypatch.setattr(binning, "run", lambda cmd, **kw: calls.append((cmd, kw)))
    # Avoid resolve_device probing torch in a conda env that is not installed here.
    monkeypatch.setattr(binning, "resolve_device", lambda device, env: device)
    # Skip the preflight so the result does not depend on the GPUs of the machine running the
    # tests; tests/test_gpu_preflight.py covers it.
    monkeypatch.setattr(binning, "gpu_preflight", lambda args, device: None)
    monkeypatch.setattr(binning, "_adopt_bins", lambda *a, **kw: None)
    return calls


def test_gpu_pin_only_fires_for_an_explicit_index():
    assert binning.gpu_pin("cuda:2") == {"CUDA_VISIBLE_DEVICES": "2"}
    assert binning.gpu_pin("cuda:0") == {"CUDA_VISIBLE_DEVICES": "0"}
    # "any GPU" is not pinned.
    assert binning.gpu_pin("cuda") is None
    assert binning.gpu_pin("cpu") is None
    # A malformed index is not pinned; the tool reports the error.
    assert binning.gpu_pin("cuda:x") is None
    assert binning.gpu_pin("cuda:") is None


def test_pin_carries_only_the_one_variable():
    """The pin holds only CUDA_VISIBLE_DEVICES, not a copy of this process's environment."""
    assert list(binning.gpu_pin("cuda:3")) == ["CUDA_VISIBLE_DEVICES"]


def test_semibin2_is_pinned_to_the_card_asked_for(tmp_path, monkeypatch):
    calls = _capture(monkeypatch)
    work = tmp_path / "work"
    work.mkdir()
    (work / "sample_0.bam").write_text("")
    binning._run_semibin2(_Args(str(tmp_path / "out"), "cuda:1"), "e", str(work), "asm.fa")
    cmd, kwargs = calls[0]
    assert kwargs["extra_env"] == {"CUDA_VISIBLE_DEVICES": "1"}
    # The flag only selects GPU use; the device comes from the pin.
    assert "--engine gpu" in cmd


def test_vamb_is_pinned_to_the_card_asked_for(tmp_path, monkeypatch):
    calls = _capture(monkeypatch)
    work = tmp_path / "work"
    work.mkdir()
    (work / "sample_0.bam").write_text("")
    binning._run_vamb(_Args(str(tmp_path / "out"), "cuda:2"), "e", str(work), "asm.fa")
    cmd, kwargs = calls[0]
    assert kwargs["extra_env"] == {"CUDA_VISIBLE_DEVICES": "2"}
    assert "--cuda" in cmd


def test_comebin_is_pinned_and_then_told_the_card_it_can_see(tmp_path, monkeypatch):
    """COMEBin is pinned and given `-d cuda` (its default), not the physical index.

    With only one device visible, the physical index would be out of range.
    """
    calls = _capture(monkeypatch)
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.setattr(binning, "_comebin_assembly", lambda args, work, assembly: assembly)
    binning._run_comebin(_Args(str(tmp_path / "out"), "cuda:3"), "e", str(work), "asm.fa")
    cmd, kwargs = calls[0]
    assert kwargs["extra_env"] == {"CUDA_VISIBLE_DEVICES": "3"}
    assert "-d cuda " in cmd + " " and "cuda:3" not in cmd


def test_unindexed_cuda_passes_through_to_comebin_unpinned(tmp_path, monkeypatch):
    calls = _capture(monkeypatch)
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.setattr(binning, "_comebin_assembly", lambda args, work, assembly: assembly)
    binning._run_comebin(_Args(str(tmp_path / "out"), "cuda"), "e", str(work), "asm.fa")
    cmd, kwargs = calls[0]
    assert kwargs["extra_env"] is None
    assert "-d cuda " in cmd + " "


def test_a_cpu_run_sets_no_cuda_variable(tmp_path, monkeypatch):
    """A CPU run sets no CUDA_VISIBLE_DEVICES."""
    calls = _capture(monkeypatch)
    work = tmp_path / "work"
    work.mkdir()
    (work / "sample_0.bam").write_text("")
    binning._run_semibin2(_Args(str(tmp_path / "out"), "cpu"), "e", str(work), "asm.fa")
    assert calls[0][1]["extra_env"] is None
