"""`--device cuda:N` has to put the training on card N, and two of three binners could not.

The defect, found on a multi-GPU host while running GPU arms: SemiBin2 is told an accelerator
*kind* (`--engine gpu`) and VAMB a bare `--cuda`, so both collapsed `cuda:1` to "use a GPU" and
torch took device 0. A 1.5-hour SemiBin2 arm trained on cuda:0 - contending with an unrelated job
there - while the sampler watching cuda:1 recorded 3 MiB and 0% utilisation for the whole run.

It escaped because it is invisible unless a multi-GPU host has something else on device 0: on one
GPU, and on an idle machine, collapsing `cuda:1` to `cuda:0` gives correct bins and plausible
timings. So these assert on the environment each command is given, which needs no GPU at all.
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
    # resolve_device would otherwise probe torch in a conda env that is not installed here.
    monkeypatch.setattr(binning, "resolve_device", lambda device, env: device)
    monkeypatch.setattr(binning, "_adopt_bins", lambda *a, **kw: None)
    return calls


def test_gpu_pin_only_fires_for_an_explicit_index():
    assert binning.gpu_pin("cuda:7") == {"CUDA_VISIBLE_DEVICES": "7"}
    assert binning.gpu_pin("cuda:0") == {"CUDA_VISIBLE_DEVICES": "0"}
    # "any GPU" needs no pin, and pinning one would be a different instruction.
    assert binning.gpu_pin("cuda") is None
    assert binning.gpu_pin("cpu") is None
    # Nothing is guessed from a malformed index; the tool reports it in its own words.
    assert binning.gpu_pin("cuda:x") is None
    assert binning.gpu_pin("cuda:") is None


def test_pin_carries_only_the_one_variable():
    """Not a copy of this process's environment, which `extra_env` would then override with."""
    assert list(binning.gpu_pin("cuda:3")) == ["CUDA_VISIBLE_DEVICES"]


def test_semibin2_is_pinned_to_the_card_asked_for(tmp_path, monkeypatch):
    calls = _capture(monkeypatch)
    work = tmp_path / "work"
    work.mkdir()
    (work / "sample_0.bam").write_text("")
    binning._run_semibin2(_Args(str(tmp_path / "out"), "cuda:5"), "e", str(work), "asm.fa")
    cmd, kwargs = calls[0]
    assert kwargs["extra_env"] == {"CUDA_VISIBLE_DEVICES": "5"}
    # The flag itself still only says *whether*, which is the whole reason the pin is needed.
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
    """COMEBin does accept an index, and must not be given the physical one once pinned.

    Inside a one-device view `cuda:7` is out of range, so the pin and the flag would contradict
    each other. All three binners go through the same mechanism; COMEBin gets `-d cuda`, which is
    its own default.
    """
    calls = _capture(monkeypatch)
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.setattr(binning, "_comebin_assembly", lambda args, work, assembly: assembly)
    binning._run_comebin(_Args(str(tmp_path / "out"), "cuda:7"), "e", str(work), "asm.fa")
    cmd, kwargs = calls[0]
    assert kwargs["extra_env"] == {"CUDA_VISIBLE_DEVICES": "7"}
    assert "-d cuda " in cmd + " " and "cuda:7" not in cmd


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
    """Pinning a CPU run would be a claim about hardware the run was told not to use."""
    calls = _capture(monkeypatch)
    work = tmp_path / "work"
    work.mkdir()
    (work / "sample_0.bam").write_text("")
    binning._run_semibin2(_Args(str(tmp_path / "out"), "cpu"), "e", str(work), "asm.fa")
    assert calls[0][1]["extra_env"] is None
