"""`--preset gpu` (formerly `slow`) selects the three-binner set and trains on a GPU by default."""

import argparse

from metawrap2.commands.run import selected_binners
from metawrap2.modules import binning


def _args(preset=None, device=None, **named):
    values = {b: False for b in binning.BINNERS}
    values.update(named)
    return argparse.Namespace(preset=preset, device=device, **values)


def test_gpu_selects_the_three_binners_and_cuda():
    args = _args("gpu")
    binning._apply_preset(args)
    assert [b for b in binning.BINNERS if getattr(args, b)] == ["metabat2", "semibin2", "comebin"]
    assert args.device == "cuda"


def test_slow_is_accepted_as_the_old_name(capsys):
    args = _args("slow")
    binning._apply_preset(args)
    assert args.preset == "gpu" and args.device == "cuda"
    assert "old name" in capsys.readouterr().out


def test_an_explicit_device_wins():
    args = _args("gpu", device="cpu")
    binning._apply_preset(args)
    assert args.device == "cpu"


def test_other_presets_keep_the_cpu_default():
    for preset in ("very-fast", "fast", "default"):
        args = _args(preset)
        binning._apply_preset(args)
        assert args.device == binning.DEFAULT_LEARNED_DEVICE


def test_named_binners_ignore_the_preset_and_its_device():
    args = _args("gpu", metabat2=True)
    binning._apply_preset(args)
    assert args.device == binning.DEFAULT_LEARNED_DEVICE and not args.comebin


def test_study_sheets_resolve_the_old_name_too():
    assert selected_binners(["--preset", "slow"]) == selected_binners(["--preset", "gpu"])
    assert selected_binners(["--preset", "gpu"]) == ["metabat2", "semibin2", "comebin"]
