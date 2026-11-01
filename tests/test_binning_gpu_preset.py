"""Binning presets: `default` and `deep` are listed; older and hidden names still resolve.

`deep` (called `gpu` in 2.1.4 and `slow` before) selects the three-binner set and trains on a GPU
by default. `fast` and `very-fast` are accepted but not listed in --help or in errors.
"""

import argparse
import re

import pytest

from metawrap2.commands.run import selected_binners
from metawrap2.modules import binning

TRIO = ["metabat2", "semibin2", "comebin"]


def _args(preset=None, device=None, **named):
    values = {b: False for b in binning.BINNERS}
    values.update(named)
    return argparse.Namespace(preset=preset, device=device, **values)


def _selected(args):
    return [b for b in binning.BINNERS if getattr(args, b)]


def test_deep_selects_the_three_binners_and_cuda():
    args = _args("deep")
    binning._apply_preset(args)
    assert _selected(args) == TRIO
    assert args.device == "cuda"


@pytest.mark.parametrize("alias", ["gpu", "slow"])
def test_gpu_and_slow_are_silent_aliases_of_deep(alias, capsys):
    args = _args(alias)
    binning._apply_preset(args)
    assert args.preset == "deep" and args.device == "cuda"
    assert _selected(args) == TRIO
    out = capsys.readouterr().out
    assert "old name" not in out and alias not in out.replace("--preset deep", "")


@pytest.mark.parametrize(
    "preset, expected",
    [
        ("default", ["metabat2", "semibin2"]),
        ("fast", ["metabat2", "maxbin2"]),
        ("very-fast", ["metabat2"]),
    ],
)
def test_cpu_presets_select_their_binners(preset, expected):
    args = _args(preset)
    binning._apply_preset(args)
    assert _selected(args) == expected
    assert args.device == binning.DEFAULT_LEARNED_DEVICE


def test_an_explicit_device_wins():
    args = _args("deep", device="cpu")
    binning._apply_preset(args)
    assert args.device == "cpu"


def test_named_binners_ignore_the_preset_and_its_device():
    args = _args("deep", metabat2=True)
    binning._apply_preset(args)
    assert args.device == binning.DEFAULT_LEARNED_DEVICE and not args.comebin


def test_study_sheets_resolve_every_accepted_name():
    for name in ("deep", "gpu", "slow"):
        assert selected_binners(["--preset", name]) == TRIO
    assert selected_binners(["--preset", "fast"]) == ["metabat2", "maxbin2"]
    assert selected_binners(["--preset", "very-fast"]) == ["metabat2"]


def test_the_parser_accepts_hidden_and_alias_names():
    base = ["-a", "a.fa", "-o", "out", "r_1.fastq", "r_2.fastq", "--preset"]
    for name in ("default", "deep", "gpu", "slow", "fast", "very-fast"):
        assert binning._parse_args(base + [name]).preset == name


def test_help_lists_only_default_and_deep(capsys):
    with pytest.raises(SystemExit):
        binning._parse_args(["-h"])
    out = capsys.readouterr().out
    assert "--preset {default,deep}" in out
    start = out.index("--preset {default,deep}")
    block = out[start : out.index("--metabat2", start)]
    for hidden in ("very-fast", "fast", "gpu", "slow"):
        assert not re.search(r"(?<![\w-])%s(?![\w-])" % re.escape(hidden), block), hidden
    assert "very-fast" not in out


def test_an_unknown_preset_lists_only_default_and_deep(capsys):
    with pytest.raises(SystemExit):
        binning._parse_args(["-a", "a.fa", "-o", "out", "r_1.fastq", "--preset", "nope"])
    error_line = capsys.readouterr().err.strip().splitlines()[-1]
    assert "(choose from default, deep)" in error_line
    for hidden in ("very-fast", "fast", "gpu", "slow"):
        assert hidden not in error_line


def test_no_binner_error_lists_only_default_and_deep(capsys):
    with pytest.raises(SystemExit):
        binning._apply_preset(_args())
    text = "".join(capsys.readouterr())
    assert "--preset {default,deep}" in text
    assert "very-fast" not in text and "gpu}" not in text
