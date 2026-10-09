"""Choosing binners by a compute tier instead of naming them.

The tiers exist because the ratio that used to pick defaults - MAGs per CPU-hour - was rejected: it is
dominated by whichever arm is cheapest rather than by which recovers the most genomes, and following it
produced a default that recovers fewer MAGs than metaWRAP 1 does. A tier puts the cost decision where it
belongs, with the person who knows their compute budget.
"""

from __future__ import annotations

import pytest

from metawrap2.modules import binning


def _args(argv):
    return binning._parse_args(argv)


BASE = ["-a", "a.fa", "-o", "out", "r_1.fastq", "r_2.fastq"]


def test_the_tiers_are_ordered_cheapest_first_and_are_not_supersets():
    """The tiers used to be nested, on the theory that stepping up can only add genomes.

    The theory was tested and withdrawn. `slow` was metaBAT2 + SemiBin2 + COMEBin on a measurement that
    turned out to be comparing against a diluted figure: COMEBin's clean-MAG count had been estimated as
    178 HQ x an 88.8% GUNC pass rate measured over all 376 bins it emitted, 198 of which were never MAG
    candidates. Counted per bin, 176 of its 178 HQ bins pass GUNC, so the three-binner ensemble's real
    margin was one bin for 33% more cpu. `slow` is COMEBin alone.

    What survives is the ordering: each tier costs more than the one below and recovers more clean MAGs
    (35, 59, 114, 176). Nesting does not.
    """
    tiers = ["very-fast", "fast", "default", "slow"]
    assert list(binning.PRESETS) == tiers, "declaration order is the order --help shows"
    assert binning.PRESETS["very-fast"] == ("metabat2",)
    assert binning.PRESETS["slow"] == ("comebin",), "one binner, and the cleanest set measured"
    assert not set(binning.PRESETS["default"]) <= set(
        binning.PRESETS["slow"]
    ), "no longer nested, and a test that asserted it was would hide the measurement that changed it"


def test_every_preset_names_binners_the_module_can_actually_run():
    """A tier is a promise the parser has to be able to keep."""
    for tier, binners in binning.PRESETS.items():
        args = binning._parse_args(BASE + ["--preset", tier])
        binning._apply_preset(args)
        for name in binners:
            assert getattr(args, name.replace("-", "_")), "%s: --%s never got set" % (tier, name)


def test_fast_to_default_swaps_maxbin2_for_semibin2_rather_than_adding():
    """The one break in the chain, pinned so it is a decision rather than an oversight.

    MaxBin2 is the dirtiest input measured - 37.7% of its bins pass GUNC - so carrying it into the
    higher tiers would add chimeras for consolidation to sort out. SemiBin2 replaces it instead.
    """
    assert "maxbin2" in binning.PRESETS["fast"]
    assert "maxbin2" not in binning.PRESETS["default"]
    assert "maxbin2" not in binning.PRESETS["slow"]


def test_every_preset_names_only_real_binners():
    for tier, binners in binning.PRESETS.items():
        for b in binners:
            assert b in binning.BINNERS, "%s lists unknown binner %r" % (tier, b)


def test_a_preset_selects_its_binners():
    args = _args(BASE + ["--preset", "slow"])
    for b in binning.PRESETS["slow"]:
        assert getattr(args, b) is False, "the parser must not pre-apply it"


def test_naming_binners_yourself_wins_over_a_preset(caplog, monkeypatch, tmp_path):
    """An explicit choice must not be silently widened - a benchmark arm depends on getting exactly it."""
    said = []
    monkeypatch.setattr(binning, "comm", lambda msg, *a, **k: said.append(msg))
    monkeypatch.setattr(
        binning, "error", lambda msg, *a, **k: pytest.fail("should not error: %s" % msg)
    )

    args = _args(BASE + ["--preset", "slow", "--metabat2"])
    binning._apply_preset(args)
    assert args.metabat2 is True
    assert args.comebin is False, "the preset must not add binners to an explicit selection"
    assert any("ignored" in m for m in said)


def test_a_preset_applies_when_no_binner_is_named(monkeypatch):
    said = []
    monkeypatch.setattr(binning, "comm", lambda msg, *a, **k: said.append(msg))
    args = _args(BASE + ["--preset", "default"])
    binning._apply_preset(args)
    assert args.metabat2 is True and args.semibin2 is True
    assert args.comebin is False and args.maxbin2 is False
    assert any("selects" in m for m in said)


def test_naming_nothing_at_all_is_still_an_error_that_mentions_the_tiers(monkeypatch):
    raised = []
    monkeypatch.setattr(binning, "error", lambda msg, *a, **k: raised.append(msg))
    monkeypatch.setattr(binning, "comm", lambda *a, **k: None)
    args = _args(BASE)
    binning._apply_preset(args)
    assert raised and "--preset" in raised[0]


def test_binning_and_coverage_agree_on_the_identity_floor(tmp_path):
    """Two modules write the same table; a floor that differs between them is a silent divergence."""
    from metawrap2 import depth as _depth
    from metawrap2.modules import coverage as _coverage

    b = binning._parse_args(["-a", "a.fa", "-o", str(tmp_path), "--metabat2", "r_1.fastq"])
    c = _coverage._parse_args(["-a", "a.fa", "-o", str(tmp_path), "r_1.fastq", "r_2.fastq"])
    assert b.min_identity == c.min_identity == _depth.PERCENT_IDENTITY


def test_the_identity_floor_can_be_put_back_to_jgis(tmp_path):
    from metawrap2 import depth as _depth

    args = binning._parse_args(
        ["-a", "a.fa", "-o", str(tmp_path), "--metabat2", "--min-identity", "97", "r_1.fastq"]
    )
    assert args.min_identity == _depth.JGI_PERCENT_IDENTITY
