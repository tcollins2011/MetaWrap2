"""Tests for `metawrap2 status` and the manifest audit it reports.

The "what would --resume do?" column must agree with what --resume does, so it comes from the
same ``can_skip`` call.
"""

from __future__ import annotations

import pytest

from metawrap2.commands import status_cmd
from metawrap2.manifest import COMPLETED, FAILED, Manifest


@pytest.fixture
def study(tmp_path):
    """A study directory with a manifest: two intact steps, one with a deleted output, one failed."""
    root = tmp_path / "study"
    (root / ".metawrap2" / "run_steps").mkdir(parents=True)

    src = root / "asm.fa"
    src.write_text(">c1\nACGT\n")
    bins = root / "BINS"
    bins.mkdir()
    (bins / "bin.1.fa").write_text(">c1\nACGT\n")
    gone = root / "kraken.out"
    gone.write_text("x")

    man = Manifest(status_cmd.manifest_path(str(root)))
    man.finish(
        "assembly",
        man.start("assembly", ["metawrap2", "assembly"], []),
        COMPLETED,
        outputs=[str(src)],
    )
    man.finish(
        "binning",
        man.start("binning", ["metawrap2", "binning"], [str(src)]),
        COMPLETED,
        outputs=[str(bins)],
        log=str(root / "LOGS" / "binning.log"),
    )
    man.finish(
        "kraken2",
        man.start("kraken2", ["metawrap2", "kraken2"], []),
        COMPLETED,
        outputs=[str(gone)],
    )
    man.finish("quant_bins", man.start("quant_bins", ["metawrap2", "quant_bins"], []), FAILED)
    gone.unlink()
    return root


# --- the table ----------------------------------------------------------------------------


def test_intact_steps_would_be_skipped_and_broken_ones_redone(study, capsys):
    """The resume column must be the decision --resume makes, with the same reason attached."""
    assert status_cmd.main([str(study)]) == 0
    out = capsys.readouterr().out
    assert "3 of 4 step(s) completed" in out
    rows = {line.split()[0]: line for line in out.splitlines() if line.startswith("  ")}
    assert "would skip" in rows["assembly"]
    assert "would skip" in rows["binning"]
    assert "would re-run" in rows["kraken2"]  # its output was deleted
    assert "is gone" in rows["kraken2"]  # and it says why
    assert "would re-run" in rows["quant_bins"]
    assert "last attempt failed" in rows["quant_bins"]

    # A step reported success but its output is gone, so --audit exits 1.
    assert status_cmd.main([str(study), "--audit"]) == 1
    out = capsys.readouterr().out
    assert "MISSING" in out and "kraken.out" in out
