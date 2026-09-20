"""Tests for `metawrap2 status` and the manifest audit it reports.

Both exist to read back information the manifest already had but nothing surfaced. The column
worth testing is the last one: "what would --resume do?" must agree with what --resume actually
does, which means it has to come from the same ``can_skip`` call and not a reimplementation.
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


# --- the table -----------------------------------------------------------------------------


def test_the_table_lists_every_step_with_its_status(study, capsys):
    assert status_cmd.main([str(study)]) == 0
    out = capsys.readouterr().out
    for name in ("assembly", "binning", "kraken2", "quant_bins"):
        assert name in out
    assert "FAILED" in out
    assert "3 of 4 step(s) completed" in out


def test_intact_steps_would_be_skipped_and_broken_ones_redone(study, capsys):
    """The resume column must be the decision --resume makes, with the same reason attached."""
    status_cmd.main([str(study)])
    out = capsys.readouterr().out
    rows = {line.split()[0]: line for line in out.splitlines() if line.startswith("  ")}
    assert "would skip" in rows["assembly"]
    assert "would skip" in rows["binning"]
    assert "would re-run" in rows["kraken2"]  # its output was deleted
    assert "is gone" in rows["kraken2"]  # and it says why
    assert "would re-run" in rows["quant_bins"]
    assert "last attempt failed" in rows["quant_bins"]


# --- the audit -----------------------------------------------------------------------------


def test_audit_reports_a_missing_output_and_exits_nonzero(study, capsys):
    """A step reported success and its output is gone: that is a real finding, so exit 1."""
    assert status_cmd.main([str(study), "--audit"]) == 1
    out = capsys.readouterr().out
    assert "MISSING" in out and "kraken.out" in out
