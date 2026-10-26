"""Tests for metawrap2.manifest and --resume.

A step re-runs when its outputs are gone or its inputs have changed.
"""

from __future__ import annotations

import os

from metawrap2 import manifest as mf


def write(path, text="hello"):
    with open(path, "w") as fh:
        fh.write(text)
    return str(path)


# --- can_skip ------------------------------------------------------------------------------


def test_can_skip_an_intact_completed_step(tmp_path):
    src = write(tmp_path / "in.fa", ">c\nACGT\n")
    out = write(tmp_path / "out.fa", ">c\nACGT\n")
    man = mf.Manifest(str(tmp_path / "m.json"))
    man.finish("binning", man.start("binning", ["c"], [src]), mf.COMPLETED, outputs=[out])
    skip, reason = man.can_skip("binning")
    assert skip is True
    assert "outputs intact" in reason and "inputs unchanged" in reason

    # Changed input: the step must re-run so results from two inputs are not mixed.
    write(tmp_path / "in.fa", ">c1\nACGTACGTACGT\n")  # re-assembled
    skip, reason = man.can_skip("binning")
    assert skip is False and "has changed since this step ran" in reason

    # Missing output: the step must re-run.
    man.finish("binning", man.start("binning", ["c"], [src]), mf.COMPLETED, outputs=[out])
    assert man.can_skip("binning")[0] is True
    os.remove(out)
    skip, reason = man.can_skip("binning")
    assert skip is False and "is gone" in reason
