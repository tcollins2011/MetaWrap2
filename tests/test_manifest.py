"""Tests for metawrap2.manifest: --resume verifying rather than assuming.

The old marker-file resume was wrong in two directions, and both are pinned here: a step whose
outputs are gone must re-run, and a step whose inputs have changed must re-run. Everything else
in this file supports those two.
"""

from __future__ import annotations

import os

from metawrap2 import manifest as mf


def write(path, text="hello"):
    with open(path, "w") as fh:
        fh.write(text)
    return str(path)


# --- can_skip: the whole point ------------------------------------------------------------


def test_can_skip_an_intact_completed_step(tmp_path):
    src = write(tmp_path / "in.fa", ">c\nACGT\n")
    out = write(tmp_path / "out.fa", ">c\nACGT\n")
    man = mf.Manifest(str(tmp_path / "m.json"))
    man.finish("binning", man.start("binning", ["c"], [src]), mf.COMPLETED, outputs=[out])
    skip, reason = man.can_skip("binning")
    assert skip is True
    assert "outputs intact" in reason and "inputs unchanged" in reason

    # The second thing marker-file resume got wrong: silently mixing two inputs' results.
    write(tmp_path / "in.fa", ">c1\nACGTACGTACGT\n")  # re-assembled
    skip, reason = man.can_skip("binning")
    assert skip is False and "has changed since this step ran" in reason

    # The first thing marker-file resume got wrong: a step whose output is gone must re-run.
    man.finish("binning", man.start("binning", ["c"], [src]), mf.COMPLETED, outputs=[out])
    assert man.can_skip("binning")[0] is True
    os.remove(out)
    skip, reason = man.can_skip("binning")
    assert skip is False and "is gone" in reason
