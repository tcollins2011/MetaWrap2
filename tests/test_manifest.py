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


# --- recording ----------------------------------------------------------------------------


def test_start_and_finish_record_a_step(tmp_path):
    src = write(tmp_path / "in.fa", ">c\nACGT\n")
    out = write(tmp_path / "out.fa", ">c\nACGT\n")
    man = mf.Manifest(str(tmp_path / "manifest.json"))
    started = man.start("assembly", ["metawrap2", "assembly"], [src])
    man.finish("assembly", started, mf.COMPLETED, outputs=[out], log="assembly.log")

    record = man.step("assembly")
    assert record is not None
    assert record.status == mf.COMPLETED
    assert record.data["command"] == ["metawrap2", "assembly"]
    assert record.data["log"] == "assembly.log"
    assert "duration_seconds" in record.data
    assert [e["path"] for e in record.inputs] == [src]
    assert [e["path"] for e in record.outputs] == [out]


# --- can_skip: the whole point ------------------------------------------------------------


def test_can_skip_an_intact_completed_step(tmp_path):
    src = write(tmp_path / "in.fa", ">c\nACGT\n")
    out = write(tmp_path / "out.fa", ">c\nACGT\n")
    man = mf.Manifest(str(tmp_path / "m.json"))
    man.finish("binning", man.start("binning", ["c"], [src]), mf.COMPLETED, outputs=[out])
    skip, reason = man.can_skip("binning")
    assert skip is True
    assert "outputs intact" in reason and "inputs unchanged" in reason


def test_cannot_skip_when_an_output_was_deleted(tmp_path):
    """The first thing marker-file resume got wrong."""
    out = write(tmp_path / "out.fa", ">c\nACGT\n")
    man = mf.Manifest(str(tmp_path / "m.json"))
    man.finish("binning", man.start("binning", ["c"], []), mf.COMPLETED, outputs=[out])
    os.remove(out)
    skip, reason = man.can_skip("binning")
    assert skip is False and "is gone" in reason


def test_cannot_skip_when_an_input_changed(tmp_path):
    """The second thing marker-file resume got wrong: silently mixing two inputs' results."""
    src = write(tmp_path / "assembly.fa", ">c1\nACGT\n")
    out = write(tmp_path / "out.fa", ">c\nACGT\n")
    man = mf.Manifest(str(tmp_path / "m.json"))
    man.finish("binning", man.start("binning", ["c"], [src]), mf.COMPLETED, outputs=[out])
    assert man.can_skip("binning")[0] is True
    write(tmp_path / "assembly.fa", ">c1\nACGTACGTACGT\n")  # re-assembled
    skip, reason = man.can_skip("binning")
    assert skip is False and "has changed since this step ran" in reason
