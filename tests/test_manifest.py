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


def test_a_slurm_job_id_is_recorded_without_overwriting_step_state(tmp_path):
    """Recording a scheduler job preserves the step's command and inputs.

    A Slurm job may outlive its controller. The job ID must survive a manifest reload so
    --resume can reconnect instead of submitting a duplicate job.
    """
    from metawrap2.manifest import Manifest

    path = tmp_path / "manifest.json"
    reads = tmp_path / "reads.fastq"
    reads.write_text("@read\nACGT\n+\nIIII\n")

    manifest = Manifest(str(path))
    command = ["python", "-m", "metawrap2.cli", "read_qc"]

    started = manifest.start(
        "read_qc:sample1",
        command,
        [str(reads)],
        grant={"threads": 8, "executor": "slurm"},
    )

    manifest.record_scheduler_job(
        "read_qc:sample1",
        "slurm",
        "123456",
    )

    reloaded = Manifest(str(path))
    record = reloaded.step("read_qc:sample1")

    assert record is not None
    assert record.status == "running"
    assert record.data["executor"] == "slurm"
    assert record.data["scheduler_job_id"] == "123456"
    assert record.data["command"] == command
    assert record.inputs[0]["path"] == str(reads)

    manifest.finish(
        "read_qc:sample1",
        started,
        "completed",
        outputs=[str(reads)],
    )

    finished = Manifest(str(path)).step("read_qc:sample1")
    assert finished is not None
    assert finished.status == "completed"
    assert finished.data["scheduler_job_id"] == "123456"
