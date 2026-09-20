"""Tests for metawrap2.manifest: --resume verifying rather than assuming.

The old marker-file resume was wrong in two directions, and both are pinned here: a step whose
outputs are gone must re-run, and a step whose inputs have changed must re-run. Everything else
in this file supports those two.
"""

from __future__ import annotations

import json
import os

from metawrap2 import manifest as mf


def write(path, text="hello"):
    with open(path, "w") as fh:
        fh.write(text)
    return str(path)


# --- fingerprint --------------------------------------------------------------------------


def test_fingerprint_is_stable_and_content_sensitive(tmp_path):
    path = write(tmp_path / "a.txt", "abc")
    first = mf.fingerprint(path)
    assert first == mf.fingerprint(path)
    write(tmp_path / "a.txt", "abd")
    assert mf.fingerprint(path) != first


def test_fingerprint_of_a_missing_file_is_none(tmp_path):
    assert mf.fingerprint(str(tmp_path / "nope")) is None


def test_fingerprint_notices_truncation_of_a_large_file(tmp_path):
    """The cheap fingerprint samples both ends, so a truncated file never looks unchanged."""
    path = str(tmp_path / "big.bin")
    with open(path, "wb") as fh:
        fh.write(os.urandom(mf.SAMPLE_BYTES * 3))
    before = mf.fingerprint(path)
    with open(path, "r+b") as fh:
        fh.truncate(mf.SAMPLE_BYTES * 2)
    assert mf.fingerprint(path) != before


def test_whole_file_fingerprint_sees_a_middle_edit_that_sampling_misses(tmp_path):
    """Why --strict-fingerprints exists: an edit between the two sampled ends."""
    path = str(tmp_path / "big.bin")
    body = bytearray(b"\0" * (mf.SAMPLE_BYTES * 3))
    with open(path, "wb") as fh:
        fh.write(body)
    cheap_before = mf.fingerprint(path)
    strict_before = mf.fingerprint(path, whole_file=True)
    with open(path, "r+b") as fh:  # same size, changed in the middle only
        fh.seek(mf.SAMPLE_BYTES + 100)
        fh.write(b"CHANGED")
    assert mf.fingerprint(path) == cheap_before
    assert mf.fingerprint(path, whole_file=True) != strict_before


# --- describe / unchanged -----------------------------------------------------------------


def test_describe_a_file_records_size_and_fingerprint(tmp_path):
    entry = mf.describe(write(tmp_path / "a.txt", "abc"))
    assert entry["kind"] == "file"
    assert entry["size"] == 3
    assert entry["fingerprint"]


def test_describe_a_missing_path(tmp_path):
    assert mf.describe(str(tmp_path / "nope"))["kind"] == "missing"


def test_describe_a_directory_captures_its_listing(tmp_path):
    bins = tmp_path / "bins"
    bins.mkdir()
    write(bins / "bin.1.fa", ">c\nACGT\n")
    before = mf.describe(str(bins))
    assert before["kind"] == "directory" and before["entries"] == 1
    # A re-run of a binner producing a different number of bins must not look unchanged.
    write(bins / "bin.2.fa", ">c\nACGT\n")
    assert mf.describe(str(bins))["fingerprint"] != before["fingerprint"]


def test_unchanged_is_true_only_while_content_is_identical(tmp_path):
    path = write(tmp_path / "a.txt", "abc")
    recorded = mf.describe(path)
    assert mf.unchanged(recorded) is True
    write(tmp_path / "a.txt", "different")
    assert mf.unchanged(recorded) is False


def test_unchanged_is_false_when_the_file_is_deleted(tmp_path):
    path = write(tmp_path / "a.txt", "abc")
    recorded = mf.describe(path)
    os.remove(path)
    assert mf.unchanged(recorded) is False


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


def test_the_manifest_is_reloaded_from_disk(tmp_path):
    path = str(tmp_path / "manifest.json")
    src = write(tmp_path / "in.fa", ">c\nACGT\n")
    man = mf.Manifest(path)
    man.finish("binning", man.start("binning", ["c"], [src]), mf.COMPLETED, outputs=[src])
    assert mf.Manifest(path).step("binning").status == mf.COMPLETED


def test_a_corrupt_manifest_is_ignored_rather_than_fatal(tmp_path):
    """A half-written manifest must cost skipping, not the run."""
    path = str(tmp_path / "manifest.json")
    write(tmp_path / "manifest.json", "{not json")
    man = mf.Manifest(path)
    assert man.step("binning") is None
    assert man.can_skip("binning")[0] is False


def test_a_manifest_from_another_schema_version_is_ignored(tmp_path):
    path = str(tmp_path / "manifest.json")
    with open(path, "w") as fh:
        json.dump(
            {"schema": mf.SCHEMA_VERSION + 99, "steps": {"binning": {"status": "completed"}}}, fh
        )
    assert mf.Manifest(path).step("binning") is None


# --- can_skip: the whole point ------------------------------------------------------------


def test_can_skip_an_intact_completed_step(tmp_path):
    src = write(tmp_path / "in.fa", ">c\nACGT\n")
    out = write(tmp_path / "out.fa", ">c\nACGT\n")
    man = mf.Manifest(str(tmp_path / "m.json"))
    man.finish("binning", man.start("binning", ["c"], [src]), mf.COMPLETED, outputs=[out])
    skip, reason = man.can_skip("binning")
    assert skip is True
    assert "outputs intact" in reason and "inputs unchanged" in reason


def test_cannot_skip_a_step_never_recorded(tmp_path):
    skip, reason = mf.Manifest(str(tmp_path / "m.json")).can_skip("binning")
    assert skip is False and "no record of it" in reason


def test_cannot_skip_a_failed_or_interrupted_step(tmp_path):
    man = mf.Manifest(str(tmp_path / "m.json"))
    for status in (mf.FAILED, mf.INTERRUPTED):
        man.finish(status, man.start(status, ["c"], []), status)
        skip, reason = man.can_skip(status)
        assert skip is False and reason == "last attempt %s" % status


def test_cannot_skip_when_an_output_was_deleted(tmp_path):
    """The first thing marker-file resume got wrong."""
    out = write(tmp_path / "out.fa", ">c\nACGT\n")
    man = mf.Manifest(str(tmp_path / "m.json"))
    man.finish("binning", man.start("binning", ["c"], []), mf.COMPLETED, outputs=[out])
    os.remove(out)
    skip, reason = man.can_skip("binning")
    assert skip is False and "is gone" in reason


def test_cannot_skip_when_an_output_was_emptied(tmp_path):
    """A tool that exits 0 having written nothing must not count as done."""
    out = write(tmp_path / "out.fa", ">c\nACGT\n")
    man = mf.Manifest(str(tmp_path / "m.json"))
    man.finish("binning", man.start("binning", ["c"], []), mf.COMPLETED, outputs=[out])
    open(out, "w").close()
    skip, reason = man.can_skip("binning")
    assert skip is False and "is now empty" in reason


def test_cannot_skip_when_an_output_directory_was_emptied(tmp_path):
    bins = tmp_path / "bins"
    bins.mkdir()
    write(bins / "bin.1.fa", ">c\nACGT\n")
    man = mf.Manifest(str(tmp_path / "m.json"))
    man.finish("binning", man.start("binning", ["c"], []), mf.COMPLETED, outputs=[str(bins)])
    os.remove(bins / "bin.1.fa")
    skip, reason = man.can_skip("binning")
    assert skip is False and "is now empty" in reason


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


def test_summary_lists_every_step(tmp_path):
    man = mf.Manifest(str(tmp_path / "m.json"))
    for name in ("read_qc", "assembly", "binning"):
        man.finish(name, man.start(name, ["c"], []), mf.COMPLETED)
    assert [s["name"] for s in man.summary()] == ["read_qc", "assembly", "binning"]


def test_save_leaves_no_temporary_file_behind(tmp_path):
    man = mf.Manifest(str(tmp_path / "m.json"))
    man.finish("binning", man.start("binning", ["c"], []), mf.COMPLETED)
    assert sorted(os.listdir(tmp_path)) == ["m.json"]
