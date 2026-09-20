"""Tests for `metawrap2 status` and the manifest audit it reports.

Both exist to read back information the manifest already had but nothing surfaced. The column
worth testing is the last one: "what would --resume do?" must agree with what --resume actually
does, which means it has to come from the same ``can_skip`` call and not a reimplementation.
"""

from __future__ import annotations

import json
import os

import pytest

from metawrap2.commands import status_cmd
from metawrap2.manifest import COMPLETED, FAILED, INTERRUPTED, Manifest


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


# --- finding the manifest ------------------------------------------------------------------


def test_manifest_path_from_a_study_directory(tmp_path):
    path = status_cmd.manifest_path(str(tmp_path / "study"))
    assert path.endswith(os.path.join(".metawrap2", "run_steps", "manifest.json"))


def test_manifest_path_accepts_a_manifest_file_directly(tmp_path):
    direct = tmp_path / "manifest.json"
    direct.write_text("{}")
    assert status_cmd.manifest_path(str(direct)) == str(direct)


def test_it_explains_itself_when_there_is_no_manifest(tmp_path, capsys):
    assert status_cmd.main([str(tmp_path)]) == 1
    out = capsys.readouterr().out
    assert "No MetaWrap2 run manifest" in out
    assert "before manifests existed" in out  # the other plausible cause


def test_an_empty_manifest_is_not_an_error(tmp_path, capsys):
    root = tmp_path / "study"
    (root / ".metawrap2" / "run_steps").mkdir(parents=True)
    Manifest(status_cmd.manifest_path(str(root))).save()
    assert status_cmd.main([str(root)]) == 0
    assert "records no steps yet" in capsys.readouterr().out


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


def test_an_interrupted_step_would_be_redone(tmp_path, capsys):
    root = tmp_path / "study"
    (root / ".metawrap2" / "run_steps").mkdir(parents=True)
    man = Manifest(status_cmd.manifest_path(str(root)))
    man.finish("assembly", man.start("assembly", ["c"], []), INTERRUPTED)
    status_cmd.main([str(root)])
    out = capsys.readouterr().out
    assert "interrupted" in out and "would re-run" in out


# --- one step in full ----------------------------------------------------------------------


def test_step_shows_inputs_outputs_and_the_command(study, capsys):
    assert status_cmd.main([str(study), "--step", "binning"]) == 0
    out = capsys.readouterr().out
    assert "metawrap2 binning" in out
    assert "asm.fa" in out  # its input
    assert "BINS" in out  # its output
    assert "binning.log" in out
    assert "assembly" not in out.replace("metawrap2 binning", "")  # only the asked-for step


def test_an_unknown_step_lists_the_known_ones(study, capsys):
    assert status_cmd.main([str(study), "--step", "nope"]) == 1
    out = capsys.readouterr().out
    assert "No step named 'nope'" in out
    assert "binning" in out


def test_verbose_shows_every_step_in_full(study, capsys):
    assert status_cmd.main([str(study), "--verbose"]) == 0
    out = capsys.readouterr().out
    assert out.count("status     ") == 4


def test_a_directory_output_is_summarised_by_its_entry_count(study, capsys):
    status_cmd.main([str(study), "--step", "binning"])
    assert "directory, 1 entry" in capsys.readouterr().out


# --- json ----------------------------------------------------------------------------------


def test_json_prints_the_manifest(study, capsys):
    assert status_cmd.main([str(study), "--json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert set(data["steps"]) == {"assembly", "binning", "kraken2", "quant_bins"}


def test_json_for_one_step(study, capsys):
    assert status_cmd.main([str(study), "--step", "binning", "--json"]) == 0
    rows = json.loads(capsys.readouterr().out)
    assert [r["name"] for r in rows] == ["binning"]


# --- the audit -----------------------------------------------------------------------------


def test_audit_reports_a_missing_output_and_exits_nonzero(study, capsys):
    """A step reported success and its output is gone: that is a real finding, so exit 1."""
    assert status_cmd.main([str(study), "--audit"]) == 1
    out = capsys.readouterr().out
    assert "MISSING" in out and "kraken.out" in out


def test_audit_is_clean_when_everything_is_present(tmp_path, capsys):
    root = tmp_path / "study"
    (root / ".metawrap2" / "run_steps").mkdir(parents=True)
    out_file = root / "final_assembly.fasta"
    out_file.write_text(">c1\nACGT\n")
    man = Manifest(status_cmd.manifest_path(str(root)))
    man.finish("assembly", man.start("assembly", ["c"], []), COMPLETED, outputs=[str(out_file)])
    assert status_cmd.main([str(root), "--audit"]) == 0
    assert "Every recorded output is present and non-empty" in capsys.readouterr().out


def test_audit_reports_an_emptied_output(tmp_path, capsys):
    root = tmp_path / "study"
    (root / ".metawrap2" / "run_steps").mkdir(parents=True)
    out_file = root / "final_assembly.fasta"
    out_file.write_text(">c1\nACGT\n")
    man = Manifest(status_cmd.manifest_path(str(root)))
    man.finish("assembly", man.start("assembly", ["c"], []), COMPLETED, outputs=[str(out_file)])
    out_file.write_text("")
    assert status_cmd.main([str(root), "--audit"]) == 1
    assert "EMPTY" in capsys.readouterr().out


def test_audit_flags_a_step_that_declared_no_outputs(tmp_path, capsys):
    """--resume can never verify such a step, so it silently degrades to a marker file."""
    root = tmp_path / "study"
    (root / ".metawrap2" / "run_steps").mkdir(parents=True)
    man = Manifest(status_cmd.manifest_path(str(root)))
    man.finish("blobology", man.start("blobology", ["c"], []), COMPLETED)
    assert status_cmd.main([str(root), "--audit"]) == 1
    out = capsys.readouterr().out
    assert "UNVERIFIED" in out and "blobology" in out


def test_audit_notes_an_undeclared_file_without_failing(tmp_path, capsys):
    """Worth naming - it is invisible to --resume - but not a failure on its own."""
    root = tmp_path / "study"
    (root / ".metawrap2" / "run_steps").mkdir(parents=True)
    declared = root / "final_assembly.fasta"
    declared.write_text(">c1\nACGT\n")
    (root / "something_nobody_claimed.txt").write_text("hello")
    man = Manifest(status_cmd.manifest_path(str(root)))
    man.finish("assembly", man.start("assembly", ["c"], []), COMPLETED, outputs=[str(declared)])
    assert status_cmd.main([str(root), "--audit"]) == 0
    assert "something_nobody_claimed.txt" in capsys.readouterr().out


def test_audit_does_not_report_the_drivers_own_directories(tmp_path, capsys):
    root = tmp_path / "study"
    (root / ".metawrap2" / "run_steps").mkdir(parents=True)
    (root / "LOGS").mkdir()
    (root / "CLEAN_READS").mkdir()
    declared = root / "final_assembly.fasta"
    declared.write_text(">c1\nACGT\n")
    man = Manifest(status_cmd.manifest_path(str(root)))
    man.finish("assembly", man.start("assembly", ["c"], []), COMPLETED, outputs=[str(declared)])
    status_cmd.main([str(root), "--audit"])
    out = capsys.readouterr().out
    assert "LOGS" not in out and "CLEAN_READS" not in out


# --- the audit API -------------------------------------------------------------------------


def test_declared_paths_collects_every_steps_outputs(study):
    man = Manifest(status_cmd.manifest_path(str(study)))
    paths = man.declared_paths()
    assert any(p.endswith("BINS") for p in paths)
    assert any(p.endswith("asm.fa") for p in paths)


def test_audit_ignores_steps_that_did_not_complete(tmp_path):
    """A failed step's missing outputs are expected, not a finding."""
    root = tmp_path / "study"
    (root / ".metawrap2" / "run_steps").mkdir(parents=True)
    man = Manifest(status_cmd.manifest_path(str(root)))
    man.finish(
        "assembly",
        man.start("assembly", ["c"], []),
        FAILED,
        outputs=[str(root / "never_written.fa")],
    )
    problems = man.audit()
    assert problems["missing"] == []
    assert problems["no_outputs"] == []
