"""Tests for `annotate` failing when every bin fails.

Checked:

* a missing Bakta database stops the step before anything is launched;
* a per-bin loop in which every bin failed fails the step;
* one failing bin out of several does not fail the step.

This guards against a run where, with no Bakta database, ``annotate --bakta`` failed on all 26 bins
and still exited 0 with its success banner.

Bakta and dbCAN are not run; ``run`` is replaced.
"""

from __future__ import annotations

import os
from typing import List

import pytest

from metawrap2.modules import annotate

#: A contig with one gene, so gene calling (Pyrodigal, run in-process) writes the per-genome files
#: the later steps need.
CONTIG = (
    "ATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGTAA"
) * 5


def _bins(directory: str, count: int) -> str:
    """A folder of *count* bins, named the way bin_refinement names them."""
    os.makedirs(directory, exist_ok=True)
    for index in range(1, count + 1):
        with open(os.path.join(directory, "bin_%03d.strict.fa" % index), "w") as handle:
            handle.write(">contig_%d\n%s\n" % (index, CONTIG))
    return directory


def _config(path: str, **databases: str) -> str:
    with open(path, "w") as handle:
        handle.write("[settings]\nuse_conda_envs = false\n\n[databases]\n")
        handle.writelines('%s = "%s"\n' % (key, value) for key, value in databases.items())
    return path


@pytest.fixture
def commands(monkeypatch) -> List[str]:
    """Record every command the module would run, and let each one succeed."""
    calls: List[str] = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return 0

    monkeypatch.setattr(annotate, "run", fake_run)
    return calls


def test_bakta_without_a_database_fails_before_launching_anything(tmp_path, commands, capsys):
    """With BAKTA_DB unset, the step exits 1 before running Bakta and says how to get the
    database."""
    config = _config(str(tmp_path / "metawrap2.toml"))
    argv = [
        "-b",
        _bins(str(tmp_path / "bins"), 3),
        "-o",
        str(tmp_path / "out"),
        "--bakta",
        "--config",
        config,
    ]
    with pytest.raises(SystemExit) as stopped:
        annotate.main(argv)
    assert stopped.value.code == 1
    out = capsys.readouterr().out
    assert "--bakta needs BAKTA_DB" in out
    assert "bakta_db download" in out, "the message has to name the way out"
    assert commands == [], "nothing may be launched before the database is checked"


def test_bakta_failing_on_every_bin_fails_the_step(tmp_path, monkeypatch, capsys):
    """When Bakta fails on every bin, the step exits 1 without the success banner."""
    db = tmp_path / "bakta_db"
    db.mkdir()
    monkeypatch.setattr(annotate, "run", lambda cmd, **kwargs: 1)
    argv = [
        "-b",
        _bins(str(tmp_path / "bins"), 3),
        "-o",
        str(tmp_path / "out"),
        "--bakta",
        "--config",
        _config(str(tmp_path / "metawrap2.toml"), BAKTA_DB=str(db)),
    ]
    with pytest.raises(SystemExit) as stopped:
        annotate.main(argv)
    assert stopped.value.code == 1
    out = capsys.readouterr().out
    assert "Bakta failed on 3 of 3 bins" in out
    assert "no Bakta output was produced" in out
    assert "SUCCESSFULLY FINISHED" not in out


def test_one_bin_failing_out_of_several_still_finishes(tmp_path, monkeypatch, capsys):
    """One failing bin out of three is reported, and the other bins are still annotated."""
    db = tmp_path / "bakta_db"
    db.mkdir()
    attempts: List[str] = []

    def fail_the_first(cmd, **kwargs):
        attempts.append(cmd)
        return 1 if len(attempts) == 1 else 0

    monkeypatch.setattr(annotate, "run", fail_the_first)
    argv = [
        "-b",
        _bins(str(tmp_path / "bins"), 3),
        "-o",
        str(tmp_path / "out"),
        "--bakta",
        "--config",
        _config(str(tmp_path / "metawrap2.toml"), BAKTA_DB=str(db)),
    ]
    assert annotate.main(argv) == 0
    out = capsys.readouterr().out
    assert "Bakta failed on" in out, "the one failure is still reported"
    assert "of 3 bins" not in out, "one failure of three is not a failed step"
    assert "SUCCESSFULLY FINISHED" in out
    assert len(attempts) == 3, "the remaining bins are still annotated"


def test_dbcan_failing_on_every_genome_fails_the_step(tmp_path, monkeypatch, capsys):
    """When dbCAN fails on every genome, `--cazy` exits 1 instead of reporting 0 CAZymes."""
    db = tmp_path / "dbcan"
    db.mkdir()
    monkeypatch.setattr(annotate, "run", lambda cmd, **kwargs: 1)
    argv = [
        "-b",
        _bins(str(tmp_path / "bins"), 2),
        "-o",
        str(tmp_path / "out"),
        "--cazy",
        "--config",
        _config(str(tmp_path / "metawrap2.toml"), DBCAN_DB=str(db)),
    ]
    with pytest.raises(SystemExit) as stopped:
        annotate.main(argv)
    assert stopped.value.code == 1
    out = capsys.readouterr().out
    assert "dbCAN failed on 2 of 2 bins" in out
    assert "SUCCESSFULLY FINISHED" not in out
