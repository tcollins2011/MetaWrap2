"""`annotate` reported success after failing on every bin. Both halves of that, and the policy.

From one real run: with no Bakta database configured, ``annotate --bakta`` launched Bakta 26 times,
failed 26 times, printed 26 warnings, produced no Bakta output and exited **0** with
``ANNOTATE PIPELINE SUCCESSFULLY FINISHED!!!``. Anything scripting on exit codes recorded that as a
finished annotation.

So two behaviours are pinned here - the configurable reason every bin would fail is checked before
anything is launched, and a per-item loop whose every item failed fails the step - plus the policy
neither may break: one pathological bin out of several still costs only that bin.

Nothing here runs Bakta or dbCAN. ``run`` is replaced, which is the only thing between this module
and a subprocess.
"""

from __future__ import annotations

import os
from typing import List

import pytest

from metawrap2.modules import annotate

#: A contig with one findable gene on it, so gene calling (which is real here - Pyrodigal runs
#: in-process) writes the per-genome files the later steps look for.
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
    """BAKTA_DB unset must stop the step, not produce one failed Bakta per bin.

    If this regresses, `--bakta` goes back to launching Bakta once per genome against no database -
    26 launches and 26 warnings on the report this comes from - for a reason that was knowable
    before the first one.
    """
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
    """Every bin failing is one cause, no output, and must not be reported as a finished run.

    This is the exact defect: 26 of 26 bins failed, the warnings were printed, and the module still
    exited 0 with its success banner.
    """
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
    """Proceeding with the rest is the right policy for a single pathological bin - keep it.

    The all-failed rule must not turn one bad genome into a failed annotation: that would cost a
    user the other 25 genomes' output for a bin that was never going to work.
    """
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
    """The CAZyme loop had the same blind spot, and the same rule applies to it.

    `--cazy` wrote a header-only cazymes.tsv and reported "0 CAZyme gene(s)" as a success when
    every dbCAN call failed - a result indistinguishable from a genuinely CAZyme-free bin set.
    """
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
