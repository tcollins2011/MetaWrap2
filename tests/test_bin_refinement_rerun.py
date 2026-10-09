"""Refining twice into the same output folder.

This is the ordinary thing to do when a run was wrong for a reason outside refinement - a binner was
re-run, a bin set was corrected - and it failed in the worst available way: the second run did all of
the work, wrote every final bin, plotted the figures, and *then* raised FileExistsError while tidying
its working sets into ``work_files``, because ``os.replace`` will overwrite a file but refuses a
directory that already has anything in it. It was found on a real 8 cpu-hour re-refinement that
returned 70 complete bins and exit code 1, with nothing in the log after "MOVING OVER TEMPORARY
FILES" - a successful run reported as a failure, having also left its remaining bin sets half-moved
for the next ``--resume`` to trip over.

CheckM is faked the way the other bin_refinement test modules fake it: with a deterministic
stand-in derived from each bin's contig set, so consolidation does real work and the run reaches the
finalisation step under test.
"""

from __future__ import annotations

import hashlib
import os

import pytest

from metawrap2 import checkm as checkm_mod
from metawrap2 import command
from metawrap2.constants import BIN_EXTENSION
from metawrap2.io.seqio import contig_id, iter_fasta
from metawrap2.modules import _common, bin_refinement

CONTIGS = {"c%02d" % i: ("ACGT" * (50000 + i * 2500)) for i in range(1, 13)}

BIN_SETS = {
    "A": {
        "bin.1": ["c01", "c02", "c03", "c04"],
        "bin.2": ["c05", "c06", "c07", "c08"],
        "bin.3": ["c09", "c10", "c11", "c12"],
    },
    "B": {
        "bin.1": ["c01", "c02", "c03", "c04"],
        "bin.2": ["c05", "c06", "c07", "c08", "c09", "c10"],
        "bin.3": ["c11", "c12"],
    },
    "C": {
        "bin.1": ["c01", "c02", "c03"],
        "bin.2": ["c05", "c06", "c07"],
        "bin.3": ["c09", "c10", "c11", "c12"],
    },
}


def _fake_run_checkm(bins_dir, **kwargs):
    rows = []
    for name in sorted(os.listdir(bins_dir)):
        if not name.endswith(BIN_EXTENSION):
            continue
        path = os.path.join(bins_dir, name)
        contigs = sorted(contig_id(h) for h, _ in iter_fasta(path))
        digest = hashlib.sha256("|".join(contigs).encode()).digest()
        rows.append(
            [
                os.path.splitext(name)[0],
                str(round(55.0 + (digest[0] / 255.0) * 45.0, 2))[:5],
                str(round((digest[1] / 255.0) * 9.0, 3))[:5],
                "0.500",
                "Bacteria",
                str(1000 + digest[2] * 37),
                str(sum(len(s) for _, s in iter_fasta(path))),
            ]
        )
    stats_path = bins_dir + checkm_mod.STATS_SUFFIX
    checkm_mod.write_stats(rows, stats_path)
    return stats_path


@pytest.fixture(autouse=True)
def fake_checkm(monkeypatch):
    monkeypatch.setattr(bin_refinement._checkm, "run_checkm", _fake_run_checkm)
    monkeypatch.setattr(checkm_mod, "run_checkm", _fake_run_checkm)
    # Nothing enters a conda env here, but env_for checks the env manager is findable first and
    # mamba is not on PATH in a non-login shell.
    monkeypatch.setattr(_common, "_require_env_manager", lambda env: None)


def _bin_sets(tmp_path, letters, sets=None):
    dirs = {}
    sets = sets or BIN_SETS
    for letter in letters:
        directory = tmp_path / ("binner_%s" % letter)
        directory.mkdir(exist_ok=True)
        for name in os.listdir(directory):
            os.remove(directory / name)
        for bin_name, contigs in sets[letter].items():
            with open(directory / ("%s.fa" % bin_name), "w") as fh:
                fh.writelines(">%s\n%s\n" % (c, CONTIGS[c]) for c in contigs)
        dirs[letter] = str(directory)
    return dirs


def _run(tmp_path, letters, force=False, sets=None, resume=False):
    """One full refinement into tmp_path/refined. Returns the output directory.

    ``force`` on the rerun, because an output folder holding a finished run is refused without it -
    which is the right default, and is also why the failure below was not caught by anything: the
    only supported way to redo a run in place is the one path nothing exercised end to end.

    Set through the runner rather than passed in argv: ``--force`` is a global flag that ``cli.main``
    strips before the module ever sees it, so a module parser handed it exits 2. Reset in a finally
    because it is module-global state that would otherwise leak into the next test.
    """
    dirs = _bin_sets(tmp_path, letters, sets)
    out = tmp_path / "refined"
    argv = ["-o", str(out), "-t", "1", "-m", "8", "-c", "50", "-x", "10", "--skip-gunc"]
    for flag, letter in zip(("-A", "-B", "-C"), letters):
        argv += [flag, dirs[letter]]
    command.set_force(force)
    command.set_resume(resume)
    try:
        assert bin_refinement.main(argv) == 0
    finally:
        command.set_force(False)
        command.set_resume(False)
    return str(out)


def _final(out):
    final = os.path.join(out, "metawrap_50_10_bins")
    return sorted(f for f in os.listdir(final) if f.endswith(BIN_EXTENSION))


def test_a_second_run_into_the_same_folder_succeeds(tmp_path):
    """The bug: the second run raised at the end, after writing everything it was asked for."""
    first = _run(tmp_path, "AB")
    expected = _final(first)
    assert expected, "the first run has to produce bins for the second to be a rerun at all"

    second = _run(tmp_path, "AB", force=True)
    assert _final(second) == expected


def test_the_second_run_leaves_no_bin_set_behind_in_the_output_root(tmp_path):
    """Half-moved working sets are the state the failure left, and what breaks the next --resume.

    Asserted on the output root rather than on work_files: "every bins* got moved" is the property
    that was violated, and it holds whatever the combination names are for this number of inputs.
    """
    _run(tmp_path, "AB")
    out = _run(tmp_path, "AB", force=True)
    assert [f for f in os.listdir(out) if f.startswith("bins")] == []
    work = os.path.join(out, "work_files")
    assert "binsA" in os.listdir(work) and "binsB" in os.listdir(work)


def test_resuming_does_not_clean_what_it_is_resuming_from(tmp_path):
    """--resume reaches the same already-exists output folder as --force and means the opposite.

    The clean was unconditional, so a resume deleted the bin sets whose checkpoint it was about to
    honour. Checked on the final bins rather than on the intermediates: a resume that destroyed its
    inputs and silently rebuilt nothing would still "succeed", and the only statement worth making is
    that the answer is unchanged.
    """
    first = _run(tmp_path, "AB")
    expected = _final(first)
    out = _run(tmp_path, "AB", resume=True)
    assert _final(out) == expected
    work = os.listdir(os.path.join(out, "work_files"))
    assert "binsA" in work and "binsB" in work, work


def test_a_rerun_with_a_corrected_input_does_not_present_the_old_one(tmp_path):
    """The reason to rerun in place is that an input changed, so the restored copy has to change too.

    Refinement copies each input set back into the output folder under the input folder's own name,
    which is where anything downstream looks for "what went in". That copy was only ever written when
    the destination did not already exist, so on a rerun it was never refreshed: a real refined/
    directory ended up advertising a 188-bin MaxBin2 set beside the results of refining the corrected
    129-bin one, with a .stats.tsv for the superseded set to go with it. Scoring the inputs from the
    output folder - the obvious thing to do - would have scored the wrong set and blamed refinement.
    """
    _run(tmp_path, "AB")
    corrected = dict(BIN_SETS)
    corrected["B"] = {"bin.1": ["c01", "c02", "c03"], "bin.2": ["c07", "c08", "c09", "c10"]}
    out = _run(tmp_path, "AB", force=True, sets=corrected)

    restored = os.path.join(out, "binner_B")
    names = sorted(
        os.path.splitext(f)[0] for f in os.listdir(restored) if f.endswith(BIN_EXTENSION)
    )
    assert names == ["bin.1", "bin.2"], names
    # The content too, not just the count: a two-bin set that still held the old bin.1 would pass a
    # count check while being the wrong sequence.
    got = sorted(
        contig_id(h) for h, _ in iter_fasta(os.path.join(restored, "bin.1" + BIN_EXTENSION))
    )
    assert got == ["c01", "c02", "c03"], got


def test_a_rerun_with_fewer_inputs_does_not_keep_the_dropped_set(tmp_path):
    """binsC from a 3-set run must not survive into a 2-set rerun of the same folder.

    The move only overwrites the sets the new run produces, so a set the new run never builds would
    otherwise sit in work_files looking like part of it - with a .stats.tsv and a contig table, and
    no way for a reader to tell it is three binners old.
    """
    _run(tmp_path, "ABC")
    out = _run(tmp_path, "AB", force=True)
    work = sorted(f for f in os.listdir(os.path.join(out, "work_files")) if f.startswith("bins"))
    assert not any(f.startswith(("binsC", "binsABC")) for f in work), work
    assert "binsA" in work and "binsB" in work and "binsAB" in work
