"""Refining twice into the same output folder.

A rerun (e.g. after a binner was re-run or a bin set corrected) used to finish all its work and
then raise FileExistsError while moving its working sets into ``work_files``, because
``os.replace`` refuses a non-empty directory. An 8 CPU-hour re-refinement wrote 70 bins and
exited 1, with the remaining bin sets half-moved.

CheckM is faked as in the other bin_refinement tests, with a deterministic stand-in derived from
each bin's contig set, so consolidation does real work and the run reaches finalisation.
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
    # No conda env is used, but env_for checks the env manager is on PATH first, and mamba is not
    # on PATH in a non-login shell.
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

    ``force`` is needed on the rerun, since an output folder holding a finished run is refused
    without it.

    Set through the runner, not argv: ``--force`` is a global flag that ``cli.main`` strips before
    the module sees it, so the module parser would exit 2. Reset in a finally because it is
    module-global state.
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
    """A second run into the same folder with --force succeeds and gives the same bins."""
    first = _run(tmp_path, "AB")
    expected = _final(first)
    assert expected, "the first run has to produce bins for the second to be a rerun at all"

    second = _run(tmp_path, "AB", force=True)
    assert _final(second) == expected
    # No half-moved working sets are left behind.
    assert [f for f in os.listdir(second) if f.startswith("bins")] == []
    work = os.path.join(second, "work_files")
    assert "binsA" in os.listdir(work) and "binsB" in os.listdir(work)


def test_resuming_does_not_clean_what_it_is_resuming_from(tmp_path):
    """--resume does not delete the bin sets its checkpoints skip rebuilding.

    The clean used to be unconditional, so a resume deleted those bin sets. Checked on the final
    bins, since a resume that deleted its inputs and rebuilt nothing could still exit 0.
    """
    first = _run(tmp_path, "AB")
    expected = _final(first)
    out = _run(tmp_path, "AB", resume=True)
    assert _final(out) == expected
    work = os.listdir(os.path.join(out, "work_files"))
    assert "binsA" in work and "binsB" in work, work
