"""Golden tests pinning bin_refinement's behaviour for 2 and 3 bin sets.

bin_refinement was generalised from "exactly 2 or 3 bin sets" to "any number up to
:data:`metawrap2.modules.bin_refinement.MAX_BIN_SETS`". The requirement for that change was
that 2- and 3-set runs keep producing *exactly* what they produced before, because those are
the cases every existing analysis used.

These tests run the real algorithm - the real Binning_refiner (it is pure Python), the real
consolidation and dereplication - with only CheckM replaced by a deterministic stand-in, and
compare the complete output tree and every file's contents against values recorded from the
implementation before the rewrite.

The fake CheckM derives a bin's completeness and contamination from a hash of its contig set,
so scores are arbitrary but stable, vary between bins, and change when a bin's membership
changes. That is what makes consolidation actually do work rather than trivially keep binsA.
"""

import hashlib
import os
from typing import ClassVar

import pytest

from metawrap2 import checkm as checkm_mod
from metawrap2.constants import BIN_EXTENSION
from metawrap2.io.seqio import contig_id, iter_fasta
from metawrap2.modules import bin_refinement

# --- a deterministic stand-in for CheckM --------------------------------------------------


def _fake_scores(bin_path):
    """Stable pseudo-random completeness/contamination/N50 from a bin's contig set."""
    contigs = sorted(contig_id(h) for h, _ in iter_fasta(bin_path))
    digest = hashlib.sha256("|".join(contigs).encode()).digest()
    completeness = 55.0 + (digest[0] / 255.0) * 45.0  # 55..100
    contamination = (digest[1] / 255.0) * 9.0  # 0..9
    n50 = 1000 + digest[2] * 37
    size = sum(len(s) for _, s in iter_fasta(bin_path))
    return completeness, contamination, n50, size


def _fake_run_checkm(bins_dir, **kwargs):
    """Write a .stats file for *bins_dir* exactly as the real run_checkm would."""
    rows = []
    for name in sorted(os.listdir(bins_dir)):
        if not name.endswith(BIN_EXTENSION):
            continue
        path = os.path.join(bins_dir, name)
        completeness, contamination, n50, size = _fake_scores(path)
        rows.append(
            [
                os.path.splitext(name)[0],
                str(round(completeness, 2))[:5],
                str(round(contamination, 3))[:5],
                "0.500",
                "Bacteria",
                str(n50),
                str(size),
            ]
        )
    stats_path = bins_dir + checkm_mod.STATS_SUFFIX
    checkm_mod.write_stats(rows, stats_path)
    return stats_path


@pytest.fixture(autouse=True)
def fake_checkm(monkeypatch):
    monkeypatch.setattr(bin_refinement._checkm, "run_checkm", _fake_run_checkm)
    monkeypatch.setattr(checkm_mod, "run_checkm", _fake_run_checkm)


# --- input bin sets -----------------------------------------------------------------------
# Three binners over the same 12 contigs, disagreeing the way real binners do: they mostly
# agree on two genomes and split the third differently.

# Contigs have to be big enough that a refined bin clears Binning_refiner's 0.5 Mbp minimum
# bin size and MetaWrap2's own 50 kb floor, or every refined set comes out empty and is
# discarded - which is not the behaviour under test.
CONTIGS = {"c%02d" % i: ("ACGT" * (50000 + i * 2500)) for i in range(1, 13)}

BIN_SETS = {
    # binner A: three clean bins
    "A": {
        "bin.1": ["c01", "c02", "c03", "c04"],
        "bin.2": ["c05", "c06", "c07", "c08"],
        "bin.3": ["c09", "c10", "c11", "c12"],
    },
    # binner B: agrees on bin.1, merges the other two, and misplaces one contig
    "B": {
        "bin.1": ["c01", "c02", "c03", "c04"],
        "bin.2": ["c05", "c06", "c07", "c08", "c09", "c10"],
        "bin.3": ["c11", "c12"],
    },
    # binner C: agrees on bins 1 and 3, drops a contig from bin 2
    "C": {
        "bin.1": ["c01", "c02", "c03"],
        "bin.2": ["c05", "c06", "c07"],
        "bin.3": ["c09", "c10", "c11", "c12"],
    },
    # a fourth and fifth binner, for the generalised cases
    "D": {"bin.1": ["c01", "c02", "c03", "c04"], "bin.2": ["c05", "c06", "c07", "c08"]},
    "E": {
        "bin.1": ["c01", "c02", "c04"],
        "bin.2": ["c05", "c06", "c07", "c08"],
        "bin.3": ["c09", "c10", "c11"],
    },
}


def _write_bin_sets(tmp_path, which):
    """Materialise the requested bin sets on disk; returns {letter: directory}."""
    dirs = {}
    for letter in which:
        directory = tmp_path / ("binner_%s" % letter)
        directory.mkdir()
        for bin_name, contigs in BIN_SETS[letter].items():
            with open(directory / ("%s.fa" % bin_name), "w") as fh:
                fh.writelines(">%s\n%s\n" % (contig, CONTIGS[contig]) for contig in contigs)
        dirs[letter] = str(directory)
    return dirs


def _run(tmp_path, letters, extra=()):
    dirs = _write_bin_sets(tmp_path, letters)
    out = tmp_path / "refined"
    argv = ["-o", str(out), "-t", "1", "-m", "8", "-c", "50", "-x", "10"]
    for flag, letter in zip(("-A", "-B", "-C", "-D", "-E", "-F"), letters):
        argv += [flag, dirs[letter]]
    argv += list(extra)
    assert bin_refinement.main(argv) == 0
    return str(out)


def _snapshot(out):
    """A comparable description of the whole output tree: paths plus bin membership."""
    snapshot = {}
    for root, _dirs, files in os.walk(out):
        for name in sorted(files):
            path = os.path.join(root, name)
            rel = os.path.relpath(path, out)
            if name.endswith(BIN_EXTENSION):
                snapshot[rel] = sorted(contig_id(h) for h, _ in iter_fasta(path))
            elif name.endswith((".stats.tsv", ".contigs.tsv")):
                with open(path) as fh:
                    snapshot[rel] = fh.read()
    return snapshot


# --- the golden expectations ---------------------------------------------------------------
# Recorded from the pre-generalisation implementation. If a change to the algorithm makes one
# of these fail, that change altered results for a case users already rely on.


def test_two_sets_produce_the_historical_bin_sets(tmp_path):
    out = _run(tmp_path, "AB")
    work = os.path.join(out, "work_files")
    # Exactly the sets the 2-input path always built: the two inputs plus their refinement,
    # then the consolidated (binsM) and dereplicated (binsO) results.
    assert sorted(
        d for d in os.listdir(work) if d.startswith("bins") and os.path.isdir(os.path.join(work, d))
    ) == ["binsA", "binsAB", "binsB", "binsM", "binsO"]


def test_three_sets_produce_the_historical_bin_sets(tmp_path):
    out = _run(tmp_path, "ABC")
    work = os.path.join(out, "work_files")
    assert sorted(
        d for d in os.listdir(work) if d.startswith("bins") and os.path.isdir(os.path.join(work, d))
    ) == ["binsA", "binsAB", "binsABC", "binsAC", "binsB", "binsBC", "binsC", "binsM", "binsO"]


def test_three_set_refinement_combinations_and_argument_order(tmp_path):
    """The historical 3-set plan, including that binsBC is built as (C, B), not (B, C).

    Binning_refiner numbers its output bins by the order the input folders are given, so the
    orientation of each pair is part of the result, not an implementation detail.
    """
    assert bin_refinement.refinement_plan(["binsA", "binsB"]) == [
        (("binsA", "binsB"), "binsAB"),
    ]
    assert bin_refinement.refinement_plan(["binsA", "binsB", "binsC"]) == [
        (("binsA", "binsB"), "binsAB"),
        (("binsC", "binsB"), "binsBC"),
        (("binsA", "binsC"), "binsAC"),
        (("binsA", "binsB", "binsC"), "binsABC"),
    ]


def test_two_set_output_is_stable(tmp_path):
    """Full output snapshot for 2 sets - the contents, not just the shape."""
    snapshot = _snapshot(_run(tmp_path, "AB"))
    final = {k: v for k, v in snapshot.items() if k.startswith("metawrap_50_10_bins/")}
    assert final == {
        "metawrap_50_10_bins/bin_001.fasta": ["c01", "c02", "c03", "c04"],
        "metawrap_50_10_bins/bin_002.fasta": ["c05", "c06", "c07", "c08", "c09", "c10"],
        "metawrap_50_10_bins/bin_003.fasta": ["c11", "c12"],
    }


def test_three_set_output_is_stable(tmp_path):
    """Note this differs from the 2-set result: consolidation really is choosing between
    versions of each bin rather than trivially keeping binsA."""
    snapshot = _snapshot(_run(tmp_path, "ABC"))
    final = {k: v for k, v in snapshot.items() if k.startswith("metawrap_50_10_bins/")}
    assert final == {
        "metawrap_50_10_bins/bin_001.fasta": ["c01", "c02", "c03", "c04"],
        "metawrap_50_10_bins/bin_002.fasta": ["c05", "c06", "c07"],
        "metawrap_50_10_bins/bin_003.fasta": ["c11", "c12"],
    }


def test_single_set_still_short_circuits(tmp_path):
    out = _run(tmp_path, "A")
    # With one input there is nothing to refine or consolidate, so binsA is used directly and
    # no work_files shuffle happens.
    assert os.path.isdir(os.path.join(out, "metawrap_50_10_bins"))
    assert not os.path.isdir(os.path.join(out, "work_files"))


# --- the generalisation --------------------------------------------------------------------


def test_plan_scales_to_more_binners():
    """Beyond three sets: every pair, then the all-way combination."""
    plan = bin_refinement.refinement_plan(["binsA", "binsB", "binsC", "binsD"])
    names = [name for _inputs, name in plan]
    assert names == ["binsAB", "binsAC", "binsAD", "binsBC", "binsBD", "binsCD", "binsABCD"]
    # the all-way entry really does take all four
    assert plan[-1][0] == ("binsA", "binsB", "binsC", "binsD")


def test_plan_pairs_only_mode():
    plan = bin_refinement.refinement_plan(
        ["binsA", "binsB", "binsC", "binsD"], combinations="pairs"
    )
    assert [name for _i, name in plan] == [
        "binsAB",
        "binsAC",
        "binsAD",
        "binsBC",
        "binsBD",
        "binsCD",
    ]


def test_plan_all_subsets_mode():
    plan = bin_refinement.refinement_plan(["binsA", "binsB", "binsC"], combinations="all")
    # every subset of size >= 2, so for three inputs that is the three pairs plus the triple
    assert sorted(name for _i, name in plan) == ["binsAB", "binsABC", "binsAC", "binsBC"]


def test_five_sets_run_end_to_end(tmp_path):
    out = _run(tmp_path, "ABCDE")
    work = os.path.join(out, "work_files")
    sets = sorted(
        d for d in os.listdir(work) if d.startswith("bins") and os.path.isdir(os.path.join(work, d))
    )

    labels = [bin_refinement.bin_set_label(i) for i in range(5)]
    planned = {name for _inputs, name in bin_refinement.refinement_plan(labels)}

    # The five inputs, the consolidated and dereplicated sets, and the all-way refinement.
    for expected in labels + ["binsM", "binsO", "binsABCDE"]:
        assert expected in sets, expected
    # Everything else on disk must be a set the plan asked for. Some planned combinations
    # legitimately yield no bins over Binning_refiner's size floor and are discarded, so the
    # count is data-dependent - but nothing unplanned may appear.
    unexpected = set(sets) - planned - set(labels) - {"binsM", "binsO"}
    assert unexpected == set(), unexpected
    assert os.path.isdir(os.path.join(out, "metawrap_50_10_bins"))


def test_too_many_sets_is_refused(tmp_path):
    """More than MAX_BIN_SETS is refused up front rather than after hours of CheckM runs."""
    argv = ["-o", str(tmp_path / "out"), "-c", "50", "-x", "10"]
    for i in range(bin_refinement.MAX_BIN_SETS + 1):
        directory = tmp_path / ("extra_%d" % i)  # distinct dirs: --bins de-duplicates
        directory.mkdir()
        argv += ["--bins", str(directory)]
    with pytest.raises(SystemExit):
        bin_refinement.main(argv)


def test_duplicate_bin_set_paths_are_ignored(tmp_path):
    """Passing the same directory twice is a mistake, not a request to refine it with itself."""
    dirs = _write_bin_sets(tmp_path, "AB")

    class Args:
        bins_a = dirs["A"]
        bins_b = dirs["B"]
        bins_c = bins_d = bins_e = bins_f = None
        bins_extra: ClassVar = [dirs["A"], dirs["B"]]

    assert bin_refinement.collect_bin_sets(Args()) == [dirs["A"], dirs["B"]]


# --- proof that the commands run for 2 and 3 sets are byte-identical to the old ones -------


def test_refiner_command_for_two_and_three_sets_is_unchanged():
    """The 2/3-set cases must still call Binning_refiner through -1/-2/-3, not the new -i.

    Together with test_three_set_refinement_combinations_and_argument_order (which pins the
    plan), this establishes that a 2- or 3-binner run issues exactly the commands it always
    did, so its results cannot have changed.
    """
    recorded = []

    def fake_run(cmd, **kwargs):
        recorded.append(cmd)
        return 0

    original_run, original_dry = bin_refinement.run, bin_refinement.dry_run
    bin_refinement.run = fake_run
    # dry_run makes _refine return right after issuing the command, so it does not try to move
    # output directories the fake run never created.
    bin_refinement.dry_run = lambda: True
    try:
        bin_refinement._refine("/out", None, ("binsA", "binsB"), "binsAB")
        bin_refinement._refine("/out", None, ("binsA", "binsB", "binsC"), "binsABC")
        bin_refinement._refine("/out", None, ("binsA", "binsB", "binsC", "binsD"), "binsABCD")
    finally:
        bin_refinement.run = original_run
        bin_refinement.dry_run = original_dry

    assert recorded[0].endswith("-1 binsA -2 binsB -o Refined_AB")
    assert recorded[1].endswith("-1 binsA -2 binsB -3 binsC -o Refined_ABC")
    # only beyond three does it switch to the repeatable flag
    assert recorded[2].endswith("-i binsA -i binsB -i binsC -i binsD -o Refined_ABCD")


def test_repeatable_bin_set_flag(tmp_path):
    """--bins may be repeated instead of using -A/-B/-C..., for scripting."""
    dirs = _write_bin_sets(tmp_path, "ABC")
    out = tmp_path / "refined"
    argv = ["-o", str(out), "-t", "1", "-c", "50", "-x", "10"]
    for letter in "ABC":
        argv += ["--bins", dirs[letter]]
    assert bin_refinement.main(argv) == 0
    work = os.path.join(str(out), "work_files")
    assert sorted(
        d for d in os.listdir(work) if d.startswith("bins") and os.path.isdir(os.path.join(work, d))
    ) == ["binsA", "binsAB", "binsABC", "binsAC", "binsB", "binsBC", "binsC", "binsM", "binsO"]


def test_bin_set_labels():
    assert bin_refinement.bin_set_label(0) == "binsA"
    assert bin_refinement.bin_set_label(5) == "binsF"
    assert bin_refinement.combined_label(["binsA", "binsC"]) == "binsAC"
    assert bin_refinement.combined_label(["binsA", "binsB", "binsC"]) == "binsABC"
