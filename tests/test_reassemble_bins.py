import os

from metawrap2.modules import reassemble_bins


def _finished_run(out, winners=("bin.1.strict", "bin.2.orig")):
    """The layout reassemble_bins leaves behind: chosen bins promoted, candidates in work_files."""
    os.makedirs(os.path.join(out, "reassembled_bins"))
    for name in winners:
        with open(os.path.join(out, "reassembled_bins", name + ".fasta"), "w") as fh:
            fh.write(">old\nACGT\n")
    with open(os.path.join(out, "reassembled_bins" + ".stats.tsv"), "w") as fh:
        fh.write("bin\tcompleteness\n")
    for item in (
        "reassembled_bins.checkm",
        "reads_for_reassembly",
        "binned_assembly",
        "reassemblies",
    ):
        os.makedirs(os.path.join(out, "work_files", item))
        with open(os.path.join(out, "work_files", item, "leftover"), "w") as fh:
            fh.write("from the first run\n")


def test_the_shuffle_survives_a_destination_that_already_exists(tmp_path):
    """The second reassembly into a folder died here, with every reassembly already done.

    ``os.replace`` overwrites a file but refuses a directory with anything in it, and every
    destination in ``work_files`` is the previous run's copy of the same intermediate. Reassembly is
    the most expensive module in the suite, so raising at the promote step is the worst place to do it.
    """
    out = str(tmp_path / "out")
    _finished_run(out)
    # this run's candidates and winners, sitting where the step expects to find them
    for item in ("reassembled_bins", "reads_for_reassembly", "binned_assembly", "reassemblies"):
        os.makedirs(os.path.join(out, item), exist_ok=True)
        with open(os.path.join(out, item, "new"), "w") as fh:
            fh.write("from this run\n")
    best = os.path.join(out, "reassembled_best_bins")
    os.makedirs(best)
    with open(os.path.join(best, "bin.1.permissive.fasta"), "w") as fh:
        fh.write(">new\nACGTACGT\n")

    reassembled = os.path.join(out, "reassembled_bins")
    reassemble_bins._shuffle_intermediates(out, best, reassembled)

    assert sorted(os.listdir(reassembled)) == ["bin.1.permissive.fasta"]
    # and work_files holds this run's candidates, not the previous run's leftovers
    moved = os.path.join(out, "work_files", "reassemblies")
    assert sorted(os.listdir(moved)) == ["new"], os.listdir(moved)
    assert not os.path.exists(best)


def test_force_clears_a_finished_reassembly(tmp_path):
    """Otherwise the previous run's winners stay in the pool this run picks its winners from.

    reassembled_bins is added to, not replaced, so the first run's promoted bins - and its .orig.fa
    copies of inputs that may since have changed - would be CheckM'd alongside this run's and could
    be chosen again. The output would look entirely normal.
    """
    out = str(tmp_path / "out")
    _finished_run(out)
    # The run record, the log and the original bins are not the previous run's results: provenance
    # is the one thing a rerun should not destroy.
    os.makedirs(os.path.join(out, "original_bins"))
    for name in ("run_config.json", "metawrap2.log", "provenance.txt"):
        with open(os.path.join(out, name), "w") as fh:
            fh.write("{}\n")
    reassemble_bins._clean_previous_run(out)
    assert not os.path.exists(os.path.join(out, "reassembled_bins"))
    assert not os.path.exists(os.path.join(out, "work_files"))
    assert not os.path.exists(os.path.join(out, "reassembled_bins.stats.tsv"))
    kept = sorted(os.listdir(out))
    assert kept == ["metawrap2.log", "original_bins", "provenance.txt", "run_config.json"], kept
