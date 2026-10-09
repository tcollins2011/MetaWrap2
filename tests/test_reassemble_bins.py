import os

from metawrap2.modules import reassemble_bins


def test_command_templates_format_cleanly():
    idx = reassemble_bins.BWA_INDEX.format(assembly="asm.fa")
    assert idx == "bwa index asm.fa"
    sp = reassemble_bins.SPADES.format(
        threads=4, mem=40, tmp="/t", contigs="c.fa", r1="a_1.fq", r2="a_2.fq", out="/o"
    )
    assert "{" not in sp and "--careful" in sp and "-t 4" in sp and "-o /o" in sp
    spn = reassemble_bins.SPADES_NANOPORE.format(
        threads=1, mem=8, tmp="/t", contigs="c.fa", r1="a_1.fq", r2="a_2.fq", nano="n.fq", out="/o"
    )
    assert "{" not in spn and "--nanopore n.fq" in spn


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
    reassemble_bins._clean_previous_run(out)
    assert not os.path.exists(os.path.join(out, "reassembled_bins"))
    assert not os.path.exists(os.path.join(out, "work_files"))
    assert not os.path.exists(os.path.join(out, "reassembled_bins.stats.tsv"))


def test_the_clean_leaves_the_rest_of_the_output_folder_alone(tmp_path):
    """It must not take the run record, the log, or the original bins with it.

    original_bins is refreshed from the input every run, and provenance is the one thing a rerun
    should not destroy - if the clean removed it there would be no record that the first run happened.
    """
    out = str(tmp_path / "out")
    _finished_run(out)
    os.makedirs(os.path.join(out, "original_bins"))
    for name in ("run_config.json", "metawrap2.log", "provenance.txt"):
        with open(os.path.join(out, name), "w") as fh:
            fh.write("{}\n")
    reassemble_bins._clean_previous_run(out)
    kept = sorted(os.listdir(out))
    assert kept == ["metawrap2.log", "original_bins", "provenance.txt", "run_config.json"], kept


def test_the_clean_is_fine_on_a_folder_that_has_no_previous_run(tmp_path):
    """--force on a fresh or half-written folder is ordinary, not an error."""
    out = str(tmp_path / "out")
    os.makedirs(out)
    reassemble_bins._clean_previous_run(out)
    assert os.listdir(out) == []
