import os

from metawrap2 import refinement


def _write_bin(path, contigs):
    """contigs: list of (name, length) -> writes a FASTA with A*length sequences."""
    with open(path, "w") as fh:
        for name, length in contigs:
            fh.write(">%s\n%s\n" % (name, "A" * length))


def _write_stats(path, rows):
    """rows: list of (bin_name, completeness, contamination, extra)."""
    with open(path, "w") as fh:
        fh.write("bin\tcompleteness\tcontamination\tstuff\tstuff2\tsize\n")
        for name, comp, cont, extra in rows:
            fh.write("%s\t%s\t%s\t0\t0\t%s\n" % (name, comp, cont, extra))


def test_overlap_percent_identical():
    a = {"c1": 100, "c2": 100}
    assert refinement._overlap_percent(a, dict(a)) == 100.0


def test_overlap_percent_disjoint():
    assert refinement._overlap_percent({"c1": 100}, {"c2": 100}) == 0.0


def test_consolidate_picks_higher_scoring_overlapping_bin(tmp_path):
    f1 = tmp_path / "binsA"
    f2 = tmp_path / "binsB"
    f1.mkdir()
    f2.mkdir()
    # Same two contigs in both -> 100% overlap. B is higher quality, so B wins.
    _write_bin(f1 / "a1.fa", [("contig_1", 1000), ("contig_2", 1000)])
    _write_bin(f2 / "b1.fa", [("contig_1", 1000), ("contig_2", 1000)])
    s1 = tmp_path / "binsA.stats"
    s2 = tmp_path / "binsB.stats"
    _write_stats(s1, [("a1", 80.0, 5.0, 2000)])
    _write_stats(s2, [("b1", 95.0, 2.0, 2000)])
    out = tmp_path / "binsM"
    n = refinement.consolidate(str(f1), str(f2), str(s1), str(s2), str(out), 70.0, 10.0)
    # one consolidated bin, and it should be the B version (higher score)
    assert n == 2  # bin_ct ends at 2 after writing bin.1
    assert (out / "bin.1.fa").exists()
    # B's contigs are identical here; assert stats line carried B's completeness (95)
    stats_out = (tmp_path / "binsM.stats").read_text()
    assert "95.0" in stats_out


def test_consolidate_keeps_unmatched_second_set_bin(tmp_path):
    f1 = tmp_path / "binsA"
    f2 = tmp_path / "binsB"
    f1.mkdir()
    f2.mkdir()
    _write_bin(f1 / "a1.fa", [("x1", 1000)])
    _write_bin(f2 / "b1.fa", [("y1", 1000)])  # disjoint -> no overlap
    s1 = tmp_path / "binsA.stats"
    s2 = tmp_path / "binsB.stats"
    _write_stats(s1, [("a1", 90.0, 1.0, 1000)])
    _write_stats(s2, [("b1", 85.0, 1.0, 1000)])
    out = tmp_path / "binsM"
    refinement.consolidate(str(f1), str(f2), str(s1), str(s2), str(out), 70.0, 10.0)
    # both bins survive since they don't overlap
    assert sorted(os.listdir(out)) == ["bin.1.fa", "bin.2.fa"]


def test_consolidate_excludes_low_quality_bins(tmp_path):
    f1 = tmp_path / "binsA"
    f2 = tmp_path / "binsB"
    f1.mkdir()
    f2.mkdir()
    _write_bin(f1 / "a1.fa", [("x1", 1000)])
    _write_bin(f2 / "b1.fa", [("y1", 1000)])
    s1 = tmp_path / "binsA.stats"
    s2 = tmp_path / "binsB.stats"
    _write_stats(s1, [("a1", 90.0, 1.0, 1000)])
    _write_stats(s2, [("b1", 40.0, 1.0, 1000)])  # below 70% completion -> dropped
    out = tmp_path / "binsM"
    refinement.consolidate(str(f1), str(f2), str(s1), str(s2), str(out), 70.0, 10.0)
    assert os.listdir(out) == ["bin.1.fa"]


def test_dereplicate_keeps_shared_contig_in_best_bin(tmp_path):
    folder = tmp_path / "bins"
    folder.mkdir()
    # contig "shared" is in both bins; bin2 has higher score so it should keep it.
    _write_bin(folder / "bin.1.fa", [("shared", 1000), ("only1", 500)])
    _write_bin(folder / "bin.2.fa", [("shared", 1000), ("only2", 500)])
    stats = tmp_path / "bins.stats"
    _write_stats(stats, [("bin.1", 80.0, 2.0, 1500), ("bin.2", 95.0, 1.0, 1500)])
    out = tmp_path / "out"
    refinement.dereplicate(str(stats), str(folder), str(out), "best")
    b1 = dict(refinement.iter_fasta(str(out / "bin.1.fa")))
    b2 = dict(refinement.iter_fasta(str(out / "bin.2.fa")))
    assert "shared" not in b1 and "only1" in b1
    assert "shared" in b2 and "only2" in b2


def test_dereplicate_remove_mode_drops_shared_from_all(tmp_path):
    folder = tmp_path / "bins"
    folder.mkdir()
    _write_bin(folder / "bin.1.fa", [("shared", 1000), ("only1", 500)])
    _write_bin(folder / "bin.2.fa", [("shared", 1000), ("only2", 500)])
    stats = tmp_path / "bins.stats"
    _write_stats(stats, [("bin.1", 80.0, 2.0, 1500), ("bin.2", 95.0, 1.0, 1500)])
    out = tmp_path / "out"
    refinement.dereplicate(str(stats), str(folder), str(out), "remove")
    b1 = dict(refinement.iter_fasta(str(out / "bin.1.fa")))
    b2 = dict(refinement.iter_fasta(str(out / "bin.2.fa")))
    assert "shared" not in b1 and "shared" not in b2
