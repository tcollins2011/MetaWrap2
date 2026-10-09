from metawrap2 import refinement


def _write_bin(path, contigs):
    """contigs: list of (name, length) -> writes a FASTA with A*length sequences."""
    with open(path, "w") as fh:
        fh.writelines(">%s\n%s\n" % (name, "A" * length) for name, length in contigs)


def _write_stats(path, rows):
    """rows: list of (bin_name, completeness, contamination, extra)."""
    with open(path, "w") as fh:
        fh.write("bin\tcompleteness\tcontamination\tstuff\tstuff2\tsize\n")
        fh.writelines(
            "%s\t%s\t%s\t0\t0\t%s\n" % (name, comp, cont, extra) for name, comp, cont, extra in rows
        )


def test_overlap_percent_identical():
    a = {"c1": 100, "c2": 100}
    assert refinement._overlap_percent(a, dict(a)) == 100.0


def test_consolidate_picks_higher_scoring_overlapping_bin(tmp_path):
    f1 = tmp_path / "binsA"
    f2 = tmp_path / "binsB"
    f1.mkdir()
    f2.mkdir()
    # Same two contigs in both -> 100% overlap. B is higher quality, so B wins.
    _write_bin(f1 / "a1.fasta", [("contig_1", 1000), ("contig_2", 1000)])
    _write_bin(f2 / "b1.fasta", [("contig_1", 1000), ("contig_2", 1000)])
    s1 = tmp_path / "binsA.stats.tsv"
    s2 = tmp_path / "binsB.stats.tsv"
    _write_stats(s1, [("a1", 80.0, 5.0, 2000)])
    _write_stats(s2, [("b1", 95.0, 2.0, 2000)])
    out = tmp_path / "binsM"
    n = refinement.consolidate(str(f1), str(f2), str(s1), str(s2), str(out), 70.0, 10.0)
    # one consolidated bin, and it should be the B version (higher score)
    assert n == 2  # bin_ct ends at 2 after writing bin.1
    assert (out / "bin_001.fasta").exists()
    # B's contigs are identical here; assert stats line carried B's completeness (95)
    stats_out = (tmp_path / "binsM.stats.tsv").read_text()
    assert "95.0" in stats_out


# --- the single-bin-set path -----------------------------------------------------------------
#
# With two or more bin sets the thresholds are applied as a side effect of consolidating. With exactly
# one there is nothing to consolidate, and the final folder - named metawrap_70_10_bins - was being
# filled with every input bin. On a COMEBin run that meant 376 bins in a folder claiming 70/10, of
# which 178 actually passed.


def test_select_passing_keeps_only_bins_over_the_threshold(tmp_path):
    folder = tmp_path / "binsA"
    folder.mkdir()
    for name in ("good", "incomplete", "contaminated", "junk"):
        _write_bin(folder / ("%s.fasta" % name), [("c_%s" % name, 1000)])
    stats = tmp_path / "binsA.stats.tsv"
    _write_stats(
        stats,
        [
            ("good", 95.0, 2.0, 1000),
            ("incomplete", 40.0, 1.0, 1000),
            ("contaminated", 98.0, 55.0, 1000),
            ("junk", 0.0, 100.0, 1000),
        ],
    )
    out = tmp_path / "metawrap_70_10_bins"
    kept = refinement.select_passing(str(folder), str(stats), str(out), 70.0, 10.0)

    assert kept == 1
    assert sorted(p.name for p in out.iterdir()) == ["bin_001.fasta"]
    assert (out / "bin_001.fasta").read_text() == (folder / "good.fasta").read_text()


def test_select_passing_renumbers_and_carries_the_scores_across(tmp_path):
    folder = tmp_path / "binsA"
    folder.mkdir()
    for name in ("z_best", "a_worst", "mid"):
        _write_bin(folder / ("%s.fasta" % name), [("c_%s" % name, 1000)])
    stats = tmp_path / "binsA.stats.tsv"
    _write_stats(
        stats,
        [
            ("z_best", 99.0, 1.0, 1000),
            ("a_worst", 71.0, 9.0, 1000),
            ("mid", 85.0, 4.0, 1000),
        ],
    )
    out = tmp_path / "final"
    assert refinement.select_passing(str(folder), str(stats), str(out), 70.0, 10.0) == 3

    lines = (tmp_path / "final.stats.tsv").read_text().splitlines()
    assert lines[0].startswith("bin\t"), "the header is carried over"
    names = [ln.split("\t")[0] for ln in lines[1:]]
    assert names == ["bin_001", "bin_002", "bin_003"], "renumbered, as consolidation does"
    # a_worst sorts first, so its row must be the one under bin_001
    assert lines[1].split("\t")[1:3] == ["71.0", "9.0"]


def test_select_passing_uses_the_same_boundary_as_consolidation(tmp_path):
    """A bin exactly on a threshold has to be treated the same by both paths, whichever way that is."""
    folder = tmp_path / "binsA"
    folder.mkdir()
    _write_bin(folder / "edge.fasta", [("c1", 1000)])
    stats = tmp_path / "binsA.stats.tsv"
    _write_stats(stats, [("edge", 70.0, 10.0, 1000)])
    out = tmp_path / "final"
    kept = refinement.select_passing(str(folder), str(stats), str(out), 70.0, 10.0)
    assert kept == (1 if refinement.passes(70.0, 10.0, 70.0, 10.0) else 0)


def test_select_passing_ignores_a_stats_row_with_no_bin_behind_it(tmp_path):
    """A stats table merged from a cache can name a bin that is not in this folder."""
    folder = tmp_path / "binsA"
    folder.mkdir()
    _write_bin(folder / "here.fasta", [("c1", 1000)])
    stats = tmp_path / "binsA.stats.tsv"
    _write_stats(stats, [("here", 90.0, 2.0, 1000), ("elsewhere", 95.0, 1.0, 1000)])
    out = tmp_path / "final"
    assert refinement.select_passing(str(folder), str(stats), str(out), 70.0, 10.0) == 1


def test_select_passing_writes_an_empty_set_rather_than_failing(tmp_path):
    folder = tmp_path / "binsA"
    folder.mkdir()
    _write_bin(folder / "junk.fasta", [("c1", 1000)])
    stats = tmp_path / "binsA.stats.tsv"
    _write_stats(stats, [("junk", 5.0, 90.0, 1000)])
    out = tmp_path / "final"
    assert refinement.select_passing(str(folder), str(stats), str(out), 70.0, 10.0) == 0
    assert out.is_dir() and not list(out.iterdir())
    assert (tmp_path / "final.stats.tsv").read_text().startswith("bin\t")
