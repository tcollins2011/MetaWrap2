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
# With two or more bin sets, consolidation applies the thresholds. With one set there is nothing
# to consolidate, so select_passing applies them. Without it, a COMEBin run put all 376 bins into
# metawrap_70_10_bins, of which 178 passed 70/10.


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
