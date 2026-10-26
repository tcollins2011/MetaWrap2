"""An empty input bin set is dropped with a warning instead of stopping bin_refinement.

A binner that finds nothing writes an empty folder and exits 0. The other sets are still refined.
With every set empty the run finishes with an empty final set and a warning.
"""

from metawrap2.modules import bin_refinement


def _bin_set(path, n):
    path.mkdir()
    for i in range(n):
        (path / ("bin.%d.fa" % i)).write_text(">c%d\nACGT\n" % i)
    return str(path)


def test_an_empty_bin_set_is_dropped_and_the_rest_are_kept(tmp_path, capsys):
    full = _bin_set(tmp_path / "semibin2_bins", 2)
    empty = _bin_set(tmp_path / "metabat2_bins", 0)
    assert bin_refinement.drop_empty_bin_sets([empty, full]) == [full]
    # The banner wraps long lines and pads them with asterisks, so a long tmp_path can split the
    # phrase. Remove the decoration before matching.
    out = capsys.readouterr()
    text = " ".join((out.out + out.err).replace("*", " ").split())
    assert "metabat2_bins" in text and "has no bins in it" in text


def test_a_folder_of_other_files_counts_as_empty(tmp_path):
    full = _bin_set(tmp_path / "a", 1)
    other = tmp_path / "b"
    other.mkdir()
    (other / "log.txt").write_text("nothing binned\n")
    assert bin_refinement.drop_empty_bin_sets([full, str(other)]) == [full]


def test_a_missing_path_is_left_for_the_input_checks(tmp_path):
    """A nonexistent path is passed through so the input checks report it as missing."""
    full = _bin_set(tmp_path / "a", 1)
    missing = str(tmp_path / "does_not_exist")
    assert bin_refinement.drop_empty_bin_sets([full, missing]) == [full, missing]


def test_every_set_empty_leaves_nothing_to_refine_and_warns(tmp_path, capsys):
    empties = [_bin_set(tmp_path / name, 0) for name in ("a", "b")]
    assert bin_refinement.drop_empty_bin_sets(empties) == []
    out = capsys.readouterr()
    assert "empty final bin set" in out.out + out.err


def test_a_run_with_only_empty_sets_finishes_with_an_empty_final_set(tmp_path):
    a = _bin_set(tmp_path / "a", 0)
    b = _bin_set(tmp_path / "b", 0)
    out = tmp_path / "REFINED"
    assert bin_refinement.main(["-o", str(out), "-A", a, "-B", b, "--skip-gunc"]) == 0
    final = out / "metawrap_70_10_bins"
    assert final.is_dir() and not list(final.iterdir())
    header = (out / "metawrap_70_10_bins.stats.tsv").read_text().splitlines()
    assert header == ["\t".join(bin_refinement._checkm.STATS_HEADER)]


def test_nothing_changes_when_every_set_has_bins(tmp_path, capsys):
    sets = [_bin_set(tmp_path / name, 1) for name in ("a", "b", "c")]
    assert bin_refinement.drop_empty_bin_sets(sets) == sets
    assert "relabelled" not in capsys.readouterr().out
