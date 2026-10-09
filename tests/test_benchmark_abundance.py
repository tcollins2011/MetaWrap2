"""Scoring an abundance table: the truth a bin is held to, and the metrics that survive three units.

The arithmetic worth pinning here is the second factor in the expected abundance - a bin is part of a
genome, so it can only account for its own share of that genome's reads - and the fact that every
metric is scale-free. A quantifier reporting coverage and one reporting a percentage have to be
comparable, or the benchmark cannot choose between them at all.
"""

import gzip
import json
import math
import os

import pytest

from metawrap2.benchmark import abundance


def write_mapping(path, rows, gz=False):
    """A CAMI2 per-sample read mapping: rows are ``(read id, genome)``."""
    text = "#anonymous_read_id\tgenome_id\ttax_id\tread_id\n" + "".join(
        "%s\t%s\t1\tsrc\n" % row for row in rows
    )
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        with open(path, "w") as fh:
            fh.write(text)
    return str(path)


def write_bin(path, contigs):
    with open(path, "w") as fh:
        fh.writelines(">%s\n%s\n" % (name, "A" * length) for name, length in contigs)


def write_table(path, samples, rows):
    with open(path, "w") as fh:
        fh.write("Genomic bins\t%s\n" % "\t".join(samples))
        fh.writelines(
            "%s\t%s\n" % (name, "\t".join(str(v) for v in values)) for name, values in rows
        )
    return str(path)


# ── the truth ─────────────────────────────────────────────────────────────────────────────────


def test_read_fractions_counts_pairs_not_mates(tmp_path):
    """Both mates come from the same genome, so counting both would just double every count.

    It would also weight a pair whose mate is missing from the truth differently from a complete one,
    and "how much of this sample is this organism" is a question about pairs.
    """
    mapping = write_mapping(
        tmp_path / "m.tsv",
        [("S0R0/1", "g_A"), ("S0R0/2", "g_A"), ("S0R1/1", "g_B"), ("S0R1/2", "g_B")],
    )
    assert abundance.read_fractions(mapping) == {"g_A": 0.5, "g_B": 0.5}


def test_read_fractions_reads_the_gzipped_truth_cami2_ships(tmp_path):
    mapping = write_mapping(
        tmp_path / "m.tsv.gz", [("S0R0/1", "g_A"), ("S0R1/1", "g_A"), ("S0R2/1", "g_B")], gz=True
    )
    fractions = abundance.read_fractions(mapping)
    assert math.isclose(fractions["g_A"], 2 / 3)
    assert math.isclose(fractions["g_B"], 1 / 3)


def test_a_truth_file_with_no_forward_reads_says_so(tmp_path):
    mapping = write_mapping(tmp_path / "m.tsv", [("S0R0/2", "g_A")])
    with pytest.raises(abundance.AbundanceError, match="no abundance truth"):
        abundance.read_fractions(mapping)


def test_a_bins_composition_is_measured_with_the_truths_own_lengths(tmp_path):
    bins = tmp_path / "bins"
    bins.mkdir()
    write_bin(bins / "bin_1.fasta", [("c1", 10), ("c2", 10), ("c9", 10)])
    composition = abundance.bin_composition(
        str(bins),
        {"c1": "g_A", "c2": "g_A", "c3": "g_B"},
        {"c1": 1000, "c2": 500, "c3": 100},
    )
    # c9 belongs to no genome in the truth, so it is counted against none - not against g_A.
    assert composition == {"bin_1": {"g_A": 1500}}


def test_a_bin_expects_its_genomes_abundance_scaled_by_how_much_of_it_the_bin_holds(tmp_path):
    """The factor that stops a quantifier being marked down for the binner's incompleteness.

    Half of a genome that is 10% of the sample can account for about 5% of the sample's sequence, and
    no more: the reads from the other half cannot be attributed to this bin by any quantifier. Holding
    the bin to the genome's full 10% would score every incomplete bin as an underestimate.
    """
    composition = {"bin_1": {"g_A": 500}}
    expected = abundance.expected_fractions(composition, {"g_A": 1000}, {"g_A": 0.10})
    assert math.isclose(expected["bin_1"], 0.05)


def test_a_bin_is_credited_with_its_majority_genome_only(tmp_path):
    """A little of a second genome is contamination, not a second abundance to be right about."""
    composition = {"bin_1": {"g_A": 900, "g_B": 100}}
    expected = abundance.expected_fractions(
        composition, {"g_A": 1000, "g_B": 1000}, {"g_A": 0.2, "g_B": 0.8}
    )
    assert math.isclose(expected["bin_1"], 0.2 * 0.9)


def test_a_bin_whose_genome_is_absent_from_the_sample_expects_zero():
    """Zero is a real expectation: reporting a genome that is not there is a quantifier's error."""
    expected = abundance.expected_fractions({"bin_1": {"g_A": 500}}, {"g_A": 1000}, {"g_B": 1.0})
    assert expected == {"bin_1": 0.0}


def test_a_bin_holding_nothing_the_truth_knows_expects_zero_rather_than_crashing():
    assert abundance.expected_fractions({"bin_1": {}}, {}, {"g_A": 1.0}) == {"bin_1": 0.0}


# ── reading the table ─────────────────────────────────────────────────────────────────────────


def test_salmons_nan_is_read_as_no_abundance_not_as_a_nan(tmp_path):
    """One NaN poisons a correlation for the whole sample, turning a rough cell into no result."""
    table = write_table(tmp_path / "t.tsv", ["s1"], [("bin_1", ["nan"]), ("bin_2", [3.0])])
    samples, rows = abundance.read_abundance_table(table)
    assert samples == ["s1"]
    assert rows == {"bin_1": [0.0], "bin_2": [3.0]}


def test_a_table_with_no_sample_columns_is_an_error(tmp_path):
    path = tmp_path / "t.tsv"
    path.write_text("Genomic bins\nbin_1\n")
    with pytest.raises(abundance.AbundanceError, match="no sample columns"):
        abundance.read_abundance_table(str(path))


# ── the metrics ───────────────────────────────────────────────────────────────────────────────


def test_a_perfect_estimate_on_a_different_scale_still_scores_perfectly(tmp_path):
    """The whole point: salmon's coverage and CoverM's percentage must be comparable.

    An estimate that is the truth times a constant is a perfect answer in a different unit, and a
    metric that called it imperfect would rank the quantifiers by which units they happen to use.
    """
    expected = {"s1": {"bin_1": 0.1, "bin_2": 0.2, "bin_3": 0.7}}
    table = write_table(
        tmp_path / "t.tsv", ["s1"], [("bin_1", [50.0]), ("bin_2", [100.0]), ("bin_3", [350.0])]
    )
    metrics = abundance.score(table, expected)
    assert math.isclose(metrics["abundance_spearman"], 1.0)
    assert math.isclose(metrics["abundance_pearson"], 1.0)
    assert math.isclose(metrics["abundance_l1"], 0.0, abs_tol=1e-12)
    assert math.isclose(metrics["abundance_l1_similarity"], 1.0)


def test_an_estimate_in_the_wrong_order_loses_the_rank_correlation(tmp_path):
    expected = {"s1": {"bin_1": 0.1, "bin_2": 0.2, "bin_3": 0.7}}
    table = write_table(
        tmp_path / "t.tsv", ["s1"], [("bin_1", [0.7]), ("bin_2", [0.2]), ("bin_3", [0.1])]
    )
    metrics = abundance.score(table, expected)
    assert math.isclose(metrics["abundance_spearman"], -1.0)


def test_the_l1_similarity_runs_the_right_way_for_the_ranker(tmp_path):
    """value_against_cost takes bigger as better, and L1 is a distance."""
    expected = {"s1": {"bin_1": 1.0, "bin_2": 0.0}}
    table = write_table(tmp_path / "t.tsv", ["s1"], [("bin_1", [0.0]), ("bin_2", [1.0])])
    metrics = abundance.score(table, expected)
    assert math.isclose(metrics["abundance_l1"], 2.0), "opposite compositions are the L1 maximum"
    assert math.isclose(metrics["abundance_l1_similarity"], 0.0)


def test_several_samples_are_averaged_not_pooled(tmp_path):
    """Pooling puts a deep sample's bins and a shallow one's on one axis, where the correlation
    measures the difference between the samples rather than the quantifier's accuracy."""
    expected = {
        "s1": {"bin_1": 0.1, "bin_2": 0.9},
        "s2": {"bin_1": 0.9, "bin_2": 0.1},
    }
    table = write_table(
        tmp_path / "t.tsv", ["s1", "s2"], [("bin_1", [1.0, 9.0]), ("bin_2", [9.0, 1.0])]
    )
    metrics = abundance.score(table, expected)
    assert math.isclose(metrics["abundance_spearman"], 1.0)


def test_a_sample_the_truth_says_nothing_about_is_skipped_not_scored_as_zero(tmp_path):
    expected = {"s1": {"bin_1": 0.1, "bin_2": 0.9}}
    table = write_table(
        tmp_path / "t.tsv", ["s1", "s_unknown"], [("bin_1", [1.0, 5.0]), ("bin_2", [9.0, 5.0])]
    )
    metrics = abundance.score(table, expected)
    assert math.isclose(metrics["abundance_spearman"], 1.0)


def test_a_table_whose_samples_are_all_unknown_to_the_truth_is_an_error(tmp_path):
    table = write_table(tmp_path / "t.tsv", ["other"], [("bin_1", [1.0]), ("bin_2", [2.0])])
    with pytest.raises(abundance.AbundanceError, match="appear in the truth"):
        abundance.score(table, {"s1": {"bin_1": 0.5, "bin_2": 0.5}})


def test_a_flat_estimate_is_no_correlation_rather_than_a_crash(tmp_path):
    """Every bin reported identically has no variance, so the correlation is undefined - but the L1
    is not, and losing the whole sample would hide a quantifier that reported nothing useful."""
    expected = {"s1": {"bin_1": 0.1, "bin_2": 0.9}}
    table = write_table(tmp_path / "t.tsv", ["s1"], [("bin_1", [5.0]), ("bin_2", [5.0])])
    metrics = abundance.score(table, expected)
    assert "abundance_spearman" not in metrics
    assert math.isclose(metrics["abundance_l1"], 0.8)


def test_ties_are_ranked_by_their_average(tmp_path):
    assert abundance._ranks([3.0, 1.0, 1.0, 5.0]) == [3.0, 1.5, 1.5, 4.0]


def test_the_real_cami2_truth_and_a_real_bin_set_line_up():
    """A guard on the join, not on the arithmetic: label and sample keys agreeing is the whole thing.

    Skipped where the benchmark data is not downloaded, which is most machines.
    """
    built = "runs/abund_bench/expected.json"
    if not os.path.isfile(built):
        pytest.skip("the cami2_gi abundance truth has not been built on this machine")
    with open(built) as fh:
        expected = json.load(fh)
    assert set(expected) == {"sample_0", "sample_1", "sample_10"}
    for by_bin in expected.values():
        assert len(by_bin) == 76
        assert all(label.startswith("bin_") for label in by_bin)


def test_a_bin_the_quantifier_left_out_is_scored_as_zero_not_dropped(tmp_path):
    """Scoring used to intersect the table with the truth, which rewarded declining to answer.

    salmon writes no row at all for a bin it attributed nothing to, while the CoverM path pre-seeds
    a row per bin - so on the same bins salmon was scored over the subset it answered for and
    CoverM over all of them, and a correlation over 2 bins is not comparable with one over 4.
    A bin missing from the table means "no abundance detected", which is the zero its own `nan`
    cells are already read as.
    """
    table = tmp_path / "bin_abundance_table.tsv"
    # bin.4 is absent: the quantifier attributed nothing to it. The truth says it is the dominant
    # genome in the sample, so omitting it should cost the arm, not be excused.
    table.write_text("Genomic bins\tsample_0\n" "bin.1\t0.1\n" "bin.2\t0.2\n" "bin.3\t0.3\n")
    truth = {"sample_0": {"bin.1": 0.1, "bin.2": 0.2, "bin.3": 0.3, "bin.4": 0.9}}
    scored = abundance.score(str(table), truth)
    # The omission is reported rather than absorbed...
    assert scored["abundance_missing_bins"] == 1.0
    # ...and it costs the arm: scored over the intersection the three rows correlate perfectly.
    assert scored["abundance_spearman"] < 1.0

    complete = tmp_path / "complete.tsv"
    complete.write_text(
        "Genomic bins\tsample_0\n" "bin.1\t0.1\n" "bin.2\t0.2\n" "bin.3\t0.3\n" "bin.4\t0.9\n"
    )
    full = abundance.score(str(complete), truth)
    assert "abundance_missing_bins" not in full
    assert full["abundance_spearman"] == 1.0
