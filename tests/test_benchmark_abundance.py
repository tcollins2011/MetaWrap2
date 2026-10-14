"""Scoring an abundance table: the truth a bin is held to, and the metrics that survive three units.

The arithmetic worth pinning here is the second factor in the expected abundance - a bin is part of a
genome, so it can only account for its own share of that genome's reads - and the fact that every
metric is scale-free. A quantifier reporting coverage and one reporting a percentage have to be
comparable, or the benchmark cannot choose between them at all.
"""

import gzip
import math

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


def write_table(path, samples, rows):
    with open(path, "w") as fh:
        fh.write("Genomic bins\t%s\n" % "\t".join(samples))
        fh.writelines(
            "%s\t%s\n" % (name, "\t".join(str(v) for v in values)) for name, values in rows
        )
    return str(path)


# ── the truth ─────────────────────────────────────────────────────────────────────────────────


def test_a_bin_expects_its_genomes_abundance_scaled_by_how_much_of_it_the_bin_holds(tmp_path):
    """The factor that stops a quantifier being marked down for the binner's incompleteness.

    Half of a genome that is 10% of the sample can account for about 5% of the sample's sequence, and
    no more: the reads from the other half cannot be attributed to this bin by any quantifier. Holding
    the bin to the genome's full 10% would score every incomplete bin as an underestimate.
    """
    # Pairs, not mates: both mates come from the same genome, so counting both would double every
    # count and weight a pair with a missing mate differently from a complete one.
    mapping = write_mapping(
        tmp_path / "m.tsv",
        [("S0R0/1", "g_A"), ("S0R0/2", "g_A"), ("S0R1/1", "g_B"), ("S0R1/2", "g_B")],
    )
    assert abundance.read_fractions(mapping) == {"g_A": 0.5, "g_B": 0.5}

    composition = {"bin_1": {"g_A": 500}}
    expected = abundance.expected_fractions(composition, {"g_A": 1000}, {"g_A": 0.10})
    assert math.isclose(expected["bin_1"], 0.05)

    # A little of a second genome is contamination, not a second abundance to be right about.
    composition = {"bin_1": {"g_A": 900, "g_B": 100}}
    expected = abundance.expected_fractions(
        composition, {"g_A": 1000, "g_B": 1000}, {"g_A": 0.2, "g_B": 0.8}
    )
    assert math.isclose(expected["bin_1"], 0.2 * 0.9)

    # Zero is a real expectation: reporting a genome that is not there is a quantifier's error.
    assert abundance.expected_fractions({"bin_1": {"g_A": 500}}, {"g_A": 1000}, {"g_B": 1.0}) == {
        "bin_1": 0.0
    }


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
