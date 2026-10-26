"""Scoring an abundance table: the expected abundance of a bin, and scale-free metrics.

Checks that a bin's expected abundance is scaled by its share of its genome's reads, and that every
metric is scale-free, so quantifiers reporting coverage and percentages can be compared.
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


# ── expected abundance ────────────────────────────────────────────────────────────────────────


def test_a_bin_expects_its_genomes_abundance_scaled_by_how_much_of_it_the_bin_holds(tmp_path):
    """A bin holding half of a genome at 10% of the sample expects about 5%.

    Without this scaling every incomplete bin would score as an underestimate.
    """
    # Counted per pair: both mates come from the same genome, and counting mates would weight a
    # pair with a missing mate differently.
    mapping = write_mapping(
        tmp_path / "m.tsv",
        [("S0R0/1", "g_A"), ("S0R0/2", "g_A"), ("S0R1/1", "g_B"), ("S0R1/2", "g_B")],
    )
    assert abundance.read_fractions(mapping) == {"g_A": 0.5, "g_B": 0.5}

    composition = {"bin_1": {"g_A": 500}}
    expected = abundance.expected_fractions(composition, {"g_A": 1000}, {"g_A": 0.10})
    assert math.isclose(expected["bin_1"], 0.05)

    # A small amount of a second genome is treated as contamination.
    composition = {"bin_1": {"g_A": 900, "g_B": 100}}
    expected = abundance.expected_fractions(
        composition, {"g_A": 1000, "g_B": 1000}, {"g_A": 0.2, "g_B": 0.8}
    )
    assert math.isclose(expected["bin_1"], 0.2 * 0.9)

    # Zero is a valid expectation; reporting an absent genome is a quantifier error.
    assert abundance.expected_fractions({"bin_1": {"g_A": 500}}, {"g_A": 1000}, {"g_B": 1.0}) == {
        "bin_1": 0.0
    }


# ── the metrics ───────────────────────────────────────────────────────────────────────────────


def test_a_perfect_estimate_on_a_different_scale_still_scores_perfectly(tmp_path):
    """An estimate equal to the truth times a constant scores as perfect.

    This keeps salmon's coverage and CoverM's percentage comparable.
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
    """A bin missing from the table is scored as zero abundance, like a `nan` cell.

    salmon writes no row for a bin it attributed nothing to, while the CoverM path writes a row per
    bin. Scoring only the intersection made the two incomparable and favoured omitted bins.
    """
    table = tmp_path / "bin_abundance_table.tsv"
    # bin.4 is absent from the table but is the dominant genome in the truth.
    table.write_text("Genomic bins\tsample_0\n" "bin.1\t0.1\n" "bin.2\t0.2\n" "bin.3\t0.3\n")
    truth = {"sample_0": {"bin.1": 0.1, "bin.2": 0.2, "bin.3": 0.3, "bin.4": 0.9}}
    scored = abundance.score(str(table), truth)
    # The omission is reported...
    assert scored["abundance_missing_bins"] == 1.0
    # ...and lowers the score (over the intersection the three rows correlate perfectly).
    assert scored["abundance_spearman"] < 1.0

    complete = tmp_path / "complete.tsv"
    complete.write_text(
        "Genomic bins\tsample_0\n" "bin.1\t0.1\n" "bin.2\t0.2\n" "bin.3\t0.3\n" "bin.4\t0.9\n"
    )
    full = abundance.score(str(complete), truth)
    assert "abundance_missing_bins" not in full
    assert full["abundance_spearman"] == 1.0
