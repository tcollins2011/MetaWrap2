"""Assessing a variant's bins: the format conversions, and the rule that stops double counting.

The tests that matter most here are the ones about arithmetic that would otherwise be wrong in a
plausible-looking way. A variant that ran three binners leaves three overlapping bin sets, and the
ranker sums every row it finds for a variant - so recording all three as the variant's own numbers
would report the same genome three times over, and an F1 score above 1. Nothing about the resulting
table looks broken, which is precisely why it is worth a test.
"""

import os

from metawrap2.benchmark import assess, scoring
from metawrap2.commands import benchmark as cmd


def write_bin(path, *contigs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.writelines(
            ">%s some length=500 description\n%s\n" % (name, "ACGT" * 10) for name in contigs
        )


def bin_set(root, name, bins):
    """A directory of bins: ``bins`` maps bin file stem to its contig names."""
    for stem, contigs in bins.items():
        write_bin(os.path.join(root, name, stem + ".fasta"), *contigs)
    return os.path.join(root, name)


def text_of(path):
    with open(path) as fh:
        return fh.read()


def stats_table(path, rows):
    """A MetaWrap2 ``.stats.tsv`` - what either CheckM's output is normalised into."""
    with open(path, "w") as fh:
        fh.write("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n")
        fh.writelines(
            "%s\t%s\t%s\t0.5\tk__Bacteria\t1000\t2000\n" % (name, completeness, contamination)
            for name, completeness, contamination in rows
        )


# ── bins to bioboxes and the gold standard ─────────────────────────────────────────────────────


def test_a_bin_directory_becomes_a_bioboxes_binning_file(tmp_path):
    """A contig id is cut at whitespace, and unbinned contigs are left out.

    The gold standard names the contig alone; a bin's header carries the assembler's description.
    Leaving the description on scores every contig as one the truth has never heard of, which AMBER
    reports as a purity of zero rather than as an error - a wrong number that looks like a result.
    """
    root = str(tmp_path)
    bins = bin_set(root, "metabat2_bins", {"bin_001": ["c1", "c2"], "bin_002": ["c3"]})
    write_bin(os.path.join(bins, "unbinned.fasta"), "c9")
    out = assess.bins_to_bioboxes(bins, os.path.join(root, "q.binning"))

    text = text_of(out)
    assert "@@SEQUENCEID\tBINID" in text
    pairs = [line.split("\t") for line in text.splitlines() if line and not line.startswith("@")]
    assert pairs == [["c1", "bin_001"], ["c2", "bin_001"], ["c3", "bin_002"]]
    assert "description" not in text
    assert "c9" not in text


def test_cami2s_anonymous_contig_ids_are_matched_to_the_names_its_assembly_uses(tmp_path):
    """A CAMI2 gold standard and its own gold assembly name the same contig differently.

    The mapping leads with an anonymous ``PC0`` that appears nowhere in ``gsa.fa``, which names the
    contig after the region it was cut from and keeps the pieces of that name in the mapping's
    ``contig_id``/``start_position``/``end_position`` columns. On the gastrointestinal challenge not
    one of 210,618 anonymous ids matched the 210,633-contig assembly, so every binning comparison
    would have died with "yielded no contig-to-genome pairs" - the dataset appearing to disagree with
    itself. The total is inclusive of both ends, so it is ``end - start + 1``.
    """
    gold = tmp_path / "mapping.tsv"
    gold.write_text(
        "#anonymous_contig_id\tgenome_id\ttax_id\tcontig_id\tnumber_reads\t"
        "start_position\tend_position\n"
        "PC0\tOTU_97.2576.0\t1501\tCP013019.1\t29013.0\t3\t4351861\n"
        "PC1\tOTU_97.326.0\t1492\tCP013252.1\t25500.0\t1\t3824894\n"
    )
    out = assess.gold_to_bioboxes(
        str(gold),
        str(tmp_path / "gold.binning"),
        lengths={
            "CP013019.1_from_3_to_4351861_total_4351859": 4351859,
            "CP013252.1_from_1_to_3824894_total_3824894": 3824894,
        },
    )
    pairs = [ln for ln in text_of(out).splitlines() if ln and not ln.startswith("@")]
    assert pairs == [
        "CP013019.1_from_3_to_4351861_total_4351859\tOTU_97.2576.0\t4351859",
        "CP013252.1_from_1_to_3824894_total_3824894\tOTU_97.326.0\t3824894",
    ]


def test_two_bin_sets_with_the_same_directory_name_are_scored_as_two_queries(tmp_path, monkeypatch):
    """Scoring two refinement arms gave two identical rows, because the queries were one file.

    bin_refinement names its final directory for the thresholds it was cut at, so every arm refined at
    ``-c 70 -x 10`` produces a ``metawrap_70_10_bins``. The converted query took that basename, so the
    second conversion wrote over the first and AMBER was handed the same path twice under two labels -
    reporting F1 0.3612 for a 74-bin set and a 60-bin set alike, with nothing to say they were the same
    file. The label is what tells the rows apart, so it names the query too.
    """
    gold = tmp_path / "gold.binning"
    gold.write_text("@Version:0.9.1\n@SampleID:tiny\n\n@@SEQUENCEID\tBINID\nc1\tgenome_A\n")
    first, second = tmp_path / "arm1", tmp_path / "arm2"
    bin_set(str(first), "metawrap_70_10_bins", {"bin_001": ["c1"]})
    bin_set(str(second), "metawrap_70_10_bins", {"bin_002": ["c1"]})

    calls = []
    monkeypatch.setattr(assess, "_run", lambda cmd, **kw: calls.append(cmd) or 0)
    assess.run_amber(
        [str(first / "metawrap_70_10_bins"), str(second / "metawrap_70_10_bins")],
        str(gold),
        str(tmp_path / "amber"),
        labels=["metabat2_maxbin2", "metabat2_concoct"],
        env=None,
    )
    queries = [arg for arg in calls[0] if arg.endswith(".binning") and "gold" not in arg]
    assert len(set(queries)) == 2, queries
    # and each query holds its own arm's bins, not the other's
    assert "bin_001" in text_of(queries[0]) and "bin_002" in text_of(queries[1])
    # AMBER pairs query to truth by sample id and exits on anything it cannot pair, so each query
    # declares the gold standard's sample id, not its bin directory's name.
    for query in queries:
        assert "@SampleID:tiny" in text_of(query)


# ── the rule that stops double counting ───────────────────────────────────────────────────────


def test_three_binners_bins_are_not_summed_into_one_variants_count(tmp_path):
    """The regression this whole rule exists for.

    Each of three binners recovered the same two good genomes. The variant recovered two, not six.
    """
    root = str(tmp_path / "metabat2_maxbin2_concoct")
    for name in ("metabat2_bins", "maxbin2_bins", "concoct_bins"):
        bin_set(root, name, {"bin_001": ["c1"], "bin_002": ["c2"]})
        stats_table(
            os.path.join(root, name + ".stats.tsv"), [("bin_001", 95, 1), ("bin_002", 93, 2)]
        )

    rows = cmd._collect_scores(
        "cami2_gi",
        "binner_set",
        "metabat2_maxbin2_concoct",
        root,
        ranked_bin_set=cmd._ranked_bin_set(assess.bin_dirs(root)),
    )
    plain = [r for r in rows if r.metric == "bins_near_complete"]
    assert plain == [], "with no consolidated set, no bin set may claim the variant's metric"
    suffixed = {r.metric for r in rows if r.metric.startswith("bins_near_complete__")}
    assert suffixed == {
        "bins_near_complete__metabat2_bins",
        "bins_near_complete__maxbin2_bins",
        "bins_near_complete__concoct_bins",
    }
    ratios = scoring.value_against_cost(
        rows
        + [scoring.Row("cami2_gi", "binner_set", "metabat2_maxbin2_concoct", "cpu_hours", 10.0)],
        "bins_near_complete",
    )
    assert ratios == [] or all(r.value == 0 for r in ratios)

    # Once a consolidated set is present, it - and only it - carries the variant's metric.
    bin_set(root, "metawrap_70_10_bins", {"bin_001": ["c1"], "bin_002": ["c2"], "bin_003": ["c3"]})
    stats_table(
        os.path.join(root, "metawrap_70_10_bins.stats.tsv"),
        [("bin_001", 95, 1), ("bin_002", 91, 3), ("bin_003", 92, 1)],
    )
    rows = cmd._collect_scores(
        "cami2_gi",
        "binner_set",
        "metabat2_maxbin2_concoct",
        root,
        ranked_bin_set=cmd._ranked_bin_set(assess.bin_dirs(root)),
    )
    plain = [r for r in rows if r.metric == "bins_near_complete"]
    assert [r.value for r in plain] == [3.0], "the refined set is the variant's number"
