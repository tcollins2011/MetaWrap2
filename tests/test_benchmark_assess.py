"""Assessing a variant's bins: the format conversions, and the rule that stops double counting.

The tests that matter most here are the ones about arithmetic that would otherwise be wrong in a
plausible-looking way. A variant that ran three binners leaves three overlapping bin sets, and the
ranker sums every row it finds for a variant - so recording all three as the variant's own numbers
would report the same genome three times over, and an F1 score above 1. Nothing about the resulting
table looks broken, which is precisely why it is worth a test.
"""

import os

import pytest

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


# ── bins to bioboxes ──────────────────────────────────────────────────────────────────────────


def test_a_bin_directory_becomes_a_bioboxes_binning_file(tmp_path):
    root = str(tmp_path)
    bins = bin_set(root, "metabat2_bins", {"bin_001": ["c1", "c2"], "bin_002": ["c3"]})
    out = assess.bins_to_bioboxes(bins, os.path.join(root, "q.binning"))

    text = text_of(out)
    assert "@@SEQUENCEID\tBINID" in text
    pairs = [line.split("\t") for line in text.splitlines() if line and not line.startswith("@")]
    assert pairs == [["c1", "bin_001"], ["c2", "bin_001"], ["c3", "bin_002"]]


def test_a_contig_id_is_cut_at_whitespace_not_taken_whole(tmp_path):
    """The gold standard names the contig alone; a bin's header carries the assembler's description.

    Leaving the description on scores every contig as one the truth has never heard of, which AMBER
    reports as a purity of zero rather than as an error - a wrong number that looks like a result.
    """
    bins = bin_set(str(tmp_path), "concoct_bins", {"bin_001": ["PC7"]})
    out = assess.bins_to_bioboxes(bins, str(tmp_path / "q.binning"))
    assert "PC7\tbin_001" in text_of(out)
    assert "description" not in text_of(out)


def test_unbinned_contigs_are_left_out_rather_than_made_into_a_genome(tmp_path):
    root = str(tmp_path)
    bins = bin_set(root, "metabat2_bins", {"bin_001": ["c1"]})
    write_bin(os.path.join(bins, "unbinned.fasta"), "c9")
    out = assess.bins_to_bioboxes(bins, os.path.join(root, "q.binning"))
    text = text_of(out)
    assert "c1\tbin_001" in text
    assert "c9" not in text


def test_an_empty_bin_directory_says_so(tmp_path):
    empty = str(tmp_path / "empty_bins")
    os.makedirs(empty)
    with pytest.raises(assess.AssessError):
        assess.bins_to_bioboxes(empty, str(tmp_path / "q.binning"))


# ── the gold standard ─────────────────────────────────────────────────────────────────────────


def test_a_gold_standard_that_already_carries_lengths_is_used_as_it_is(tmp_path):
    gold = tmp_path / "gold.binning"
    gold.write_text(
        "@Version:0.9.1\n@SampleID:x\n\n@@SEQUENCEID\tBINID\t_LENGTH\nc1\tgenome_A\t40\n"
    )
    out = assess.gold_to_bioboxes(str(gold), str(tmp_path / "converted.binning"))
    assert out == str(gold), "a file that needs no conversion should not be copied"


def test_a_bioboxes_gold_standard_without_lengths_has_them_added(tmp_path):
    """AMBER takes sequence lengths from the gold standard alone, and dies in pandas without them.

    A bioboxes file is otherwise ready to use, so the fix is to add the column rather than to refuse
    the file - and the sample id it declares has to survive, because the queries are named after it.
    """
    gold = tmp_path / "gold.binning"
    gold.write_text("@Version:0.9.1\n@SampleID:tiny\n\n@@SEQUENCEID\tBINID\nc1\tgenome_A\n")
    out = assess.gold_to_bioboxes(
        str(gold), str(tmp_path / "converted.binning"), lengths={"c1": 400}
    )
    assert out != str(gold)
    text = text_of(out)
    assert "@SampleID:tiny" in text
    assert "@@SEQUENCEID\tBINID\t_LENGTH" in text
    assert "c1\tgenome_A\t400" in text


def test_a_gold_standard_with_no_lengths_anywhere_says_where_they_come_from(tmp_path):
    gold = tmp_path / "gold.binning"
    gold.write_text("@Version:0.9.1\n@SampleID:x\n\n@@SEQUENCEID\tBINID\nc1\tgenome_A\n")
    with pytest.raises(assess.AssessError, match="--assembly"):
        assess.gold_to_bioboxes(str(gold), str(tmp_path / "converted.binning"))


def test_sequence_lengths_are_read_from_the_assembly_fasta(tmp_path):
    fasta = tmp_path / "assembly.fa"
    fasta.write_text(">c1 len=8 multi=3\nACGT\nACGT\n>c2 whatever\nAC\n")
    assert assess.fasta_lengths(str(fasta)) == {"c1": 8, "c2": 2}


def test_a_cami2_mapping_table_is_converted_to_bioboxes(tmp_path):
    gold = tmp_path / "mapping.tsv"
    gold.write_text(
        "#anonymous_contig_id\tgenome_id\ttax_id\tcontig_id\n"
        "PC0\tOTU_97.2576.0\t1501\tCP013019.1\n"
        "PC1\tOTU_97.326.0\t1492\tCP013252.1\n"
    )
    out = assess.gold_to_bioboxes(
        str(gold), str(tmp_path / "gold.binning"), lengths={"PC0": 4351859, "PC1": 3824894}
    )
    text = text_of(out)
    assert "@@SEQUENCEID\tBINID\t_LENGTH" in text
    assert "PC0\tOTU_97.2576.0\t4351859" in text
    assert "1501" not in text, "taxonomy columns are not part of the binning format"


def test_a_contig_the_assembly_does_not_have_is_left_out_of_the_truth(tmp_path):
    """A contig no binner could have seen does not belong in the denominator of anybody's recall."""
    gold = tmp_path / "mapping.tsv"
    gold.write_text("#contig\tgenome\nPC0\tgenome_A\nPC9\tgenome_B\n")
    out = assess.gold_to_bioboxes(str(gold), str(tmp_path / "gold.binning"), lengths={"PC0": 500})
    pairs = [ln for ln in text_of(out).splitlines() if ln and not ln.startswith("@")]
    assert pairs == ["PC0\tgenome_A\t500"]


def test_a_contig_repeated_by_the_raw_mapping_is_not_put_in_two_bins(tmp_path):
    """CAMI2 ships a raw mapping that repeats a contig once per source region, and a fixed one.

    Pointing the harness at the raw file should give the same answer as the fixed one rather than an
    AMBER error, because which of the two is on disk is not something the person scoring chose.
    """
    gold = tmp_path / "raw.tsv"
    gold.write_text(
        "#anonymous_contig_id\tgenome_id\n"
        "PC0\tOTU_97.2576.0\n"
        "PC0\tOTU_97.2576.0\n"
        "PC1\tOTU_97.326.0\n"
    )
    out = assess.gold_to_bioboxes(
        str(gold), str(tmp_path / "gold.binning"), lengths={"PC0": 100, "PC1": 200}
    )
    pairs = [ln for ln in text_of(out).splitlines() if ln and not ln.startswith("@")]
    assert pairs == ["PC0\tOTU_97.2576.0\t100", "PC1\tOTU_97.326.0\t200"]


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


def test_an_anonymous_id_the_assembly_does_use_is_still_preferred(tmp_path):
    """The first column wins where it matches, so a dataset using anonymous ids is unaffected."""
    gold = tmp_path / "mapping.tsv"
    gold.write_text(
        "#anonymous_contig_id\tgenome_id\ttax_id\tcontig_id\tnumber_reads\t"
        "start_position\tend_position\n"
        "PC0\tgenome_A\t1\tCP1\t10.0\t1\t100\n"
    )
    out = assess.gold_to_bioboxes(
        str(gold),
        str(tmp_path / "gold.binning"),
        lengths={"PC0": 100, "CP1_from_1_to_100_total_100": 100},
    )
    pairs = [ln for ln in text_of(out).splitlines() if ln and not ln.startswith("@")]
    assert pairs == ["PC0\tgenome_A\t100"]


def test_a_mapping_from_another_dataset_is_reported_unmatched_not_renamed(tmp_path):
    """Reconstruction only applies to rows of CAMI2's shape; anything else must not be invented."""
    gold = tmp_path / "mapping.tsv"
    gold.write_text("#contig\tgenome\tnotes\nfoo\tgenome_A\tsome text\n")
    with pytest.raises(assess.AssessError, match="nothing to score against"):
        assess.gold_to_bioboxes(str(gold), str(tmp_path / "gold.binning"), lengths={"bar": 1})


def test_a_mapping_with_no_pairs_in_it_is_an_error_not_an_empty_truth(tmp_path):
    gold = tmp_path / "junk.tsv"
    gold.write_text("# nothing but a comment\n")
    with pytest.raises(assess.AssessError, match="nothing to score against"):
        assess.gold_to_bioboxes(str(gold), str(tmp_path / "gold.binning"), lengths={"PC0": 1})


def test_a_query_declares_the_gold_standards_sample_id_not_its_bin_directorys_name(
    tmp_path, monkeypatch
):
    """AMBER pairs query to truth by sample id and exits on anything it cannot pair.

    Naming a query after the directory it came from made every query unpairable, with
    ``Sample ID concoct_bins in concoct_bins not found in the genome binning gold standard``. The
    bin set is told apart by AMBER's -l label instead; the sample id belongs to the gold standard.
    """
    root = str(tmp_path)
    gold = tmp_path / "gold.binning"
    gold.write_text("@Version:0.9.1\n@SampleID:tiny\n\n@@SEQUENCEID\tBINID\nc1\tgenome_A\n")
    bin_set(root, "concoct_bins", {"bin_001": ["c1"]})
    bin_set(root, "metabat2_bins", {"bin_001": ["c1"]})

    calls = []
    monkeypatch.setattr(assess, "_run", lambda cmd, **kw: calls.append(cmd) or 0)
    assess.run_amber(assess.bin_dirs(root), str(gold), str(tmp_path / "amber"), env=None)

    queries = [arg for arg in calls[0] if arg.endswith("_bins.binning")]
    assert len(queries) == 2
    for query in queries:
        assert "@SampleID:tiny" in text_of(query)
    assert "concoct_bins,metabat2_bins" in calls[0], "the bin set is the -l label"


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


def test_duplicate_labels_are_refused_rather_than_scored(tmp_path):
    """Two rows under one name cannot be told apart, so it is a mistake worth stopping for."""
    gold = tmp_path / "gold.binning"
    gold.write_text("@Version:0.9.1\n@SampleID:tiny\n\n@@SEQUENCEID\tBINID\nc1\tgenome_A\n")
    bin_set(str(tmp_path), "a_bins", {"bin_001": ["c1"]})
    bin_set(str(tmp_path), "b_bins", {"bin_001": ["c1"]})
    with pytest.raises(assess.AssessError, match="same label more than once"):
        assess.run_amber(
            [str(tmp_path / "a_bins"), str(tmp_path / "b_bins")],
            str(gold),
            str(tmp_path / "amber"),
            labels=["arm", "arm"],
            env=None,
        )


def test_a_label_per_bin_set_or_it_is_an_error(tmp_path):
    """AMBER matches -l to queries by position, so a short list relabels every row after the gap."""
    gold = tmp_path / "gold.binning"
    gold.write_text("@Version:0.9.1\n@SampleID:tiny\n\n@@SEQUENCEID\tBINID\nc1\tgenome_A\n")
    bin_set(str(tmp_path), "a_bins", {"bin_001": ["c1"]})
    bin_set(str(tmp_path), "b_bins", {"bin_001": ["c1"]})
    with pytest.raises(assess.AssessError, match="labels for"):
        assess.run_amber(
            [str(tmp_path / "a_bins"), str(tmp_path / "b_bins")],
            str(gold),
            str(tmp_path / "amber"),
            labels=["only_one"],
            env=None,
        )


def test_a_gold_standard_with_no_sample_id_is_reported_rather_than_crashing_amber(tmp_path):
    gold = tmp_path / "gold.binning"
    gold.write_text("@Version:0.9.1\n\n@@SEQUENCEID\tBINID\nc1\tgenome_A\n")
    with pytest.raises(assess.AssessError, match="@SampleID"):
        assess.gold_sample_id(str(gold))


def test_a_multi_sample_gold_standard_is_refused_instead_of_scored_against_one_sample(tmp_path):
    """One co-assembly cannot be scored against several samples' truth: the other samples' contigs
    would count as misassigned, which reads as a bad binner rather than as the wrong gold file."""
    gold = tmp_path / "gold.binning"
    gold.write_text(
        "@Version:0.9.1\n@SampleID:sample_0\n\n@@SEQUENCEID\tBINID\nc1\tgenome_A\n\n"
        "@Version:0.9.1\n@SampleID:sample_1\n\n@@SEQUENCEID\tBINID\nc2\tgenome_B\n"
    )
    assert assess.bioboxes_sample_ids(str(gold)) == ["sample_0", "sample_1"]
    with pytest.raises(assess.AssessError, match="declares 2 samples"):
        assess.gold_sample_id(str(gold))


# ── which assessors run ───────────────────────────────────────────────────────────────────────


def test_a_trimmer_comparisons_metrics_start_neither_assessor():
    metrics = ["reads_kept", "bases_kept", "q30_bases_kept"]
    assert not assess.needs_amber(metrics)
    assert not assess.needs_checkm(metrics)


def test_a_binning_comparisons_metrics_ask_for_both():
    metrics = ["bins_near_complete", "amber_f1_bp"]
    assert assess.needs_amber(metrics)
    assert assess.needs_checkm(metrics)


def test_only_non_empty_bin_directories_are_found(tmp_path):
    root = str(tmp_path)
    bin_set(root, "metabat2_bins", {"bin_001": ["c1"]})
    os.makedirs(os.path.join(root, "maxbin2_bins"))  # ran, found nothing
    os.makedirs(os.path.join(root, "work_files"))
    assert [os.path.basename(p) for p in assess.bin_dirs(root)] == ["metabat2_bins"]


# ── the rule that stops double counting ───────────────────────────────────────────────────────


def test_a_consolidated_set_is_the_variants_answer_when_there_is_one():
    sets = ["/w/metabat2_bins", "/w/concoct_bins", "/w/metawrap_70_10_bins"]
    assert cmd._ranked_bin_set(sets) == "metawrap_70_10_bins"


def test_the_consolidated_set_is_recognised_at_whatever_thresholds_it_was_cut_at():
    """The thresholds are the user's, and the name carries them: -c 50 -x 5 is as valid.

    This used to be a check against the literal "metawrap2_bins", which bin_refinement has never
    written - so a refined run looked like several raw overlapping sets with nothing consolidated,
    and the harness declined to rank the one answer that was there. Pinning the name the module
    does not produce is what let that survive, so the fictional one is pinned as *not* matching.
    """
    assert (
        cmd._ranked_bin_set(["/w/metabat2_bins", "/w/metawrap_50_5_bins"]) == "metawrap_50_5_bins"
    )
    assert cmd._ranked_bin_set(["/w/metabat2_bins", "/w/metawrap2_bins"]) is None


def test_two_consolidated_sets_are_not_an_answer_either():
    """A run resumed at different thresholds. Picking one is the guess this rule exists to refuse."""
    sets = ["/w/metawrap_70_10_bins", "/w/metawrap_50_5_bins"]
    assert cmd._ranked_bin_set(sets) is None


def test_a_lone_bin_set_is_the_variants_answer():
    assert cmd._ranked_bin_set(["/w/metabat2_bins"]) == "metabat2_bins"


def test_three_raw_bin_sets_are_not_an_answer_yet():
    """No arbitrary winner: the sets overlap, so none of them is the variant's number."""
    sets = ["/w/metabat2_bins", "/w/maxbin2_bins", "/w/concoct_bins"]
    assert cmd._ranked_bin_set(sets) is None


def test_an_unranked_bin_sets_metrics_are_kept_but_renamed():
    assert cmd._for_bin_set("bins_near_complete", "maxbin2_bins", "metawrap_70_10_bins") == (
        "bins_near_complete__maxbin2_bins"
    )
    assert cmd._for_bin_set("bins_near_complete", "metawrap_70_10_bins", "metawrap_70_10_bins") == (
        "bins_near_complete"
    )


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


def test_the_consolidated_set_carries_the_variants_metric_and_the_raw_ones_do_not(tmp_path):
    root = str(tmp_path / "metabat2_concoct")
    for name in ("metabat2_bins", "concoct_bins"):
        bin_set(root, name, {"bin_001": ["c1"]})
        stats_table(os.path.join(root, name + ".stats.tsv"), [("bin_001", 95, 1)])
    bin_set(root, "metawrap_70_10_bins", {"bin_001": ["c1"], "bin_002": ["c2"]})
    stats_table(
        os.path.join(root, "metawrap_70_10_bins.stats.tsv"),
        [("bin_001", 95, 1), ("bin_002", 91, 3)],
    )

    rows = cmd._collect_scores(
        "cami2_gi",
        "binner_set",
        "metabat2_concoct",
        root,
        ranked_bin_set=cmd._ranked_bin_set(assess.bin_dirs(root)),
    )
    plain = [r for r in rows if r.metric == "bins_near_complete"]
    assert [r.value for r in plain] == [2.0], "the refined set is the variant's number"


def test_checkms_own_report_is_not_counted_alongside_the_normalised_one(tmp_path):
    """CheckM2 writes quality_report.tsv and MetaWrap2 writes .stats.tsv from it. One run, one set
    of rows - reading both would double every count, under two different source names."""
    root = str(tmp_path / "metabat2")
    bin_set(root, "metabat2_bins", {"bin_001": ["c1"]})
    stats_table(os.path.join(root, "metabat2_bins.stats.tsv"), [("bin_001", 95, 1)])
    os.makedirs(os.path.join(root, "metabat2_bins.checkm"))
    with open(os.path.join(root, "metabat2_bins.checkm", "quality_report.tsv"), "w") as fh:
        fh.write("Name\tCompleteness\tContamination\nbin_001\t95\t1\n")

    rows = cmd._collect_scores(
        "cami2_gi", "binner_set", "metabat2", root, ranked_bin_set="metabat2_bins"
    )
    counts = [r for r in rows if r.metric == "bins_near_complete"]
    assert [r.value for r in counts] == [1.0]


# ── AMBER's table ─────────────────────────────────────────────────────────────────────────────


def test_one_amber_run_over_three_bin_sets_is_read_as_three(tmp_path):
    """AMBER takes several queries at once and tags them. Folding them together would report the
    last tool's numbers as though they were the run's."""
    path = tmp_path / "results.tsv"
    path.write_text(
        "Tool\tF1 score (bp)\tARI (bp)\n"
        "metabat2_bins\t0.81\t0.75\n"
        "maxbin2_bins\t0.62\t0.51\n"
        "concoct_bins\t0.70\t0.64\n"
    )
    by_tool = scoring.amber_metrics_by_tool(str(path))
    assert set(by_tool) == {"metabat2_bins", "maxbin2_bins", "concoct_bins"}
    assert by_tool["metabat2_bins"]["amber_f1_bp"] == 0.81
    assert by_tool["concoct_bins"]["amber_ari_bp"] == 0.64


def test_an_unlabelled_amber_table_still_reads(tmp_path):
    path = tmp_path / "results.tsv"
    path.write_text("F1 score (bp)\tARI (bp)\n0.81\t0.75\n")
    by_tool = scoring.amber_metrics_by_tool(str(path))
    assert by_tool == {"": {"amber_f1_bp": 0.81, "amber_ari_bp": 0.75}}


def test_ambers_f1_is_not_summed_across_three_bin_sets(tmp_path):
    """An F1 above 1 is the signature of this bug, and it sorts first when it happens."""
    root = str(tmp_path / "three")
    for name in ("metabat2_bins", "maxbin2_bins", "concoct_bins"):
        bin_set(root, name, {"bin_001": ["c1"]})
    with open(os.path.join(root, "results.tsv"), "w") as fh:
        fh.write(
            "Tool\tF1 score (bp)\n" "metabat2_bins\t0.81\nmaxbin2_bins\t0.62\nconcoct_bins\t0.70\n"
        )

    rows = cmd._collect_scores(
        "cami2_gi",
        "binner_set",
        "three",
        root,
        ranked_bin_set=cmd._ranked_bin_set(assess.bin_dirs(root)),
    )
    assert [r for r in rows if r.metric == "amber_f1_bp"] == []
    assert sum(1 for r in rows if r.metric.startswith("amber_f1_bp__")) == 3


# ── GUNC at score time ────────────────────────────────────────────────────────────────────────


def test_gunc_is_asked_for_separately_from_checkm():
    """Both score bins, and they do not answer the same question - which is why the metric exists.

    CheckM counts single-copy markers, so a bin that is 95% one genome and 5% another reads as barely
    contaminated; AMBER weights by base pairs, so the same bin costs almost nothing in purity.
    """
    assert assess.needs_gunc(["gunc_pass_fraction"])
    assert not assess.needs_gunc(["bins_near_complete", "amber_f1_bp"])
    assert not assess.needs_checkm(["gunc_pass_fraction"])


def test_a_bin_set_that_already_has_a_gunc_report_is_not_scored_again(tmp_path, monkeypatch):
    """GUNC's diamond stage measured 3.2 to 5.7 CPU-hours per bin set on this project's own arms, so
    re-reading a table to add one metric must not cost an afternoon of somebody's machine."""
    root = str(tmp_path)
    scored = bin_set(root, "metabat2_bins", {"bin_001": ["c1"]})
    fresh = bin_set(root, "maxbin2_bins", {"bin_001": ["c2"]})
    os.makedirs(scored + ".gunc")
    with open(os.path.join(scored + ".gunc", "GUNC.progenomes_2.1.maxCSS_level.tsv"), "w") as fh:
        fh.write("genome\tpass.GUNC\nbin_001\tTrue\n")

    ran = []

    def fake_run_gunc(bins_dir, **kw):
        ran.append(bins_dir)
        os.makedirs(kw["out_dir"], exist_ok=True)
        path = os.path.join(kw["out_dir"], "GUNC.progenomes_2.1.maxCSS_level.tsv")
        with open(path, "w") as fh:
            fh.write("genome\tpass.GUNC\nbin_001\tFalse\n")
        return path

    monkeypatch.setattr(assess._gunc, "run_gunc", fake_run_gunc)
    written = assess.run_chimerism(
        [scored, fresh], threads=4, db="db.dmnd", env=None, say=lambda *_: None
    )
    assert ran == [fresh]
    # Both reports come back, so the skipped set is still scored from the report it already had.
    assert len(written) == 2

    # --reassess is the way to redo it when the bins themselves changed.
    ran.clear()
    assess.run_chimerism(
        [scored, fresh], threads=4, db="db.dmnd", env=None, reassess=True, say=lambda *_: None
    )
    assert ran == [scored, fresh]


def test_a_missing_gunc_database_is_reported_rather_than_run_into(tmp_path, monkeypatch):
    """A 13 GB database that was never fetched is the likeliest reason this metric is missing, and
    hearing it from the scorer beats hearing it from a diamond run that dies an hour in."""
    root = str(tmp_path)
    bin_set(root, "metabat2_bins", {"bin_001": ["c1"]})
    monkeypatch.setattr(assess, "gunc_settings", lambda: (None, ""))
    said = []
    _dirs, outcome = assess.assess_bins(
        root,
        metrics=["gunc_pass_fraction"],
        gold_binning=None,
        threads=4,
        memory_gb=8,
        say=said.append,
    )
    assert "no GUNC database" in outcome["gunc"]
    assert any("GUNC skipped" in line for line in said)


def test_a_set_with_its_checkm_table_is_still_outstanding_for_gunc(tmp_path):
    """The skip is per assessor, not per bin set.

    One file used to stand in for "assessed", so adding GUNC as a second per-set assessor would have
    meant every set that had already been through CheckM was skipped before GUNC ever ran - and the
    metric would have come back empty with nothing saying why.
    """
    path = bin_set(str(tmp_path), "metabat2_bins", {"bin_001": ["c1"]})
    with open(path + ".stats.tsv", "w") as fh:
        fh.write("bin\tcompleteness\tcontamination\nbin_001\t95.0\t1.0\n")

    assert cmd._already_assessed(path, ["bins_near_complete"])
    assert not cmd._already_assessed(path, ["bins_near_complete", "gunc_pass_fraction"])

    os.makedirs(path + ".gunc")
    with open(os.path.join(path + ".gunc", "GUNC.progenomes_2.1.maxCSS_level.tsv"), "w") as fh:
        fh.write("genome\tpass.GUNC\nbin_001\tTrue\n")
    assert cmd._already_assessed(path, ["bins_near_complete", "gunc_pass_fraction"])


# --- finding a consolidated set that lives one level down -------------------------------------


def _bin_set(path, n=1):
    os.makedirs(path, exist_ok=True)
    for i in range(n):
        with open(os.path.join(path, "bin_%03d.fasta" % (i + 1)), "w") as fh:
            fh.write(">c%d\nACGTACGTAC\n" % i)
    return path


def test_a_consolidated_set_inside_a_subdirectory_is_found(tmp_path):
    """A binning arm is refined as a second step, which writes into a subdirectory of the variant.

    With a flat scan those arms reported "several raw bin sets and no consolidated one" and went
    unranked while their consolidated set sat one level down - six arms of the binner comparison,
    every ensemble among them.
    """
    root = tmp_path / "metabat2_semibin2"
    _bin_set(str(root / "metabat2_bins"), 2)
    _bin_set(str(root / "semibin2_bins"), 3)
    _bin_set(str(root / "refined" / "metawrap_70_10_bins"), 4)

    found = [os.path.relpath(p, str(root)) for p in assess.bin_dirs(str(root))]
    assert os.path.join("refined", "metawrap_70_10_bins") in found
    assert "metabat2_bins" in found and "semibin2_bins" in found


def test_the_raw_sets_refinement_restores_beside_its_output_are_not_counted_twice(tmp_path):
    """Refinement copies its inputs back under their original names; taking those would double-count."""
    root = tmp_path / "variant"
    _bin_set(str(root / "metabat2_bins"), 2)
    _bin_set(str(root / "refined" / "metabat2_bins"), 2)
    _bin_set(str(root / "refined" / "metawrap_70_10_bins"), 3)

    found = [os.path.relpath(p, str(root)) for p in assess.bin_dirs(str(root))]
    assert found.count("metabat2_bins") == 1
    assert os.path.join("refined", "metabat2_bins") not in found


def test_an_empty_nested_set_is_not_offered_as_an_answer(tmp_path):
    root = tmp_path / "variant"
    _bin_set(str(root / "metabat2_bins"), 1)
    (root / "refined" / "metawrap_70_10_bins").mkdir(parents=True)
    found = [os.path.relpath(p, str(root)) for p in assess.bin_dirs(str(root))]
    assert found == ["metabat2_bins"]


def test_a_file_where_a_subdirectory_might_be_does_not_stop_the_scan(tmp_path):
    root = tmp_path / "variant"
    _bin_set(str(root / "metabat2_bins"), 1)
    (root / "run.stdout").write_text("log output\n")
    assert len(assess.bin_dirs(str(root))) == 1
