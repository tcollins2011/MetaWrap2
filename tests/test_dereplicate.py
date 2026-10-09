"""Tests for the dereplicate module's clustering, which is the part that decides the catalogue.

skani is a black box here on purpose: everything that turns its table into a set of representative
genomes is a pure function, so what is asserted is the thing that can be got wrong - which genomes
end up collapsed, and which copy of a duplicated genome is the one kept.
"""

import os

import pytest

from metawrap2.modules import dereplicate


def _edge(ref, query, ani, af_ref=100.0, af_query=100.0):
    return dereplicate.Edge(ref, query, ani, af_ref, af_query)


def _bins(tmp_path, name, bins):
    """A bin folder, each entry (bin name -> sequence length in bases)."""
    directory = tmp_path / name
    directory.mkdir(parents=True)
    for bin_name, length in bins.items():
        (directory / (bin_name + ".fasta")).write_text(">c1\n" + "ACGT" * (length // 4) + "\n")
    return str(directory)


def _stats(path, rows):
    """A CheckM table of the shape this project writes, beside a bin folder."""
    with open(path + ".stats.tsv", "w") as fh:
        fh.write("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n")
        fh.writelines(
            "%s\t%s\t%s\t0.5\tBacteria\t1000\t1000\n" % (name, completeness, contamination)
            for name, (completeness, contamination) in rows.items()
        )


# --- clustering ---------------------------------------------------------------------------


def test_a_pair_above_the_threshold_becomes_one_cluster():
    clusters = dereplicate.cluster_genomes(["a", "b"], [_edge("a", "b", 99.0)], ani=95.0)
    assert clusters == [["a", "b"]]


def test_a_pair_below_the_threshold_stays_apart():
    clusters = dereplicate.cluster_genomes(["a", "b"], [_edge("a", "b", 94.9)], ani=95.0)
    assert clusters == [["a"], ["b"]]


def test_a_genome_nothing_matched_is_still_in_the_output():
    """Every input has to appear somewhere, or the catalogue silently loses bins."""
    clusters = dereplicate.cluster_genomes(["a", "b", "c"], [_edge("a", "b", 99.0)], ani=95.0)
    assert clusters == [["a", "b"], ["c"]]


def test_identity_without_enough_alignment_is_not_believed():
    """99% ANI over 5% of a genome is a shared prophage, not the same organism."""
    edges = [_edge("a", "b", 99.0, af_ref=5.0, af_query=4.0)]
    assert dereplicate.cluster_genomes(["a", "b"], edges, min_aligned_fraction=10.0) == [
        ["a"],
        ["b"],
    ]


def test_the_larger_aligned_fraction_decides():
    """A half-complete bin aligns fully to its complete neighbour; the neighbour does not to it.

    Taking the smaller of the two fractions would refuse to collapse exactly the duplicates this
    module exists for, so the rule is the larger.
    """
    edges = [_edge("partial", "complete", 99.0, af_ref=98.0, af_query=51.0)]
    assert dereplicate.cluster_genomes(
        ["partial", "complete"], edges, min_aligned_fraction=90.0
    ) == [["partial", "complete"]]


def test_clusters_are_ordered_by_input_and_not_by_edge_order():
    """Two runs over the same bins must give the same catalogue, whatever order skani emitted."""
    names = ["a", "b", "c", "d"]
    forwards = dereplicate.cluster_genomes(names, [_edge("a", "c", 99.0), _edge("b", "d", 99.0)])
    backwards = dereplicate.cluster_genomes(names, [_edge("d", "b", 99.0), _edge("c", "a", 99.0)])
    assert forwards == backwards == [["a", "c"], ["b", "d"]]


def test_single_linkage_chains_which_is_the_documented_cost():
    """Pinned as a property, not as an accident: a-b and b-c merge all three even without a-c.

    The docstring on cluster_genomes says this happens and why it is accepted, so a future change to
    complete linkage should have to change this test deliberately rather than discover it.
    """
    edges = [_edge("a", "b", 96.0), _edge("b", "c", 96.0)]
    assert dereplicate.cluster_genomes(["a", "b", "c"], edges, ani=95.0) == [["a", "b", "c"]]


def test_edges_naming_a_genome_that_was_not_asked_about_are_ignored():
    edges = [_edge("a", "elsewhere", 99.0)]
    assert dereplicate.cluster_genomes(["a", "b"], edges) == [["a"], ["b"]]


# --- parsing skani ------------------------------------------------------------------------


def test_skani_columns_are_read_by_name(tmp_path):
    """--detailed and --ci insert columns, so a positional read would take the wrong one."""
    table = tmp_path / "ani.tsv"
    table.write_text(
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        "a.fa\tb.fa\t99.28\t100.00\t97.50\tc1\tc1\n"
    )
    (edge,) = dereplicate.parse_skani(str(table))
    assert (edge.ref, edge.query, edge.ani, edge.af_query) == ("a.fa", "b.fa", 99.28, 97.5)


def test_a_table_that_is_not_skani_output_is_rejected(tmp_path):
    table = tmp_path / "ani.tsv"
    table.write_text("genome1\tgenome2\tidentity\n a\tb\t99\n")
    with pytest.raises(ValueError, match="not a skani table"):
        dereplicate.parse_skani(str(table))


def test_an_empty_table_is_no_edges_rather_than_an_error(tmp_path):
    """skani writes the header and nothing else when every pair is screened out."""
    table = tmp_path / "ani.tsv"
    table.write_text(
        "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
    )
    assert dereplicate.parse_skani(str(table)) == []


# --- gathering bins -----------------------------------------------------------------------


def test_bins_are_renamed_by_their_folder_so_two_samples_can_both_have_bin_001(tmp_path):
    a = _bins(tmp_path, "sampleA_bins", {"bin_001": 400})
    b = _bins(tmp_path, "sampleB_bins", {"bin_001": 400})
    assert [g.name for g in dereplicate.collect_genomes([a, b])] == [
        "sampleA_bins__bin_001",
        "sampleB_bins__bin_001",
    ]


def test_two_folders_with_the_same_name_are_disambiguated(tmp_path):
    """The cohort layout is one metawrap_70_10_bins per sample, all called the same thing."""
    first = _bins(tmp_path / "s1", "metawrap_70_10_bins", {"bin_001": 400})
    second = _bins(tmp_path / "s2", "metawrap_70_10_bins", {"bin_001": 400})
    names = [g.name for g in dereplicate.collect_genomes([first, second])]
    assert names == ["metawrap_70_10_bins__bin_001", "metawrap_70_10_bins__bin_001_2"]


# --- choosing the representative ----------------------------------------------------------


def test_the_best_scoring_genome_represents_its_cluster():
    scores = {"a": 70.0, "b": 90.0}
    sizes = {"a": 5_000_000, "b": 3_000_000}
    assert dereplicate.choose_representative(["a", "b"], scores, sizes) == "b"


def test_a_scored_genome_beats_an_unscored_one_however_large():
    """Size says a bin is big, not that it is one organism: a 6 Mb chimera must not win."""
    scores = {"small": 85.0}
    sizes = {"small": 3_000_000, "chimera": 6_000_000}
    assert dereplicate.choose_representative(["chimera", "small"], scores, sizes) == "small"


def test_with_nothing_scored_the_largest_genome_wins():
    sizes = {"a": 1_000, "b": 2_000}
    assert dereplicate.choose_representative(["a", "b"], {}, sizes) == "b"


def test_the_choice_is_stable_when_everything_ties():
    sizes = {"a": 1_000, "b": 1_000}
    assert dereplicate.choose_representative(["a", "b"], {}, sizes) == "a"


# --- reading the CheckM tables this project already writes ---------------------------------


def test_scores_come_from_the_checkm_table_beside_the_bins(tmp_path):
    folder = _bins(tmp_path, "sampleA_bins", {"bin_001": 400, "bin_002": 400})
    _stats(folder, {"bin_001": (95.0, 1.0), "bin_002": (60.0, 10.0)})
    genomes = dereplicate.collect_genomes([folder])
    scores = dereplicate.read_quality([folder], genomes)
    # completeness - 5 * contamination, plus the size tiebreak, so 95-5 beats 60-50 by a long way.
    assert scores["sampleA_bins__bin_001"] > scores["sampleA_bins__bin_002"]


def test_a_bin_with_no_row_is_absent_rather_than_scored_zero(tmp_path):
    """ "Never scored" and "scored badly" pick different representatives, so they stay apart."""
    folder = _bins(tmp_path, "sampleA_bins", {"bin_001": 400, "bin_002": 400})
    _stats(folder, {"bin_001": (95.0, 1.0)})
    genomes = dereplicate.collect_genomes([folder])
    scores = dereplicate.read_quality([folder], genomes)
    assert "sampleA_bins__bin_001" in scores
    assert "sampleA_bins__bin_002" not in scores


def test_a_folder_with_no_checkm_table_simply_has_no_scores(tmp_path):
    folder = _bins(tmp_path, "sampleA_bins", {"bin_001": 400})
    genomes = dereplicate.collect_genomes([folder])
    assert dereplicate.read_quality([folder], genomes) == {}


def test_scores_are_not_taken_from_another_folders_table(tmp_path):
    """Both folders have a bin_001; matching on bin name alone would score the wrong genome."""
    a = _bins(tmp_path, "sampleA_bins", {"bin_001": 400})
    b = _bins(tmp_path, "sampleB_bins", {"bin_001": 400})
    _stats(a, {"bin_001": (95.0, 1.0)})
    genomes = dereplicate.collect_genomes([a, b])
    scores = dereplicate.read_quality([a, b], genomes)
    assert list(scores) == ["sampleA_bins__bin_001"]


# --- the ANI recorded against each representative -----------------------------------------


def test_each_member_records_its_identity_to_the_representative():
    edges = [_edge("a", "rep", 99.1), _edge("rep", "b", 97.3)]
    identities = dereplicate.ani_to_representative(["a", "b", "rep"], "rep", edges)
    assert identities == {"a": 99.1, "b": 97.3, "rep": 100.0}


def test_a_member_pulled_in_by_chaining_has_no_identity_to_the_representative():
    """Which is the point of recording it: this is what an over-merged cluster looks like."""
    edges = [_edge("rep", "middle", 96.0), _edge("middle", "far", 96.0)]
    identities = dereplicate.ani_to_representative(["rep", "middle", "far"], "rep", edges)
    assert identities["middle"] == 96.0
    assert "far" not in identities


# --- the whole module, with skani stubbed out ---------------------------------------------


def test_the_catalogue_keeps_one_genome_per_cluster(tmp_path, monkeypatch):
    """End to end over a fake skani table: two samples sharing one genome give three genomes."""
    a = _bins(tmp_path, "sampleA_bins", {"bin_001": 4000, "bin_002": 400})
    b = _bins(tmp_path, "sampleB_bins", {"bin_001": 4000})
    _stats(b, {"bin_001": (95.0, 1.0)})
    out = str(tmp_path / "out")
    # No conda env, so the test does not need mamba to be on PATH; skani is stubbed anyway.
    config = tmp_path / "mw2.toml"
    config.write_text("[settings]\nuse_conda_envs = false\n")

    def fake_run(cmd, **kwargs):
        # Write the table skani would have written for these three genomes.
        target = cmd.split("-o ")[1].split()[0]
        with open(target, "w") as fh:
            fh.write(
                "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\t"
                "Ref_name\tQuery_name\n"
            )
            fh.write(
                "%s\t%s\t99.5\t100.0\t100.0\tc1\tc1\n"
                % (os.path.join(a, "bin_001.fasta"), os.path.join(b, "bin_001.fasta"))
            )
        return 0

    monkeypatch.setattr(dereplicate, "run", fake_run)
    argv = ["-b", a, "-b", b, "-o", out, "-t", "2", "--config", str(config)]
    assert dereplicate.main(argv) == 0

    kept = sorted(os.listdir(os.path.join(out, dereplicate.CATALOGUE_DIR)))
    # The scored copy of the shared genome, not sample A's unscored one.
    assert kept == ["sampleA_bins__bin_002.fasta", "sampleB_bins__bin_001.fasta"]

    with open(os.path.join(out, dereplicate.CLUSTERS_TABLE)) as fh:
        rows = fh.read().splitlines()
    assert rows[0].split("\t") == list(dereplicate.CLUSTER_COLUMNS)
    # Every input bin is accounted for, not just the ones that were kept.
    assert len(rows) == 4
