"""The functional hierarchy, and the two ways of weighting it.

What has to hold for a percentage over this hierarchy to mean anything:

* **Totals are conserved.** A KO belongs to several pathways - 12,025 of 28,545 do - so the obvious
  implementation counts a gene once per pathway and produces percentages that add to more than 100.
  Every level of the output has to sum to the input total.
* **A sample is not the average of its genomes.** The per-sample distribution is weighted by
  abundance, so a rare genome with unusual metabolism contributes little. Getting this wrong by
  averaging per-genome percentages would overstate rare organisms, which is the opposite of what a
  community profile is for.
* **Nothing is silently dropped.** KOs KEGG does not place, and the branches of BRITE that are
  classification scaffolding rather than biology, stay in the table and are flagged. Dropping them
  makes every percentage beside them an overestimate.
"""

from __future__ import annotations

import pytest

from metawrap2 import hierarchy as H

#: A small BRITE ko00001 in KEGG's own htext layout: two biological branches, the two scaffolding
#: branches, a KO on two pathways, and a KO whose description has no gene symbol.
HTEXT = """+D\tKO
!
A09100 Metabolism
B
B  09101 Carbohydrate metabolism
C    00010 Glycolysis / Gluconeogenesis [PATH:ko00010]
D      K00844  HK; hexokinase [EC:2.7.1.1]
D      K00845  glk; glucokinase [EC:2.7.1.2]
C    00500 Starch and sucrose metabolism [PATH:ko00500]
D      K00844  HK; hexokinase [EC:2.7.1.1]
B  09102 Energy metabolism
C    00190 Oxidative phosphorylation [PATH:ko00190]
D      K00330  nuoA; NADH-quinone oxidoreductase subunit A [EC:7.1.1.2]
A09180 Brite Hierarchies
B  09182 Protein families: genetic information processing
C    03000 Transcription factors [BR:ko03000]
D      K23479  CCAAT/enhancer binding protein (C/EBP), other
A09190 Not Included in Pathway or Brite
B  09191 Unclassified: metabolism
C    99980 Enzymes with EC numbers
D      K01000  xyz; an enzyme [EC:1.1.1.1 1.1.1.2]
"""


@pytest.fixture(scope="module")
def brite(tmp_path_factory):
    """A KEGG_BRITE directory as `install-db kegg` writes it, from the htext above."""
    directory = tmp_path_factory.mktemp("kegg")
    htext = directory / "ko00001.keg"
    htext.write_text(HTEXT)
    rows = H.project_htext(str(htext), str(directory / H.HIERARCHY_FILE))
    H.write_snapshot(str(directory), str(htext), H.BRITE_URL, rows)
    return directory


@pytest.fixture(scope="module")
def tree(brite):
    return H.load_hierarchy(str(brite))


def _fake(paths, level_a_ids=None):
    return H.Hierarchy(
        paths=paths,
        gene_names={ko: ko for ko in paths},
        level_a_ids=level_a_ids or {"Metabolism": "09100", "Brite Hierarchies": "09180"},
    )


# ── the downloaded hierarchy ─────────────────────────────────────────────────────────────────


def test_the_projection_reads_keggs_htext(brite, tree):
    """Each D line becomes one row per pathway it sits under, with its symbol, name and EC."""
    import csv
    import gzip

    with gzip.open(str(brite / H.HIERARCHY_FILE), "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) == 6
    hexokinase = [r for r in rows if r["ko"] == "K00844"]
    assert {r["pathway"] for r in hexokinase} == {
        "Glycolysis / Gluconeogenesis",
        "Starch and sucrose metabolism",
    }, "the [PATH:] tag is stripped, and a KO is listed under every pathway it is on"
    assert hexokinase[0]["gene_symbol"] == "HK" and hexokinase[0]["ec"] == "2.7.1.1"
    no_symbol = next(r for r in rows if r["ko"] == "K23479")
    assert no_symbol["gene_symbol"] == ""
    assert no_symbol["gene_name"] == "CCAAT/enhancer binding protein (C/EBP), other"
    assert next(r for r in rows if r["ko"] == "K01000")["ec"] == "1.1.1.1 1.1.1.2"


def test_every_row_names_the_download_it_came_from(tree):
    """BRITE is revised, so a percentage is only comparable with one over the same hierarchy."""
    assert tree.snapshot.startswith("ko00001 ") and "md5:" in tree.snapshot
    rows = H.hierarchy_rows("genome", "b", H.distribute({"K00844": 1.0}, tree), tree)
    assert {r["hierarchy_snapshot"] for r in rows} == {tree.snapshot}


# ── conserving the total ─────────────────────────────────────────────────────────────────────


def test_a_kos_weight_is_split_across_its_paths_not_counted_on_each(tree):
    """The whole reason percentages here add up."""
    fake = _fake({"K1": (("Metabolism", "B1", "P1"), ("Metabolism", "B2", "P2"))})
    distribution = H.distribute({"K1": 1.0}, fake)
    assert sum(distribution.values()) == pytest.approx(1.0)
    assert distribution[("Metabolism", "B1", "P1")] == pytest.approx(0.5)
    # Every depth of the real hierarchy sums to the input total.
    weights = {ko: 1.0 for ko in tree.paths}
    distribution = H.distribute(weights, tree)
    for depth in (1, 2, 3):
        assert sum(H.collapse(distribution, depth).values()) == pytest.approx(len(weights))


# ── gene counts and abundance weighting ──────────────────────────────────────────────────────


def test_a_sample_profile_weights_each_genome_by_its_abundance():
    """Not the average of per-genome profiles: a rare genome contributes in proportion to its cells.

    Genome A is 10x as abundant as B, so the community's KO profile looks like A's.
    """
    per_genome = {"A": {"K1": 10.0}, "B": {"K2": 10.0}}
    profile = H.weighted_sample_profile(per_genome, {"A": 0.9, "B": 0.09})
    assert profile["K1"] > profile["K2"] * 9


# ── Krona input ──────────────────────────────────────────────────────────────────────────────


# ── weighting by gene abundance rather than by MAG abundance ─────────────────────────────────


# ── weighting by contig depth: the route taken when no bins are passed ───────────────────────


# ── recovering the contig from the gene id ───────────────────────────────────────────────────


def test_the_contig_is_split_off_at_the_last_underscore_not_the_first():
    """MEGAHIT names contigs with five underscores in them, so the first underscore is `k141`.

    Splitting there would put every gene in the assembly on one contig, and that contig's depth
    would then weight the entire sample - a single number where thousands were measured.
    """
    gene = "k141_73537_length_916682_cov_50.9988_12"
    assert H.contig_from_gene(gene) == "k141_73537_length_916682_cov_50.9988"
    assert H.contig_from_gene("contig_24_7") == "contig_24"
