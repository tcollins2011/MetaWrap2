"""Parsers for other tools' annotation tables: KOfam ko_list, Pfam accessions, dbCAN, GECCO, BRITE.

In these files a missing value can look like a result. A KOfam profile with no curated threshold
must parse as None, not 0.0 (which accepts every hit); dbCAN's lone ``-`` must parse as empty, not
a substrate named ``-``; GECCO's placeholder row for a genome with no clusters is not a cluster.

The dbCAN and GECCO fixtures are real output (dbCAN 5.2.9 and GECCO over CAMI2 bins; see
``tests/fixtures/cazy/SOURCE.md`` and ``tests/fixtures/bgc/SOURCE.md``). KOfam's ko_list and
KEGG's htext are written inline in their own layouts, as in ``tests/test_kofam.py`` and
``tests/test_hierarchy.py``, since the layout is what is tested.
"""

from __future__ import annotations

import csv
import dataclasses
import gzip
import os

import pytest

from metawrap2 import bgc, cazy, hierarchy, kofam, pfam

CAZY = os.path.join(os.path.dirname(__file__), "fixtures", "cazy")
BGC = os.path.join(os.path.dirname(__file__), "fixtures", "bgc")

#: KOfam's ko_list in its twelve-column layout. The third row is a profile with no cutoff, marked
#: with "-" in the threshold, score_type, profile_type and F-measure columns.
KO_LIST = (
    "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\talen\tmlen\teff_nseq"
    "\tre/pos\tdefinition\n"
    "K00001\t320.40\tdomain\tall\t0.883405\t97\t92\t560\t288\t1.26\t0.589\talcohol dehydrogenase"
    " [EC:1.1.1.1]\n"
    "K00002\t457.17\tfull\tall\t0.447532\t3471\t3318\t7322\t528\t8.57\t0.590\talcohol"
    " dehydrogenase (NADP+) [EC:1.1.1.2]\n"
    "K00492\t-\t-\t-\t-\t1\t1\t396\t396\t0.39\t0.593\t1,3,7-trimethyluric acid 5-monooxygenase\n"
)

#: A minimal BRITE ko00001 htext: one biological branch and one scaffolding branch.
HTEXT = """+D\tKO
!
A09100 Metabolism
B  09101 Carbohydrate metabolism
C    00010 Glycolysis / Gluconeogenesis [PATH:ko00010]
D      K00844  HK; hexokinase [EC:2.7.1.1]
A09180 Brite Hierarchies
B  09182 Protein families: genetic information processing
C    03000 Transcription factors [BR:ko03000]
D      K23479  CCAAT/enhancer binding protein (C/EBP), other
"""


def _header(path):
    with open(path, newline="") as fh:
        return next(csv.reader(fh, delimiter="\t"))


# ── KOfam's ko_list ───────────────────────────────────────────────────────────────────────────


def test_ko_list_columns_are_read_by_name_from_its_own_header(tmp_path):
    """``knum``, ``threshold``, ``score_type``, ``profile_type``, ``F-measure``, ``definition``.

    Columns are read by name: ``F-measure`` is hyphenated, and ``definition`` is the last of twelve
    columns and holds free text. A positional read would break if a column were added upstream.
    """
    path = tmp_path / "ko_list"
    path.write_text(KO_LIST)
    profiles = kofam.load_ko_list(str(path))
    assert set(profiles) == {"K00001", "K00002", "K00492"}
    assert profiles["K00001"].definition == "alcohol dehydrogenase [EC:1.1.1.1]"
    assert profiles["K00001"].profile_type == "all"
    assert profiles["K00001"].f_measure == pytest.approx(0.883405)


def test_a_profile_with_no_curated_cutoff_has_none_rather_than_a_threshold_of_zero(tmp_path):
    """A ``-`` threshold parses as None, not 0.0 (which would accept every hit).

    Around 40% of KOfam's profiles have no curated cutoff; parsed as 0.0, each would accept hits
    at any score that look like normal KO calls.
    """
    path = tmp_path / "ko_list"
    path.write_text(KO_LIST)
    uncurated = kofam.load_ko_list(str(path))["K00492"]
    assert kofam.NO_THRESHOLD == "-"
    assert uncurated.threshold is None
    assert uncurated.has_threshold is False
    assert uncurated.score_type is None
    assert uncurated.f_measure is None
    assert kofam.load_ko_list(str(path))["K00001"].has_threshold is True


# ── Pfam's accessions and its release ─────────────────────────────────────────────────────────


def test_a_pfam_accession_keeps_its_version_and_joins_without_it():
    """``PF00001.23`` is family ``PF00001`` at version 23; both forms are kept.

    The versioned accession is output for reproducibility. The family is used for joins, since
    Pfam bumps versions without changing the family and a version mismatch would drop rows.
    """
    hit = pfam.DomainHit(
        gene="k141_1_1",
        accession="PF00001.23",
        name="7tm_1",
        start=41,
        end=290,
        bit_score=120.5,
        evalue=1e-35,
        curated_cutoff=True,
    )
    assert hit.accession == "PF00001.23"
    assert hit.family == "PF00001"
    assert hit.length == 250
    # Both forms are in the output table.
    row = pfam.domain_row(hit, "bin_003")
    assert row["pfam_accession"] == "PF00001.23"
    assert row["pfam_family"] == "PF00001"
    # An unversioned accession is its own family; a non-Pfam accession is passed through as is.
    assert dataclasses.replace(hit, accession="PF00001").family == "PF00001"
    assert dataclasses.replace(hit, accession="custom_hmm").family == "custom_hmm"


def test_the_pfam_release_is_read_from_pfams_own_relnotes_and_not_from_a_path(tmp_path):
    """The Pfam release is read from ``relnotes.txt``, not the directory name.

    An unknown release is reported as ``"unknown"``, not guessed or left blank, since Pfam
    renumbers families between releases.
    """
    directory = tmp_path / "pfam-37.4-mislabelled-as-36"
    directory.mkdir()
    (directory / "relnotes.txt").write_text(
        "              Pfam release 37.0\n              ===================\n"
    )
    assert pfam.release_from_notes(str(directory)) == "37.0"
    assert pfam.release_from_notes(str(directory / "relnotes.txt")) == "37.0"
    assert pfam.release_from_notes(str(tmp_path / "no-notes-here")) == "unknown"


def test_dbcan_overview_columns_are_dbcans_own_spelling():
    """``Gene ID``, ``EC#``, ``dbCAN_hmm``, ``dbCAN_sub``, ``DIAMOND``, ``#ofTools``.

    Checked against the real header: a ``#`` in one column name, a space in another, and the three
    method columns named after the methods. The header also has ``Recommend Results``.
    """
    header = set(_header(os.path.join(CAZY, "overview.tsv")))
    for key, column in cazy.OVERVIEW_COLUMNS.items():
        assert column in header, "%s -> %r is not a column dbCAN writes" % (key, column)
    assert cazy.OVERVIEW_COLUMNS["tools"] == "#ofTools"
    assert cazy.OVERVIEW_COLUMNS["recommended"] == "Recommend Results"


def test_a_partly_unmapped_substrate_list_drops_only_the_dash():
    """dbCAN's ``-;peptidoglycan`` (two domains, one unmapped) gives one substrate.

    The list is cleaned member by member, since checking only for a lone dash would count a
    substrate named ``-``.
    """
    assert cazy._clean_list("-;peptidoglycan") == "peptidoglycan"
    assert cazy._clean_list("-") == ""
    assert cazy._clean_list("host glycan;host glycan") == "host glycan"

    overview = os.path.join(CAZY, "overview.tsv")
    with open(overview, newline="") as fh:
        raw = {
            row[cazy.OVERVIEW_COLUMNS["substrate"]] for row in csv.DictReader(fh, delimiter="\t")
        }
    assert "-;peptidoglycan" in raw, "the fixture still carries the case this guards"
    substrates = {call.substrate for call in cazy.parse_overview(overview, "g")}
    assert "peptidoglycan" in substrates
    assert "-;peptidoglycan" not in substrates


def test_gecco_cluster_columns_are_geccos_own():
    """``sequence_id``, ``cluster_id``, ``average_p``, ``max_p``, plus one ``*_probability`` per class.

    The per-class probabilities are read by name from the class list, not positionally, since GECCO
    adds product classes between releases.
    """
    header = set(_header(os.path.join(BGC, "clusters.tsv")))
    for key, column in bgc.CLUSTER_COLUMNS.items():
        assert column in header, "%s -> %r is not a column GECCO writes" % (key, column)
    for name in bgc.BGC_TYPES:
        assert "%s_probability" % name in header, "no probability column for %r" % name

    clusters = bgc.parse_clusters(os.path.join(BGC, "clusters.tsv"), "bin_006")
    assert len(clusters) == 4
    first = clusters[0]
    assert set(first.probabilities) == set(bgc.BGC_TYPES)
    assert first.average_p == pytest.approx(0.9763342716696883)


def test_a_cluster_on_an_unmeasured_contig_reports_none_rather_than_a_clean_flank():
    """Without a contig length, cluster completeness is None.

    A zero would mean "at the contig edge", and a large default would mark every cluster complete.
    """
    path = os.path.join(BGC, "clusters.tsv")
    unmeasured = bgc.parse_clusters(path, "bin_006")[0]
    assert unmeasured.contig_length == 0
    assert unmeasured.distance_to_contig_end is None
    assert unmeasured.is_truncated is False

    measured = bgc.parse_clusters(path, "bin_006", {unmeasured.contig: 9234})[0]
    assert measured.distance_to_contig_end == 391
    assert measured.is_truncated is True
    assert bgc.EDGE_WINDOW == 1_000


# ── KEGG's BRITE htext ────────────────────────────────────────────────────────────────────────


def test_the_brite_projection_keeps_its_column_order_and_strips_keggs_tags(tmp_path):
    """One row per (KO, pathway), ten columns in a fixed order, with ``[PATH:ko00010]`` removed.

    The long table, the wide table and the Krona input all rely on this column order. The
    ``[PATH:]`` and ``[EC:]`` tags are KEGG markup and are removed so names are stable across
    releases.
    """
    htext = tmp_path / "ko00001.keg"
    htext.write_text(HTEXT)
    out = tmp_path / hierarchy.HIERARCHY_FILE
    assert hierarchy.project_htext(str(htext), str(out)) == 2

    with gzip.open(str(out), "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        assert tuple(reader.fieldnames or ()) == hierarchy.HIERARCHY_FILE_COLUMNS
        rows = {row["ko"]: row for row in reader}
    assert rows["K00844"]["pathway"] == "Glycolysis / Gluconeogenesis"
    assert rows["K00844"]["pathway_id"] == "00010"
    assert rows["K00844"]["gene_symbol"] == "HK"
    assert rows["K00844"]["ec"] == "2.7.1.1"
    # A D line without the "symbol; name" separator is all name, no symbol.
    assert rows["K23479"]["gene_symbol"] == ""
    assert rows["K23479"]["gene_name"].startswith("CCAAT/enhancer binding protein")


def test_britess_scaffolding_branches_are_flagged_by_id_not_by_their_display_name(tmp_path):
    """Non-pathway branches such as ``09180 Brite Hierarchies`` are identified by id.

    KEGG may reword display names, and this branch is the largest in the file (22,768 of 64,956
    rows). The branches stay in the table; the flag lets a pathway view exclude them.
    """
    htext = tmp_path / "ko00001.keg"
    htext.write_text(HTEXT)
    hierarchy.project_htext(str(htext), str(tmp_path / hierarchy.HIERARCHY_FILE))
    hierarchy.write_snapshot(str(tmp_path), str(htext), hierarchy.BRITE_URL, 2)
    tree = hierarchy.load_hierarchy(str(tmp_path))

    assert hierarchy.NON_PATHWAY_BRANCHES == ("09180", "09190")
    assert tree.level_a_ids["Metabolism"] == "09100"
    assert tree.level_a_ids["Brite Hierarchies"] == "09180"
    assert tree.is_pathway_branch("Metabolism") is True
    assert tree.is_pathway_branch("Brite Hierarchies") is False
    # Both KOs are still in the hierarchy; KOs in non-pathway branches are not dropped.
    assert "K00844" in tree and "K23479" in tree
    # A KO that KEGG does not place gets one Unclassified path.
    assert tree.paths_for("K99999") == ((hierarchy.UNCLASSIFIED,) * 3,)
