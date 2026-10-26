"""Tests for the AMRFinderPlus, RGI and argNorm output parsers.

``tests/test_amr.py`` covers the context ladder; these cover the tools' file formats, where a
mismatch gives empty results without an error:

* AMRFinderPlus renamed columns between 3.x and 4.x (``Gene symbol`` to ``Element symbol``;
  ``Contig id`` keeps its lowercase ``id``). Read by the old names, every call has an empty symbol
  and drug class.
* RGI writes the ARO as a bare number; CARD, hAMRonization and argNorm write it ``ARO:``-prefixed.
  Without normalising, calls from the two backends never match.
* argNorm's first line is a ``#`` comment, which a plain DictReader takes as the header.

The RGI and argNorm fixtures are real output; see ``tests/fixtures/amr/rgi/SOURCE.md``.
"""

from __future__ import annotations

import csv
import os

import pytest

from metawrap2 import amr

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "amr")
RGI = os.path.join(FIXTURES, "rgi", "bin_007.txt")
ARGNORM = os.path.join(FIXTURES, "argnorm_output.tsv")

CONTIG = "k141_57263_length_92857_cov_80.8395"


def _amrfinder_report(path, rows):
    """Write a report in AMRFinderPlus' column order, using the module's column names."""
    header = list(amr.AMRFINDER_COLUMNS.values())
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.writelines(
            "\t".join(str(row.get(column, "")) for column in header) + "\n" for row in rows
        )
    return str(path)


# ── AMRFinderPlus ─────────────────────────────────────────────────────────────────────────────


def test_amrfinderplus_column_names_carry_their_spaces_and_percent_signs(tmp_path):
    """Column names match AMRFinderPlus 4.x, including spaces and ``%`` signs.

    A renamed column, such as 4.x's ``Gene symbol`` to ``Element symbol``, would parse as an empty
    value.
    """
    assert amr.AMRFINDER_COLUMNS["contig"] == "Contig id"
    assert amr.AMRFINDER_COLUMNS["symbol"] == "Element symbol"
    assert amr.AMRFINDER_COLUMNS["identity"] == "% Identity to reference"
    assert amr.AMRFINDER_COLUMNS["coverage"] == "% Coverage of reference"
    assert amr.AMRFINDER_COLUMNS["stop"] == "Stop"

    path = _amrfinder_report(
        tmp_path / "report.tsv",
        [
            {
                "Contig id": CONTIG,
                "Start": "14927",
                "Stop": "15718",
                "Strand": "-",
                "Element symbol": "aph(3')-IIIa",
                "Element name": "aminoglycoside O-phosphotransferase APH(3')-IIIa",
                "Type": "AMR",
                "Subtype": "AMR",
                "Class": "AMINOGLYCOSIDE",
                "Subclass": "KANAMYCIN",
                "Method": "EXACTX",
                "% Identity to reference": "100.00",
                "% Coverage of reference": "100.00",
            }
        ],
    )
    call = amr.parse_amrfinder(path, "bin_007", {CONTIG: 92857}, "2026-08-07.1")[0]
    assert call.symbol == "aph(3')-IIIa"
    assert call.drug_class == "AMINOGLYCOSIDE"
    assert call.identity == pytest.approx(100.0)
    assert call.method == "EXACTX"


def test_a_renamed_amrfinder_contig_column_is_an_error_not_a_resistome_of_nothing(tmp_path):
    """A missing contig column raises KeyError naming AMRFinderPlus.

    Without the contig, the length lookup fails and every call gets distance 0 to the contig end,
    which the context ladder reads as too close to the edge to place.
    """
    path = tmp_path / "report.tsv"
    header = [c.replace("Contig id", "Contig") for c in amr.AMRFINDER_COLUMNS.values()]
    path.write_text("\t".join(header) + "\n")
    with pytest.raises(KeyError) as caught:
        amr.parse_amrfinder(str(path), "bin_007", {}, "2026-08-07.1")
    assert "AMRFinderPlus" in str(caught.value)


def test_rgi_column_names_are_the_ones_rgi_writes():
    """RGI column names match a real RGI table header.

    RGI mixes underscores (``Best_Hit_ARO``) and spaces (``Drug Class``) in one header, and CARD's
    documentation spells some names differently, so the check uses a real file.
    """
    with open(RGI, newline="") as fh:
        header = set(next(csv.reader(fh, delimiter="\t")))
    for key, column in amr.RGI_COLUMNS.items():
        assert column in header, "%s -> %r is not a column RGI writes" % (key, column)
    assert amr.RGI_COLUMNS["aro_name"] == "Best_Hit_ARO"
    assert amr.RGI_COLUMNS["drug_class"] == "Drug Class"
    assert amr.RGI_COLUMNS["coverage"] == "Percentage Length of Reference Sequence"


def test_rgi_writes_a_bare_aro_number_and_the_join_needs_the_prefix():
    """RGI's bare ARO numbers (``3002647``) are prefixed to ``ARO:3002647`` for the join.

    Without the prefix, no RGI call matches argNorm output. An already-prefixed value is left as is.
    """
    with open(RGI, newline="") as fh:
        raw = {row["ARO"] for row in csv.DictReader(fh, delimiter="\t")}
    assert all(not value.startswith("ARO:") for value in raw)

    calls = amr.parse_rgi(RGI, "bin_007", {CONTIG: 92857}, database_version="4.0.1")
    assert all(call.aro.startswith("ARO:") for call in calls)
    assert "ARO:3002647" in {call.aro for call in calls}
    assert amr._aro_accession("ARO:3002647") == "ARO:3002647"
    assert amr._aro_accession("") == ""

    # With the prefix, a call found by both backends joins on the same accession.
    normalised = {row.aro for row in amr.parse_argnorm(ARGNORM)}
    assert "ARO:3002647" in normalised & {call.aro for call in calls}


def test_a_gene_argnorm_could_not_place_comes_back_blank_not_as_a_placeholder(tmp_path):
    """A gene argNorm could not map keeps its row with an empty ARO, name and drug list.

    A new gene may not be in the ontology yet, so the row is kept, without ``"nan"`` or zero.
    """
    with open(ARGNORM) as src:
        text = src.read()
    lines = text.splitlines()
    header = next(line for line in lines if not line.startswith("#"))
    columns = header.split("\t")
    blank = ["" for _ in columns]
    blank[columns.index("input_file_name")] = "bin_009"
    blank[columns.index("gene_symbol")] = "newly_described_arg"
    path = tmp_path / "argnorm_output.tsv"
    path.write_text(text + "\t".join(blank) + "\n")

    rows = {row.gene: row for row in amr.parse_argnorm(str(path))}
    assert "newly_described_arg" in rows, "an unplaced gene is kept, not dropped"
    assert rows["newly_described_arg"].aro == ""
    assert rows["newly_described_arg"].aro_name == ""
    assert rows["newly_described_arg"].drugs == ()
    # The comment line before the header is skipped.
    assert "cepA" in rows
