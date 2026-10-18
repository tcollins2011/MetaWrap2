"""AMRFinderPlus, RGI and argNorm: three tools, three spellings of the same resistance call.

Complements ``tests/test_amr.py``, which pins the context ladder. These pin the formats the three
tools write, which is where this module's silent failures live:

* AMRFinderPlus renamed every column between 3.x and 4.x (``Gene symbol`` to ``Element symbol``,
  ``Contig id`` keeping its lowercase ``id``). Read by the old names, every call parses with an
  empty symbol and an empty drug class - a resistome of unnamed genes.
* RGI writes the ARO as a bare number where CARD, hAMRonization and argNorm all write ``ARO:``-
  prefixed. Joined raw, the two backends never agree on a single call.
* argNorm's first line is a ``#`` comment, so a plain DictReader takes it as the header and returns
  every column under the key ``# argNorm version: 1.1.0`` - an empty result, not an error.

The RGI and argNorm fixtures are vendored real output; see ``tests/fixtures/amr/rgi/SOURCE.md``.
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
    """A report in AMRFinderPlus' own column order, built from the names the module reads by."""
    header = list(amr.AMRFINDER_COLUMNS.values())
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.writelines(
            "\t".join(str(row.get(column, "")) for column in header) + "\n" for row in rows
        )
    return str(path)


# ── AMRFinderPlus ─────────────────────────────────────────────────────────────────────────────


def test_amrfinderplus_column_names_carry_their_spaces_and_percent_signs(tmp_path):
    """``Contig id``, ``Element symbol``, ``% Identity to reference`` - as AMRFinderPlus 4.x writes them.

    The ``%`` and the spaces are the point: these are not identifiers anyone would guess, and the
    4.x rename of ``Gene symbol`` to ``Element symbol`` is exactly the kind of upstream change that
    parses as an absence. An absent symbol is a resistance call with no gene name on it.
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
    """The contig column is checked by name, because every context claim is built on it.

    Missed, the contig is empty, the contig length lookup misses, and every call reports distance 0
    to the contig end - which the context ladder reads as "too close to the edge to say". A whole
    resistome of unplaceable genes, with nothing logged.
    """
    path = tmp_path / "report.tsv"
    header = [c.replace("Contig id", "Contig") for c in amr.AMRFINDER_COLUMNS.values()]
    path.write_text("\t".join(header) + "\n")
    with pytest.raises(KeyError) as caught:
        amr.parse_amrfinder(str(path), "bin_007", {}, "2026-08-07.1")
    assert "AMRFinderPlus" in str(caught.value)


def test_rgi_column_names_are_the_ones_rgi_writes():
    """Checked against a real RGI table: ``Best_Hit_ARO`` underscored, ``Drug Class`` spaced.

    RGI mixes both conventions in one header, so neither can be derived from the other. The pin is
    on the real file's header rather than on CARD's documentation, which spells some of them
    differently again.
    """
    with open(RGI, newline="") as fh:
        header = set(next(csv.reader(fh, delimiter="\t")))
    for key, column in amr.RGI_COLUMNS.items():
        assert column in header, "%s -> %r is not a column RGI writes" % (key, column)
    assert amr.RGI_COLUMNS["aro_name"] == "Best_Hit_ARO"
    assert amr.RGI_COLUMNS["drug_class"] == "Drug Class"
    assert amr.RGI_COLUMNS["coverage"] == "Percentage Length of Reference Sequence"


def test_rgi_writes_a_bare_aro_number_and_the_join_needs_the_prefix():
    """RGI's ``ARO`` column holds ``3002647``; everything else in the ecosystem writes ``ARO:3002647``.

    Joined raw against argNorm's output, no call from either backend ever matches one from the
    other - which reads as two tools finding entirely different genes rather than as a key format
    mismatch. Pinned in both directions so an already-prefixed value is not doubly prefixed.
    """
    with open(RGI, newline="") as fh:
        raw = {row["ARO"] for row in csv.DictReader(fh, delimiter="\t")}
    assert all(not value.startswith("ARO:") for value in raw)

    calls = amr.parse_rgi(RGI, "bin_007", {CONTIG: 92857}, database_version="4.0.1")
    assert all(call.aro.startswith("ARO:") for call in calls)
    assert "ARO:3002647" in {call.aro for call in calls}
    assert amr._aro_accession("ARO:3002647") == "ARO:3002647"
    assert amr._aro_accession("") == ""

    # And the prefixed accessions are what make a call found by both backends one finding.
    normalised = {row.aro for row in amr.parse_argnorm(ARGNORM)}
    assert "ARO:3002647" in normalised & {call.aro for call in calls}


def test_a_gene_argnorm_could_not_place_comes_back_blank_not_as_a_placeholder(tmp_path):
    """Absence from the ontology is not evidence a call is wrong - a new gene may just not be in it.

    An empty ARO is therefore a legitimate value, and must not become the string ``"nan"``, a zero,
    or a row dropped from the table: the gene was still found, it just has no ontology term.
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
    # The comment line that precedes the header is still skipped rather than read as column names.
    assert "cepA" in rows
