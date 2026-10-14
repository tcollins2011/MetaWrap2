"""Resistance calls, their context, and the three claims that must stay distinct.

The arithmetic here is trivial; the distinctions are not, and every one of them corresponds to a way
published resistome results get over-read:

* **"looked and found nothing" is not "did not look."** With no mobile-element caller run, the
  honest verdict is ``context_unknown``, and reporting ``no_mge_in_window`` instead turns an unasked
  question into a negative finding.
* **"nothing within 10 kb" needs 10 kb of sequence.** A gene 2 kb from the end of its contig cannot
  support that claim however long the window is.
* **detection is not phenotype.** Every row is L1 or L2 of Cao et al.'s ladder, and the level is
  derived from the context rather than from how confident the detection was.

The fixture is a real finding from CAMI2 sample 10: ant(6)-Ia, aph(3')-IIIa and sat4 within 2.3 kb
on one 92,857 bp contig - a known aminoglycoside/streptothricin cassette.
"""

from __future__ import annotations

import os

import pytest

from metawrap2 import amr

CONTIG = "k141_57263_length_92857_cov_80.8395"
CONTIG_LENGTH = 92857


def _call(
    symbol="aph(3')-IIIa",
    start=14927,
    end=15718,
    element_type="AMR",
    method="EXACTX",
    drug_class="AMINOGLYCOSIDE",
    contig=CONTIG,
    contig_length=CONTIG_LENGTH,
):
    distance = min(start - 1, contig_length - end) if contig_length else 0
    return amr.ArgCall(
        genome="bin_007",
        contig=contig,
        start=start,
        end=end,
        strand="-",
        symbol=symbol,
        name=symbol,
        element_type=element_type,
        subtype=element_type,
        drug_class=drug_class,
        subclass=drug_class,
        method=method,
        identity=100.0,
        coverage=100.0,
        contig_length=contig_length,
        distance_to_contig_end=max(distance, 0),
        tool="amrfinderplus",
        database_version="2026-08-07.1",
    )


# ── parsing ──────────────────────────────────────────────────────────────────────────────────


def _write_report(path, rows):
    header = list(amr.AMRFINDER_COLUMNS.values())
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.writelines(
            "\t".join(str(row.get(column, "")) for column in header) + "\n" for row in rows
        )
    return str(path)


def test_contig_length_comes_from_the_assembly_not_the_report(tmp_path):
    """AMRFinderPlus reports coordinates and not contig lengths, so a wrapper passing its output
    through cannot tell a gene mid-megabase from one at the edge of what was assembled."""
    path = _write_report(
        tmp_path / "r.tsv",
        [
            {
                "Contig id": CONTIG,
                "Start": 14927,
                "Stop": 15718,
                "Strand": "-",
                "Element symbol": "aph(3')-IIIa",
                "Type": "AMR",
                "Method": "EXACTX",
                "% Identity to reference": "100.00",
            }
        ],
    )
    calls = amr.parse_amrfinder(path, "bin_007", {CONTIG: CONTIG_LENGTH}, "2026-08-07.1")
    assert len(calls) == 1
    assert calls[0].contig_length == CONTIG_LENGTH
    assert calls[0].distance_to_contig_end == 14926

    # A contig missing from the assembly yields zero rather than a wrong distance: zero is the
    # conservative answer, it forces context_unknown rather than inventing room.
    path = _write_report(
        tmp_path / "absent.tsv",
        [{"Contig id": "absent", "Start": 500, "Stop": 100, "Type": "AMR", "Method": "EXACTX"}],
    )
    call = amr.parse_amrfinder(path, "b", {}, "v")[0]
    assert (call.start, call.end) == (100, 500), "reversed coordinates are normalised"
    assert (call.contig_length, call.distance_to_contig_end) == (0, 0)
    assert amr.classify_context(call, [])[0] is amr.Context.UNKNOWN


# ── context ──────────────────────────────────────────────────────────────────────────────────


def test_a_gene_too_close_to_a_contig_end_cannot_support_an_absence_claim():
    """The measured reason this matters: MEGAHIT's ARG-bearing contigs average 284 bp.

    And no caller run is unknown, not absence: None means the question was never asked; an empty
    list means it was asked and answered.
    """
    assert amr.classify_context(_call(), None)[0] is amr.Context.UNKNOWN
    assert amr.classify_context(_call(), [])[0] is amr.Context.NO_MGE_IN_WINDOW
    near_edge = _call(start=200, end=700, contig_length=1000)
    assert near_edge.distance_to_contig_end == 199
    assert amr.classify_context(near_edge, [], window=10_000)[0] is amr.Context.UNKNOWN
    # The same gene with room either side can support it.
    roomy = _call(start=20_000, end=20_500, contig_length=92_857)
    assert amr.classify_context(roomy, [], window=10_000)[0] is amr.Context.NO_MGE_IN_WINDOW


# ── the evidence ladder ──────────────────────────────────────────────────────────────────────


def test_the_level_follows_the_context_not_the_confidence_of_detection():
    """A 100%-identity exact match with an unobserved neighbourhood is still only L1.

    L2 is "linked to a host and a genetic neighbourhood". How certain the detection was does not
    bear on whether the neighbourhood was seen.
    """
    assert amr.evidence_level(amr.Context.UNKNOWN) is amr.EvidenceLevel.PRESENCE
    assert amr.evidence_level(amr.Context.MGE_NEARBY) is amr.EvidenceLevel.CONTEXT
    assert amr.evidence_level(amr.Context.NO_MGE_IN_WINDOW) is amr.EvidenceLevel.CONTEXT
    # No level above two can be produced: L3 needs expression and L4 a phenotype assay.
    assert {level.value for level in amr.EvidenceLevel} == {"L1_presence", "L2_context"}
    assert set(amr.LEVEL_LANGUAGE) == {"L1_presence", "L2_context"}
    for phrasing in amr.LEVEL_LANGUAGE.values():
        assert "confers" not in phrasing


# ── clustering and summary ───────────────────────────────────────────────────────────────────


# ── hAMRonization ────────────────────────────────────────────────────────────────────────────


# ── argNorm: one name per gene, and the drugs it implicates ──────────────────────────────────

ARGNORM_OUTPUT = os.path.join(os.path.dirname(__file__), "fixtures", "amr", "argnorm_output.tsv")


def test_argnorm_output_is_parsed_past_its_comment_line():
    """argNorm writes `# argNorm version: ...` first, which csv would take as the header.

    Handed the file directly, DictReader keys every column under that comment and the result is
    silently empty rather than an error - so the comment handling is the thing being tested.
    """
    rows = amr.parse_argnorm(ARGNORM_OUTPUT)
    assert len(rows) == 4
    assert [row.gene for row in rows] == ["cepA", "ant(6)-Ia", "aph(3')-IIIa", "sat4"]
    assert [row.genome for row in rows] == ["bin_003", "bin_007", "bin_007", "bin_007"]


# ── mobile-element context, from geNomad ─────────────────────────────────────────────────────

GENOMAD_FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "amr", "genomad")


def test_proviruses_are_read_with_their_host_coordinates():
    """geNomad's provirus table is its only output with coordinates on the host contig.

    The `seq_name` column names the trimmed viral piece; `source_seq` names the contig an ARG call
    would carry, and using the wrong one gives elements no call can ever be near.
    """
    elements = amr.read_proviruses(os.path.join(GENOMAD_FIXTURES, "assembly_provirus.tsv"))
    assert len(elements) == 3
    first = elements[0]
    assert first.contig == "k141_71203_length_191108_cov_24.3493"
    assert (first.start, first.end) == (16390, 39498)
    assert first.source == "genomad"


# ── RGI, the second backend ──────────────────────────────────────────────────────────────────

RGI_REPORT = os.path.join(os.path.dirname(__file__), "fixtures", "amr", "rgi", "bin_007.txt")
RGI_LENGTHS = {
    "k141_18354_length_118866_cov_78.9813": 118866,
    "k141_57263_length_92857_cov_80.8395": 92857,
    "k141_11256_length_72581_cov_80.2039": 72581,
    "k141_13458_length_32395_cov_77.6371": 32395,
}


def test_rgi_output_parses_into_the_same_calls_as_amrfinder():
    """One dataclass for both backends: context, evidence and clustering are properties of a call."""
    calls = amr.parse_rgi(RGI_REPORT, "bin_007", RGI_LENGTHS, database_version="4.0.2")
    assert len(calls) == 5
    assert {c.tool for c in calls} == {"rgi"}
    assert {c.database_version for c in calls} == {"4.0.2"}
    assert all(c.element_type == "AMR" for c in calls)
    perfect = [c for c in calls if c.method == "Perfect"]
    assert len(perfect) == 1
    assert perfect[0].symbol == "APH(3')-IIIa"
    assert perfect[0].identity == pytest.approx(100.0)
    # RGI reports Start > Stop on the minus strand, and every consumer here assumes otherwise.
    for call in calls:
        assert call.contig_length == RGI_LENGTHS[call.contig]
        assert call.start <= call.end
        assert 0 <= call.distance_to_contig_end <= call.contig_length
    # RGI writes a bare number; hAMRonization, argNorm and CARD all write ARO:<number>.
    assert perfect[0].aro == "ARO:3002647"


# ── backend agreement ────────────────────────────────────────────────────────────────────────
