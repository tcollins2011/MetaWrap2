"""The bin-quality tables other tools write: GUNC's maxCSS report, and CheckM1's and CheckM2's.

These are regressions over *other people's file formats*, which is a narrower thing than a unit test
over our own arithmetic. Each upstream release can rename a column, change a unit, or retype a
boolean, and none of those arrive as an exception - they arrive as a plausible number:

* GUNC writes ``pass.GUNC`` as ``True``/``False`` when the column stays boolean and as ``1.0``/
  ``0.0`` when a ``nan`` in it forces a float dtype. A reader that accepted only the words scored a
  **166-bin set at 0% non-chimeric**, which reads as every bin being a chimera.
* CheckM1 spells the bin column ``Bin Id`` and CheckM2 spells it ``Name``; a reader that knew one
  scored nothing at all from the other, with no error.

The tables here are hand-written in each tool's own layout rather than vendored, which is the same
choice ``tests/test_gunc.py`` and ``tests/test_benchmark_scoring.py`` already make: the column names
are the thing under test, so they are spelled out in the test rather than hidden in a fixture.
"""

from __future__ import annotations

import pytest

from metawrap2 import checkm, gunc
from metawrap2.benchmark import scoring

#: GUNC's maxCSS report header, exactly as GUNC writes it.
GUNC_HEADER = (
    "genome\tn_genes_called\tn_genes_mapped\ttaxonomic_level\t"
    "clade_separation_score\tpass.GUNC\n"
)

#: The same three bins with the verdict column written both ways GUNC writes it. `bin.3` is a bin
#: GUNC could not assess - 20 of 900 genes represented in the reference - which it reports as `nan`.
GUNC_AS_FLOATS = GUNC_HEADER + (
    "bin.1\t900\t880\tphylum\t0.0200\t1.0\n"
    "bin.2\t900\t850\tkingdom\t0.9100\t0.0\n"
    "bin.3\t900\t20\tphylum\t0.0100\tnan\n"
)
GUNC_AS_WORDS = GUNC_HEADER + (
    "bin.1\t900\t880\tphylum\t0.0200\tTrue\n"
    "bin.2\t900\t850\tkingdom\t0.9100\tFalse\n"
    "bin.3\t900\t20\tphylum\t0.0100\t\n"
)

#: CheckM1's bin_stats_ext.tsv: a bin name, a tab, and a Python dict repr. Not a table at all.
CHECKM1_BIN_STATS_EXT = (
    "bin.13\t{'marker lineage': 'k__Bacteria (UID203)', 'Completeness': 95.52, "
    "'Contamination': 1.23, 'GC': 0.4123, 'N50 (contigs)': 24680, 'Genome size': 3012345}\n"
)

#: CheckM1's own summary table, which spells the bin column `Bin Id`.
CHECKM1_TABLE = (
    "Bin Id\tMarker lineage\tCompleteness\tContamination\nbin.13\tk__Bacteria\t95.52\t1.23\n"
)

#: CheckM2's quality_report.tsv, with the columns in the order CheckM2 writes them.
CHECKM2_QUALITY_REPORT = (
    "Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\t"
    "Coding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\t"
    "Total_Coding_Sequences\n"
    "bin.13\t95.52\t1.23\tNeural Network (Specific Model)\t11\t0.885\t24680\t301\t3012345\t"
    "0.41\t2800\n"
)


def _report(directory, text, name="GUNC.progenomes_2.1.maxCSS_level.tsv"):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text(text)
    return str(path)


# ── GUNC's pass column, in both of the types it has ───────────────────────────────────────────


def test_a_float_typed_pass_column_is_read_as_a_verdict(tmp_path):
    """``1.0``/``0.0`` mean the same as ``True``/``False``, and both came out of one GUNC version.

    This is the regression behind a real wrong number: a 166-bin set whose column was written as
    floats scored a GUNC pass fraction of 0, i.e. every bin chimeric. 0 is believable, which is why
    it cost a week.
    """
    assert gunc.passed_value("1.0") is True
    assert gunc.passed_value("0.0") is False
    assert gunc.passed_value("True") is True
    assert gunc.passed_value("False") is False

    floats = gunc.parse_report(_report(tmp_path / "f", GUNC_AS_FLOATS))
    words = gunc.parse_report(_report(tmp_path / "w", GUNC_AS_WORDS))
    assert floats["bin.1"].passed is True and words["bin.1"].passed is True
    assert floats["bin.2"].passed is False and words["bin.2"].passed is False


def test_a_bin_gunc_could_not_assess_is_none_rather_than_a_failure():
    """``nan`` is "too little reference representation to say", which is not a failed bin.

    Read as False it becomes a chimera call that GUNC never made. Every spelling of an absence has
    to come back None, because GUNC's own writer produces several of them.
    """
    assert gunc.passed_value("nan") is None
    assert gunc.passed_value("NA") is None
    assert gunc.passed_value("") is None
    assert gunc.passed_value(None) is None
    # And not a stray True/False from Python's own truthiness of a non-empty string.
    assert gunc.passed_value("not a verdict") is None


def test_the_gunc_report_is_found_and_read_by_guncs_own_column_names(tmp_path):
    """``genome``, ``clade_separation_score``, ``taxonomic_level``, ``pass.GUNC``.

    GUNC writes several tables per run; the maxCSS one is picked by its own prefix and suffix rather
    than by position in the directory. A renamed column is the dangerous case: `pass.GUNC` going
    missing is not an error, it falls back to the CSS threshold, so the pin is on the names.
    """
    assert gunc.REPORT_PREFIX == "GUNC." and gunc.REPORT_SUFFIX == ".maxCSS_level.tsv"
    out = tmp_path / "gunc"
    out.mkdir()
    (out / "GUNC.progenomes_2.1.all_levels.tsv").write_text(GUNC_AS_FLOATS)
    (out / "GUNC.progenomes_2.1.maxCSS_level.tsv").write_text(GUNC_AS_FLOATS)
    assert gunc.find_report(str(out)).endswith(".maxCSS_level.tsv")
    assert gunc.find_report(str(tmp_path / "nothing-here")) is None

    verdicts = gunc.parse_report(gunc.find_report(str(out)))
    assert verdicts["bin.2"].css == pytest.approx(0.91)
    assert verdicts["bin.2"].level == "kingdom"
    # Without the verdict column the threshold stands in, so the bin is still judged rather than
    # dropped - but it is judged on CSS alone, which is the more forgiving of the two rules.
    no_column = _report(tmp_path, "genome\tclade_separation_score\nbin.2\t0.9100\n", "x.tsv")
    assert gunc.parse_report(no_column)["bin.2"].passed is False


def test_the_bin_key_keeps_the_dot_in_a_metawrap_bin_name(tmp_path):
    """GUNC names a genome with ``basename.split(file_suffix)[0]``, so `bin.13` stays `bin.13`.

    ``os.path.splitext`` would make it `bin`, collapsing every bin in a set onto one key and
    annotating the stats table for one bin with whichever row happened to be last.
    """
    report = _report(tmp_path, GUNC_HEADER + "bin.13\t900\t880\tphylum\t0.0200\t1.0\n")
    assert set(gunc.parse_report(report)) == {"bin.13"}
    # The run command passes the suffix *with* its dot for the same reason; "fasta" without one
    # would name every bin "bin.13." and match nothing in the stats table.
    cmd = gunc.GUNC_RUN.format(bins="b", ext=".fasta", db="db.dmnd", out="o", threads=4)
    assert "--file_suffix .fasta" in cmd


# ── the same column, read by the benchmark ────────────────────────────────────────────────────


def test_the_benchmarks_gunc_fraction_counts_floats_as_verdicts_too(tmp_path):
    """One pass, one fail and one unassessable: 1/3, not 0/3.

    The benchmark reads the column through the same helper precisely because reading it separately
    is what produced the 0% on the 166-bin set. It is a fraction rather than a count so variants
    that produced different numbers of bins stay comparable.
    """
    path = _report(tmp_path, GUNC_AS_FLOATS)
    assert scoring.gunc_pass_fraction(path) == pytest.approx(1.0 / 3.0)
    assert scoring.gunc_pass_fraction(_report(tmp_path / "w", GUNC_AS_WORDS)) == pytest.approx(
        1.0 / 3.0
    )


def test_a_report_with_no_verdict_column_scores_none_rather_than_zero(tmp_path):
    """A missing column means "this run says nothing about chimerism", not "nothing passed".

    None travels into the scoring table as an absent metric; 0.0 would travel as a measured result
    and rank the variant last on a question that was never asked of it.
    """
    path = _report(tmp_path, "genome\tclade_separation_score\nbin.1\t0.0200\n", "no_verdict.tsv")
    assert scoring.gunc_pass_fraction(path) is None
    assert scoring.gunc_unassessed_fraction(path) is None


# ── CheckM1 and CheckM2, which share nothing but their numbers ─────────────────────────────────


def test_checkm1s_extended_stats_are_python_dict_reprs_not_a_table(tmp_path):
    """``bin_stats_ext.tsv`` is a name, a tab, and a `repr` of a dict - so it has no header at all.

    The keys are CheckM1's own, spaces and parentheses included: ``N50 (contigs)``, ``Genome size``,
    ``marker lineage``. A renamed key here is a KeyError, which is the intended outcome; what must
    never happen is a zero appearing where a genome size should be.
    """
    path = tmp_path / "bin_stats_ext.tsv"
    path.write_text(CHECKM1_BIN_STATS_EXT)
    rows = checkm.summarize(str(path))
    assert checkm.STATS_HEADER == [
        "bin",
        "completeness",
        "contamination",
        "GC",
        "lineage",
        "N50",
        "size",
    ]
    assert rows == [["bin.13", "95.52", "1.23", "0.412", "Bacteria (UID203)", "24680", "3012345"]]
    # N50 before size, matching STATS_HEADER. Swapped, a 3 Mbp genome reads as a 3 Mbp N50 and
    # every bin in the set looks like a single closed contig.
    assert int(rows[0][5]) < int(rows[0][6])


def test_checkm2s_quality_report_is_read_by_name_against_its_own_spelling(tmp_path):
    """``Name``, ``Genome_Size``, ``GC_Content``, ``Contig_N50`` - CheckM2's underscores, not CheckM1's.

    Read positionally this file is a disaster: ``Contig_N50`` sits at index 6 and ``Genome_Size`` at
    8, with ``Coding_Density`` and ``Average_Gene_Length`` interleaved, so an off-by-one reports a
    coding density as a genome size. A column CheckM2 renames raises instead of parsing to empty.
    """
    path = tmp_path / "quality_report.tsv"
    path.write_text(CHECKM2_QUALITY_REPORT)
    rows = checkm.summarize_checkm2(str(path))
    assert rows == [["bin.13", "95.52", "1.23", "0.41", "checkm2", "24680", "3012345"]]

    renamed = tmp_path / "renamed.tsv"
    renamed.write_text(CHECKM2_QUALITY_REPORT.replace("Genome_Size", "Genome_size"))
    with pytest.raises(RuntimeError) as caught:
        checkm.summarize_checkm2(str(renamed))
    assert "Unexpected CheckM2 quality_report.tsv columns" in str(caught.value)


def test_both_checkm_spellings_of_the_bin_column_reach_the_scorer(tmp_path):
    """``Bin Id`` and ``Name``, plus CheckM2's ``Completeness_General``.

    A benchmark that understood one spelling scored zero bins from the other tool and reported it as
    a variant with no bins rather than as a file it could not read. The bin name is passed through
    unchanged - not extension-stripped - so `bin.13` joins to the bin set it came from.
    """
    checkm1 = tmp_path / "checkm.tsv"
    checkm1.write_text(CHECKM1_TABLE)
    checkm2 = tmp_path / "quality_report.tsv"
    checkm2.write_text(CHECKM2_QUALITY_REPORT)
    assert scoring.bin_quality(str(checkm1)) == [("bin.13", 95.52, 1.23)]
    assert scoring.bin_quality(str(checkm2)) == [("bin.13", 95.52, 1.23)]

    general = tmp_path / "general.tsv"
    general.write_text("Name\tCompleteness_General\tContamination\nbin.13\t95.52\t1.23\n")
    assert scoring.bin_quality(str(general)) == [("bin.13", 95.52, 1.23)]
