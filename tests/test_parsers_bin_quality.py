"""The bin-quality tables other tools write: GUNC's maxCSS report, and CheckM1's and CheckM2's.

Regression tests for third-party file formats. An upstream release can rename a column, change a
unit or retype a boolean, and the result is a wrong number, not an exception:

* GUNC writes ``pass.GUNC`` as ``True``/``False``, or as ``1.0``/``0.0`` when a ``nan`` forces a
  float dtype. A reader that accepted only the words scored a 166-bin set at 0% non-chimeric.
* CheckM1 names the bin column ``Bin Id`` and CheckM2 ``Name``; a reader that knew one scored
  nothing from the other, without an error.

The tables are written in each tool's layout inside the tests, as in ``tests/test_gunc.py`` and
``tests/test_benchmark_scoring.py``, so the column names under test are visible.
"""

from __future__ import annotations

import pytest

from metawrap2 import checkm, gunc
from metawrap2.benchmark import scoring

#: GUNC's maxCSS report header, as GUNC writes it.
GUNC_HEADER = (
    "genome\tn_genes_called\tn_genes_mapped\ttaxonomic_level\t"
    "clade_separation_score\tpass.GUNC\n"
)

#: The same three bins with the verdict column in both of GUNC's forms. GUNC could not assess
#: `bin.3` (20 of 900 genes represented in the reference) and reports it as `nan`.
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

#: CheckM1's bin_stats_ext.tsv: a bin name, a tab, and a Python dict repr.
CHECKM1_BIN_STATS_EXT = (
    "bin.13\t{'marker lineage': 'k__Bacteria (UID203)', 'Completeness': 95.52, "
    "'Contamination': 1.23, 'GC': 0.4123, 'N50 (contigs)': 24680, 'Genome size': 3012345}\n"
)

#: CheckM1's summary table, with the bin column named `Bin Id`.
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


# ── GUNC's pass column, boolean and float ─────────────────────────────────────────────────────


def test_a_float_typed_pass_column_is_read_as_a_verdict(tmp_path):
    """``1.0``/``0.0`` are read as ``True``/``False``; one GUNC version writes both forms.

    A 166-bin set with a float column previously scored a GUNC pass fraction of 0.
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
    """``nan`` (too little reference representation) is read as None, not as a failed bin.

    GUNC writes absence in several forms, and each must return None.
    """
    assert gunc.passed_value("nan") is None
    assert gunc.passed_value("NA") is None
    assert gunc.passed_value("") is None
    assert gunc.passed_value(None) is None
    # Not True/False from the truthiness of a non-empty string.
    assert gunc.passed_value("not a verdict") is None


def test_the_gunc_report_is_found_and_read_by_guncs_own_column_names(tmp_path):
    """``genome``, ``clade_separation_score``, ``taxonomic_level``, ``pass.GUNC``.

    GUNC writes several tables per run; the maxCSS one is found by prefix and suffix. A missing
    `pass.GUNC` column falls back to the CSS threshold without an error, so the names are tested.
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
    # Without the verdict column the bin is judged on the CSS threshold alone, the more lenient
    # of the two rules.
    no_column = _report(tmp_path, "genome\tclade_separation_score\nbin.2\t0.9100\n", "x.tsv")
    assert gunc.parse_report(no_column)["bin.2"].passed is False


def test_the_bin_key_keeps_the_dot_in_a_metawrap_bin_name(tmp_path):
    """GUNC names a genome with ``basename.split(file_suffix)[0]``, so `bin.13` stays `bin.13`.

    ``os.path.splitext`` would give `bin` for every bin, so all rows would share one key.
    """
    report = _report(tmp_path, GUNC_HEADER + "bin.13\t900\t880\tphylum\t0.0200\t1.0\n")
    assert set(gunc.parse_report(report)) == {"bin.13"}
    # The run command passes the suffix with its dot; "fasta" alone would give "bin.13." and
    # match nothing in the stats table.
    cmd = gunc.GUNC_RUN.format(bins="b", ext=".fasta", db="db.dmnd", out="o", threads=4)
    assert "--file_suffix .fasta" in cmd


# ── the pass column in the benchmark ──────────────────────────────────────────────────────────


def test_the_benchmarks_gunc_fraction_counts_floats_as_verdicts_too(tmp_path):
    """One pass, one fail and one unassessable: 1/3, not 0/3.

    The benchmark reads the column through the same helper (a separate reader produced the 0% on
    the 166-bin set). It is a fraction so variants with different bin counts can be compared.
    """
    path = _report(tmp_path, GUNC_AS_FLOATS)
    assert scoring.gunc_pass_fraction(path) == pytest.approx(1.0 / 3.0)
    assert scoring.gunc_pass_fraction(_report(tmp_path / "w", GUNC_AS_WORDS)) == pytest.approx(
        1.0 / 3.0
    )


def test_a_report_with_no_verdict_column_scores_none_rather_than_zero(tmp_path):
    """A missing verdict column gives None (no chimerism data), not 0.0.

    None is an absent metric in the scoring table; 0.0 would rank the variant last.
    """
    path = _report(tmp_path, "genome\tclade_separation_score\nbin.1\t0.0200\n", "no_verdict.tsv")
    assert scoring.gunc_pass_fraction(path) is None
    assert scoring.gunc_unassessed_fraction(path) is None


# ── CheckM1 and CheckM2 ────────────────────────────────────────────────────────────────────────


def test_checkm1s_extended_stats_are_python_dict_reprs_not_a_table(tmp_path):
    """``bin_stats_ext.tsv`` is a name, a tab and a dict `repr`, with no header.

    The keys are CheckM1's, with spaces and parentheses: ``N50 (contigs)``, ``Genome size``,
    ``marker lineage``. A renamed key raises KeyError instead of giving a zero.
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
    # N50 before size, matching STATS_HEADER.
    assert int(rows[0][5]) < int(rows[0][6])


def test_checkm2s_quality_report_is_read_by_name_against_its_own_spelling(tmp_path):
    """``Name``, ``Genome_Size``, ``GC_Content``, ``Contig_N50``: CheckM2's column names.

    Columns are read by name: ``Contig_N50`` is at index 6 and ``Genome_Size`` at 8, with
    ``Coding_Density`` and ``Average_Gene_Length`` between them. A renamed column raises.
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

    A benchmark that knew one name scored zero bins from the other tool. The bin name is passed
    through unchanged (not extension-stripped), so `bin.13` joins to its bin set.
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
