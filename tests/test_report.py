"""The run report PDF and its "What was not run" page.

A page whose inputs are missing is left out, so the report could look complete when it is not.
Checked: every page not drawn is named on the "What was not run" page, that page is always written,
and the numbers on drawn pages match their source files.
"""

from __future__ import annotations

import os
import re
import shutil
from pathlib import Path

import pytest

from metawrap2 import report as R
from metawrap2.commands import report_cmd

FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures", "report")
#: Every module output directory in the fixture run. See fixtures/report/SOURCE.md.
ALL_PARTS = ("REFINEMENT", "CLASSIFY", "QUANT", "READ_QC", "ASSEMBLY", "PATHWAYS", "AMR", "VIRUSES")

#: How many ``/Type /Page`` objects a PDF declares. Matched with a trailing non-``s`` so ``/Pages``,
#: the single node that lists them, is not counted as a tenth page.
PDF_PAGE = re.compile(rb"/Type\s*/Page[^s]")


def _run(tmp_path, *parts: str) -> str:
    """A run directory holding only *parts* of the fixture run, so inputs can be taken away."""
    root = tmp_path / "run"
    root.mkdir(exist_ok=True)
    for part in parts or ALL_PARTS:
        shutil.copytree(os.path.join(FIXTURES, part), root / part, dirs_exist_ok=True)
    return str(root)


@pytest.fixture
def pages(monkeypatch):
    """Every page's title and body, recorded as the document is written.

    The PDF's text is glyph-indexed on disk and cannot be read back, so pages are recorded as they
    are written. Both text pages and embedded-PNG pages are recorded.
    """
    recorded: list = []
    text_page, image_page = R._Pdf.page, R._Pdf.image_page

    def page(self, title, lines, footer="", chart=None):
        recorded.append((title, list(lines), footer))
        text_page(self, title, lines, footer, chart)

    def image(self, title, images, caption=""):
        recorded.append((title, list(images), caption))
        image_page(self, title, images, caption)

    monkeypatch.setattr(R._Pdf, "page", page)
    monkeypatch.setattr(R._Pdf, "image_page", image)
    return recorded


def _titles(pages: list) -> list:
    return [title for title, _lines, _footer in pages]


def _page(pages: list, title: str) -> tuple:
    for entry in pages:
        if entry[0] == title:
            return entry
    raise AssertionError("no page titled %r; pages were %s" % (title, _titles(pages)))


# ── the document ─────────────────────────────────────────────────────────────────────────────


def test_a_pdf_is_written_with_a_page_for_every_section(tmp_path, pages):
    """One file with ten pages, and the PDF's declared page count matches the pages drawn.

    The ten are: the summary, read QC, assembly, the genome table, abundance, one functional wheel,
    the resistance and virus summaries, "what was not run", and the provenance appendix. The count
    catches a page lost to a swallowed exception.
    """
    out = str(tmp_path / "report.pdf")
    result = R.build_report([_run(tmp_path)], out)
    data = Path(out).read_bytes()
    assert result.pages == 10
    assert len(PDF_PAGE.findall(data)) == result.pages
    assert len(data) > 20_000, "a ten-page report with a wheel in it is not a 20 kB file"
    assert _titles(pages) == [
        "MetaWrap2 report",
        "Read QC",
        "Assembly",
        "Genomes (4)",
        "Abundance across 2 sample(s)",
        "Function: genome.bin_001",
        "Resistance and virulence",
        "Viruses",
        "What was not run",
        "Appendix: provenance",
    ]


def test_the_first_page_is_a_summary_and_the_first_result_is_the_second_page(tmp_path, pages):
    """Page one is a one-line-per-module summary and page two is the first result.

    The provenance dumps used to be on the title page and filled six pages on the reference run.
    """
    R.build_report([_run(tmp_path)], str(tmp_path / "r.pdf"), title="Fixture run")
    summary = "\n".join(_page(pages, "Fixture run")[1])
    assert _titles(pages)[0] == "Fixture run"
    assert _titles(pages)[1] == "Read QC"
    # Module, dates and duration, on one page.
    assert "bin_refinement" in summary and "2026-01-01T12:00:00" in summary
    assert len(_page(pages, "Fixture run")[1]) <= R.PAGE_LINES
    assert "threads" not in summary  # a resolved parameter, which belongs in the appendix


def test_the_provenance_appendix_keeps_every_parameter_tool_and_database(tmp_path, pages):
    """The appendix, read from run_config.json, keeps all parameters, tools and databases."""
    R.build_report([_run(tmp_path, "REFINEMENT")], str(tmp_path / "r.pdf"))
    body = "\n".join(_page(pages, "Appendix: provenance")[1])
    assert "bin_refinement" in body and "20260101T120000-abcd1234" in body
    assert "2026-01-01T12:00:00" in body and "2026-01-01T12:41:00" in body
    assert "GUNC 1.0.6" in body  # tool_versions
    assert "gunc_db_progenomes2.1.dmnd" in body  # databases
    assert "contamination" in body and "10.0" in body  # resolved parameters
    assert "threads" in body and "8" in body
    assert _titles(pages)[-1].startswith("Appendix: provenance")


# ── what was not run ─────────────────────────────────────────────────────────────────────────


def test_the_final_page_names_exactly_the_modules_whose_output_was_absent(tmp_path, pages):
    """A refinement-only run lists exactly the eight modules with no output."""
    result = R.build_report([_run(tmp_path, "REFINEMENT")], str(tmp_path / "r.pdf"))
    assert [entry.module for entry in result.absent] == [
        "read_qc",
        "assembly",
        "classify_bins",
        "quant_bins",
        "pathways / annotate",
        "amr",
        "viruses",
        "strains",
    ]
    body = "\n".join(_page(pages, "What was not run")[1])
    for entry in result.absent:
        assert entry.module in body and entry.looked_for in body


def test_a_module_that_was_not_run_gets_no_page_of_its_own(tmp_path, pages):
    """A module with no output gets no page, not an empty one."""
    R.build_report([_run(tmp_path, "REFINEMENT")], str(tmp_path / "r.pdf"))
    assert "Read QC" not in _titles(pages)
    assert "Assembly" not in _titles(pages)


def test_the_file_to_look_for_is_named_in_full_and_not_truncated(tmp_path, pages):
    """Filenames on the "What was not run" page are not truncated.

    The page is laid out by hand because the table formatter truncates cells at 34 characters
    ("viral_contigs.tsv / plasmid_conti…").
    """
    R.build_report([_run(tmp_path, "REFINEMENT")], str(tmp_path / "r.pdf"))
    body = "\n".join(_page(pages, "What was not run")[1])
    assert "viral_contigs.tsv / plasmid_contigs.tsv / votus.tsv" in body
    assert R.ELLIPSIS not in body


def test_the_what_was_not_run_page_is_written_even_when_nothing_was_absent(tmp_path, pages):
    """The "What was not run" page is written even when nothing is absent.

    The fixture run lacks strains output and has one empty Krona distribution, so strains tables
    are added and the empty distribution removed.
    """
    root = _run(tmp_path)
    os.remove(os.path.join(root, "PATHWAYS", "hierarchy", "krona", "genome.bin_003.all.txt"))
    strains = os.path.join(root, "STRAINS")
    os.mkdir(strains)
    with open(os.path.join(strains, "populations.tsv"), "w") as fh:
        fh.write(
            "genome\tsample\tcoverage\tbreadth\tnucl_diversity\nbin_001\tsampleA\t24.1\t0.98\t0.0012\n"
        )
    with open(os.path.join(strains, "strain_comparisons.tsv"), "w") as fh:
        fh.write("genome\tsample_1\tsample_2\tpopANI\tpercent_compared\n")

    result = R.build_report([root], str(tmp_path / "r.pdf"))
    assert result.absent == []
    _title, body, footer = _page(pages, "What was not run")
    # Last result page; only the provenance appendix follows it.
    assert _titles(pages)[-2:] == ["What was not run", "Appendix: provenance"]
    assert "produced output" in "\n".join(body)
    assert footer == R.ABSENCE_NOTE


# ── the pages that were drawn ────────────────────────────────────────────────────────────────


def test_reads_in_and_out_come_from_whichever_report_the_run_left(tmp_path):
    """Read counts come from FastQC's HTML or fastp's JSON; both samples get a row.

    FastQC counts are summed over both mates; R1 alone would be half the paired-sample total.
    """
    counts = {entry.sample: entry for entry in R.read_qc_counts([_run(tmp_path, "READ_QC")])}
    assert (counts["sampleA"].reads_in, counts["sampleA"].reads_out) == (2000, 1800)
    assert (counts["sampleB"].reads_in, counts["sampleB"].reads_out) == (5000, 4500)


def test_a_count_no_report_stated_is_blank_rather_than_zero(tmp_path, pages):
    """sampleC has only a post-QC report, so its reads-in is blank, not zero."""
    counts = {entry.sample: entry for entry in R.read_qc_counts([_run(tmp_path, "READ_QC")])}
    assert counts["sampleC"].reads_in is None
    R.build_report([_run(tmp_path, "READ_QC")], str(tmp_path / "r.pdf"))
    row = [line for line in _page(pages, "Read QC")[1] if line.startswith("sampleC")]
    assert row and "0" not in row[0].split()[1:2]


def test_assembly_numbers_are_quasts_own_rows_and_not_its_length_cuts(tmp_path, pages):
    """Assembly metrics are copied from QUAST's report.tsv, without the length-cut variants.

    QUAST writes about a dozen ``>= 5000 bp`` variants of each metric; only the base metrics are
    shown.
    """
    ((name, metrics),) = R.assembly_stats([_run(tmp_path, "ASSEMBLY")])
    assert name == "final_assembly"
    assert metrics["N50"] == "4083" and metrics["Total length"] == "55479790"
    R.build_report([_run(tmp_path, "ASSEMBLY")], str(tmp_path / "r.pdf"))
    body = "\n".join(_page(pages, "Assembly")[1])
    assert "auN" not in body and ">= 5000 bp" not in body


def test_a_wheel_is_embedded_for_every_krona_input_with_a_distribution(tmp_path, pages):
    """One sunburst wheel per Krona input with a distribution.

    Wheels are grouped two per page by unit, so a genome's two wheels (with and without BRITE's
    classification branches) are side by side. The fixture has one drawable input: one page with
    one wheel, so an odd last wheel is not lost.
    """
    assert len(R.krona_inputs([_run(tmp_path, "PATHWAYS")])) == 2
    R.build_report([_run(tmp_path, "PATHWAYS")], str(tmp_path / "r.pdf"))
    wheels = [entry for entry in pages if entry[0].startswith("Function:")]
    assert [entry[0] for entry in wheels] == ["Function: genome.bin_001"]
    assert len(wheels[0][1]) == 1


def test_two_wheels_of_one_unit_go_on_one_page(tmp_path, pages):
    """A unit's two wheels are drawn on the same page."""
    krona = tmp_path / "run" / "PATHWAYS" / R.KRONA_DIR
    _run(tmp_path, "PATHWAYS")
    shutil.copyfile(krona / "genome.bin_001.all.txt", krona / "genome.bin_001.pathways.txt")
    R.build_report([str(tmp_path / "run")], str(tmp_path / "r.pdf"))
    wheels = [entry for entry in pages if entry[0].startswith("Function:")]
    assert [entry[0] for entry in wheels] == ["Function: genome.bin_001"]
    assert len(wheels[0][1]) == R.WHEELS_PER_PAGE


def test_a_wheels_caption_names_its_source_relative_to_the_run_directory(tmp_path, pages):
    """Wheel captions give the source path relative to the run directory."""
    root = _run(tmp_path, "PATHWAYS")
    R.build_report([root], str(tmp_path / "r.pdf"))
    caption = _page(pages, "Function: genome.bin_001")[2]
    assert "PATHWAYS/hierarchy/krona/genome.bin_001.all.txt" in caption
    assert str(tmp_path) not in caption


def test_an_empty_distribution_is_named_rather_than_drawn_as_a_blank_wheel(tmp_path):
    """A genome with no KO is listed as absent instead of drawn as a blank wheel."""
    result = R.build_report([_run(tmp_path, "PATHWAYS")], str(tmp_path / "r.pdf"))
    assert [entry for entry in result.absent if "bin_003" in entry.looked_for]


def test_the_first_genomes_page_is_headed_genomes_and_not_a_tier_name(tmp_path, pages, monkeypatch):
    """The heading is the first text in the page's content stream.

    A figure draws children in ascending zorder (Axes default 0, figure text 3), so the tier chart
    was written before the heading. Extracted text (copy-paste, search, screen readers) then began
    with the chart's tick labels, and "Genomes (26)" came tenth.
    """
    captured: list = []
    add = R._Pdf.add

    def record(self, figure):
        captured.append(figure)
        add(self, figure)

    monkeypatch.setattr(R._Pdf, "add", record)
    R.build_report([_run(tmp_path, "REFINEMENT")], str(tmp_path / "r.pdf"))

    at = _titles(pages).index("Genomes (3)")
    figure = captured[at]
    heading = figure.texts[0]
    assert heading.get_text() == "Genomes (3)"
    # The chart's axes must have a higher zorder than the heading, so it is drawn after it. Compare
    # against the chart's own axes; comparing to 0 would pass on a page without a chart.
    assert figure.axes, "the first genomes page carries the tier chart"
    assert figure.axes[0].get_zorder() > heading.get_zorder()


def test_the_genome_pages_are_paginated_rather_than_shrunk(tmp_path, pages):
    """Large genome tables span several pages instead of being scaled down.

    120 genomes is the size of a dereplicated catalogue from a small cohort.
    """
    root = tmp_path / "big" / "REFINE"
    root.mkdir(parents=True)
    with open(root / "metawrap_70_10_bins.stats.tsv", "w") as fh:
        fh.write("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n")
        fh.writelines(
            "bin_%03d\t90.5\t1.5\t0.41\tClostridiales\t40000\t2500000\n" % number
            for number in range(1, 121)
        )
    R.build_report([str(tmp_path / "big")], str(tmp_path / "r.pdf"))
    genome_pages = [title for title in _titles(pages) if title.startswith("Genomes (120)")]
    assert len(genome_pages) > 1
    assert genome_pages[1].endswith("(continued)")


def test_the_abundance_heatmap_quant_bins_already_drew_is_embedded(tmp_path, pages):
    """quant_bins' existing heatmap PNG is embedded instead of re-plotted."""
    root = _run(tmp_path, "QUANT")
    heatmap = os.path.join(root, "QUANT", R.ABUNDANCE_HEATMAP)
    shutil.copyfile(
        os.path.join(FIXTURES, "PATHWAYS", "hierarchy", "krona", "genome.bin_001.all.txt"),
        os.path.join(root, "QUANT", "unused.txt"),
    )
    from metawrap2 import sunburst

    sunburst.draw_from_krona_text(os.path.join(root, "QUANT", "unused.txt"), heatmap)
    R.build_report([root], str(tmp_path / "r.pdf"))
    assert _page(pages, "Abundance heatmap")[1] == [(heatmap, "")]


def test_a_module_that_ran_and_found_nothing_is_not_reported_as_absent(tmp_path, pages):
    """A module that ran and found nothing (header-only virulence.tsv) is not listed as absent."""
    result = R.build_report([_run(tmp_path, "AMR")], str(tmp_path / "r.pdf"))
    assert "amr" not in [entry.module for entry in result.absent]
    body = "\n".join(_page(pages, "Resistance and virulence")[1])
    assert "virulence calls" in body and "resistance gene calls" in body


# ── the command ──────────────────────────────────────────────────────────────────────────────


def test_the_command_writes_the_pdf_and_the_full_table_beside_it(tmp_path, capsys):
    """The command writes the PDF and an untruncated genome table beside it."""
    out = str(tmp_path / "study.pdf")
    assert report_cmd.main([_run(tmp_path), "-o", out]) == 0
    table = str(tmp_path / ("study" + report_cmd.MAG_TABLE_SUFFIX))
    assert os.path.getsize(out) > 20_000
    header = Path(table).read_text().splitlines()[0].split("\t")
    assert header[0] == "genome" and "abundance.sampleA" in header
    assert out in capsys.readouterr().out


def test_the_command_says_what_is_not_in_the_report(tmp_path, capsys):
    """The command prints the "What was not run" list to the terminal."""
    assert report_cmd.main([_run(tmp_path, "REFINEMENT"), "-o", str(tmp_path / "r.pdf")]) == 0
    out = capsys.readouterr().out
    assert "strains" in out and "viral_contigs.tsv" in out
    assert R.ABSENCE_NOTE in out


def test_a_mistyped_run_directory_is_an_error_rather_than_an_empty_report(tmp_path, capsys):
    """A nonexistent run directory is an error and writes no PDF."""
    assert report_cmd.main([str(tmp_path / "nope"), "-o", str(tmp_path / "r.pdf")]) == 1
    assert "Not a directory" in capsys.readouterr().out
    assert not os.path.exists(str(tmp_path / "r.pdf"))


# ── studies with several bin sets ────────────────────────────────────────────────────────────


def _study(tmp_path, bins_per_sample: dict) -> str:
    """A per-sample study, laid out as ``metawrap2 run`` writes one: one bin set per sample."""
    root = tmp_path / "study" / "BIN_REFINEMENT"
    for sample, count in bins_per_sample.items():
        directory = root / sample
        directory.mkdir(parents=True)
        with open(directory / "metawrap_70_10_bins.stats.tsv", "w") as fh:
            fh.write("bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n")
            fh.writelines(
                "bin_%03d\t95.0\t1.0\t0.41\tClostridiales\t40000\t2500000\n" % number
                for number in range(1, count + 1)
            )
    return str(tmp_path / "study")


def test_the_genomes_heading_says_how_many_bin_sets_it_covers(tmp_path, pages):
    """The genomes heading gives the total and the number of bin sets.

    It used to read "Genomes (24)" for one sample of a 62-genome, three-sample study.
    """
    R.build_report(
        [_study(tmp_path, {"sample_6": 24, "sample_7": 15, "sample_13": 23})],
        str(tmp_path / "r.pdf"),
    )
    headings = [title for title in _titles(pages) if title.startswith("Genomes")]
    assert headings[0] == "Genomes (62 across 3 bin sets)"


def test_one_bin_set_keeps_the_plain_heading(tmp_path, pages):
    """A single bin set gets the plain heading."""
    R.build_report([_study(tmp_path, {"sample_6": 4})], str(tmp_path / "r.pdf"))
    headings = [title for title in _titles(pages) if title.startswith("Genomes")]
    assert headings == ["Genomes (4)"]


def test_the_footer_says_which_bin_sets_an_input_could_not_be_joined_to(tmp_path, pages):
    """The footer names the bin sets an input could not be joined to.

    Unjoined and unwritten inputs both leave empty cells; the footer distinguishes them. In a cohort
    run, taxonomy is keyed on dereplicated genome names and usually does not join per-sample sets.
    """
    study = _study(tmp_path, {"sample_6": 2, "sample_7": 2})
    classify = tmp_path / "study" / "CLASSIFY"
    classify.mkdir()
    (classify / "bin_taxonomy.tsv").write_text("derep_cluster_1.fasta\tBacillota;Clostridia\n")
    R.build_report([study], str(tmp_path / "r.pdf"))
    footer = _page(pages, "Genomes (5 across 2 bin sets)")[2]
    assert "matching no genome in the bin set named" in footer
    assert "taxonomy for sample_6, sample_7" in footer
    assert "dereplicated genome names" in footer
