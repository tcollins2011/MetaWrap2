"""The one PDF a finished run did not have, and the page that says what is not in it.

A run writes twenty-odd correct tables and nothing a collaborator can open. The risk in assembling
them into one document is not that it looks wrong - it is that it reads as complete when it is not,
because a page whose inputs were missing simply is not there. These tests pin that: every page that
was not drawn is named on the last page, that page is written whether or not anything is missing,
and the numbers on the pages that were drawn are the ones in the files they came from.
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

    The PDF's text is glyph-indexed by the time it is on disk and cannot be read back, so the only
    way to assert what a page *says* is to watch it being written. Both page kinds are recorded,
    text pages and the embedded-PNG ones, because a wheel that was not drawn is as much a gap as a
    table that was not.
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
    """One file, ten pages, and the page count the PDF declares is the one that was drawn.

    The ten are: the summary, read QC, assembly, the genome table, abundance, one functional wheel,
    the resistance and virus summaries, "what was not run", and the provenance appendix. Pinned as
    a number because a page lost to an exception swallowed somewhere would otherwise still produce
    a plausible-looking PDF.
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
    """What was run, when and where - one page, one line per module, then straight into results.

    This is a regression: the per-run provenance dumps used to be on the title page, and on the
    reference run they filled six pages, so the first actual result was page seven. A reader
    opening the report to see what the run found read six pages of configuration first.
    """
    R.build_report([_run(tmp_path)], str(tmp_path / "r.pdf"), title="Fixture run")
    summary = "\n".join(_page(pages, "Fixture run")[1])
    assert _titles(pages)[0] == "Fixture run"
    assert _titles(pages)[1] == "Read QC"
    # The one-liner: module, dates and duration, and nothing that needs a second page.
    assert "bin_refinement" in summary and "2026-01-01T12:00:00" in summary
    assert len(_page(pages, "Fixture run")[1]) <= R.PAGE_LINES
    assert "threads" not in summary  # a resolved parameter, which belongs in the appendix


def test_the_provenance_appendix_keeps_every_parameter_tool_and_database(tmp_path, pages):
    """Moved to the back, not dropped. Read out of run_config.json rather than re-derived.

    The appendix is the answer to "what exactly was this run", so losing a block on the way to the
    back of the document would be worse than the six pages of front matter it replaced.
    """
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
    """Eight modules produced nothing for a refinement-only run, and all eight are named.

    Exactly: a page missing from this list is a page a reader would take for a negative result, and
    a module on it that did run sends somebody looking for output that is already there.
    """
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
    """Skipped, not drawn empty. An empty read-QC page would read as a sample with no reads."""
    R.build_report([_run(tmp_path, "REFINEMENT")], str(tmp_path / "r.pdf"))
    assert "Read QC" not in _titles(pages)
    assert "Assembly" not in _titles(pages)


def test_the_file_to_look_for_is_named_in_full_and_not_truncated(tmp_path, pages):
    """The point of the page is to send a reader to a file, so the filename cannot be elided.

    Laid out by hand for this reason: the table formatter truncates a cell at 34 characters, which
    turns the virus module's three filenames into "viral_contigs.tsv / plasmid_conti…".
    """
    R.build_report([_run(tmp_path, "REFINEMENT")], str(tmp_path / "r.pdf"))
    body = "\n".join(_page(pages, "What was not run")[1])
    assert "viral_contigs.tsv / plasmid_contigs.tsv / votus.tsv" in body
    assert R.ELLIPSIS not in body


def test_the_what_was_not_run_page_is_written_even_when_nothing_was_absent(tmp_path, pages):
    """Its own absence would be ambiguous: a reader could not tell whether the report had checked.

    The fixture run has no strains output and one empty Krona distribution, so both are supplied
    here - the strains tables as the headers that module writes, and the empty wheel removed.
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
    # Last of the result pages: only the provenance appendix is allowed to follow it.
    assert _titles(pages)[-2:] == ["What was not run", "Appendix: provenance"]
    assert "produced output" in "\n".join(body)
    assert footer == R.ABSENCE_NOTE


# ── the pages that were drawn ────────────────────────────────────────────────────────────────


def test_reads_in_and_out_come_from_whichever_report_the_run_left(tmp_path):
    """FastQC states the count in its HTML and fastp in its JSON; both samples get a row.

    FastQC counts are summed over the mates, because a page reporting R1 only would be wrong by a
    factor of exactly two for every paired sample and look entirely plausible.
    """
    counts = {entry.sample: entry for entry in R.read_qc_counts([_run(tmp_path, "READ_QC")])}
    assert (counts["sampleA"].reads_in, counts["sampleA"].reads_out) == (2000, 1800)
    assert (counts["sampleB"].reads_in, counts["sampleB"].reads_out) == (5000, 4500)


def test_a_count_no_report_stated_is_blank_rather_than_zero(tmp_path, pages):
    """sampleC has only a post-QC report, so its reads-in is unknown, not zero.

    Zero reads in and 700 reads out is an impossible row that a reader would have to disbelieve.
    """
    counts = {entry.sample: entry for entry in R.read_qc_counts([_run(tmp_path, "READ_QC")])}
    assert counts["sampleC"].reads_in is None
    R.build_report([_run(tmp_path, "READ_QC")], str(tmp_path / "r.pdf"))
    row = [line for line in _page(pages, "Read QC")[1] if line.startswith("sampleC")]
    assert row and "0" not in row[0].split()[1:2]


def test_assembly_numbers_are_quasts_own_rows_and_not_its_length_cuts(tmp_path, pages):
    """Copied from report.tsv, and only the metrics a summary page wants.

    QUAST writes a dozen ``>= 5000 bp`` variants of each metric; putting them on this page would
    make a reader pick one, and they would not all pick the same one.
    """
    ((name, metrics),) = R.assembly_stats([_run(tmp_path, "ASSEMBLY")])
    assert name == "final_assembly"
    assert metrics["N50"] == "4083" and metrics["Total length"] == "55479790"
    R.build_report([_run(tmp_path, "ASSEMBLY")], str(tmp_path / "r.pdf"))
    body = "\n".join(_page(pages, "Assembly")[1])
    assert "auN" not in body and ">= 5000 bp" not in body


def test_a_wheel_is_embedded_for_every_krona_input_with_a_distribution(tmp_path, pages):
    """One wheel per Krona input, by contract - drawn by sunburst, not by a second renderer.

    Grouped two to a page by unit, so the two wheels of one genome (with and without BRITE's
    scaffolding branches) are side by side. The fixture has one drawable input, so one page with
    one wheel on it - a unit with an odd number of wheels must not lose the last one.
    """
    assert len(R.krona_inputs([_run(tmp_path, "PATHWAYS")])) == 2
    R.build_report([_run(tmp_path, "PATHWAYS")], str(tmp_path / "r.pdf"))
    wheels = [entry for entry in pages if entry[0].startswith("Function:")]
    assert [entry[0] for entry in wheels] == ["Function: genome.bin_001"]
    assert len(wheels[0][1]) == 1


def test_two_wheels_of_one_unit_go_on_one_page(tmp_path, pages):
    """The pair exists to be compared, and comparing them across a page turn is not possible."""
    krona = tmp_path / "run" / "PATHWAYS" / R.KRONA_DIR
    _run(tmp_path, "PATHWAYS")
    shutil.copyfile(krona / "genome.bin_001.all.txt", krona / "genome.bin_001.pathways.txt")
    R.build_report([str(tmp_path / "run")], str(tmp_path / "r.pdf"))
    wheels = [entry for entry in pages if entry[0].startswith("Function:")]
    assert [entry[0] for entry in wheels] == ["Function: genome.bin_001"]
    assert len(wheels[0][1]) == R.WHEELS_PER_PAGE


def test_a_wheels_caption_names_its_source_relative_to_the_run_directory(tmp_path, pages):
    """An absolute path on a cluster is 60 characters before the part that identifies the file."""
    root = _run(tmp_path, "PATHWAYS")
    R.build_report([root], str(tmp_path / "r.pdf"))
    caption = _page(pages, "Function: genome.bin_001")[2]
    assert "PATHWAYS/hierarchy/krona/genome.bin_001.all.txt" in caption
    assert str(tmp_path) not in caption


def test_an_empty_distribution_is_named_rather_than_drawn_as_a_blank_wheel(tmp_path):
    """A genome with no KO at all is a real result, and a missing wheel reads as a missing genome."""
    result = R.build_report([_run(tmp_path, "PATHWAYS")], str(tmp_path / "r.pdf"))
    assert [entry for entry in result.absent if "bin_003" in entry.looked_for]


def test_the_first_genomes_page_is_headed_genomes_and_not_a_tier_name(tmp_path, pages, monkeypatch):
    """The heading has to be the page's *first* text, not just the text highest up the page.

    A regression with no visual symptom. A figure draws its children in ascending zorder, an Axes
    defaults to 0 and figure-level text to 3, so the tier bar chart went into the PDF's content
    stream ahead of the heading above it. The page looked right and read wrong: the extracted text
    of this page began "near", "complete", "high", "quality" - the chart's tick labels - with
    "Genomes (26)" tenth. Anything reading the text rather than the pixels (copy-paste, a search
    index, a screen reader) saw a quality tier where the title belongs.
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
    # The mechanism, stated as the assertion: the chart is drawn after the heading, so the heading
    # reaches the content stream first. Comparing text to 0 would pass by accident on a page with
    # no chart at all, so the chart's own axes is what is compared against.
    assert figure.axes, "the first genomes page carries the tier chart"
    assert figure.axes[0].get_zorder() > heading.get_zorder()


def test_the_genome_pages_are_paginated_rather_than_shrunk(tmp_path, pages):
    """A 300-genome catalogue scaled onto one page is a grey rectangle.

    120 genomes is the size of a real dereplicated catalogue from a small cohort.
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
    """An existing PNG, not a second heatmap: a report that re-plotted it could normalise it."""
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
    """virulence.tsv is header-only. Nothing found is a result; nothing run is not.

    These are the two states the final page exists to keep apart, and here they are in one module.
    """
    result = R.build_report([_run(tmp_path, "AMR")], str(tmp_path / "r.pdf"))
    assert "amr" not in [entry.module for entry in result.absent]
    body = "\n".join(_page(pages, "Resistance and virulence")[1])
    assert "virulence calls" in body and "resistance gene calls" in body


# ── the command ──────────────────────────────────────────────────────────────────────────────


def test_the_command_writes_the_pdf_and_the_full_table_beside_it(tmp_path, capsys):
    """The PDF truncates long lineages to stay legible, so the untruncated table ships with it."""
    out = str(tmp_path / "study.pdf")
    assert report_cmd.main([_run(tmp_path), "-o", out]) == 0
    table = str(tmp_path / ("study" + report_cmd.MAG_TABLE_SUFFIX))
    assert os.path.getsize(out) > 20_000
    header = Path(table).read_text().splitlines()[0].split("\t")
    assert header[0] == "genome" and "abundance.sampleA" in header
    assert out in capsys.readouterr().out


def test_the_command_says_what_is_not_in_the_report(tmp_path, capsys):
    """The same list as the final page, on the terminal, so it is seen before the PDF is opened."""
    assert report_cmd.main([_run(tmp_path, "REFINEMENT"), "-o", str(tmp_path / "r.pdf")]) == 0
    out = capsys.readouterr().out
    assert "strains" in out and "viral_contigs.tsv" in out
    assert R.ABSENCE_NOTE in out


def test_a_mistyped_run_directory_is_an_error_rather_than_an_empty_report(tmp_path, capsys):
    """Every page would say the module was absent, which is the one wrong answer this must not give."""
    assert report_cmd.main([str(tmp_path / "nope"), "-o", str(tmp_path / "r.pdf")]) == 1
    assert "Not a directory" in capsys.readouterr().out
    assert not os.path.exists(str(tmp_path / "r.pdf"))


# ── a study is more than one bin set ─────────────────────────────────────────────────────────


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
    """ "Genomes (24)" for a study that recovered 62 was a page a reader could only misread.

    The count alone did not say it was one sample's, and the other two samples appeared nowhere.
    """
    R.build_report(
        [_study(tmp_path, {"sample_6": 24, "sample_7": 15, "sample_13": 23})],
        str(tmp_path / "r.pdf"),
    )
    headings = [title for title in _titles(pages) if title.startswith("Genomes")]
    assert headings[0] == "Genomes (62 across 3 bin sets)"


def test_one_bin_set_keeps_the_plain_heading(tmp_path, pages):
    """No "across 1 bin sets" on a per-run report, and no sample column to explain."""
    R.build_report([_study(tmp_path, {"sample_6": 4})], str(tmp_path / "r.pdf"))
    headings = [title for title in _titles(pages) if title.startswith("Genomes")]
    assert headings == ["Genomes (4)"]


def test_the_footer_says_which_bin_sets_an_input_could_not_be_joined_to(tmp_path, pages):
    """Blank because the file does not describe these genomes, not because it was never written.

    Both leave identical empty cells, so only this sentence tells a reader which happened - and in
    a cohort run, where taxonomy is keyed on dereplicated names, the second is the normal case.
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
