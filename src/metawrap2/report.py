"""One PDF per run, assembled from the files the run already wrote.

A finished study is twenty-odd correct TSVs, a handful of PNGs and a Krona HTML per genome. There
is no single place to look at it and nothing a collaborator can open: the person who ran it knows
which file answers which question, and nobody else does. This builds the thing that can be attached
to a message - a one-page summary of what was run, read QC, assembly, the per-MAG table,
abundance, the functional wheels, and the resistance/virus/strain summaries - in that order, in one
PDF, with the full per-run provenance as an appendix at the back.

**Every page is a picture of a file.** Nothing here measures, estimates, normalises or re-scores
anything. Counts are row counts of tables that already exist, the quality tiers come from
:mod:`metawrap2.magtable`, and the functional wheels are drawn by :mod:`metawrap2.sunburst` from the
same Krona input file Krona is given, so the wheel in the PDF and the wheel in the browser cannot
disagree. If a page here ever needs a number that is not already on disk, the module that owns that
number should write it and this should read it.

**Absence is stated, never implied.** A page whose inputs are missing is skipped, and every skipped
page is named on a "what was not run" page that is not optional - the last page before the
appendix. A report that simply had no resistance page
would be read as a run with no resistance genes, which is the one mistake a document of this shape
can make that matters: "not run" and "ran and found nothing" are different results, and a reader
cannot tell them apart from a gap.

``matplotlib``'s PDF backend is the whole toolchain on purpose. It is already a hard dependency and
needs no network, no browser and no system libraries beyond the ones the bundled plotting scripts
already require - so this renders the same on a cluster node with no outbound route as on a laptop.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
from typing import TYPE_CHECKING, Callable, Dict, List, NamedTuple, Optional, Sequence, Tuple

from . import __version__
from . import magtable as _magtable
from . import sunburst as _sunburst
from .constants import STATS_SUFFIX
from .provenance import CONFIG_NAME, ENVIRONMENT_NAME, PROVENANCE_NAME

if TYPE_CHECKING:  # import-time cost of matplotlib is paid by the command, not by `metawrap2 -h`
    from matplotlib.axes import Axes
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.figure import Figure

#: A4 landscape, in inches. Landscape because the per-MAG table is the reason this document exists
#: and it is twelve columns wide before a single sample is quantified; in portrait it either wraps
#: or shrinks to an unreadable font.
PAGE_SIZE = (11.69, 8.27)

#: Monospaced body text at this size, and the width in characters that fits a page at it. Tables
#: are laid out as text rather than with ``Axes.table``, which cannot paginate and sizes its cells
#: from the figure rather than from the content.
BODY_FONT = 7.5
MIN_FONT = 4.0
PAGE_CHARS = 158
#: Body lines that fit under a title, and under a title plus a chart.
PAGE_LINES = 44
PAGE_LINES_WITH_CHART = 26

#: zorder for an axes that shares a page with the page heading, so the heading is drawn - and
#: therefore written into the PDF's content stream - *first*.
#:
#: Not cosmetic. A figure's children are drawn in ascending zorder, an Axes defaults to 0 and
#: figure-level text to 3, so a page carrying a chart put the chart's artists into the stream ahead
#: of its own heading. The page looked right and read wrong: the extracted text of the first genomes
#: page began "near", "complete", "high", "quality" - the bar chart's tick labels - and its heading
#: came tenth. Anything that reads a PDF's text rather than its pixels (copy-paste, a search index,
#: a screen reader, `pymupdf.Page.get_text`) saw a tier name where the title belongs.
ARTIST_OVER_TEXT = 4

#: A cell wider than this is truncated. GTDB lineages reach 120 characters and one of them would
#: set the width of a whole column, which at twelve columns shrinks the font past legibility. The
#: untruncated value is in the TSV written beside the PDF, and the page footer says so.
MAX_CELL = 34
ELLIPSIS = "…"

#: How deep below a run directory to look. Shared with :mod:`metawrap2.magtable` so the two agree
#: about what "under this run directory" means - a report that found the stats table but not the
#: assembly beside it would be a confusing half-report.
MAX_DEPTH = _magtable.MAX_DEPTH

#: QUAST metrics worth a page, in the order a reader wants them. QUAST writes 40-odd rows and most
#: of the rest are its ``>= 5000 bp`` cuts, which belong in QUAST's own report.
QUAST_METRICS = ("# contigs", "Total length", "N50", "N90", "Largest contig", "GC (%)")

#: FastQC states the read count in its HTML report and in ``fastqc_data.txt``. Matched loosely
#: because the markup has gained whitespace and attributes between FastQC versions, and a pattern
#: that silently stopped matching would leave the read-QC page blank rather than absent - which is
#: the exact confusion the "what was not run" page exists to prevent.
TOTAL_SEQUENCES_HTML = re.compile(r"Total Sequences\s*</td>\s*<td>\s*([0-9,]+)")
TOTAL_SEQUENCES_TEXT = re.compile(r"^Total Sequences\t([0-9,]+)", re.MULTILINE)

#: The two FastQC report directories ``read_qc`` writes, before and after cleaning, and fastp's
#: JSON - preferred where the run used ``--fastp``, because it states both counts in one file.
PRE_QC_DIR = "pre-QC_report"
POST_QC_DIR = "post-QC_report"
FASTP_JSON = "fastp.json"

#: Where ``pathways`` leaves the Krona inputs, and what ``quant_bins`` calls the heatmap it draws.
KRONA_DIR = os.path.join("hierarchy", "krona")
ABUNDANCE_HEATMAP = "bin_abundance_heatmap.png"

#: Functional wheels to a page. Two, because ``pathways`` writes two per unit and they are the same
#: distribution with and without BRITE's scaffolding branches - a pair that exists to be compared.
#: At one per page the reference run's 18 wheels were 18 pages of a 31-page report.
WHEELS_PER_PAGE = 2

#: Printed on the "what was not run" page. The wording is deliberate: a reader who takes a missing page for a
#: negative result has been misled by the document rather than by the data.
ABSENCE_NOTE = (
    "A module listed here produced no output under the directories this report was built from. "
    "That is not a negative result - nothing looked. Re-run the module, or point this report at "
    "the directory that holds its output."
)

__all__ = [
    "ABSENCE_NOTE",
    "PAGE_SIZE",
    "Absence",
    "ReadCounts",
    "Report",
    "assembly_stats",
    "build_report",
    "krona_inputs",
    "read_qc_counts",
]


class Absence(NamedTuple):
    """A page that was not drawn, and what was looked for."""

    module: str
    looked_for: str


class ReadCounts(NamedTuple):
    """Reads in and out for one sample. Either may be None, meaning no report stated it."""

    sample: str
    reads_in: Optional[int]
    reads_out: Optional[int]


class Report(NamedTuple):
    """What :func:`build_report` wrote."""

    path: str
    pages: int
    absent: List[Absence]
    #: The per-MAG table, written beside the PDF because it is the only place the untruncated
    #: lineages and full-precision abundances exist.
    mag_table: Optional[str]


# ── finding things under a run directory ─────────────────────────────────────────────────────


def _dirs(run_dirs: Sequence[str]) -> List[str]:
    """Each run directory and its subdirectories, at most :data:`MAX_DEPTH` down, deduplicated.

    Two depths matter in practice: a caller may point this at one module's output directory
    (``READ_QC/ERR011347``) or at the study directory that holds all of them, and the same marker
    file is one level further down in the second case.
    """
    seen: List[str] = []
    for root in run_dirs:
        if not os.path.isdir(root):
            continue
        root = os.path.abspath(root)
        base_depth = root.count(os.sep)
        for parent, subdirs, _names in os.walk(root):
            if parent.count(os.sep) - base_depth >= MAX_DEPTH:
                subdirs[:] = []
            subdirs.sort()
            if parent not in seen:
                seen.append(parent)
    return seen


def _with(run_dirs: Sequence[str], *relative: str) -> List[str]:
    """Every existing ``<dir>/<relative>`` under *run_dirs*, in directory order."""
    found: List[str] = []
    for directory in _dirs(run_dirs):
        path = os.path.join(directory, *relative)
        if os.path.exists(path):
            found.append(path)
    return found


def _relative(path: str, run_dirs: Sequence[str]) -> str:
    """*path* as ``<run directory>/<the rest>``, or unchanged if it is under none of them.

    A caption is the one line on a page with no room to spare, and an absolute path can spend most
    of it on directories before the part that identifies the file. The run directory's own name is kept at the front so the caption still says *which*
    of several inputs a wheel came from.
    """
    target = os.path.abspath(path)
    for root in run_dirs:
        root = os.path.abspath(root)
        if target == root or target.startswith(root + os.sep):
            return os.path.join(os.path.basename(root), os.path.relpath(target, root))
    return path


def _table(path: str) -> Tuple[List[str], List[List[str]]]:
    """A TSV with a header row as ``(header, rows)``; ``([], [])`` for an empty file."""
    rows = list(_magtable.read_tsv(path))
    return (rows[0], rows[1:]) if rows else ([], [])


def _row_count(path: str) -> int:
    _header, rows = _table(path)
    return len(rows)


def _column_counts(path: str, column: str) -> List[Tuple[str, int]]:
    """How many rows hold each distinct value of *column*, commonest first.

    A row count of an existing table, grouped by one of that table's own columns. Nothing is
    weighted, deduplicated or normalised, so this cannot disagree with the table it counts. An
    unrecognised column name gives nothing rather than a guess at a similar one.
    """
    header, rows = _table(path)
    lowered = [name.strip().lower() for name in header]
    if column.lower() not in lowered:
        return []
    at = lowered.index(column.lower())
    counts: Dict[str, int] = {}
    for row in rows:
        value = (row[at].strip() if at < len(row) else "") or _magtable.NOT_MEASURED
        counts[value] = counts.get(value, 0) + 1
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))


# ── reading each module's output ─────────────────────────────────────────────────────────────


def _int(text: str) -> Optional[int]:
    try:
        return int(text.replace(",", ""))
    except ValueError:
        return None


def _search(pattern: re.Pattern[str], path: str) -> Optional[int]:
    with open(path, errors="replace") as fh:
        match = pattern.search(fh.read())
    return _int(match.group(1)) if match else None


def _fastqc_total(report_dir: str) -> Optional[int]:
    """Reads in a FastQC report directory, summed over the mates it reported on.

    Summed rather than taken from R1, because a page that said "12.6 M reads" for a paired sample
    holding 25.2 M would be wrong by exactly a factor of two and look entirely plausible. The page
    footer says the number is a sum over mates.

    None, not 0, when no report states a count: FastQC's HTML has changed shape between versions,
    and a zero here would read as a sample with no reads in it.
    """
    if not os.path.isdir(report_dir):
        return None
    total: Optional[int] = None
    for name in sorted(os.listdir(report_dir)):
        path = os.path.join(report_dir, name)
        count: Optional[int] = None
        if name.endswith("_fastqc.html") and os.path.isfile(path):
            count = _search(TOTAL_SEQUENCES_HTML, path)
        elif name.endswith("_fastqc") and os.path.isdir(path):
            data = os.path.join(path, "fastqc_data.txt")
            count = _search(TOTAL_SEQUENCES_TEXT, data) if os.path.isfile(data) else None
        if count is not None:
            total = count if total is None else total + count
    return total


def _fastp_totals(path: str) -> Tuple[Optional[int], Optional[int]]:
    """fastp's before/after read counts, which it writes into one JSON."""
    try:
        with open(path) as fh:
            loaded = json.load(fh)
    except (OSError, ValueError):
        return None, None
    summary = loaded.get("summary") or {} if isinstance(loaded, dict) else {}
    before = (summary.get("before_filtering") or {}).get("total_reads")
    after = (summary.get("after_filtering") or {}).get("total_reads")
    return (before if isinstance(before, int) else None, after if isinstance(after, int) else None)


def read_qc_counts(run_dirs: Sequence[str]) -> List[ReadCounts]:
    """Reads in and out per sample, from whichever report ``read_qc`` left behind.

    ``read_qc`` writes one output directory per sample, named for the sample, holding FastQC
    reports before and after cleaning (and ``fastp.json`` when fastp did the trimming). The sample
    name is that directory's, which is where ``read_qc`` put it - not re-derived from a FASTQ
    filename, which is the step that goes wrong on ``_1`` versus ``_R1`` naming.
    """
    counts: List[ReadCounts] = []
    for directory in _dirs(run_dirs):
        fastp = os.path.join(directory, FASTP_JSON)
        pre, post = os.path.join(directory, PRE_QC_DIR), os.path.join(directory, POST_QC_DIR)
        if os.path.isfile(fastp):
            reads_in, reads_out = _fastp_totals(fastp)
        elif os.path.isdir(pre) or os.path.isdir(post):
            reads_in, reads_out = _fastqc_total(pre), _fastqc_total(post)
        else:
            continue
        counts.append(ReadCounts(os.path.basename(directory), reads_in, reads_out))
    return counts


def assembly_stats(run_dirs: Sequence[str]) -> List[Tuple[str, Dict[str, str]]]:
    """QUAST's metrics per assembly, as ``(assembly name, metric -> value)``.

    Read from ``QUAST_out/report.tsv`` rather than measured off the FASTA. ``assembly`` already
    runs QUAST, and the numbers on this page have to be the numbers in that report: a second N50
    implementation with a different definition of the contig set would be a page-long disagreement
    with the file sitting next to it.
    """
    stats: List[Tuple[str, Dict[str, str]]] = []
    for path in _with(run_dirs, "QUAST_out", "report.tsv"):
        header, rows = _table(path)
        for index, name in enumerate(header[1:]):
            metrics = {
                row[0].strip(): row[index + 1].strip()
                for row in rows
                if row[0].strip() in QUAST_METRICS and index + 1 < len(row)
            }
            if metrics:
                stats.append((name, metrics))
    return stats


def krona_inputs(run_dirs: Sequence[str]) -> List[str]:
    """Every Krona input file ``pathways`` wrote, sorted. One wheel per file, by contract."""
    found: List[str] = []
    for directory in _with(run_dirs, KRONA_DIR):
        found.extend(
            os.path.join(directory, name)
            for name in sorted(os.listdir(directory))
            if name.endswith(".txt")
        )
    return found


def _wheel(krona_txt: str, work_dir: str) -> Optional[str]:
    """The wheel PNG for a Krona input, reusing the one beside it when it is current.

    ``pathways`` already draws these into the Krona directory, and redrawing is not free: the
    reference run's hierarchy directory holds 18 input files at about 1.3 s a wheel, so a report
    over it would spend 24 s re-rendering pictures that are already on disk. A PNG older than its
    input is redrawn, because a stale wheel is worse than a slow report.
    """
    beside = os.path.splitext(krona_txt)[0] + ".png"
    if os.path.isfile(beside) and os.path.getmtime(beside) >= os.path.getmtime(krona_txt):
        return beside
    name = os.path.splitext(os.path.basename(krona_txt))[0] + ".png"
    return _sunburst.draw_from_krona_text(krona_txt, os.path.join(work_dir, name))


def _run_records(run_dirs: Sequence[str]) -> List[Tuple[str, Dict[str, object]]]:
    """Each ``run_config.json`` under *run_dirs* as ``(directory, record)``."""
    records: List[Tuple[str, Dict[str, object]]] = []
    for path in _with(run_dirs, CONFIG_NAME):
        try:
            with open(path) as fh:
                loaded = json.load(fh)
        except (OSError, ValueError):
            # A half-written run_config.json is not worth failing a report over - the directory
            # still appears on the summary page, without its provenance.
            continue
        if isinstance(loaded, dict):
            records.append((os.path.dirname(path), loaded))
    return records


# ── laying out a page ────────────────────────────────────────────────────────────────────────


def _truncate(text: str, width: int) -> str:
    return text if len(text) <= width else text[: width - 1] + ELLIPSIS


def _cell(row: Sequence[str], index: int) -> str:
    return str(row[index]) if index < len(row) else ""


def _unscored_note(table: _magtable.Table) -> str:
    """How many rows carry no completeness, and where those genomes came from.

    Without this the page is ambiguous in a way that matters. The table is the union of its inputs,
    so a genome that appears only in `quant_bins`' abundance table - a different, usually larger bin
    set than the refined one - gets a row with every quality cell blank. On the reference run that
    is 13 of 26 rows, and a reader seeing half the table empty will reasonably suspect the report
    rather than conclude that two inputs describe two different bin sets.
    """
    if "completeness" not in table.columns:
        return ""
    column = table.columns.index("completeness")
    unscored = sum(1 for row in table.rows if not str(row[column]).strip())
    if not unscored:
        return ""
    return (
        " %d of %d genome(s) have no quality row in this bin set's stats table and are listed "
        "because another input - usually an abundance table over a different bin set - names them."
        % (unscored, len(table.rows))
    )


def _table_lines(columns: Sequence[str], rows: Sequence[Sequence[str]]) -> List[str]:
    """A monospaced table: header, rule, then one line per row, every column its own width."""
    widths = [
        min(MAX_CELL, max([len(str(columns[i]))] + [len(_cell(row, i)) for row in rows]))
        for i in range(len(columns))
    ]
    header = "  ".join(_truncate(str(name), w).ljust(w) for name, w in zip(columns, widths))
    lines = [header, "-" * len(header)]
    for row in rows:
        lines.append("  ".join(_truncate(_cell(row, i), w).ljust(w) for i, w in enumerate(widths)))
    return lines


def _font_for(lines: Sequence[str]) -> float:
    """Shrink the body font until the widest line fits the page, down to :data:`MIN_FONT`.

    Shrunk rather than wrapped or clipped: a wrapped table row stops being a row, and a clipped one
    silently loses its last columns - which for the MAG table are the per-sample abundances, the
    columns a reader came for.
    """
    widest = max((len(line) for line in lines), default=0)
    if widest <= PAGE_CHARS:
        return BODY_FONT
    return max(MIN_FONT, BODY_FONT * PAGE_CHARS / widest)


class _Pdf:
    """The document being written: its pages, and what was not drawn."""

    def __init__(self, pdf: PdfPages, run_dirs: Sequence[str]) -> None:
        self._pdf = pdf
        self.run_dirs = list(run_dirs)
        self.pages = 0
        self.absent: List[Absence] = []

    def missing(self, module: str, looked_for: str) -> None:
        """Record a page that was not drawn. Repeats are dropped, not listed twice.

        Two pages can depend on one file - the MAG table and the abundance page both need
        ``quant_bins``' matrix - and that page is a list of what is missing, not a list of
        pages that wanted it.
        """
        entry = Absence(module, looked_for)
        if entry not in self.absent:
            self.absent.append(entry)

    def add(self, figure: Figure) -> None:
        import matplotlib.pyplot as plt

        self._pdf.savefig(figure)
        plt.close(figure)
        self.pages += 1

    def page(
        self,
        title: str,
        lines: Sequence[str],
        footer: str = "",
        chart: Optional[Callable[[Axes], None]] = None,
    ) -> None:
        """One page: a title, a monospaced body, an optional chart above it, an optional footer."""
        import matplotlib.pyplot as plt

        figure = plt.figure(figsize=PAGE_SIZE)
        figure.text(0.04, 0.955, title, fontsize=13, weight="bold", va="top")
        top = 0.90
        if chart is not None:
            axes = figure.add_axes((0.06, 0.62, 0.40, 0.26), zorder=ARTIST_OVER_TEXT)
            chart(axes)
            top = 0.56
        figure.text(
            0.04,
            top,
            "\n".join(lines),
            fontsize=_font_for(lines),
            family="monospace",
            va="top",
            linespacing=1.5,
        )
        if footer:
            figure.text(0.04, 0.025, footer, fontsize=6.5, color="#4a5058", va="bottom", wrap=True)
        self.add(figure)

    def paginate(
        self,
        title: str,
        head: Sequence[str],
        body: Sequence[str],
        footer: str = "",
        chart: Optional[Callable[[Axes], None]] = None,
    ) -> None:
        """*head* repeated on every page, *body* split over as many pages as it needs.

        Paginated rather than scaled to fit, because the MAG table is the point of the document and
        a 300-genome catalogue shrunk onto one page is a grey rectangle. At least one page is always
        written, so an empty body still produces the header and the footer that explains it.
        """
        remaining = list(body)
        first = True
        while True:
            room = (PAGE_LINES_WITH_CHART if (first and chart) else PAGE_LINES) - len(head)
            chunk, remaining = remaining[: max(1, room)], remaining[max(1, room) :]
            self.page(
                title if first else "%s (continued)" % title,
                list(head) + chunk,
                footer,
                chart if first else None,
            )
            first = False
            if not remaining:
                return

    def image_page(self, title: str, images: Sequence[Tuple[str, str]], caption: str = "") -> None:
        """A page of one or two existing PNGs, each scaled to fit its half and not cropped.

        Two per page because the two wheels of one unit are the same distribution with and without
        BRITE's scaffolding branches, and the only thing anyone does with that pair is compare them
        - which is impossible with one of them on the previous page. It also halves the page count:
        the reference run's hierarchy directory holds 18 wheels, which was 18 of a 31-page report.
        """
        import matplotlib.image as mpimg
        import matplotlib.pyplot as plt

        figure = plt.figure(figsize=PAGE_SIZE)
        figure.text(0.04, 0.965, title, fontsize=12, weight="bold", va="top")
        width = 0.92 / max(1, len(images))
        for index, (png, label) in enumerate(images):
            left = 0.04 + index * width
            axes = figure.add_axes((left, 0.07, width, 0.82), zorder=ARTIST_OVER_TEXT + index)
            axes.imshow(mpimg.imread(png))
            axes.set_axis_off()
            if label:
                figure.text(left + width / 2, 0.925, label, fontsize=9, ha="center", va="top")
        if caption:
            figure.text(0.04, 0.025, caption, fontsize=6.5, color="#4a5058", va="bottom")
        self.add(figure)


def _tier_chart(counts: Sequence[Tuple[str, int]]) -> Callable[[Axes], None]:
    """A bar per quality tier, including the empty tiers and the unmeasured genomes."""

    def draw(axes: Axes) -> None:
        labels = [label for label, _n in counts]
        values = [count for _label, count in counts]
        axes.bar(range(len(values)), values, color="#1f6fb4", width=0.68)
        axes.set_xticks(range(len(labels)))
        axes.set_xticklabels([label.replace("_", "\n") for label in labels], fontsize=6)
        axes.set_ylabel("genomes", fontsize=7)
        axes.tick_params(axis="y", labelsize=6)
        axes.set_title("quality tier distribution", fontsize=8)
        for spine in ("top", "right"):
            axes.spines[spine].set_visible(False)
        for index, value in enumerate(values):
            axes.text(index, value, str(value), ha="center", va="bottom", fontsize=6)

    return draw


# ── the pages ────────────────────────────────────────────────────────────────────────────────


def _title_page(doc: _Pdf, title: str, records: Sequence[Tuple[str, Dict[str, object]]]) -> None:
    """One page: what was run, when, which MetaWrap2, and where. Nothing a reader has to scroll past.

    One page, and one line per module. The full per-run dumps - resolved parameters, tool versions,
    databases - go to :func:`_provenance_appendix` at the end instead. They used to be here, and on
    the reference run they filled **six pages before the first result**: somebody opening the report
    to see what the run found read six pages of configuration first. Provenance has to be in the
    document, but it is reference material, and reference material goes at the back.
    """
    body = ["report written by MetaWrap2 %s" % __version__, "", "run directories"]
    for directory in doc.run_dirs:
        body.append("  %s%s" % (directory, "" if os.path.isdir(directory) else "   (not found)"))
    body.append("")

    if not records:
        doc.missing("provenance", "%s under the run directories" % CONFIG_NAME)
        doc.paginate(title, [], body, footer=_NO_CONFIG_FOOTER)
        return

    rows = [
        [
            str(record.get("module") or "?"),
            str(record.get("started") or ""),
            str(record.get("finished") or ""),
            _duration(record.get("duration_seconds")),
            str(record.get("metawrap2_version") or ""),
            os.path.basename(directory.rstrip(os.sep)),
        ]
        for directory, record in records
    ]
    body.append("what was run")
    body.extend(
        "  " + line
        for line in _table_lines(
            ("module", "started", "finished", "duration", "MetaWrap2", "directory"), rows
        )
    )
    doc.paginate(
        title,
        [],
        body,
        footer=(
            "Resolved parameters, tool versions and databases for each of these runs are in the "
            "appendix at the end of this report, and in full in each directory's %s and %s."
            % (CONFIG_NAME, ENVIRONMENT_NAME)
        ),
    )


def _duration(seconds: object) -> str:
    """A run's wall time, or blank when it did not record one."""
    if not isinstance(seconds, (int, float)):
        return _magtable.BLANK
    return "%.0f s" % seconds if seconds < 600 else "%.1f h" % (seconds / 3600.0)


#: What a run calls itself, in the order a reader wants it. Read out of ``run_config.json`` rather
#: than reconstructed, so the dates in the appendix are the dates the run recorded for itself.
_RUN_IDENTITY = (
    ("module", "module"),
    ("run id", "run_id"),
    ("MetaWrap2", "metawrap2_version"),
    ("started", "started"),
    ("finished", "finished"),
    ("conda env", "conda_env"),
)

_NO_CONFIG_FOOTER = (
    "No %s found under these directories, so resolved parameters, tool versions and databases are "
    "unavailable; any %s found is quoted in the appendix." % (CONFIG_NAME, PROVENANCE_NAME)
)


def _provenance_appendix(doc: _Pdf, records: Sequence[Tuple[str, Dict[str, object]]]) -> None:
    """The full per-run provenance, after the results and after "what was not run".

    Deliberately last. A reader who wants the exact thread count or the GUNC database build will go
    looking for it and find it; a reader who wants to know what the run recovered never has to read
    past it. Nothing is dropped relative to what the title page used to show.
    """
    body: List[str] = []
    for directory, record in records:
        body.append(os.path.basename(directory.rstrip(os.sep)))
        for label, key in _RUN_IDENTITY:
            value = record.get(key)
            if value:
                body.append("  %-12s %s" % (label, value))
        duration = record.get("duration_seconds")
        if isinstance(duration, (int, float)):
            body.append("  %-12s %.0f s" % ("duration", duration))
        body.extend(_provenance_lines(directory, record))
        body.append("")
    if not records:
        body = _provenance_text(doc.run_dirs)
    if not body:
        return
    doc.paginate(
        "Appendix: provenance",
        [],
        body,
        footer="Read back from each run's %s; nothing here is re-derived." % CONFIG_NAME,
    )


def _provenance_text(run_dirs: Sequence[str], limit: int = 30) -> List[str]:
    """The first *limit* lines of any ``provenance.txt``, for a run with no machine-readable record.

    Worth quoting rather than ignoring: it is the human-readable history of which steps produced
    this output, which is most of what the summary page is for, and an appendix that could only
    show a path would be worse than one that shows a step list it cannot parse.
    """
    lines: List[str] = []
    for path in _with(run_dirs, PROVENANCE_NAME):
        with open(path, errors="replace") as fh:
            quoted = fh.read().splitlines()[:limit]
        lines.append("%s" % path)
        lines.extend("  %s" % line.rstrip() for line in quoted)
        lines.append("")
    return lines


def _provenance_lines(directory: str, record: Dict[str, object]) -> List[str]:
    """One module's provenance block: resolved parameters, tool versions, databases."""
    lines: List[str] = []
    for label, key in (
        ("parameters", "parameters"),
        ("tool versions", "tool_versions"),
        ("databases", "databases"),
    ):
        block = record.get(key)
        if not isinstance(block, dict) or not block:
            continue
        lines.append("  %s" % label)
        lines.extend("    %-22s %s" % (name, block[name]) for name in sorted(block))
    if os.path.isfile(os.path.join(directory, ENVIRONMENT_NAME)):
        lines.append("  %s recorded (every conda package and build)" % ENVIRONMENT_NAME)
    return lines


def _read_qc_page(doc: _Pdf) -> None:
    counts = read_qc_counts(doc.run_dirs)
    if not counts:
        doc.missing("read_qc", "%s / %s / %s per sample" % (PRE_QC_DIR, POST_QC_DIR, FASTP_JSON))
        return
    rows = [[entry.sample, _count(entry.reads_in), _count(entry.reads_out)] for entry in counts]
    lines = _table_lines(("sample", "reads in", "reads out"), rows)
    doc.paginate(
        "Read QC",
        lines[:2],
        lines[2:],
        footer=(
            "Reads summed over both mates of a pair, as FastQC and fastp report them. A blank is a "
            "count no report stated, which is not a count of zero."
        ),
    )


def _count(value: Optional[int]) -> str:
    """A read count, or :data:`metawrap2.magtable.BLANK` where none was stated."""
    return _magtable.BLANK if value is None else "%d" % value


def _assembly_page(doc: _Pdf) -> None:
    stats = assembly_stats(doc.run_dirs)
    if not stats:
        doc.missing("assembly", os.path.join("QUAST_out", "report.tsv"))
        return
    rows = [
        [name] + [metrics.get(key, _magtable.BLANK) for key in QUAST_METRICS]
        for name, metrics in stats
    ]
    lines = _table_lines(("assembly",) + QUAST_METRICS, rows)
    doc.paginate(
        "Assembly",
        lines[:2],
        lines[2:],
        footer="Copied from QUAST's own report.tsv, not re-measured from the FASTA.",
    )


#: Which module writes each of the MAG table's optional inputs, so "what was not run" can name the
#: module a reader has to re-run rather than the file they would have to produce by hand.
_MAG_INPUTS = {
    "taxonomy": ("classify_bins", _magtable.TAXONOMY_NAME),
    "abundance": ("quant_bins", _magtable.ABUNDANCE_NAME),
    "gunc": ("bin_refinement --gunc", "GUNC.*.maxCSS_level.tsv"),
}


def _unjoined_note(table: _magtable.Table) -> str:
    """Which inputs exist but matched none of a bin set's genomes, and why those cells are blank.

    The page would otherwise be ambiguous in the opposite direction from :func:`_unscored_note`: a
    study-wide table that is present and simply does not describe these genomes leaves exactly the
    same blank cells as one that was never written, and only this sentence tells them apart.

    It is the normal state of a cohort run. ``classify_bins`` and ``quant_bins`` work on the
    *dereplicated* catalogue, whose genome names need not resemble any one sample's ``bin_001``, so
    there is no join to make and nothing to put in the cell but a guess.
    """
    if not table.unjoined:
        return ""
    by_kind: Dict[str, List[str]] = {}
    for kind, sample in table.unjoined:
        by_kind.setdefault(kind, []).append(sample)
    return (
        " Present, but matching no genome in the bin set named, so blank there rather than guessed:"
        " %s - a cohort's taxonomy and abundance are keyed on dereplicated genome names, which need"
        " not match a per-sample bin name."
        % "; ".join(
            "%s for %s" % (kind, ", ".join(samples)) for kind, samples in sorted(by_kind.items())
        )
    )


def _genomes_heading(table: _magtable.Table) -> str:
    """``Genomes (62 across 3 bin sets)``: a heading that states what it covers.

    The count alone was a trap on a per-sample study. The page said "Genomes (24)" for a study that
    recovered 62 across three samples, and a reader had nothing to tell them that the number was
    one sample's.
    """
    if table.bin_sets > 1:
        return "Genomes (%d across %d bin sets)" % (len(table.rows), table.bin_sets)
    return "Genomes (%d)" % len(table.rows)


def _mag_pages(doc: _Pdf, mag_table_path: Optional[str]) -> Optional[str]:
    """The per-MAG table, paginated, with the quality tier distribution above the first page."""
    inputs = _magtable.discover(doc.run_dirs)
    table = _magtable.build(inputs)
    for kind, (module, filename) in _MAG_INPUTS.items():
        if kind in table.absent:
            doc.missing(module, filename)
    if not table.rows:
        doc.missing("bin_refinement / binning", "a *%s table with genomes in it" % STATS_SUFFIX)
        return None

    written = _magtable.write(table, mag_table_path) if mag_table_path else None
    lines = _table_lines(table.columns, table.rows)
    doc.paginate(
        _genomes_heading(table),
        lines[:2],
        lines[2:],
        footer=(
            "A blank cell is a measurement that was not made, never a zero. Inputs with no file at "
            "all, hence blank for every genome: %s.%s%s Cells are truncated to %d characters "
            "here%s."
            % (
                ", ".join(sorted(table.absent)) or "none",
                _unscored_note(table),
                _unjoined_note(table),
                MAX_CELL,
                "; the full table is %s" % written if written else "",
            )
        ),
        chart=_tier_chart(_magtable.tier_counts(table)),
    )
    return written


def _abundance_pages(doc: _Pdf) -> None:
    """``quant_bins``' matrix as it was written, plus the heatmap that module already drew."""
    path = _magtable.discover(doc.run_dirs).abundance
    if path is None:
        return  # named by _mag_pages, which looks for the same file
    samples, values = _magtable.read_abundance(path)
    rows = [
        [genome] + [values[genome].get(sample, _magtable.BLANK) for sample in samples]
        for genome in values
    ]
    lines = _table_lines(["genome"] + samples, rows)
    doc.paginate(
        "Abundance across %d sample(s)" % len(samples),
        lines[:2],
        lines[2:],
        footer=(
            "Copied from %s. The unit is whichever quantifier ran - coverage for salmon and the "
            "depth route, relative abundance for sylph - which that run's %s records."
            % (os.path.basename(path), CONFIG_NAME)
        ),
    )
    heatmap = os.path.join(os.path.dirname(path), ABUNDANCE_HEATMAP)
    if os.path.isfile(heatmap):
        doc.image_page(
            "Abundance heatmap",
            [(heatmap, "")],
            caption="Drawn by quant_bins: %s" % _relative(heatmap, doc.run_dirs),
        )


def _krona_unit(krona_txt: str) -> Tuple[str, str]:
    """``genome.bin_003.pathways.txt`` -> ``("genome.bin_003", "pathways")``.

    ``pathways`` writes two Krona inputs per unit - one with BRITE's scaffolding branches and one
    without - distinguished by that last dotted component. Split rather than matched against a list
    of known variants, so a third variant lands on a page of its own instead of disappearing.
    """
    stem = os.path.splitext(os.path.basename(krona_txt))[0]
    unit, _dot, variant = stem.rpartition(".")
    return (unit, variant) if unit else (stem, "")


def _function_pages(doc: _Pdf, work_dir: str) -> None:
    """The functional wheels, two to a page, drawn by :mod:`metawrap2.sunburst`.

    Grouped by unit, so a genome's scaffolding-included and pathway-only wheels sit side by side -
    the only comparison anybody makes between them.
    """
    inputs = krona_inputs(doc.run_dirs)
    if not inputs:
        doc.missing("pathways / annotate", os.path.join(KRONA_DIR, "*.txt"))
        return
    drawn: Dict[str, List[Tuple[str, str, str]]] = {}
    for krona_txt in inputs:
        unit, variant = _krona_unit(krona_txt)
        png = _wheel(krona_txt, work_dir)
        if png is None:
            # An empty distribution is a real result - a genome with no KO at all - and sunburst
            # returns None rather than drawing a blank wheel. Named on "what was not run" rather
            # skipped in silence, because a missing wheel otherwise reads as a missing genome.
            doc.missing("pathways", "a non-empty distribution in %s" % os.path.basename(krona_txt))
            continue
        drawn.setdefault(unit, []).append((png, variant, krona_txt))
    for unit, wheels in drawn.items():
        for start in range(0, len(wheels), WHEELS_PER_PAGE):
            pair = wheels[start : start + WHEELS_PER_PAGE]
            doc.image_page(
                "Function: %s" % unit,
                [(png, variant) for png, variant, _source in pair],
                caption="Drawn from %s, the same file(s) Krona is given."
                % ", ".join(_relative(source, doc.run_dirs) for _png, _v, source in pair),
            )


def _summary_page(
    doc: _Pdf,
    title: str,
    module: str,
    tables: Sequence[Tuple[str, str, str]],
    footer: str,
) -> None:
    """A page of row counts and per-column breakdowns for one module's tables.

    *tables* is ``(filename, what one row is, the column to break down by)``, with an empty column
    meaning "count the rows and nothing else". Row counts only: how many resistance gene calls, how
    many viral contigs in each tier. Anything beyond counting a table's own rows belongs in the
    module that wrote the table.
    """
    lines: List[str] = []
    for filename, unit, column in tables:
        for path in _with(doc.run_dirs, filename):
            lines.append(path)
            lines.append("  %-30s %d" % (unit, _row_count(path)))
            for value, count in _column_counts(path, column)[:12] if column else ():
                lines.append("    %-28s %d" % (_truncate(value, 28), count))
            lines.append("")
    if not lines:
        doc.missing(module, " / ".join(name for name, _unit, _column in tables))
        return
    doc.paginate(title, [], lines, footer=footer)


def _absence_page(doc: _Pdf) -> None:
    """The page that is not optional: every module whose output was absent, named.

    Drawn even when nothing is missing, because its own absence would be ambiguous - a reader could
    not tell whether the report had checked.
    """
    if doc.absent:
        # Laid out by hand rather than through _table_lines, which truncates a cell at MAX_CELL.
        # Truncation is right for a lineage and wrong here: "viral_contigs.tsv / plasmid_conti…"
        # does not tell a reader which file to go and look for.
        lines = ["%-26s %s" % ("module / page", "looked for"), ""]
        lines += ["%-26s %s" % (entry.module, entry.looked_for) for entry in doc.absent]
    else:
        lines = ["Every module this report covers produced output."]
    lines += ["", "searched, at most %d directories deep:" % MAX_DEPTH]
    lines += ["  %s" % directory for directory in doc.run_dirs]
    doc.paginate("What was not run", [], lines, footer=ABSENCE_NOTE)


# ── the document ─────────────────────────────────────────────────────────────────────────────

#: ``amr``, ``viruses`` and ``strains`` each write a few flat tables and nothing else worth a plot,
#: so each gets one page of row counts. Declared as data so adding a fourth module is a row here.
_SUMMARIES = (
    (
        "Resistance and virulence",
        "amr",
        (("amr.tsv", "resistance gene calls", "class"), ("virulence.tsv", "virulence calls", "")),
        "Row counts of amr's own tables, broken down by their `class` column.",
    ),
    (
        "Viruses",
        "viruses",
        (
            ("viral_contigs.tsv", "viral contigs", "tier"),
            ("plasmid_contigs.tsv", "plasmid contigs", ""),
            ("votus.tsv", "vOTUs", ""),
        ),
        "Tiers are MIUViG's, as viruses wrote them - not re-derived here.",
    ),
    (
        "Strains",
        "strains",
        (
            ("populations.tsv", "genome x sample profiles", "genome"),
            ("strain_comparisons.tsv", "pairwise comparisons", ""),
        ),
        "Row counts of strains' own tables; that run's %s holds the thresholds." % CONFIG_NAME,
    ),
)


def build_report(
    run_dirs: Sequence[str],
    out_path: str,
    *,
    title: str = "",
    mag_table_path: Optional[str] = None,
) -> Report:
    """Write one PDF for *run_dirs* and return what was written and what was absent.

    Pages are added in a fixed order whatever the run contains, so two reports of two runs can be
    read the same way, and a page whose inputs are missing is recorded rather than quietly dropped.
    """
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.backends.backend_pdf import PdfPages

    directory = os.path.dirname(os.path.abspath(out_path))
    if directory:
        os.makedirs(directory, exist_ok=True)
    heading = title or "MetaWrap2 report"

    with (
        tempfile.TemporaryDirectory(prefix="metawrap2-report-") as work_dir,
        PdfPages(out_path) as pdf,
    ):
        doc = _Pdf(pdf, run_dirs)
        pdf.infodict()["Title"] = heading
        pdf.infodict()["Creator"] = "MetaWrap2 %s" % __version__

        # Read once and passed to both the summary page and the appendix, so the one-liner at the
        # front and the full dump at the back cannot describe different sets of runs.
        records = _run_records(run_dirs)

        _title_page(doc, heading, records)
        _read_qc_page(doc)
        _assembly_page(doc)
        written = _mag_pages(doc, mag_table_path)
        _abundance_pages(doc)
        _function_pages(doc, work_dir)
        for summary_title, module, tables, footer in _SUMMARIES:
            _summary_page(doc, summary_title, module, tables, footer)
        _absence_page(doc)
        _provenance_appendix(doc, records)

    return Report(path=out_path, pages=doc.pages, absent=doc.absent, mag_table=written)
