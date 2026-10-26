"""GUNC helpers: chimerism detection for refined bins.

CheckM scores marker-gene completeness and redundancy. A bin made of halves of two different
genomes can still show a full, uncontaminated-looking marker set. GUNC maps the bin's genes against
a reference taxonomy and checks whether they agree on where the bin sits, so it complements CheckM
and both scores go in the same ``.stats.tsv``.

This module runs GUNC and merges its result into that table. It does not drop bins; filtering on
the GUNC columns is left to the caller.
"""

from __future__ import annotations

import csv
import os
from typing import Dict, List, NamedTuple, Optional

from .command import run, runner
from .constants import BIN_EXTENSION

#: GUNC's clade-separation-score threshold; above it a bin is called chimeric. 0.45 is GUNC's
#: published cut-off (Orakov et al., Genome Biology 2021;22:178), kept for comparability.
CSS_THRESHOLD = 0.45

#: The report GUNC writes: one row per bin, at the taxonomic level that scored highest. Matched by
#: prefix and suffix because the middle carries the reference database version
#: (``GUNC.progenomes_2.1.maxCSS_level.tsv``).
REPORT_PREFIX = "GUNC."
REPORT_SUFFIX = ".maxCSS_level.tsv"

#: ``--file_suffix`` includes the leading dot. GUNC derives each genome's name with
#: ``basename.split(file_suffix)[0]``, so a suffix of "fasta" yields "bin_001.", which matches no
#: row in the stats table. ``--temp_dir`` is set because GUNC otherwise puts diamond's temporary
#: files in the working directory.
GUNC_RUN = (
    "gunc run --input_dir {bins} --file_suffix {ext} --db_file {db} "
    "--out_dir {out} --temp_dir {out} --threads {threads}"
)

#: Shown whenever GUNC cannot run: install the database, or skip the check.
DB_HINT = (
    "GUNC needs its reference database; install it with 'metawrap2 install-db gunc' and point "
    "GUNC_DB at the .dmnd file, or pass --skip-gunc."
)

__all__ = [
    "CSS_THRESHOLD",
    "DB_HINT",
    "STATS_COLUMNS",
    "Verdict",
    "annotate_stats",
    "find_report",
    "parse_report",
    "run_gunc",
]

#: The columns :func:`annotate_stats` appends to a ``.stats.tsv``. Exported, like
#: ``checkm.STATS_HEADER`` and ``mimag.STATS_COLUMNS``, for code that builds the full table header
#: (`metawrap2 score`).
STATS_COLUMNS: tuple = ("gunc_css", "gunc_pass")


class Verdict(NamedTuple):
    """GUNC's result for one bin."""

    css: float
    passed: bool
    level: str


def find_report(out_dir: str) -> Optional[str]:
    """GUNC's maxCSS table under *out_dir*, or None if it is not there."""
    for base, _dirs, names in os.walk(out_dir):
        for name in sorted(names):
            if name.startswith(REPORT_PREFIX) and name.endswith(REPORT_SUFFIX):
                return os.path.join(base, name)
    return None


def passed_value(raw: object) -> Optional[bool]:
    """GUNC's ``pass.GUNC`` cell as a verdict, or None when it does not state one.

    GUNC writes this column as ``True``/``False`` when pandas kept it boolean, and as
    ``1.0``/``0.0`` when a ``nan`` forced a float dtype. The same GUNC version produced both in one
    benchmark run. Accepting only the words gave a 166-bin set written as floats a GUNC pass
    fraction of 0.

    ``nan`` is None, not False: GUNC writes it when a bin has too little reference representation
    to be assessed.
    """
    if raw is None:
        return None
    text = str(raw).strip().lower()
    if text in ("", "nan", "na", "none"):
        return None
    if text in ("true", "yes", "pass", "passed"):
        return True
    if text in ("false", "no", "fail", "failed"):
        return False
    try:
        return float(text) != 0.0
    except ValueError:
        return None


def parse_report(path: str) -> Dict[str, Verdict]:
    """GUNC's verdict per bin, keyed by bin name.

    ``pass.GUNC`` is used where present, because it combines clade separation with whether the bin
    is represented well enough in the reference. Without it, :data:`CSS_THRESHOLD` is applied,
    which agrees with GUNC for well-represented bins and is more lenient for poorly represented
    ones.
    """
    verdicts: Dict[str, Verdict] = {}
    with open(path, newline="") as fh:
        for record in csv.DictReader(fh, delimiter="\t"):
            columns = {key.strip().lower(): key for key in record if key}
            name = (record.get(columns.get("genome", ""), "") or "").strip()
            raw_css = record.get(columns.get("clade_separation_score", ""), "")
            if not name or raw_css in (None, ""):
                continue
            try:
                css = float(raw_css)
            except ValueError:
                # A bin GUNC could not score, or a header repeated mid-file. The bin is left
                # unannotated ("not assessed").
                continue
            stated = passed_value(record.get(columns.get("pass.gunc", ""), ""))
            passed = stated if stated is not None else css <= CSS_THRESHOLD
            level = (record.get(columns.get("taxonomic_level", ""), "") or "").strip()
            verdicts[name] = Verdict(css=css, passed=passed, level=level)
    return verdicts


def run_gunc(
    bins_dir: str,
    *,
    threads: int,
    db: str,
    env: Optional[str],
    out_dir: Optional[str] = None,
    log_path: Optional[str] = None,
) -> Optional[str]:
    """Run GUNC over *bins_dir* and return its maxCSS table, or None under ``--dry-run``."""
    out = out_dir or (bins_dir + ".gunc")
    if not runner.dry_run:
        os.makedirs(out, exist_ok=True)  # gunc does not create its own --out_dir
    run(
        GUNC_RUN.format(
            bins=bins_dir,
            ext=BIN_EXTENSION,
            db=db,
            out=out,
            threads=threads,
        ),
        env=env,
        tool="gunc",
        log_path=log_path,
        hint=DB_HINT,
    )
    if runner.dry_run:
        return None
    report = find_report(out)
    if report is None:
        raise RuntimeError(
            "GUNC did not produce a %s*%s report in %s" % (REPORT_PREFIX, REPORT_SUFFIX, out)
        )
    return report


def annotate_stats(stats_path: str, verdicts: Dict[str, Verdict]) -> int:
    """Add ``gunc_css`` and ``gunc_pass`` to a ``.stats.tsv``; return how many rows were annotated.

    A bin GUNC did not score gets empty cells, not a pass or a fail.
    """
    if not os.path.isfile(stats_path):
        return 0
    with open(stats_path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    if not rows:
        return 0

    css_column, pass_column = STATS_COLUMNS
    header = rows[0]
    if css_column in header:
        # Already annotated (a --resume re-ran this step): overwrite the columns instead of
        # appending them a second time.
        at = {name: i for i, name in enumerate(header)}
        css_at, pass_at = at[css_column], at[pass_column]
    else:
        css_at, pass_at = len(header), len(header) + 1
        header = header + list(STATS_COLUMNS)

    width = len(header)
    out: List[List[str]] = [header]
    annotated = 0
    for row in rows[1:]:
        row = row + [""] * (width - len(row))
        verdict = verdicts.get(row[0])
        if verdict is not None:
            row[css_at] = "%.4f" % verdict.css
            row[pass_at] = "yes" if verdict.passed else "no"
            annotated += 1
        out.append(row)
    with open(stats_path, "w") as fh:
        fh.writelines("\t".join(row) + "\n" for row in out)
    return annotated
