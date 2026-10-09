"""GUNC helpers: chimerism detection for refined bins.

CheckM answers "does this bin have one copy of every marker gene it should have". A bin assembled
from two halves of two different genomes answers yes - it has a full, uncontaminated-looking marker
set, and marker-gene scoring cannot see that the halves are unrelated. GUNC can: it maps the bin's
genes against a reference taxonomy and asks whether they agree about where the bin sits. That makes
it a complement to CheckM rather than a replacement, which is why both scores end up in the same
``.stats.tsv`` instead of one overruling the other.

What this module does is run GUNC and merge its verdict into that table. It deliberately does not
drop bins: a chimera is worth knowing about, and a run that silently deleted bins on a second
tool's opinion would be hard to audit. Filtering is the caller's decision, and the column is what
makes it possible.
"""

from __future__ import annotations

import csv
import os
from typing import Dict, List, NamedTuple, Optional

from .command import run, runner
from .constants import BIN_EXTENSION

#: GUNC's clade-separation-score threshold. Above it, the genes in a bin disagree about its taxonomy
#: by more than chance explains, and the bin is called chimeric. 0.45 is GUNC's own published
#: cut-off (Orakov et al., Genome Biology 2021;22:178), kept rather than re-tuned so that
#: "chimeric" here means what it means in the literature.
CSS_THRESHOLD = 0.45

#: The report GUNC writes: one row per bin, at whichever taxonomic level scored highest. Matched by
#: prefix and suffix because the middle carries the reference database version
#: (``GUNC.progenomes_2.1.maxCSS_level.tsv``), and a hardcoded version would find nothing the first
#: time that database is updated.
REPORT_PREFIX = "GUNC."
REPORT_SUFFIX = ".maxCSS_level.tsv"

#: ``--file_suffix`` carries the leading dot on purpose. GUNC derives each genome's name with
#: ``basename.split(file_suffix)[0]``, so a suffix of "fasta" yields "bin_001." - a name that matches
#: nothing in the stats table, which would silently annotate zero bins. ``--temp_dir`` is set
#: because GUNC defaults it to the working directory, and diamond's temporaries do not belong there.
GUNC_RUN = (
    "gunc run --input_dir {bins} --file_suffix {ext} --db_file {db} "
    "--out_dir {out} --temp_dir {out} --threads {threads}"
)

#: Said whenever GUNC cannot run, because both ways out are worth naming: get the database, or
#: decide the check is not wanted on this run.
DB_HINT = (
    "GUNC needs its reference database; install it with 'metawrap2 install-db gunc' and point "
    "GUNC_DB at the .dmnd file, or pass --skip-gunc."
)

__all__ = [
    "CSS_THRESHOLD",
    "DB_HINT",
    "Verdict",
    "annotate_stats",
    "find_report",
    "parse_report",
    "run_gunc",
]


class Verdict(NamedTuple):
    """GUNC's answer for one bin."""

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


def parse_report(path: str) -> Dict[str, Verdict]:
    """GUNC's verdict per bin, keyed by bin name.

    ``pass.GUNC`` is used where GUNC wrote it, because that column is its whole judgement: clade
    separation *and* whether the bin was represented well enough in the reference for the score to
    mean anything. Where the column is absent the threshold stands in, which agrees with GUNC for a
    well-represented bin and is more forgiving for a poorly represented one.
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
                # A bin GUNC could not score, or a header repeated mid-file. Skipping it leaves the
                # bin unannotated, which reads as "not assessed" rather than as "passed".
                continue
            raw_pass = record.get(columns.get("pass.gunc", ""), "")
            if raw_pass not in (None, ""):
                passed = str(raw_pass).strip().lower() in ("true", "1", "yes", "pass")
            else:
                passed = css <= CSS_THRESHOLD
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

    A bin GUNC did not score gets empty cells rather than a default. Writing an unassessed bin as
    passing would be the exact mistake this step exists to catch, and writing it as failing would
    condemn a bin over a measurement that was never made.
    """
    if not os.path.isfile(stats_path):
        return 0
    with open(stats_path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    if not rows:
        return 0

    header = rows[0]
    if "gunc_css" in header:
        # Already annotated - a --resume that re-runs this step would otherwise append the two
        # columns a second time and leave a table nothing can parse.
        at = {name: i for i, name in enumerate(header)}
        css_at, pass_at = at["gunc_css"], at["gunc_pass"]
    else:
        css_at, pass_at = len(header), len(header) + 1
        header = header + ["gunc_css", "gunc_pass"]

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
