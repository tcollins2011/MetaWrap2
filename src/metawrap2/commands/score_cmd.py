"""``metawrap2 score -o OUT BIN_DIR [BIN_DIR ...]``: score any directory of genomes.

Runs the scoring steps of ``bin_refinement`` (CheckM for completeness and contamination, GUNC for
chimerism, MIMAG for the tier) on directories of genomes from any source, and writes the same
``<bin set>.stats.tsv``. Every bin set in a run is scored with the same tool, env and database.
CheckM1 and CheckM2 disagree most in the high-completeness tail, where the 90/5 near-complete line
sits, so bin sets scored by different estimators cannot be compared on near-complete counts.

    metawrap2 score -o SCORED refined/metawrap_70_10_bins
    metawrap2 score -o SCORED -x .fa metawrap1_run/metawrap_70_10_bins refined/metawrap_70_10_bins
    metawrap2 score -o SCORED --checkm1 published_mags/      # against CheckM1-scored literature

CheckM comes from :mod:`metawrap2.checkm`, GUNC from :mod:`metawrap2.gunc` and the MIMAG tier from
:mod:`metawrap2.mimag`; the table's columns are the ones those modules write (see
:func:`stats_columns`). This module runs them and prints the comparison table. It does not modify
or remove any bin.

Author: Gherman Uritskiy.
"""

from __future__ import annotations

import argparse
import os
import shutil
from collections import Counter
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

from .. import checkm as _checkm
from .. import gunc as _gunc
from .. import magtable as _magtable
from .. import mimag as _mimag
from ..benchmark.assess import CHECKM_MODULE
from ..benchmark.scoring import QUALITY_TIERS, bin_quality, quality_counts
from ..cliformat import make_parser, percent
from ..config import load_settings
from ..constants import (
    BIN_EXTENSION,
    DEFAULT_MAX_CONTAMINATION,
    DEFAULT_MIN_COMPLETION,
    STATS_SUFFIX,
)
from ..io.seqio import fasta_files
from ..modules import bin_refinement as _refinement_module
from ..modules._common import (
    absolutize_paths,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    memory_arg,
    optional_env_for,
    resolve_memory_gb,
    resolve_threads,
    start_run,
    threads_arg,
    warning,
)
from ..utils import strip_trailing_slash

# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: The two columns :func:`metawrap2.gunc.annotate_stats` appends, in its order. Imported from
#: :mod:`metawrap2.gunc`, like ``checkm.STATS_HEADER`` and ``mimag.STATS_COLUMNS``.
GUNC_STATS_COLUMNS: Tuple[str, ...] = _gunc.STATS_COLUMNS

#: Cache file for per-bin CheckM scores, with the estimator in the filename. Scores are keyed on a
#: hash of the bin's contigs (see :class:`metawrap2.checkm.BinScoreCache`), so a ``--resume`` that
#: toggles ``--checkm1`` must not read the other estimator's scores.
SCORE_CACHE = "checkm_scores.%s.json"

#: The near-complete tier (90/5) shown in the comparison table, and its thresholds. Taken from
#: :data:`metawrap2.benchmark.scoring.QUALITY_TIERS`, so the count matches ``docs/benchmarks.md``.
NEAR_COMPLETE_TIER, NEAR_COMPLETE_COMPLETENESS, NEAR_COMPLETE_CONTAMINATION = QUALITY_TIERS[0]

#: Printed under the comparison table.
COMPARISON_NOTE = (
    "All rows were scored with the same tool and database. CheckM1 and CheckM2 disagree most for "
    "near-complete bins, so compare pipelines only on scores from one tool. 'clean' means "
    "near-complete and not chimeric by GUNC, counted per genome."
)
# ───────────────────────────────────────────────────────────────────────────────────────────


def stats_columns(with_gunc: bool = True) -> List[str]:
    """The columns a scored bin set's table holds, in order.

    Built from the modules that write them, so the table matches the one ``bin_refinement``
    writes. The order is the order the passes run in: :func:`metawrap2.checkm.write_stats` writes
    the base header, :func:`metawrap2.gunc.annotate_stats` appends two columns and
    :func:`metawrap2.mimag.annotate_stats` appends three.

    *with_gunc* is False under ``--no-gunc``, which matches ``bin_refinement --skip-gunc``: the
    GUNC columns are left out when GUNC did not run.
    """
    return (
        list(_checkm.STATS_HEADER)
        + (list(GUNC_STATS_COLUMNS) if with_gunc else [])
        + list(_mimag.STATS_COLUMNS)
    )


# ── -c and -x ────────────────────────────────────────────────────────────────────────────────
# `bin_refinement` uses -c/--completeness and -x/--contamination; CheckM, GUNC and GTDB-Tk use -x
# for the FASTA extension. Here -c keeps bin_refinement's meaning, -x is the extension (metaWRAP 1
# writes .fa and MetaWrap2 writes .fasta), and the contamination ceiling is --contamination only.


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 score",
        usage="metawrap2 score [options] -o output_dir BIN_DIR [BIN_DIR ...]",
        description="Score one or more directories of genome FASTAs with one scorer, writing the "
        "same per-genome table bin_refinement writes. Several directories are compared side by "
        "side.",
    )
    p.add_argument(
        "bin_dirs",
        nargs="+",
        metavar="BIN_DIR",
        help="a directory of genome FASTAs; repeat to score and compare several bin sets",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    threads_arg(p)
    memory_arg(p)
    p.add_argument(
        "-c",
        "--completeness",
        type=percent(),
        default=DEFAULT_MIN_COMPLETION,
        help="completeness floor for the 'passing' count (same flag as bin_refinement)",
    )
    p.add_argument(
        "--contamination",
        type=percent(),
        default=DEFAULT_MAX_CONTAMINATION,
        help="contamination ceiling for the 'passing' count (bin_refinement's -x; -x is the "
        "FASTA extension here)",
    )
    p.add_argument(
        "-x",
        "--extension",
        default=BIN_EXTENSION,
        help="extension of the genome FASTAs, with or without the leading dot (MetaWrap2 writes "
        "%s, metaWRAP 1 writes .fa)" % BIN_EXTENSION,
    )
    p.add_argument(
        "--checkm1",
        action="store_true",
        help="score with CheckM1 instead of CheckM2, to compare with published CheckM1 numbers. "
        "Slower, and pplacer needs ~40GB of memory per thread",
    )
    p.add_argument(
        "--no-gunc",
        dest="run_gunc",
        action="store_false",
        help="don't run the GUNC chimerism check, so its 6.7 GB reference database is not "
        "needed",
    )
    p.add_argument(
        "--score-everything",
        action="store_true",
        help="score every genome with CheckM, including ones too small to pass. Without it, those "
        "are recorded as 0/100 completeness/contamination without being scored",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def normalise_extension(extension: str) -> str:
    """Return *extension* with a leading dot (``fa`` becomes ``.fa``); empty means the default.

    CheckM takes ``-x fasta`` and GUNC ``--file_suffix .fasta``, so users may type either form.
    """
    text = extension.strip()
    if not text:
        return BIN_EXTENSION
    return text if text.startswith(".") else "." + text


def set_names(bin_dirs: Sequence[str]) -> List[str]:
    """One output name per input directory: its basename, unless two of them share one.

    The name is used for the bin set's output files (``<name>.stats.tsv``, ``<name>.gunc/``) and
    its row in the comparison table. Both pipelines name a consolidated bin set after its
    thresholds (``metawrap_70_10_bins``), so collisions are common. A colliding name is prefixed
    with its parent directory's name, and one that still collides gets a counter.
    """
    # Strip trailing slashes first: `os.path.dirname("/runs/a/bins/")` is "/runs/a/bins", which
    # would prefix the name with itself instead of its parent.
    paths = [strip_trailing_slash(path) for path in bin_dirs]
    bases = [os.path.basename(path) or "bins" for path in paths]
    shared = {name for name, n in Counter(bases).items() if n > 1}
    names: List[str] = []
    taken: Dict[str, int] = {}
    for path, base in zip(paths, bases):
        wanted = base
        if base in shared:
            parent = os.path.basename(os.path.dirname(path))
            if parent:
                wanted = "%s.%s" % (parent, base)
        candidate, suffix = wanted, 1
        while candidate in taken:
            suffix += 1
            candidate = "%s.%d" % (wanted, suffix)
        taken[candidate] = 1
        names.append(candidate)
    return names


def _stage(bin_dir: str, staged: str, extension: str) -> int:
    """Link *bin_dir*'s genomes into *staged* under the canonical extension. Returns how many.

    :mod:`metawrap2.checkm` and :mod:`metawrap2.gunc` pass
    :data:`metawrap2.constants.BIN_EXTENSION` to every tool, so genomes with another extension are
    linked under that one (``bin_refinement`` normalises extensions in its working copy the same
    way). The scorers write beside their input directory, so staging under ``-o`` also keeps the
    input directory unchanged.

    Symlinks are used to avoid copying large bin sets, as in
    :func:`metawrap2.checkm.run_checkm_cached`; if the filesystem refuses a link, the file is
    copied. The unbinned catch-all (``unbinned.fasta``) is skipped by
    :func:`metawrap2.io.seqio.fasta_files`, since it is not a genome.
    """
    ensure_dir(staged)
    found = fasta_files(bin_dir, (extension,))
    for path in found:
        link = os.path.join(staged, os.path.splitext(os.path.basename(path))[0] + BIN_EXTENSION)
        if os.path.lexists(link):
            os.remove(link)  # a --force or --resume into the same output directory
        try:
            os.symlink(os.path.abspath(path), link)
        except OSError:
            shutil.copy2(path, link)
    return len(found)


class Scored(NamedTuple):
    """What one bin set's scoring produced, and the counts read back out of its own table."""

    #: The name this bin set's files are under, from :func:`set_names`.
    name: str
    #: The directory the genomes came from, so a surprising row can be traced to its input.
    source: str
    #: The ``.stats.tsv`` written for it.
    stats: str
    #: How many genomes were staged and scored.
    genomes: int
    #: Genomes clearing the run's ``-c``/``--contamination``.
    passing: int
    #: Genomes clearing each :data:`metawrap2.benchmark.scoring.QUALITY_TIERS` tier, plus the total.
    tiers: Dict[str, int]
    #: Genomes GUNC called not chimeric, how many it scored, and how many near-complete genomes
    #: it passed. Counted per genome from the table.
    gunc_pass: Optional[int]
    gunc_scored: Optional[int]
    clean_near_complete: Optional[int]


def summarise_table(
    name: str,
    source: str,
    stats: str,
    *,
    completeness: float,
    contamination: float,
    with_gunc: bool,
) -> Scored:
    """Read a finished ``.stats.tsv`` back and count what is in it.

    The counts come from the file, so they can be checked against the table. Completeness and
    contamination are read with :func:`metawrap2.benchmark.scoring.bin_quality`, which drops rows
    it cannot parse; the other columns with :func:`metawrap2.magtable.read_stats`, which treats a
    blank cell as not measured.

    ``clean_near_complete`` counts genomes that are near-complete and not called chimeric, decided
    from each genome's own row.
    """
    if not os.path.isfile(stats):
        # --dry-run: nothing ran and there is no table, so all counts are zero.
        return Scored(name, source, stats, 0, 0, quality_counts([]), None, None, None)
    rows = _magtable.read_stats(stats)
    quality = bin_quality(stats)

    def passed_gunc(genome: str) -> str:
        return (rows.get(genome, {}).get("gunc_pass", "") or "").strip().lower()

    # Same count as bin_refinement (bounds inclusive).
    passing = len(_refinement_module.good_bins(stats, completeness, contamination))
    near_complete = [
        genome
        for genome, complete, contaminated in quality
        if complete >= NEAR_COMPLETE_COMPLETENESS and contaminated <= NEAR_COMPLETE_CONTAMINATION
    ]
    gunc_scored = sum(1 for genome in rows if passed_gunc(genome) in ("yes", "no"))
    gunc_pass = sum(1 for genome in rows if passed_gunc(genome) == "yes")
    clean = sum(1 for genome in near_complete if passed_gunc(genome) == "yes")
    tiers = quality_counts(quality)
    return Scored(
        name=name,
        source=source,
        stats=stats,
        genomes=len(rows),
        passing=passing,
        tiers=tiers,
        gunc_pass=gunc_pass if with_gunc else None,
        gunc_scored=gunc_scored if with_gunc else None,
        clean_near_complete=clean if with_gunc else None,
    )


def _write_empty_stats(stats: str, with_gunc: bool) -> None:
    """A header-only table for a bin set with no genomes in it.

    Writes the full header, where ``bin_refinement.write_empty_result`` writes CheckM's columns
    only, so it has the same columns as the other bin sets' tables.
    """
    with open(stats, "w") as fh:
        fh.write("\t".join(stats_columns(with_gunc)) + "\n")


def score_one(
    bin_dir: str,
    name: str,
    args,
    *,
    env: Optional[str],
    checkm2_env: Optional[str],
    gunc_db: str,
    cache: _checkm.BinScoreCache,
) -> Scored:
    """Score one directory of genomes, leaving ``<name>.stats.tsv`` in the output directory."""
    staged = os.path.join(args.output, name)
    genomes = _stage(bin_dir, staged, args.extension)
    stats = staged + STATS_SUFFIX

    if genomes == 0 and not dry_run():
        # A warning, not an error: the other bin sets are still scored, and this one gets an empty
        # table. The message names the extension because a wrong -x is the usual cause.
        warning(
            "No %s files in %s, so there is nothing to score in it.\n"
            "If the genomes are there, they have a different extension - pass -x with theirs "
            "(metaWRAP 1 writes .fa, MetaWrap2 writes %s)."
            % (args.extension, bin_dir, BIN_EXTENSION)
        )
        _write_empty_stats(stats, args.run_gunc)
        return summarise_table(
            name,
            bin_dir,
            stats,
            completeness=args.completeness,
            contamination=args.contamination,
            with_gunc=args.run_gunc,
        )

    comm(
        "scoring %d genome(s) from %s with %s"
        % (genomes, bin_dir, "CheckM1" if args.checkm1 else "CheckM2")
    )
    _checkm.run_checkm_cached(
        staged,
        cache=cache,
        prefilter=not args.score_everything,
        say=comm,
        threads=args.threads,
        mem_gb=args.memory,
        env=env,
        checkm2=not args.checkm1,
        checkm2_env=checkm2_env,
        log_path=staged + ".checkm.log",
    )

    if args.run_gunc:
        report = _gunc.run_gunc(
            staged,
            threads=args.threads,
            # Placeholder for --dry-run when no database is configured.
            db=gunc_db or "<GUNC_DB>",
            env=env,
            out_dir=staged + ".gunc",
            log_path=staged + ".gunc.log",
        )
        verdicts = _gunc.parse_report(str(report)) if report else {}
        # Called even when GUNC scored nothing, so every table written with GUNC on has its two
        # columns. Genomes GUNC could not assess get empty cells.
        if not dry_run():
            annotated = _gunc.annotate_stats(stats, verdicts)
            chimeric = sorted(genome for genome, v in verdicts.items() if not v.passed)
            comm("GUNC scored %d genome(s) in %s; %d annotated" % (len(verdicts), name, annotated))
            if chimeric:
                warning(
                    "GUNC calls %d of %d genomes in %s chimeric (clade separation score above "
                    "%s): %s\nThey are kept and counted, with gunc_pass=no in the table - a "
                    "genome two organisms wide can look complete and clean to CheckM."
                    % (len(chimeric), len(verdicts), name, _gunc.CSS_THRESHOLD, ", ".join(chimeric))
                )

    if not dry_run():
        # There is no rRNA/tRNA evidence for these genomes, so genomes clearing 90/5 are tagged
        # high_quality_unverified_rrna.
        assessments = _checkm.tag_mimag(stats)
        if assessments:
            comm(_mimag.summarise(assessments))
            if _mimag.tally(assessments)[_mimag.MimagTag.HIGH_UNVERIFIED]:
                comm(_mimag.NO_EVIDENCE_HINT)

    return summarise_table(
        name,
        bin_dir,
        stats,
        completeness=args.completeness,
        contamination=args.contamination,
        with_gunc=args.run_gunc,
    )


def _threshold_label(completeness: float, contamination: float) -> str:
    """The run's thresholds as a column heading, e.g. ``>=70/<=10``."""
    return ">=%g/<=%g" % (completeness, contamination)


def comparison_table(results: Sequence[Scored], completeness: float, contamination: float) -> str:
    """One row per bin set, so two of them can be read side by side.

    Printed when more than one directory is scored.
    """
    scorer_counts = any(r.gunc_pass is not None for r in results)
    headings = [
        "bin set",
        "genomes",
        _threshold_label(completeness, contamination),
        "near-complete",
    ]
    if scorer_counts:
        headings += ["gunc pass", "clean near-complete"]
    rows = []
    for result in results:
        row = [
            result.name,
            str(result.genomes),
            str(result.passing),
            str(result.tiers.get(NEAR_COMPLETE_TIER, 0)),
        ]
        if scorer_counts:
            row += [
                "%d/%d" % (result.gunc_pass or 0, result.gunc_scored or 0),
                str(result.clean_near_complete or 0),
            ]
        rows.append(row)
    widths = [max(len(cell) for cell in column) for column in zip(headings, *rows)]

    def line(cells: Sequence[str]) -> str:
        # Bin set names left-aligned, counts right-aligned.
        laid_out = [cells[0].ljust(widths[0])] + [
            cell.rjust(width) for cell, width in zip(cells[1:], widths[1:])
        ]
        return "  " + "  ".join(laid_out).rstrip()

    return "\n".join(
        [line(headings), "  " + "  ".join("-" * width for width in widths)]
        + [line(row) for row in rows]
    )


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    resolve_threads(args, settings)
    resolve_memory_gb(args, settings)
    absolutize_paths(args)
    # absolutize_paths only handles the list arguments it knows (`reads`, `seqs`).
    args.bin_dirs = [os.path.abspath(os.path.expanduser(path)) for path in args.bin_dirs]
    args.extension = normalise_extension(args.extension)

    missing = [path for path in args.bin_dirs if not os.path.isdir(path)]
    if missing:
        # Checked before any env or tool, so a mistyped path fails at once.
        error(
            "Not a directory: %s\nEach BIN_DIR must be a directory of genome FASTAs."
            % ", ".join(missing)
        )

    # CheckM1 and GUNC are in bin_refinement's env, which metawrap2.benchmark.assess also uses.
    env = env_for(CHECKM_MODULE, settings)
    # CheckM2 is the default, as in bin_refinement, so its env is checked up front.
    checkm2_env = (
        None
        if args.checkm1
        else optional_env_for("checkm2", settings, "CheckM2, which scores unless --checkm1")
    )
    gunc_db = settings.db("GUNC_DB") if args.run_gunc else ""
    if args.run_gunc and not (gunc_db and os.path.exists(gunc_db)):
        # A dry run warns and continues printing the plan. A real run stops before CheckM runs.
        (warning if dry_run() else error)(
            "GUNC runs by default and its database is not configured (GUNC_DB=%s). %s\n"
            "Or pass --no-gunc to score completeness and contamination only."
            % (gunc_db or "<not set>", _gunc.DB_HINT)
        )

    rec = start_run("score", args, env, settings, inputs=list(args.bin_dirs))
    try:
        names = set_names(args.bin_dirs)
        # One cache for all bin sets: bin sets refined from one assembly often share identical
        # bins, and a score depends only on the contigs. Kept under -o for --resume.
        cache = _checkm.BinScoreCache(
            os.path.join(args.output, SCORE_CACHE % ("checkm1" if args.checkm1 else "checkm2"))
        )
        announcement("SCORING %d BIN SET(S)" % len(args.bin_dirs))
        results = [
            score_one(
                bin_dir,
                name,
                args,
                env=env,
                checkm2_env=checkm2_env,
                gunc_db=gunc_db,
                cache=cache,
            )
            for bin_dir, name in zip(args.bin_dirs, names)
        ]
        cache.save()

        print()
        for result in results:
            print("  %s  (%s)" % (result.name, result.source))
            print("    %-26s %d" % ("genomes", result.genomes))
            print(
                "    %-26s %d"
                % (
                    _threshold_label(args.completeness, args.contamination),
                    result.passing,
                )
            )
            print(
                "    %-26s %d"
                % (
                    "near-complete (%g/%g)"
                    % (NEAR_COMPLETE_COMPLETENESS, NEAR_COMPLETE_CONTAMINATION),
                    result.tiers.get(NEAR_COMPLETE_TIER, 0),
                )
            )
            if result.gunc_pass is not None:
                print(
                    "    %-26s %d of %d scored"
                    % ("not chimeric (GUNC)", result.gunc_pass, result.gunc_scored or 0)
                )
                print("    %-26s %d" % ("clean near-complete", result.clean_near_complete or 0))
            print("    %-26s %s" % ("table", result.stats))
            print()

        if len(results) > 1:
            print("  Scored with %s:\n" % ("CheckM1" if args.checkm1 else "CheckM2"))
            print(comparison_table(results, args.completeness, args.contamination))
            print("\n  %s\n" % COMPARISON_NOTE)
        announcement("SCORING FINISHED SUCCESSFULLY!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
