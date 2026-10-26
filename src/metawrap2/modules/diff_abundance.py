"""MetaWrap2 diff_abundance: features that differ between sample groups, tested compositionally.

Relative abundance tables are compositional: values are sum-constrained, so one feature rising
pushes the others down, and a t-test or Wilcoxon on them tests the constraint. The pipeline
previously stopped at feature tables and left this analysis to users.

The test is ANCOM-BC from scikit-bio, which avoids a dependency on R (other implementations are R
packages).

Input is the long tables `quantify_function` writes and a two-column sample sheet of groups. The
design is ``feature ~ group`` only; covariates, repeated measures and interactions depend on the
experiment and are not guessed.
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import Dict, List

from .. import biomout as _biom
from .. import diffabund as _da
from ..cliformat import make_parser
from ..config import load_settings
from ..constants import TSV_SUFFIX
from ._common import (
    absolutize_paths,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    error,
    finish_run,
    require_file,
    start_run,
    validate_inputs,
    warning,
)

# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: One row per feature tested.
RESULTS_TABLE = "differential_abundance" + TSV_SUFFIX
#: The settings used (correction method, prevalence filter, ...), needed to compare adjusted
#: p-values between runs.
SETTINGS_FILE = "differential_abundance_settings" + TSV_SUFFIX
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 diff_abundance",
        usage="metawrap2 diff_abundance [options] -o out -i table.tsv -g groups.tsv",
        description="Differential abundance between sample groups, by ANCOM-BC.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument(
        "-i",
        "--table",
        required=True,
        help="a long feature table from `metawrap2 quantify_function`: columns sample, feature, "
        "value",
    )
    p.add_argument(
        "-g",
        "--groups",
        required=True,
        help="TSV of sample<TAB>group, with a header. Samples absent from it are not tested",
    )
    p.add_argument(
        "--term",
        help="which model term to report (default: the first non-intercept term, which for a "
        "two-group comparison is the only one)",
    )
    p.add_argument(
        "--alpha",
        type=float,
        default=_da.ALPHA,
        metavar="P",
        help="significance threshold for the adjusted p-value (default: %(default)s)",
    )
    p.add_argument(
        "--p-adjust",
        default=_da.P_ADJUST,
        help="multiple-testing correction (default: %(default)s), needed when testing e.g. 573 "
        "modules or 28,000 KOs at once. The method chosen changes which results are significant",
    )
    p.add_argument(
        "--min-prevalence",
        type=float,
        default=_da.MIN_PREVALENCE,
        metavar="FRACTION",
        help="only test features present in at least this fraction of samples (default: "
        "%(default)s, test everything). Dropping a feature changes every other feature's "
        "multiple-testing correction",
    )
    p.add_argument(
        "--biom",
        action="store_true",
        help="also write the input table as BIOM beside the results",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _read_long_table(path: str) -> Dict[str, Dict[str, float]]:
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        for required in ("sample", "feature", "value"):
            if required not in fields:
                error(
                    "%s has no %r column; this takes the long tables quantify_function writes "
                    "(found %r)" % (path, required, fields)
                )
        return _biom.long_table_to_values(list(reader))


def _read_groups(path: str) -> Dict[str, str]:
    """sample -> group, from a two-column TSV.

    A header with ``sample`` as its first column is required, so a data row is never taken for a
    header and dropped.
    """
    groups: Dict[str, str] = {}
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        error("%s is empty" % path)
    header = [cell.strip().lower() for cell in rows[0]]
    if header[:1] != ["sample"]:
        error(
            "%s needs a header whose first column is 'sample' (found %r). The header is required "
            "so that a data row is not read as a header" % (path, rows[0][:2])
        )
    for row in rows[1:]:
        if len(row) < 2 or not row[0].strip():
            continue
        groups[row[0].strip()] = row[1].strip()
    return groups


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    absolutize_paths(args)
    validate_inputs(
        "diff_abundance",
        [
            (require_file, args.table, "feature table (-i)"),
            (require_file, args.groups, "group sheet (-g)"),
        ],
    )

    rec = start_run("diff_abundance", args, None, settings, inputs=[args.table, args.groups])
    try:
        announcement("BEGIN DIFF_ABUNDANCE PIPELINE!")
        ensure_dir(args.output)
        if dry_run():
            comm("would test every feature for a difference between groups with ANCOM-BC")
            announcement("DIFF_ABUNDANCE PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        values = _read_long_table(args.table)
        groups = _read_groups(args.groups)
        shared = sorted(set(values) & set(groups))
        missing = sorted(set(values) - set(groups))
        if missing:
            warning(
                "%d sample(s) in the table have no group and are not tested: %s"
                % (len(missing), ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else ""))
            )
        comm(
            "%d sample(s) with a group, %d feature(s) in the table"
            % (len(shared), len({f for row in values.values() for f in row}))
        )

        try:
            results, used = _da.test_groups(
                values,
                groups,
                term=args.term,
                alpha=args.alpha,
                p_adjust=args.p_adjust,
                min_prevalence=args.min_prevalence,
            )
        except _da.StatsUnavailable as problem:
            error(str(problem))
            return 1
        except ValueError as problem:
            error(str(problem))
            return 1

        table = os.path.join(args.output, RESULTS_TABLE)
        with open(table, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=_da.RESULT_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for result in sorted(results, key=lambda r: (r.q_value, r.feature)):
                writer.writerow(_da.result_row(result, used))

        with open(os.path.join(args.output, SETTINGS_FILE), "w", newline="") as fh:
            settings_writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            settings_writer.writerow(["setting", "value"])
            for key, value in sorted(used.items()):
                settings_writer.writerow([key, value])

        significant = [r for r in results if r.significant]
        # The settings dict is untyped (it is written to a provenance file), so narrow this field.
        compared = used.get("groups")
        group_names = [str(g) for g in compared] if isinstance(compared, list) else []
        comm(
            "tested %d feature(s) across %s; %d significant at alpha %s after %s correction"
            % (
                len(results),
                " vs ".join(group_names),
                len(significant),
                args.alpha,
                args.p_adjust,
            )
        )
        for result in sorted(significant, key=lambda r: r.q_value)[:5]:
            comm(
                "  %s  log2FC %+.2f  q %.3g"
                % (result.feature, result.log_fold_change, result.q_value)
            )
        comm(_da.COMPOSITIONAL_CAVEAT)
        comm("-> %s" % table)

        if args.biom:
            biom_path = os.path.join(
                args.output, os.path.basename(args.table)[: -len(TSV_SUFFIX)] + ".biom"
            )
            try:
                n = _biom.write_biom(values, biom_path)
                comm("%d features -> %s" % (n, biom_path))
            except (_biom.BiomUnavailable, ValueError) as problem:
                warning("could not write BIOM: %s" % problem)

        announcement("DIFF_ABUNDANCE PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
