"""MetaWrap2 quantify_function: sample-level gene and pathway abundance.

Combines per-genome annotation (`annotate`, `pathways`) with genome abundances to give the
abundance of each feature per sample, in the units of :mod:`metawrap2.funcabund`. Also reports the
fraction of each sample's abundance that the annotated genomes account for.

* Coverage-based abundances (CoverM) are not divided by genome size; read counts are. The kind is
  set with ``--abundance-kind``, not inferred from the column name.
* Each unit (relative, per genome equivalent) is written to its own file, with what it may be
  compared across.
* A genome missing from the abundance table contributes nothing to that sample.
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from typing import Dict, List, Optional

from .. import biomout as _biom
from .. import funcabund as _fa
from .. import hierarchy as _hierarchy
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
    resolve_threads,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_THREADS = 1
#: Output file names, one file per unit.
KO_TABLE = "ko_abundance_%s" + TSV_SUFFIX
MODULE_TABLE = "pathway_abundance_%s" + TSV_SUFFIX
#: Fraction of each sample's abundance covered by annotated genomes.
ASSIGNABLE_TABLE = "assignable_fraction" + TSV_SUFFIX
#: Minimum module completeness for a module to count toward a sample's pathway abundance. 0.8 is
#: the usual threshold for "present" in the literature; it can be changed on the command line.
MIN_MODULE_COMPLETENESS = 0.8
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 quantify_function",
        usage="metawrap2 quantify_function [options] -o out -k kos.tsv -a abundance.tsv",
        description="Sample-level gene and pathway abundance from per-genome annotation.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument(
        "-k", "--kos", required=True, help="per-gene KO table from `metawrap2 annotate --kofam`"
    )
    p.add_argument(
        "-a",
        "--abundance",
        required=True,
        help="bin_abundance_table.tsv from `metawrap2 quant_bins`: genomes as rows, samples as "
        "columns",
    )
    p.add_argument(
        "-m",
        "--modules",
        help="modules.tsv from `metawrap2 pathways`, to also quantify pathways",
    )
    threads_arg(p)
    p.add_argument(
        "--abundance-kind",
        default=_fa.AbundanceKind.COVERAGE.value,
        choices=[kind.value for kind in _fa.AbundanceKind],
        help="whether the abundance table holds coverage-proportional numbers (CoverM, sylph, any "
        "depth-derived figure) or read counts (default: %(default)s). Read counts are divided by "
        "genome length once; coverage is already per cell and is not divided",
    )
    p.add_argument(
        "--genome-lengths",
        help="TSV of genome<TAB>length_bp, required only for --abundance-kind read_count",
    )
    p.add_argument(
        "--microbial-fraction",
        help="output of `singlem microbial_fraction`, which supplies the denominator for the "
        "per-genome-equivalent unit. Without it, that unit is not reported",
    )
    p.add_argument(
        "--min-module-completeness",
        type=float,
        default=MIN_MODULE_COMPLETENESS,
        metavar="FRACTION",
        help="a genome's module counts toward the sample only above this completeness "
        "(default: %(default)s)",
    )
    p.add_argument(
        "--biom",
        action="store_true",
        help="also write each table as BIOM 2.1 beside the TSV, for QIIME 2 / phyloseq / "
        "scikit-bio. Needs the biom extra: pip install 'metawrap2[biom]'",
    )
    p.add_argument(
        "--evidence",
        default="curated",
        choices=("curated", "all"),
        help="which KO calls to use (default: %(default)s)",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _read_ko_counts(path: str, evidence: str) -> Dict[str, Dict[str, float]]:
    """genome -> KO -> copy number, from the per-gene KO table.

    A gene with several KOs splits its weight, so a genome's total is its number of annotated
    genes (as in :mod:`metawrap2.hierarchy`).
    """
    per_gene: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if evidence == "curated" and row.get("evidence") != "curated":
                continue
            per_gene[row["genome"]][row["gene_id"]].append(row["ko"])
    return {
        genome: _hierarchy.weights_from_gene_kos(genes.values())
        for genome, genes in per_gene.items()
    }


def _read_module_presence(path: str, minimum: float) -> Dict[str, Dict[str, float]]:
    """genome -> module -> 1.0 for modules complete enough to count.

    Binary, because summing fractional completeness across genomes would credit the community
    with pathways no single genome has.
    """
    out: Dict[str, Dict[str, float]] = defaultdict(dict)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                completeness = float(row["completeness_pathwise"])
            except (KeyError, ValueError):
                continue
            if completeness >= minimum:
                out[row["genome"]][row["module"]] = 1.0
    return dict(out)


def _read_abundance(
    path: str, kind: _fa.AbundanceKind, lengths: Optional[Dict[str, int]]
) -> List[_fa.GenomeAbundance]:
    entries: List[_fa.GenomeAbundance] = []
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 2:
            error("%s has no sample columns" % path)
            return entries
        samples = header[1:]
        for row in reader:
            if len(row) < 2:
                continue
            genome = row[0]
            for extension in (".fasta", ".fa", ".fna", ".fas"):
                if genome.endswith(extension):
                    genome = genome[: -len(extension)]
                    break
            for sample, raw in zip(samples, row[1:]):
                try:
                    value = float(raw)
                except ValueError:
                    continue
                entries.append(
                    _fa.GenomeAbundance(
                        genome=genome,
                        sample=sample,
                        value=value,
                        kind=kind,
                        genome_length_bp=(lengths or {}).get(genome),
                    )
                )
    return entries


def _read_lengths(path: Optional[str]) -> Optional[Dict[str, int]]:
    if not path:
        return None
    out: Dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                out[parts[0]] = int(float(parts[1]))
            except ValueError:
                continue
    return out


def _write_units(
    output: str,
    template: str,
    feature_type: str,
    per_sample: Dict[str, Dict[str, float]],
    per_genome: Dict[str, Dict[str, float]],
    abundances: List[_fa.GenomeAbundance],
    kind: _fa.AbundanceKind,
    fractions: Dict[str, float],
    biom: bool = False,
) -> List[str]:
    """Write one file per unit and return the paths written."""
    written: List[str] = []

    for unit in (_fa.Unit.RELATIVE, _fa.Unit.PER_GENOME_EQUIVALENT):
        rows: List[Dict[str, object]] = []
        for sample, features in sorted(per_sample.items()):
            if unit is _fa.Unit.RELATIVE:
                values = _fa.to_relative(features)
            else:
                fraction = fractions.get(sample)
                if not fraction:
                    continue
                values = _fa.to_genome_equivalents(features, fraction)
            rows.extend(
                _fa.abundance_rows(
                    sample,
                    values,
                    unit,
                    feature_type,
                    kind,
                    _fa.contributors_per_feature(per_genome, abundances, sample),
                )
            )
        if not rows:
            continue
        path = os.path.join(output, template % unit.value)
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=_fa.ABUNDANCE_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        written.append(path)
        if biom:
            # Built from the rows written to the TSV, so the two files match.
            biom_path = path[: -len(TSV_SUFFIX)] + ".biom"
            try:
                n = _biom.write_biom(
                    _biom.long_table_to_values(rows),
                    biom_path,
                    sample_metadata=_biom.unit_metadata(rows),
                    table_id=os.path.basename(biom_path),
                )
                written.append(biom_path)
                comm("%d features -> %s" % (n, biom_path))
            except _biom.BiomUnavailable as problem:
                warning(str(problem))
            except ValueError as problem:
                warning("could not write %s: %s" % (biom_path, problem))
    return written


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)

    checks = [
        (require_file, args.kos, "KO table (-k)"),
        (require_file, args.abundance, "abundance table (-a)"),
    ]
    if args.modules:
        checks.append((require_file, args.modules, "modules table (-m)"))
    if args.microbial_fraction:
        checks.append((require_file, args.microbial_fraction, "microbial fraction"))
    validate_inputs("quantify_function", checks)

    kind = _fa.AbundanceKind(args.abundance_kind)
    if kind is _fa.AbundanceKind.READ_COUNT and not args.genome_lengths:
        error(
            "--abundance-kind read_count needs --genome-lengths: read counts are divided "
            "by genome length to make them proportional to cell numbers."
        )

    rec = start_run("quantify_function", args, None, settings, inputs=[args.kos, args.abundance])
    # No checkpoint: this module only does fast table arithmetic over existing files.
    try:
        announcement("BEGIN QUANTIFY_FUNCTION PIPELINE!")
        ensure_dir(args.output)
        if dry_run():
            comm("would quantify gene and pathway abundance per sample")
            announcement("QUANTIFY_FUNCTION PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        lengths = _read_lengths(args.genome_lengths)
        ko_counts = _read_ko_counts(args.kos, args.evidence)
        abundances = _read_abundance(args.abundance, kind, lengths)
        if not ko_counts:
            error("no KO calls in %s matching --evidence %s" % (args.kos, args.evidence))
        if not abundances:
            error("no usable abundances in %s" % args.abundance)

        samples = sorted({entry.sample for entry in abundances})
        quantified = {entry.genome for entry in abundances}
        missing = sorted(set(ko_counts) - quantified)
        if missing:
            warning(
                "%d annotated genome(s) have no abundance and contribute nothing: %s"
                % (len(missing), ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else ""))
            )
        comm(
            "%d sample(s), %d genome(s) with abundance, %d with KO annotation"
            % (len(samples), len(quantified), len(ko_counts))
        )

        fractions: Dict[str, float] = {}
        if args.microbial_fraction:
            fractions = _fa.parse_singlem_fraction(args.microbial_fraction)
            comm("microbial fraction for %d sample(s)" % len(fractions))
        else:
            comm(
                "no --microbial-fraction, so the per-genome-equivalent unit is not reported. "
                "It is the only unit that is comparable across samples"
            )

        ko_per_sample = _fa.feature_abundance(ko_counts, abundances)
        written = _write_units(
            args.output,
            KO_TABLE,
            "ko",
            ko_per_sample,
            ko_counts,
            abundances,
            kind,
            fractions,
            args.biom,
        )
        for path in written:
            comm("-> %s" % path)

        if args.modules:
            presence = _read_module_presence(args.modules, args.min_module_completeness)
            module_per_sample = _fa.feature_abundance(presence, abundances)
            comm(
                "%d genome(s) have at least one module at >=%.0f%% complete"
                % (len(presence), 100 * args.min_module_completeness)
            )
            for path in _write_units(
                args.output,
                MODULE_TABLE,
                "kegg_module",
                module_per_sample,
                presence,
                abundances,
                kind,
                fractions,
                args.biom,
            ):
                comm("-> %s" % path)

        entries = _fa.assignable(abundances, ko_counts)
        path = os.path.join(args.output, ASSIGNABLE_TABLE)
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=_fa.ASSIGNABLE_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for entry in entries:
                writer.writerow(_fa.assignable_row(entry))
        comm("-> %s" % path)
        for entry in entries:
            comm(
                "%s: %.1f%% of quantified abundance is annotated"
                % (entry.sample, 100 * entry.assignable_fraction)
            )
        comm(
            "this is a fraction of the genomes in the abundance table, not of the community. "
            "`singlem appraise` estimates the community share, which is usually much smaller"
        )
        announcement("QUANTIFY_FUNCTION PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
