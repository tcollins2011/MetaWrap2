"""MetaWrap2 pathways: KEGG module completeness and functional profiles of genomes and samples.

Takes the per-gene KO table written by ``annotate --kofam`` and produces:

* **Module completeness per genome**, with both the pathwise and stepwise definitions, since
  published tools disagree on the algorithm. See :mod:`metawrap2.kegg`.
* **The functional hierarchy** (gene -> chemical process -> pathway) as percentages: per genome,
  with every gene weighted equally, and per sample, with each genome weighted by its abundance.
  The per-sample profile is not the average of the per-genome ones.
* **Krona plots** of both.

Completeness is computed from vendored definitions and the hierarchy from a vendored BRITE
projection. The only external program is Krona's ``ktImportText``, for the plots.
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from typing import Dict, List, Tuple

from .. import depth as _depth
from .. import hierarchy as _hierarchy
from .. import kegg as _kegg
from .. import sunburst as _sunburst
from ..cliformat import make_parser
from ..config import MODULE_ENVS, load_settings
from ..constants import TSV_SUFFIX
from ..progress import bar
from ._common import (
    absolutize_paths,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    make_checkpoint,
    require_file,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

#: Krona's importer is in the kraken2 env.
CONDA_ENV = MODULE_ENVS["kraken2"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# One plot per input file. -n names the root of the wheel; without it Krona labels it "all".
#
# {root} is quoted: Krona takes input files positionally, so an unquoted name with a space
# ("genome bin_006") becomes a second, nonexistent input path.
KRONA = "ktImportText {inputs} -o {out} -n '{root}'"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_THREADS = 1
#: genome x module, both completeness scores and every flag that qualifies them.
MODULES_TABLE = "modules" + TSV_SUFFIX
#: Long hierarchy tables, one for genomes and one for samples.
HIERARCHY_DIR = "hierarchy"
GENOME_HIERARCHY = "genome_hierarchy" + TSV_SUFFIX
SAMPLE_HIERARCHY = "sample_hierarchy" + TSV_SUFFIX
#: The same shape, weighted by per-gene abundance instead of per-MAG abundance. A separate file
#: so the two profiles are not summed together.
GENE_SAMPLE_HIERARCHY = "sample_hierarchy_by_gene_abundance" + TSV_SUFFIX
#: The same shape, weighted by per-contig depth. Needs no MAGs and no extra read mapping, since
#: `binning` already wrote the depth table.
CONTIG_SAMPLE_HIERARCHY = "sample_hierarchy_by_contig_depth" + TSV_SUFFIX
KRONA_DIR = os.path.join(HIERARCHY_DIR, "krona")
#: Modules below this pathwise completeness are left out of the table. The default 0 keeps all
#: 573 rows per genome, so absent modules are distinguishable from unreported ones.
MIN_COMPLETENESS = 0.0
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 pathways",
        usage="metawrap2 pathways [options] -o output_dir -k kos.tsv",
        description="KEGG module completeness and the functional hierarchy, per genome and "
        "per sample.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument(
        "-k",
        "--kos",
        required=True,
        help="per-gene KO table from `metawrap2 annotate --kofam` (kos.tsv)",
    )
    p.add_argument(
        "-C",
        "--contig-depth",
        help="a per-contig depth table from `metawrap2 binning`/`coverage` (contigName, contigLen, "
        "totalAvgDepth, then one column per sample). Weights every gene by its contig's depth. "
        "Use it when the analysis is not about MAGs: it covers the whole "
        "assembly, binned or not, and needs no extra mapping. Less precise per gene than -G, since "
        "a contig's depth is assumed uniform across its genes",
    )
    p.add_argument(
        "-G",
        "--gene-abundance",
        help="gene abundance in long format from `metawrap2 catalogue --quantify` (sample, "
        "feature, value). Weights the sample profile by the abundance of each gene instead of "
        "each MAG, giving the whole-sample potential: a gene "
        "catalogue built from the assembly includes unbinned sequence. Can be "
        "given with -a; both profiles are then written and each row names its "
        "weighting",
    )
    p.add_argument(
        "-a",
        "--abundance",
        help="bin_abundance_table.tsv from `metawrap2 quant_bins`; without it only per-genome "
        "distributions are written, because a community profile needs the abundance of each "
        "genome",
    )
    threads_arg(p)
    p.add_argument(
        "--min-completeness",
        type=float,
        default=MIN_COMPLETENESS,
        metavar="FRACTION",
        help="omit modules below this pathwise completeness (default: %(default)s, keep all)",
    )
    p.add_argument(
        "--brite-db",
        metavar="DIR",
        help="directory (or file) holding the BRITE ko00001 projection `install-db kegg` writes "
        "(default: KEGG_BRITE from metawrap2.toml). Without it, module completeness is written "
        "and the hierarchy and Krona plots are skipped",
    )
    p.add_argument(
        "--skip-krona",
        action="store_true",
        help="write the hierarchy tables but not the Krona plots",
    )
    p.add_argument(
        "--evidence",
        default="curated",
        choices=("curated", "all"),
        help="which KO calls to use: 'curated' keeps those that passed a curated KOfam "
        "threshold, 'all' keeps every call including relative and thresholdless ones "
        "(default: %(default)s)",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _read_kos(path: str, evidence: str) -> Tuple[Dict[str, Dict[str, float]], Dict[str, int]]:
    """(genome -> KO -> gene-count weight, genome -> genes with any KO).

    A gene with several KOs contributes a fraction to each (see
    :func:`metawrap2.hierarchy.weights_from_gene_kos`), so a genome's weights sum to its number of
    annotated genes and the percentages are shares of its genes.
    """
    per_gene: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if evidence == "curated" and row.get("evidence") != "curated":
                continue
            per_gene[row["genome"]][row["gene_id"]].append(row["ko"])

    weights: Dict[str, Dict[str, float]] = {}
    annotated: Dict[str, int] = {}
    for genome, genes in per_gene.items():
        weights[genome] = _hierarchy.weights_from_gene_kos(genes.values())
        annotated[genome] = len(genes)
    return weights, annotated


def _read_gene_abundance(path: str) -> Dict[str, Dict[str, float]]:
    """sample -> gene -> abundance, from `catalogue --quantify`'s long table.

    Long format, as used by `diff_abundance` and the BIOM writer.
    """
    out: Dict[str, Dict[str, float]] = defaultdict(dict)
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        for required in ("sample", "feature", "value"):
            if required not in fields:
                error(
                    "%s has no %r column; this takes the long gene-abundance table "
                    "`catalogue --quantify` writes (found %r)" % (path, required, fields)
                )
        for row in reader:
            try:
                value = float(row["value"])
            except (TypeError, ValueError):
                continue
            if value:
                out[row["sample"]][row["feature"]] = value
    return dict(out)


def _read_contig_depth(path: str) -> Dict[str, Dict[str, float]]:
    """sample -> contig -> depth, from the per-contig depth table.

    The table has three leading columns, then a mean-depth and a ``-var`` column per sample. The
    variance columns are skipped by name.
    """
    out: Dict[str, Dict[str, float]] = defaultdict(dict)
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header or header[0] != "contigName":
            error(
                "%s does not look like a contig depth table: its first column is %r, not "
                "'contigName'" % (path, header[0] if header else None)
            )
        samples = [
            (index, name)
            for index, name in enumerate(header)
            if index >= len(_depth.LEADING_COLUMNS) and not name.endswith("-var")
        ]
        if not samples:
            error("%s has no per-sample depth columns" % path)
        for row in reader:
            if not row:
                continue
            contig = row[0]
            for index, name in samples:
                if index >= len(row):
                    continue
                try:
                    value = float(row[index])
                except ValueError:
                    continue
                if value:
                    out[name][contig] = value
    return dict(out)


def _gene_kos(path: str) -> Dict[str, List[str]]:
    """gene -> its KOs, ignoring the genome column.

    Quantified catalogue genes are named by source, not genome, so only the gene id is matched.
    The ids agree when the catalogue itself was annotated; the caller warns when they do not.
    """
    out: Dict[str, List[str]] = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = (row.get("gene_id") or "").strip()
            ko = (row.get("ko") or "").strip()
            if gene and ko:
                out[gene].append(ko)
    return dict(out)


def _read_abundance(path: str) -> Dict[str, Dict[str, float]]:
    """sample -> genome -> abundance, from quant_bins' table (bins as rows, samples as columns).

    Row labels have the bin's file extension (``bin.1.fa``) and the KO table's genome names do
    not, so the extension is stripped.
    """
    out: Dict[str, Dict[str, float]] = defaultdict(dict)
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header or len(header) < 2:
            error("%s does not look like an abundance table (no sample columns)" % path)
            return {}
        samples = header[1:]
        # Add every sample up front, so main() can warn about a sample whose column did not parse.
        for sample in samples:
            out[sample] = {}
        for row in reader:
            if len(row) < 2:
                continue
            genome = row[0]
            stem = genome
            for extension in (".fasta", ".fa", ".fna", ".fas"):
                if stem.endswith(extension):
                    stem = stem[: -len(extension)]
                    break
            for sample, value in zip(samples, row[1:]):
                try:
                    out[sample][stem] = float(value)
                except ValueError:
                    continue
    return dict(out)


def _write_modules(path: str, weights: Dict[str, Dict[str, float]], min_completeness: float) -> int:
    modules = _kegg.load_modules()
    ko_counts = _kegg.ko_module_counts(modules)
    written = 0
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=_kegg.MODULE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for genome in bar(
            sorted(weights), desc="scoring modules", unit=" genome", total=len(weights)
        ):
            present = set(weights[genome])
            for score in _kegg.score_all(present, modules, ko_counts):
                if score.completeness_pathwise < min_completeness:
                    continue
                writer.writerow(_kegg.module_row(genome, score, modules[score.module].is_signature))
                written += 1
    return written


def _write_hierarchy(
    path: str,
    unit: str,
    per_unit: Dict[str, Dict[str, float]],
    tree: _hierarchy.Hierarchy,
    weighting: str = _hierarchy.WEIGHTING_GENES,
) -> Dict[str, Dict[Tuple[str, str, str], float]]:
    """Write the long table and return each unit's distribution, for the Krona step."""
    distributions: Dict[str, Dict[Tuple[str, str, str], float]] = {}
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=_hierarchy.HIERARCHY_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for name in sorted(per_unit):
            distribution = _hierarchy.distribute(per_unit[name], tree)
            distributions[name] = distribution
            writer.writerows(
                _hierarchy.hierarchy_rows(unit, name, distribution, tree, weighting=weighting)
            )
    return distributions


def _draw_krona(
    args,
    env,
    distributions: Dict[str, Dict[Tuple[str, str, str], float]],
    tree: _hierarchy.Hierarchy,
    unit: str,
) -> None:
    """One plot per unit, twice: the whole hierarchy and the pathway branches alone.

    The full wheel shows how much of a genome KEGG cannot place in a pathway; that branch is the
    largest in BRITE and dominates the plot, so a pathway-only wheel is also drawn.
    """
    krona_dir = os.path.join(args.output, KRONA_DIR)
    ensure_dir(krona_dir)
    for name, distribution in bar(
        sorted(distributions.items()), desc="krona (%s)" % unit, unit=" plot"
    ):
        for suffix, pathways_only in (("all", False), ("pathways", True)):
            text = os.path.join(krona_dir, "%s.%s.%s.txt" % (unit, name, suffix))
            rows = _hierarchy.write_krona_input(
                distribution, text, hierarchy=tree, pathways_only=pathways_only
            )
            if not rows:
                continue
            html = os.path.join(krona_dir, "%s.%s.%s.html" % (unit, name, suffix))
            rc = run(
                KRONA.format(inputs=text, out=html, root="%s %s" % (unit, name)),
                env=env,
                tool="ktImportText",
                check=False,
            )
            if rc != 0:
                warning(
                    "ktImportText failed for %s %s; the input TSV is kept at %s"
                    % (unit, name, text)
                )
            # A static image of the same wheel, for reports; Krona writes only interactive HTML.
            try:
                _sunburst.draw_from_krona_text(
                    text, title="%s %s (%s)" % (unit.replace("_", " "), name, suffix)
                )
            except Exception as problem:  # noqa: BLE001 - a picture must not fail the analysis
                warning("could not draw the static wheel for %s %s: %s" % (unit, name, problem))


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)

    checks = [(require_file, args.kos, "KO table (-k)")]
    if args.abundance:
        checks.append((require_file, args.abundance, "abundance table (-a)"))
    validate_inputs("pathways", checks)

    env = env_for("kraken2", settings) if not args.skip_krona else None
    rec = start_run("pathways", args, env, settings, inputs=[args.kos])

    ckpt = make_checkpoint(args.output)
    try:
        announcement("BEGIN PATHWAYS PIPELINE!")
        ensure_dir(args.output)

        if dry_run():
            comm("would score 573 KEGG modules and build the functional hierarchy")
            announcement("PATHWAYS PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        weights, annotated = _read_kos(args.kos, args.evidence)
        if not weights:
            error(
                "no KO calls in %s matching --evidence %s. Nothing to do."
                % (args.kos, args.evidence)
            )
        comm(
            "%d genome(s), %d genes with a KO, %d distinct KOs"
            % (
                len(weights),
                sum(annotated.values()),
                len({ko for kos in weights.values() for ko in kos}),
            )
        )

        modules_path = os.path.join(args.output, MODULES_TABLE)
        if ckpt.todo("modules"):
            rows = _write_modules(modules_path, weights, args.min_completeness)
            comm("%d module rows -> %s" % (rows, modules_path))
            ckpt.done("modules")
        else:
            comm("skipping module completeness (already done; --resume)")

        # The flag overrides the config key, as for every other database flag.
        brite = args.brite_db or settings.db("KEGG_BRITE")
        if not brite or not os.path.isfile(_hierarchy.hierarchy_path(brite)):
            warning(
                "Skipping the functional hierarchy and Krona plots: no BRITE hierarchy. Install it "
                "with `metawrap2 install-db kegg`, set KEGG_BRITE in metawrap2.toml, or pass "
                "--brite-db. Module completeness is written."
            )
            announcement("PATHWAYS PIPELINE SUCCESSFULLY FINISHED!!!")
            finish_run(rec, args.output)
            return 0
        tree = _hierarchy.load_hierarchy(brite)
        hierarchy_dir = os.path.join(args.output, HIERARCHY_DIR)
        ensure_dir(hierarchy_dir)

        genome_distributions: Dict[str, Dict[Tuple[str, str, str], float]] = {}
        if ckpt.todo("genome_hierarchy"):
            genome_distributions = _write_hierarchy(
                os.path.join(hierarchy_dir, GENOME_HIERARCHY),
                "genome",
                weights,
                tree,
                weighting=_hierarchy.WEIGHTING_GENES,
            )
            comm("per-genome hierarchy -> %s" % os.path.join(hierarchy_dir, GENOME_HIERARCHY))
            ckpt.done("genome_hierarchy")
        else:
            comm("skipping per-genome hierarchy (already done; --resume)")

        sample_distributions: Dict[str, Dict[Tuple[str, str, str], float]] = {}
        if args.abundance:
            abundance = _read_abundance(args.abundance)
            per_sample: Dict[str, Dict[str, float]] = {}
            for sample, genomes in abundance.items():
                matched = set(genomes) & set(weights)
                if not matched:
                    warning(
                        "no genome in sample %s matches a genome in the KO table; skipping the "
                        "sample (its profile would be empty)" % sample
                    )
                    continue
                per_sample[sample] = _hierarchy.weighted_sample_profile(weights, genomes)
            if per_sample:
                sample_distributions = _write_hierarchy(
                    os.path.join(hierarchy_dir, SAMPLE_HIERARCHY),
                    "sample",
                    per_sample,
                    tree,
                    weighting=_hierarchy.WEIGHTING_MAG_ABUNDANCE,
                )
                comm(
                    "per-sample hierarchy (abundance-weighted) -> %s"
                    % os.path.join(hierarchy_dir, SAMPLE_HIERARCHY)
                )
        gene_distributions: Dict[str, Dict[Tuple[str, str, str], float]] = {}
        if args.gene_abundance:
            per_sample_genes = _read_gene_abundance(args.gene_abundance)
            gene_kos = _gene_kos(args.kos)
            profiles: Dict[str, Dict[str, float]] = {}
            # Named differently from `abundance`, which holds the per-genome table above.
            for sample, gene_values in sorted(per_sample_genes.items()):
                matched = set(gene_values) & set(gene_kos)
                if not matched:
                    warning(
                        "no gene in sample %s matches a gene in the KO table. The usual cause is "
                        "that the catalogue was annotated separately from the genomes, so the gene "
                        "ids differ (`catalogue` prefixes its genes with their source)" % sample
                    )
                    continue
                profiles[sample] = _hierarchy.weights_from_gene_abundance(gene_kos, gene_values)
                comm(
                    "%s: %d of %d quantified gene(s) have a KO"
                    % (sample, len(matched), len(gene_values))
                )
            if profiles:
                gene_distributions = _write_hierarchy(
                    os.path.join(hierarchy_dir, GENE_SAMPLE_HIERARCHY),
                    "sample",
                    profiles,
                    tree,
                    weighting=_hierarchy.WEIGHTING_GENE_ABUNDANCE,
                )
                comm(
                    "per-sample hierarchy (gene-abundance-weighted) -> %s"
                    % os.path.join(hierarchy_dir, GENE_SAMPLE_HIERARCHY)
                )
                comm(
                    "this profile covers the whole sample, including unbinned "
                    "sequence; the MAG-weighted one covers the binned part. The two measure "
                    "different fractions of the sample; each row names its "
                    "weighting"
                )

        contig_distributions: Dict[str, Dict[Tuple[str, str, str], float]] = {}
        if args.contig_depth:
            per_sample_contigs = _read_contig_depth(args.contig_depth)
            contig_kos = _gene_kos(args.kos)
            contig_profiles: Dict[str, Dict[str, float]] = {}
            contigs = {_hierarchy.contig_from_gene(gene) for gene in contig_kos}
            for sample, depths in sorted(per_sample_contigs.items()):
                matched = contigs & set(depths)
                if not matched:
                    warning(
                        "no contig in sample %s matches a contig the KO table's genes came from. "
                        "The usual cause is a depth table built over a different assembly than the "
                        "one annotated" % sample
                    )
                    continue
                contig_profiles[sample] = _hierarchy.weights_from_contig_depth(contig_kos, depths)
                comm(
                    "%s: %d of %d annotated contig(s) have depth"
                    % (sample, len(matched), len(contigs))
                )
            if contig_profiles:
                contig_distributions = _write_hierarchy(
                    os.path.join(hierarchy_dir, CONTIG_SAMPLE_HIERARCHY),
                    "sample",
                    contig_profiles,
                    tree,
                    weighting=_hierarchy.WEIGHTING_CONTIG_DEPTH,
                )
                comm(
                    "per-sample hierarchy (contig-depth-weighted) -> %s"
                    % os.path.join(hierarchy_dir, CONTIG_SAMPLE_HIERARCHY)
                )
                comm("this profile covers the whole assembly, binned and unbinned contigs")

        if not args.abundance and not args.gene_abundance and not args.contig_depth:
            comm(
                "no abundance table (-a): per-genome distributions only. A community profile needs "
                "the abundance of each genome."
            )

        if not args.skip_krona:
            if genome_distributions:
                _draw_krona(args, env, genome_distributions, tree, "genome")
            if sample_distributions:
                _draw_krona(args, env, sample_distributions, tree, "sample")
            if gene_distributions:
                _draw_krona(args, env, gene_distributions, tree, "sample_by_gene_abundance")
            if contig_distributions:
                _draw_krona(args, env, contig_distributions, tree, "sample_by_contig_depth")
            comm("krona plots -> %s" % os.path.join(args.output, KRONA_DIR))

        announcement("PATHWAYS PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
