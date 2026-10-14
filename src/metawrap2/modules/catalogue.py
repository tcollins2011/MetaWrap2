"""MetaWrap2 catalogue: a non-redundant gene catalogue, for the functional answers that need no bins.

Most assembled sequence never reaches a MAG - across 92,187 metagenomes about 10% of assembled bases
end up in a medium- or high-quality one - so a pipeline that can only answer functional questions
about bins is answering them about a tenth of the data, and a biased tenth. This is the other path:
call genes on all the contigs, collapse them into one catalogue across however many samples, annotate
the representatives once, and quantify by mapping reads back.

It is a first-class path rather than a degraded mode, and deliberately joinable with the binned one:
same identifiers, same tables, so "what fraction of the community's function did my MAGs capture" is
a question the two outputs can answer together.

Where the cost goes, and why clustering is worth it at all: annotating 5 million genes from ten
samples is ten times the work of annotating the 500,000 distinct ones among them. The catalogue is
the thing that makes the bin-free path affordable, not a tidying step afterwards.

**Gene ids are prefixed with their source before anything else happens.** :mod:`metawrap2.genecall`
names a gene ``<contig>_<n>``, and every assembler numbers contigs from one, so two samples routinely
produce the same gene id. ``mmseqs createdb`` keys on the name, so without prefixing it merges those
genes *before clustering runs* - two identical inputs of 1,772 genes each went in and 1,772 came out,
not 3,544. That is not a tidying detail: it silently drops genes in exactly the multi-sample case
this module exists for. The catalogue therefore carries ``<source>~<gene_id>`` names, the same
convention and separator :mod:`metawrap2.modules.coverage` uses for its contig catalogue, and
``members.tsv`` carries the original gene id and source beside the prefixed one so the join back to
the per-sample tables still works.

**Nucleotide clustering, with the proteins carried along.** The catalogue has two jobs: reads are
mapped against it (which needs nucleotides) and it is annotated (which wants proteins). Clustering
the nucleotide genes and then pulling out the representatives' proteins gives one catalogue that
does both, with a single definition of what a cluster is. Clustering proteins instead would give a
catalogue whose members differ in nucleotide sequence by more than the mapping step assumes.
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import Dict, Iterable, List, Optional, Set, Tuple

from ..cliformat import make_parser
from ..config import MODULE_ENVS, load_settings
from ..constants import TSV_SUFFIX
from ..io.seqio import contig_id, iter_fasta
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

#: MMseqs2 comes with the annotate env already (Bakta pulls it in), so the bin-free path costs no
#: extra install.
CONDA_ENV = MODULE_ENVS["annotate"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# linclust, not cluster: linear-time, which is what makes a ten-sample catalogue tractable. The
# quality difference against `mmseqs cluster` is small at 95% identity, where the clusters are tight
# enough that the greedy incremental algorithm has little room to differ.
#
# --cov-mode 1 is coverage of the *target*, which with --cluster-mode 2 (set-cover, longest first)
# means "the member is covered by the representative". Mode 0, coverage of both, would refuse to
# cluster a gene fragment with the full-length gene it is part of - and a metagenome is full of gene
# fragments, so mode 0 produces a catalogue inflated with pieces of genes it already contains.
#
# --cluster-mode 2 makes the longest sequence the representative, which is the convention the
# published gene catalogues use and the right choice when members may be fragments.
MMSEQS_CREATEDB = "mmseqs createdb {inputs} {db}"
MMSEQS_LINCLUST = (
    "mmseqs linclust {db} {clusters} {tmp} --min-seq-id {identity} -c {coverage} "
    "--cov-mode 1 --cluster-mode 2 --threads {threads}"
)
MMSEQS_TSV = "mmseqs createtsv {db} {db} {clusters} {out}"
MMSEQS_REPSEQ = "mmseqs createsubdb {clusters} {db} {rep_db}"
MMSEQS_CONVERT = "mmseqs convert2fasta {rep_db} {out}"
# Reads against the catalogue, which is what turns a gene list into gene abundance - and is the half
# of the bin-free path that makes it comparable with the binned one.
#
# --no-clip matters here and nowhere else in the pipeline. Salmon trims a trailing run of >=10 A's
# from any reference, because in a transcriptome that is a poly-A tail; in a gene catalogue it is
# sequence. A gene really ending in ten A's would otherwise be indexed shorter than it is, and
# length is the denominator of every unit salmon reports.
SALMON_INDEX = "salmon index -p {threads} -t {catalogue} -i {index} --no-clip"
# --libType A rather than IU: a catalogue is quantified against whatever reads a user has, including
# single-end and stranded libraries, and salmon's auto-detection is better than our assumption.
# --meta switches to plain EM with uniform initialisation, which is the setting for a reference that
# is a mixture of organisms rather than one transcriptome.
SALMON_QUANT = (
    "salmon quant -i {index} --libType A {reads} -o {out} --meta -p {threads} --no-version-check"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_THREADS = 1
#: The conventional thresholds for a metagenomic gene catalogue: 95% nucleotide identity over 90%
#: of the shorter sequence. Stated rather than inherited so a gene count here means the same thing
#: as one in a paper that used them.
IDENTITY = 0.95
COVERAGE = 0.90
#: Separates a source's label from the original gene id in the catalogue. ``~`` for the same reasons
#: coverage.SEPARATOR uses it: legal in a fasta header, meaningless to every aligner, and absent
#: from the gene ids genecall produces. A gene id containing it would make the split cut in the wrong
#: place, so the inputs are checked rather than trusted.
SEPARATOR = "~"
#: Representative genes, nucleotide and protein.
CATALOGUE_FFN = "catalogue.ffn"
CATALOGUE_FAA = "catalogue.faa"
#: gene -> representative, which is what lets per-sample counts be aggregated onto the catalogue.
#: Both the prefixed name and the original (source, gene id) are written, because the catalogue is
#: keyed on the first and every per-sample table is keyed on the second.
MEMBERS_TABLE = "members" + TSV_SUFFIX
MEMBERS_COLUMNS = (
    "representative",
    "member",
    "member_source",
    "member_gene_id",
    "is_representative",
)
#: Gene abundance per sample, long format - the shape `quantify_function`, `diff_abundance` and the
#: BIOM writer all take.
ABUNDANCE_TABLE = "gene_abundance" + TSV_SUFFIX
#: What fraction of each sample's reads the catalogue accounted for. Written beside the abundances
#: because it is the denominator question: a catalogue built from three samples and quantified
#: against a fourth can map very little of it, and the abundances alone do not say so.
MAPPING_TABLE = "gene_abundance_mapping" + TSV_SUFFIX
#: One row per source, for the "was clustering worth it" question.
SUMMARY_TABLE = "catalogue_summary" + TSV_SUFFIX
SUMMARY_COLUMNS = ("source", "genes_in", "representatives_from_this_source")
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 catalogue",
        usage="metawrap2 catalogue [options] -o output_dir (-d annotate_dir ... | -g genes.ffn ...)",
        description="Collapse genes from any number of samples into one non-redundant catalogue.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "-d",
        "--annotate-dir",
        nargs="+",
        help="one or more `metawrap2 annotate` output directories; their genes/*.ffn and the "
        "matching .faa are used",
    )
    source.add_argument(
        "-g", "--genes", nargs="+", help="nucleotide gene files (.ffn) to cluster directly"
    )
    p.add_argument(
        "-p",
        "--proteins",
        nargs="+",
        help="protein files (.faa) matching --genes, so the catalogue carries proteins too",
    )
    threads_arg(p)
    p.add_argument(
        "--identity",
        type=float,
        default=IDENTITY,
        metavar="FRACTION",
        help="minimum sequence identity to cluster, 0-1 (default: %(default)s)",
    )
    p.add_argument(
        "--coverage",
        type=float,
        default=COVERAGE,
        metavar="FRACTION",
        help="minimum alignment coverage of the member, 0-1 (default: %(default)s)",
    )
    p.add_argument(
        "--separator",
        default=SEPARATOR,
        help="separates a source label from the gene id in the catalogue (default %s). Only worth "
        "changing if a gene id contains it, which the run checks for" % SEPARATOR,
    )
    p.add_argument(
        "--quantify",
        nargs="+",
        metavar="SAMPLE=R1[,R2]",
        help="also map reads against the finished catalogue with salmon, one argument per sample: "
        "`gut1=gut1_1.fq,gut1_2.fq`, or `gut1=gut1.fq` for single-end. Writes gene abundance in "
        "the long format quantify_function and diff_abundance read",
    )
    p.add_argument(
        "--unit",
        choices=("tpm", "reads"),
        default="tpm",
        help="which salmon column to report (default: %(default)s). tpm is length-normalised and "
        "comparable between genes; reads is the raw count and is not",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _collect(args) -> Tuple[List[str], List[str]]:
    """(nucleotide files, protein files), from annotate directories or explicit paths."""
    if args.genes:
        return list(args.genes), list(args.proteins or [])
    genes: List[str] = []
    proteins: List[str] = []
    for directory in args.annotate_dir:
        genes_dir = os.path.join(directory, "genes")
        if not os.path.isdir(genes_dir):
            warning("%s has no genes/ directory; skipping it" % directory)
            continue
        for name in sorted(os.listdir(genes_dir)):
            if not name.endswith(".ffn"):
                continue
            genes.append(os.path.join(genes_dir, name))
            protein = os.path.join(genes_dir, name[: -len(".ffn")] + ".faa")
            if os.path.isfile(protein):
                proteins.append(protein)
    return genes, proteins


def source_label(path: str) -> str:
    """The label a file's genes are prefixed with: its basename without the extension."""
    return os.path.splitext(os.path.basename(path))[0]


def write_prefixed(source: str, label: str, out_path: str, separator: str = SEPARATOR) -> int:
    """Copy *source* with every gene id prefixed ``<label><sep>``, returning the count.

    Raises if a gene id already contains the separator, because then the prefix could not be split
    off again unambiguously and every downstream join would be silently wrong for those genes.
    """
    written = 0
    with open(out_path, "w") as out:
        for header, sequence in iter_fasta(source):
            gene = contig_id(header)
            if separator in gene:
                raise ValueError(
                    "gene id %r in %s already contains the separator %r, so a prefix could not be "
                    "removed again. Pass a different --separator." % (gene, source, separator)
                )
            out.write(">%s%s%s\n%s\n" % (label, separator, gene, sequence))
            written += 1
    return written


def split_prefixed(name: str, separator: str = SEPARATOR) -> Tuple[str, str]:
    """(source label, original gene id) from a prefixed catalogue name.

    Splits on the *first* separator: a label cannot contain one (it is a filename stem) and a gene
    id cannot either (write_prefixed refuses), so the first is the right one.
    """
    label, found, gene = name.partition(separator)
    return (label, gene) if found else ("", name)


def parse_members(path: str) -> Dict[str, str]:
    """gene -> representative, from mmseqs' cluster TSV (representative first, member second)."""
    members: Dict[str, str] = {}
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 2:
                continue
            members[row[1]] = row[0]
    return members


def write_subset(source: str, keep: Set[str], out_path: str) -> int:
    """Write the records of *source* whose id is in *keep*, in the source's order."""
    written = 0
    with open(out_path, "w") as out:
        for header, sequence in iter_fasta(source):
            if contig_id(header) in keep:
                out.write(">%s\n%s\n" % (header, sequence))
                written += 1
    return written


def _append_subset(sources: Iterable[str], keep: Set[str], out_path: str) -> int:
    written = 0
    with open(out_path, "w") as out:
        for source in sources:
            for header, sequence in iter_fasta(source):
                if contig_id(header) in keep:
                    out.write(">%s\n%s\n" % (header, sequence))
                    written += 1
    return written


def _parse_quantify(specs: List[str]) -> List[Tuple[str, List[str]]]:
    """``sample=r1[,r2]`` -> (sample, [reads]). Validated here so a typo fails before salmon runs.

    One argument per sample rather than a sample sheet, because a catalogue is usually quantified
    against the same handful of read files it was built from, and a second file format for three
    paths is a cost with no return. A sample name that repeats is refused: it would silently
    overwrite a column of the abundance table.
    """
    out: List[Tuple[str, List[str]]] = []
    seen: Set[str] = set()
    for spec in specs:
        sample, _, paths = spec.partition("=")
        sample = sample.strip()
        reads = [path.strip() for path in paths.split(",") if path.strip()]
        if not sample or not reads:
            error("--quantify takes sample=reads, e.g. gut1=gut1_1.fq,gut1_2.fq (got %r)" % spec)
        if len(reads) > 2:
            error(
                "--quantify takes one or two read files per sample (got %d in %r)"
                % (len(reads), spec)
            )
        if sample in seen:
            error(
                "sample %r is given twice to --quantify; each sample is one column of the "
                "abundance table and two columns cannot share a name" % sample
            )
        seen.add(sample)
        for path in reads:
            require_file(path, "reads for %s (--quantify)" % sample)
        out.append((sample, reads))
    return out


def _quantify(args, settings, catalogue_ffn: str, samples: List[Tuple[str, List[str]]]) -> None:
    """Map each sample's reads against the catalogue and write gene abundance, long format.

    This is the half of the bin-free path that makes it answer the same questions as the binned
    one: `quantify_function` takes exactly this table, so a KO's abundance in the whole sample can
    be computed without any bin existing.

    salmon's mapping rate is written out per sample beside the abundances. It is the bin-free
    equivalent of the binned fraction `bin_refinement` reports, and it is the number that says
    whether the catalogue describes the sample it was quantified against or only part of it.
    """
    # salmon lives with the binning tools, not with the annotate env this module otherwise uses -
    # the same arrangement `viruses` has for skani, and it saves a second copy of the tool.
    env = env_for("quant_bins", settings)
    work = os.path.join(args.output, "work_files", "salmon")
    ensure_dir(work)
    index = os.path.join(work, "index")
    rc = run(
        SALMON_INDEX.format(threads=args.threads, catalogue=catalogue_ffn, index=index),
        env=env,
        tool="salmon index",
        check=False,
    )
    if rc != 0:
        warning("salmon index failed, so the catalogue is not quantified")
        return

    values: Dict[str, Dict[str, float]] = {}
    mapping: Dict[str, Optional[float]] = {}
    for sample, reads in samples:
        out_dir = os.path.join(work, sample)
        flags = "-1 %s -2 %s" % (reads[0], reads[1]) if len(reads) == 2 else "-r %s" % reads[0]
        rc = run(
            SALMON_QUANT.format(index=index, reads=flags, out=out_dir, threads=args.threads),
            env=env,
            tool="salmon quant",
            check=False,
        )
        quant = os.path.join(out_dir, "quant.sf")
        if rc != 0 or not os.path.isfile(quant):
            warning(
                "salmon failed on %s; that sample is absent from %s rather than present as zeroes"
                % (sample, ABUNDANCE_TABLE)
            )
            continue
        values[sample] = _read_quant(quant, args.unit)
        mapping[sample] = _mapping_rate(out_dir)

    if not values:
        warning("no sample quantified; no abundance table written")
        return

    table = os.path.join(args.output, ABUNDANCE_TABLE)
    with open(table, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "feature", "value", "unit"])
        for sample in sorted(values):
            for gene, value in sorted(values[sample].items()):
                if value:
                    writer.writerow([sample, gene, "%.6g" % value, args.unit])
    comm("gene abundance (%s) -> %s" % (args.unit, table))

    rates = os.path.join(args.output, MAPPING_TABLE)
    with open(rates, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "percent_mapped", "genes_detected", "genes_in_catalogue"])
        for sample in sorted(values):
            rate = mapping.get(sample)
            writer.writerow(
                [
                    sample,
                    # Blank, not zero, when salmon did not report one: a missing mapping rate and a
                    # rate of nothing mapping are different findings.
                    "" if rate is None else "%.2f" % rate,
                    sum(1 for value in values[sample].values() if value),
                    len(values[sample]),
                ]
            )
    comm("mapping rate per sample -> %s" % rates)
    for sample in sorted(mapping):
        rate = mapping[sample]
        if rate is not None and rate < 50:
            warning(
                "%s: only %.1f%% of its reads mapped to the catalogue, so these abundances "
                "describe a minority of that sample. A catalogue describes the samples it was "
                "built from; quantifying a different one against it measures the overlap"
                % (sample, rate)
            )


def _read_quant(path: str, unit: str) -> Dict[str, float]:
    """gene -> TPM or read count, from salmon's quant.sf."""
    column = "TPM" if unit == "tpm" else "NumReads"
    out: Dict[str, float] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                out[row["Name"]] = float(row[column])
            except (KeyError, TypeError, ValueError):
                continue
    return out


def _mapping_rate(out_dir: str) -> Optional[float]:
    """Percent of a sample's reads salmon mapped, from its meta_info.json. None if unavailable.

    Read out of salmon's own record rather than computed from the table: a read mapping to several
    genes contributes to several rows, so summing the table would count it more than once and give
    a rate above 100.
    """
    import json

    path = os.path.join(out_dir, "aux_info", "meta_info.json")
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as fh:
            info = json.load(fh)
    except (OSError, ValueError):
        return None
    rate = info.get("percent_mapped")
    return float(rate) if isinstance(rate, (int, float)) else None


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)

    genes, proteins = _collect(args)
    if not genes:
        error("no gene files found. Pass -g, or -d pointing at `annotate` output directories")
    validate_inputs("catalogue", [])
    if not 0 < args.identity <= 1 or not 0 < args.coverage <= 1:
        error("--identity and --coverage are fractions between 0 and 1")

    env = env_for("annotate", settings)
    rec = start_run("catalogue", args, env, settings, inputs=list(genes))

    ckpt = make_checkpoint(args.output)
    try:
        announcement("BEGIN CATALOGUE PIPELINE!")
        ensure_dir(args.output)
        work = os.path.join(args.output, "work_files")
        ensure_dir(work)
        comm(
            "%d gene file(s), %d protein file(s); clustering at %.0f%% identity over %.0f%% of the "
            "member" % (len(genes), len(proteins), 100 * args.identity, 100 * args.coverage)
        )

        db = os.path.join(work, "genes_db")
        clusters = os.path.join(work, "clusters")
        tmp = os.path.join(work, "tmp")
        raw_members = os.path.join(work, "clusters" + TSV_SUFFIX)
        members_tsv = os.path.join(args.output, MEMBERS_TABLE)

        # Prefix before anything touches mmseqs: gene ids collide across assemblies, and createdb
        # keys on the name, so unprefixed inputs lose genes before clustering even starts.
        prefixed_dir = os.path.join(work, "prefixed")
        ensure_dir(prefixed_dir)
        prefixed_genes: List[str] = []
        labels: Dict[str, str] = {}
        if not dry_run():
            total_in = 0
            for path in genes:
                label = source_label(path)
                if label in labels:
                    error(
                        "two inputs are both called %r (%s and %s); their genes would be "
                        "indistinguishable in the catalogue" % (label, labels[label], path)
                    )
                labels[label] = path
                out_path = os.path.join(prefixed_dir, label + ".ffn")
                try:
                    total_in += write_prefixed(path, label, out_path, args.separator)
                except ValueError as problem:
                    error(str(problem))
                prefixed_genes.append(out_path)
            comm("%d gene(s) in, prefixed with %d source label(s)" % (total_in, len(labels)))
            prefixed_proteins: List[str] = []
            for path in proteins:
                label = source_label(path)
                out_path = os.path.join(prefixed_dir, label + ".faa")
                try:
                    write_prefixed(path, label, out_path, args.separator)
                except ValueError as problem:
                    error(str(problem))
                prefixed_proteins.append(out_path)
        else:
            prefixed_genes, prefixed_proteins = list(genes), list(proteins)

        if ckpt.todo("cluster"):
            for command, tool in (
                (
                    MMSEQS_CREATEDB.format(inputs=" ".join(prefixed_genes), db=db),
                    "mmseqs createdb",
                ),
                (
                    MMSEQS_LINCLUST.format(
                        db=db,
                        clusters=clusters,
                        tmp=tmp,
                        identity=args.identity,
                        coverage=args.coverage,
                        threads=args.threads,
                    ),
                    "mmseqs linclust",
                ),
                (MMSEQS_TSV.format(db=db, clusters=clusters, out=raw_members), "mmseqs createtsv"),
            ):
                if run(command, env=env, tool=tool, check=False) != 0:
                    error("%s failed; the catalogue cannot be built without it" % tool)
            ckpt.done("cluster")
        else:
            comm("skipping clustering (already done; --resume)")

        if dry_run():
            announcement("CATALOGUE PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        members = parse_members(raw_members)
        representatives = set(members.values())
        if not representatives:
            error("clustering produced no representatives; see %s" % raw_members)

        with open(members_tsv, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=MEMBERS_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for member, representative in sorted(members.items()):
                label, gene = split_prefixed(member, args.separator)
                writer.writerow(
                    {
                        "representative": representative,
                        "member": member,
                        "member_source": label,
                        "member_gene_id": gene,
                        "is_representative": int(member == representative),
                    }
                )
        comm(
            "%d genes -> %d representatives (%.1f%% of the redundancy removed)"
            % (
                len(members),
                len(representatives),
                100.0 * (1 - len(representatives) / float(len(members))) if members else 0.0,
            )
        )

        ffn = os.path.join(args.output, CATALOGUE_FFN)
        written = _append_subset(prefixed_genes, representatives, ffn)
        comm("%d representative gene(s) -> %s" % (written, ffn))
        if written != len(representatives):
            warning(
                "wrote %d of %d representatives: %d id(s) in the cluster table were not found in "
                "the input files, which usually means --genes and the cluster table are from "
                "different runs" % (written, len(representatives), len(representatives) - written)
            )

        if proteins:
            faa = os.path.join(args.output, CATALOGUE_FAA)
            n_proteins = _append_subset(prefixed_proteins, representatives, faa)
            comm("%d representative protein(s) -> %s" % (n_proteins, faa))
        else:
            comm(
                "no protein files given (-p, or -d which finds them): the catalogue has nucleotides "
                "only, so annotating it means calling genes on it again"
            )

        # From the prefix rather than by re-reading the inputs: the prefix is unambiguous, which is
        # the whole reason it exists, so no gene can be credited to the wrong source.
        summary = os.path.join(args.output, SUMMARY_TABLE)
        per_source_in: Dict[str, int] = {}
        per_source_rep: Dict[str, int] = {}
        for member in members:
            label, _gene = split_prefixed(member, args.separator)
            per_source_in[label] = per_source_in.get(label, 0) + 1
            if member in representatives:
                per_source_rep[label] = per_source_rep.get(label, 0) + 1
        with open(summary, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for label in sorted(per_source_in):
                writer.writerow(
                    {
                        "source": label,
                        "genes_in": per_source_in[label],
                        "representatives_from_this_source": per_source_rep.get(label, 0),
                    }
                )
        comm("per-source summary -> %s" % summary)
        if args.quantify:
            _quantify(args, settings, ffn, _parse_quantify(args.quantify))
        else:
            comm(
                "annotate the catalogue once (`metawrap2 annotate -a %s` or --kofam on the .faa); "
                "--quantify sample=reads maps reads against it for gene abundance" % ffn
            )
        announcement("CATALOGUE PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
