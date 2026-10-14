"""MetaWrap2 annotate: call genes once over a bin set or a whole assembly, then annotate them.

This replaces ``annotate_bins``, which shelled out to one Prokka run per bin and reorganised the
output into three folders. Two things were wrong with that, and neither is fixed by swapping the
annotator:

* **It produced nothing joinable.** A GFF per bin with product names in it answers "what is this
  gene called" and no other question. Pathway completeness needs KEGG orthologs; AMR needs a gene's
  neighbourhood; abundance needs a stable gene identifier that survives into a read-mapping step.
  None of those can be recovered from a folder of per-bin GFFs.
* **It had no answer for a user who did not bin.** Most assembled sequence never reaches a MAG, and
  ``-b bins`` was the only input the module took.

So gene calling happens once, here, in :mod:`metawrap2.genecall`, and it takes either a bin folder
(``-b``) or an assembly (``-a``). Everything after it - orthologs, domains, pathways, resistance -
reads the gene table this writes rather than calling genes again. The cost argument is the smaller
half of that: the reason is that two callers disagreeing by one gene produce two answers that cannot
be joined.

This file is meant to be read and edited: the commands each external tool runs live in the COMMANDS
block below, and the tunable numbers in the CONSTANTS block.
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import Dict, List, Optional

from .. import bgc, cazy, genecall, kofam, pfam
from ..cliformat import make_parser
from ..config import MODULE_ENVS, load_settings
from ..constants import FASTA_EXTENSIONS, TSV_SUFFIX
from ..io.seqio import contig_id, iter_fasta
from ..progress import bar
from ._common import (
    absolutize_paths,
    announcement,
    check_fasta,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    list_bins,
    make_checkpoint,
    require_nonempty_dir,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["annotate_bins"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# Bakta is given the genes this module already called (--skip-cds with a --regions file is not a
# stable interface across Bakta versions), so for now it re-predicts them and is used for what it is
# good at: naming products and finding the non-coding features - tRNAs via tRNAscan-SE, rRNAs via
# Infernal/Rfam - that Prodigal does not look for. Those are what MIMAG's rRNA/tRNA criteria need.
BAKTA = (
    "bakta --threads {threads} --output {outdir} --prefix {prefix} --force "
    "--keep-contig-headers {db} {input}"
)
# CAZyme families, over the proteins this module has already called. --mode protein for that reason:
# dbCAN would otherwise call its own genes with Prodigal, and two gene sets over one genome cannot be
# joined - a CAZyme call would then belong to a gene that is in no other table here.
#
# All three methods run (the default), because the disagreement between them is the output's most
# useful column: dbCAN's own guidance is to trust a family found by at least two, and a third of the
# calls on one real bin were found by one.
DBCAN = (
    "run_dbcan CAZyme_annotation --mode protein --input_raw_data {proteins}"
    " --output_dir {outdir} --db_dir {db} --threads {threads}"
)
# Biosynthetic gene clusters, over the genes this module has already called - hence --gff. GECCO
# calls its own genes otherwise, and while its Pyrodigal names happen to match ours, "happens to
# match" is not a join: a BGC's member genes have to be the genes in genes.tsv.
#
# --force-tsv because a genome with no clusters otherwise writes no table at all, and "no clusters"
# would be indistinguishable from "the run failed" - which for most genomes is the common case.
GECCO = (
    "gecco run -g {genome} --gff {gff} -o {outdir} -j {threads} --force-tsv --no-color --no-markup"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_THREADS = 1
#: Where per-genome gene files go, one set per genome, named for the genome.
GENES_DIR = "genes"
#: The per-gene table every downstream module reads. Long, one row per gene.
GENES_TABLE = "genes" + TSV_SUFFIX
#: The per-gene KO table, written by --kofam. One row per accepted assignment, carrying how it
#: qualified - see metawrap2.kofam on why "how" is not an implementation detail.
KOS_TABLE = "kos" + TSV_SUFFIX
#: One row per Pfam domain found, written by --pfam. Per domain rather than per gene: a protein can
#: carry several and the same one twice, and an architecture is only recoverable if each is a row.
PFAM_TABLE = "pfam_domains" + TSV_SUFFIX
#: One row per gene with its ordered domain architecture.
PFAM_ARCHITECTURE_TABLE = "pfam_architectures" + TSV_SUFFIX
#: The per-gene CAZyme table, written by --cazy. One row per gene with a CAZyme call, carrying which
#: of dbCAN's three methods found it.
CAZY_TABLE = "cazymes" + TSV_SUFFIX
#: Per-class counts per genome: the summary a reader actually quotes.
CAZY_CLASS_TABLE = "cazyme_classes" + TSV_SUFFIX
#: Where dbCAN's own output is kept, one directory per genome.
CAZY_DIR = "dbcan"
#: One row per predicted biosynthetic gene cluster, written by --bgc.
BGC_TABLE = "gene_clusters" + TSV_SUFFIX
#: Per-genome cluster counts, split by whether a cluster runs into unassembled sequence.
BGC_SUMMARY_TABLE = "gene_cluster_summary" + TSV_SUFFIX
#: Where GECCO's own output is kept, one directory per genome - including the per-cluster GenBank
#: files, which are what an antiSMASH sideload or a manual look needs.
BGC_DIR = "gecco"
#: Per-genome summary: how the genes were called and what they look like in aggregate. This is the
#: file to look at when a genome's functional annotation seems wrong, because the usual cause is
#: visible here - a coding density far off 0.85-0.90, or a bin that fell back to meta mode.
SUMMARY_TABLE = "genome_summary" + TSV_SUFFIX
#: Columns of SUMMARY_TABLE.
SUMMARY_COLUMNS = (
    "genome",
    "n_contigs",
    "total_bp",
    "n_genes",
    "n_partial_genes",
    "coding_density",
    "calling_mode",
)
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 annotate",
        usage="metawrap2 annotate [options] -o output_dir (-b bin_folder | -a assembly.fasta)",
        description="Call genes once over a bin set or an assembly, and annotate them.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "-b", "--bins", help="folder of bins in fasta format; each is annotated as one genome"
    )
    source.add_argument(
        "-a",
        "--assembly",
        help="a whole assembly, annotated without bins (genes are called in metagenomic mode)",
    )
    threads_arg(p)
    p.add_argument(
        "--bakta",
        action="store_true",
        help="also run Bakta for product names and non-coding features (needs a Bakta database)",
    )
    p.add_argument(
        "--bakta-db",
        help="path to the Bakta database directory (default: BAKTA_DB from metawrap2.toml)",
    )
    p.add_argument(
        "--kofam",
        action="store_true",
        help="assign KEGG orthologs to the called genes with PyHMMER over the KOfam profiles "
        "(needs KOFAM_DB)",
    )
    p.add_argument(
        "--kofam-db",
        help="path to the KOfam HMM profile directory (default: KOFAM_DB from metawrap2.toml)",
    )
    p.add_argument(
        "--kofam-relative-fraction",
        type=float,
        default=None,
        metavar="FRACTION",
        help="accept an otherwise-unannotated gene's best hit when it reaches this fraction of the "
        "profile's curated threshold, flagged as evidence=relative (default: off, i.e. curated "
        "thresholds only, which reproduces KofamScan)",
    )
    p.add_argument(
        "--kofam-thresholdless-evalue",
        type=float,
        default=None,
        metavar="E",
        help="also accept hits to the 1866 profiles KOfam publishes no threshold for, when their "
        "E-value is at most E, flagged as evidence=no_threshold (default: off - see the module "
        "docstring on why that default is strict)",
    )
    p.add_argument(
        "--pfam",
        action="store_true",
        help="also assign Pfam domains to every called protein, at each family's own gathering "
        "threshold. This is the annotation for the genes KEGG cannot place - on real bins only 41%% "
        "of genes take a curated KO - and it needs PFAM_DB (`metawrap2 install-db pfam`, 399 MB)",
    )
    p.add_argument(
        "--pfam-db",
        metavar="DIR",
        help="directory holding Pfam-A.hmm (default: PFAM_DB from metawrap2.toml)",
    )
    p.add_argument(
        "--cazy",
        action="store_true",
        help="also assign CAZyme families to the called proteins with dbCAN (needs DBCAN_DB). "
        "CAZy answers which glycosidic bond an enzyme acts on, which a KO table does not",
    )
    p.add_argument(
        "--cazy-db",
        help="path to the dbCAN database directory (default: DBCAN_DB from metawrap2.toml)",
    )
    p.add_argument(
        "--cazy-methods",
        default=None,
        metavar="LIST",
        help="which of dbCAN's methods to run, comma-separated from diamond,hmm,dbCANsub "
        "(default: all three). Dropping one makes every call look less supported, since the "
        "number of agreeing methods is reported per gene",
    )
    p.add_argument(
        "--bgc",
        action="store_true",
        help="also predict biosynthetic gene clusters with GECCO over the called genes. Needs no "
        "database. A cluster is defined by gene adjacency, so the table says how close each one "
        "sits to a contig edge",
    )
    p.add_argument(
        "--meta",
        action="store_true",
        help="force metagenomic gene calling even for bins (default: train a model per bin when it "
        "has at least %d bp of sequence)" % genecall.MIN_TRAIN_BP,
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _inputs(args) -> List[tuple]:
    """[(genome name, fasta path, meta mode)], in the order they will be annotated.

    One entry per bin, or a single entry for an assembly. An assembly is called in metagenomic mode
    whatever ``--meta`` says, because training one coding model on many organisms trains it on an
    average of them - see genecall for why that is the one choice not left to the user.
    """
    if args.assembly:
        name = os.path.splitext(os.path.basename(args.assembly))[0]
        return [(name, args.assembly, True)]
    names = list_bins(args.bins)
    return [(os.path.splitext(name)[0], os.path.join(args.bins, name), args.meta) for name in names]


def _call_and_write(name: str, path: str, meta: bool, genes_dir: str) -> genecall.GenomeGenes:
    called = genecall.call_genome(path, name=name, meta=meta)
    genecall.write_faa(called, os.path.join(genes_dir, name + ".faa"))
    genecall.write_ffn(called, os.path.join(genes_dir, name + ".ffn"))
    genecall.write_gff(called, os.path.join(genes_dir, name + ".gff"))
    return called


def _bakta_db(args, settings) -> Optional[str]:
    db = args.bakta_db or settings.db("BAKTA_DB")
    if args.bakta and not db:
        comm("No Bakta database given (--bakta-db / BAKTA_DB); relying on Bakta's own default.")
    return db


def _kofam_profiles(args, settings) -> List[str]:
    """The HMM files to search, from a directory of them or a single concatenated file.

    KOfam ships one ``.hmm`` per KO - 28,453 of them - which is a directory, not a file. Both
    layouts are accepted because a user who has concatenated them into one file (which hmmsearch
    and PyHMMER both read happily, and which is faster to open) should not have to undo that.
    """
    db = args.kofam_db or settings.db("KOFAM_DB")
    if not db:
        error(
            "--kofam needs the KOfam profiles and thresholds: run `metawrap2 install-db kofam`, "
            "or pass --kofam-db or set KOFAM_DB in metawrap2.toml"
        )
    if os.path.isfile(db):
        return [db]
    if not os.path.isdir(db):
        error("KOfam profile path is neither a file nor a directory: %s" % db)
    # KOfam's tarball extracts to profiles/, sometimes nested one level down.
    for candidate in (db, os.path.join(db, "profiles")):
        if os.path.isdir(candidate):
            found = sorted(
                os.path.join(candidate, name)
                for name in os.listdir(candidate)
                if name.endswith((".hmm", ".hmm.gz"))
            )
            if found:
                return found
    error("no .hmm profiles found under %s" % db)
    return []


def _ko_genes(args) -> Dict[str, set]:
    """genome -> the genes that already have a KEGG ortholog, from this run's own KO table.

    Read off disk rather than threaded through main, which also makes it work on a `--resume` where
    the KOfam pass ran in an earlier invocation. Empty when --kofam was not used, and the Pfam
    report then omits the increment rather than claiming every domain is new information.
    """
    path = os.path.join(args.output, KOS_TABLE)
    if not os.path.isfile(path):
        return {}
    out: Dict[str, set] = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            genome = (row.get("genome") or "").strip()
            gene = (row.get("gene_id") or "").strip()
            if genome and gene:
                out.setdefault(genome, set()).add(gene)
    return out


def _run_pfam(args, genes_dir: str, todo: List[tuple], ko_genes_by_genome=None) -> None:
    """Search every genome's proteins against Pfam-A and write the domain tables.

    In-process through PyHMMER, like the KOfam search, so there is no second way of searching HMMs
    in this package and no subprocess to wrap. Per genome for the reason the CAZyme pass is: every
    table here is keyed on (genome, gene), and two genomes sharing a gene id would merge.
    """
    settings = load_settings(args.config)
    db = args.pfam_db or settings.db("PFAM_DB")
    if not db:
        error(
            "--pfam needs Pfam-A: pass --pfam-db or set PFAM_DB in metawrap2.toml "
            "(`metawrap2 install-db pfam`, 399 MB)"
        )
    hmm = db if db.endswith((".hmm", ".hmm.gz")) else os.path.join(db, "Pfam-A.hmm")
    if not os.path.isfile(hmm):
        error("no Pfam-A.hmm under %s" % db)
    release = pfam.release_from_notes(db)
    comm("Pfam release %s, calling domains at each family's gathering threshold" % release)
    if dry_run():
        comm("would assign Pfam domains to %d genome(s)" % len(todo))
        return

    domain_path = os.path.join(args.output, PFAM_TABLE)
    arch_path = os.path.join(args.output, PFAM_ARCHITECTURE_TABLE)
    totals = {"genes": 0, "annotated": 0, "domains": 0, "rescued": 0}
    with (
        open(domain_path, "w", newline="") as domain_fh,
        open(arch_path, "w", newline="") as arch_fh,
    ):
        domain_writer = csv.DictWriter(
            domain_fh, fieldnames=pfam.DOMAIN_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        domain_writer.writeheader()
        arch_writer = csv.DictWriter(
            arch_fh, fieldnames=pfam.ARCHITECTURE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        arch_writer.writeheader()
        for name, _path, _meta in bar(todo, desc="Pfam", unit=" genome", total=len(todo)):
            faa = os.path.join(genes_dir, name + ".faa")
            if not os.path.isfile(faa):
                warning("no proteins for %s; skipping Pfam annotation" % name)
                continue
            proteins = [(contig_id(header), sequence) for header, sequence in iter_fasta(faa)]
            if not proteins:
                continue
            try:
                hits = pfam.search(proteins, hmm, threads=args.threads)
            except Exception as problem:  # noqa: BLE001 - one genome failing must not lose the rest
                warning("Pfam search failed on %s: %s" % (name, problem))
                continue
            for hit in hits:
                domain_writer.writerow(pfam.domain_row(hit, name))
            for row in pfam.architecture_rows(hits, name):
                arch_writer.writerow(row)
            known = (ko_genes_by_genome or {}).get(name)
            found = pfam.coverage(hits, [gene for gene, _ in proteins], known)
            totals["genes"] += found.genes
            totals["annotated"] += found.genes_with_a_domain
            totals["domains"] += found.domains
            if found.genes_rescued is not None:
                totals["rescued"] += found.genes_rescued

    share = 100.0 * totals["annotated"] / totals["genes"] if totals["genes"] else 0.0
    comm(
        "%d domain(s) over %d of %d gene(s) (%.1f%%) -> %s"
        % (totals["domains"], totals["annotated"], totals["genes"], share, domain_path)
    )
    if ko_genes_by_genome:
        comm(
            "%d gene(s) have a Pfam domain and no KEGG ortholog - the increment the extra search "
            "bought" % totals["rescued"]
        )
        comm(pfam.INCREMENT_CAVEAT)
    comm("-> %s" % arch_path)
    comm(pfam.DOMAIN_CAVEAT)


def _run_cazy(args, genes_dir: str, todo: List[tuple]) -> None:
    """Run dbCAN over each genome's proteins and write the CAZyme tables.

    Per genome rather than over a concatenation, because every table in this module is keyed on
    (genome, gene) and dbCAN's output is keyed on gene alone - two genomes sharing a gene id would
    merge silently, which is the same defect the catalogue module prefixes its ids to avoid.
    """
    settings = load_settings(args.config)
    # Its own lookup rather than main's `env`, which is None unless --bakta was given: dbCAN lives
    # in the annotate env whether or not Bakta is being run, and passing None here would reach for
    # run_dbcan on the user's PATH - the mistake `check` used to make.
    env = env_for("annotate", settings)
    db = args.cazy_db or settings.db("DBCAN_DB")
    if not db:
        error(
            "--cazy needs the dbCAN database: pass --cazy-db or set DBCAN_DB in metawrap2.toml "
            "(`metawrap2 install-db dbcan`, 7.4 GB)"
        )
    if dry_run():
        comm("would assign CAZyme families to %d genome(s)" % len(todo))
        return

    root = os.path.join(args.output, CAZY_DIR)
    ensure_dir(root)
    gene_path = os.path.join(args.output, CAZY_TABLE)
    class_path = os.path.join(args.output, CAZY_CLASS_TABLE)
    methods = " --methods %s" % args.cazy_methods if args.cazy_methods else ""
    total = confident = 0
    with open(gene_path, "w", newline="") as gene_fh, open(class_path, "w", newline="") as class_fh:
        gene_writer = csv.DictWriter(
            gene_fh, fieldnames=cazy.CAZYME_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        gene_writer.writeheader()
        class_writer = csv.DictWriter(
            class_fh, fieldnames=cazy.CLASS_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        class_writer.writeheader()
        for name, _path, _meta in bar(todo, desc="CAZymes", unit=" genome", total=len(todo)):
            faa = os.path.join(genes_dir, name + ".faa")
            if not os.path.isfile(faa):
                warning("no proteins for %s; skipping CAZyme annotation" % name)
                continue
            outdir = os.path.join(root, name)
            rc = run(
                DBCAN.format(proteins=faa, outdir=outdir, db=db, threads=args.threads) + methods,
                env=env,
                tool="run_dbcan",
                check=False,
            )
            overview = os.path.join(outdir, "overview.tsv")
            if rc != 0 or not os.path.isfile(overview):
                warning("dbCAN failed on %s; it contributes no CAZyme calls" % name)
                continue
            calls = cazy.parse_overview(overview, name)
            for call in calls:
                gene_writer.writerow(cazy.cazyme_row(call))
            for summary in cazy.summarise_classes(calls, name):
                class_writer.writerow(cazy.class_row(summary))
            total += len(calls)
            confident += sum(1 for call in calls if call.is_confident)

    comm(
        "%d CAZyme gene(s), %d of them found by at least %d of dbCAN's three methods -> %s"
        % (total, confident, cazy.MIN_TOOLS, gene_path)
    )
    comm(
        "a CBM is a binding module, not an enzyme: %s separates it from the catalytic classes "
        "rather than summing them" % os.path.basename(class_path)
    )


def _run_bgc(args, genes_dir: str, todo: List[tuple]) -> None:
    """Predict gene clusters per genome and write both BGC tables.

    Per genome, like the CAZyme pass and for the same reason: every table here is keyed on (genome,
    gene), and GECCO's output is keyed on contig and gene alone.
    """
    settings = load_settings(args.config)
    env = env_for("annotate", settings)
    if dry_run():
        comm("would predict gene clusters in %d genome(s)" % len(todo))
        return

    root = os.path.join(args.output, BGC_DIR)
    ensure_dir(root)
    cluster_path = os.path.join(args.output, BGC_TABLE)
    summary_path = os.path.join(args.output, BGC_SUMMARY_TABLE)
    total = truncated = 0
    with (
        open(cluster_path, "w", newline="") as cluster_fh,
        open(summary_path, "w", newline="") as summary_fh,
    ):
        cluster_writer = csv.DictWriter(
            cluster_fh, fieldnames=bgc.CLUSTER_TABLE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        cluster_writer.writeheader()
        summary_writer = csv.DictWriter(
            summary_fh, fieldnames=bgc.SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        summary_writer.writeheader()
        for name, path, _meta in bar(todo, desc="gene clusters", unit=" genome", total=len(todo)):
            gff = os.path.join(genes_dir, name + ".gff")
            if not os.path.isfile(gff):
                warning("no GFF for %s; skipping cluster prediction" % name)
                continue
            outdir = os.path.join(root, name)
            rc = run(
                GECCO.format(genome=path, gff=gff, outdir=outdir, threads=args.threads),
                env=env,
                tool="gecco run",
                check=False,
            )
            table = os.path.join(outdir, "%s.clusters.tsv" % os.path.splitext(name)[0])
            if not os.path.isfile(table):
                # GECCO names its outputs after the input file's stem, which is not always the
                # genome name (a bin called `bin.1.fa` has stem `bin.1`). Fall back to whatever
                # cluster table it wrote rather than guessing again.
                found = sorted(
                    os.path.join(outdir, entry)
                    for entry in (os.listdir(outdir) if os.path.isdir(outdir) else [])
                    if entry.endswith(".clusters.tsv")
                )
                table = found[0] if found else ""
            if rc != 0 or not table:
                warning("gecco failed on %s; it contributes no clusters" % name)
                continue
            lengths = {contig_id(header): len(seq) for header, seq in iter_fasta(path)}
            clusters = bgc.parse_clusters(table, name, lengths)
            for cluster in clusters:
                cluster_writer.writerow(bgc.cluster_row(cluster))
            summary_writer.writerow(bgc.summary_row(bgc.summarise(clusters, name)))
            total += len(clusters)
            truncated += sum(1 for cluster in clusters if cluster.is_truncated)

    comm("%d gene cluster(s), %d of them at a contig edge -> %s" % (total, truncated, cluster_path))
    if total:
        comm(bgc.TRUNCATION_CAVEAT)


def _run_kofam(args, genes_dir: str, todo: List[tuple]) -> None:
    """Search every genome's proteins against the KOfam profiles and write the KO table."""
    settings = load_settings(args.config)
    profile_files = _kofam_profiles(args, settings)
    policy = kofam.ThresholdPolicy(
        relative_fraction=args.kofam_relative_fraction,
        thresholdless_max_evalue=args.kofam_thresholdless_evalue,
    )
    comm("%d KOfam profile file(s); policy: %s" % (len(profile_files), policy.describe()))
    if dry_run():
        comm("would assign KOs to %d genome(s)" % len(todo))
        return

    db = args.kofam_db or settings.db("KOFAM_DB")
    ko_list = kofam.find_ko_list(db)
    if ko_list is None:
        error(
            "KOfam's thresholds (ko_list) were not found next to %s. Install the profiles and "
            "thresholds together with `metawrap2 install-db kofam`." % db
        )
    if kofam.ko_list_md5(ko_list) != kofam.KO_LIST_MD5:
        warning(
            "%s differs from the KOfam threshold list this release was validated on. KO calls "
            "will follow its cutoffs; each row records which list was used." % ko_list
        )
    profiles = kofam.load_ko_list(ko_list)
    ko_list_tag = kofam.ko_list_label(ko_list)
    out_path = os.path.join(args.output, KOS_TABLE)
    counts = {evidence: 0 for evidence in kofam.Evidence}
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=kofam.KO_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for name, _path, _meta in bar(todo, desc="assigning KOs", unit=" genome", total=len(todo)):
            faa = os.path.join(genes_dir, name + ".faa")
            if not os.path.isfile(faa):
                warning("no proteins for %s; skipping KO assignment" % name)
                continue
            hits = kofam.search(faa, profile_files, threads=args.threads)
            for call in kofam.assign(hits, profiles, policy):
                writer.writerow(kofam.ko_row(name, call, policy, ko_list_tag))
                counts[call.evidence] += 1
    total = sum(counts.values())
    comm(
        "assigned %d KOs (%s)"
        % (
            total,
            ", ".join("%d %s" % (n, e.value) for e, n in counts.items() if n),
        )
    )
    comm("per-gene KO table: %s" % out_path)


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)

    if args.bins:
        checks = [
            (
                lambda path, what: require_nonempty_dir(path, what, FASTA_EXTENSIONS),
                args.bins,
                "bin folder (-b)",
            )
        ]
        source = args.bins
    else:
        checks = [(check_fasta, args.assembly, "assembly (-a)")]
        source = args.assembly
    validate_inputs("annotate", checks)

    env = env_for("annotate_bins", settings) if args.bakta else None
    bakta_db = _bakta_db(args, settings)
    rec = start_run("annotate", args, env, settings, inputs=[source])

    ckpt = make_checkpoint(args.output)
    try:
        announcement("BEGIN ANNOTATION PIPELINE!")
        ensure_dir(args.output)
        genes_dir = os.path.join(args.output, GENES_DIR)
        ensure_dir(genes_dir)

        todo = _inputs(args)
        if not todo:
            error("No fasta input found (%s). Exiting." % "/".join(FASTA_EXTENSIONS))
        comm("%d genome(s) to annotate" % len(todo))

        if ckpt.todo("genecall"):
            if dry_run():
                comm("would call genes in %d genome(s) with Pyrodigal" % len(todo))
            else:
                genes_table = os.path.join(args.output, GENES_TABLE)
                summary_table = os.path.join(args.output, SUMMARY_TABLE)
                with (
                    open(genes_table, "w", newline="") as gfh,
                    open(summary_table, "w", newline="") as sfh,
                ):
                    genes_out = csv.DictWriter(
                        gfh, fieldnames=genecall.GENE_COLUMNS, delimiter="\t", lineterminator="\n"
                    )
                    genes_out.writeheader()
                    summary_out = csv.DictWriter(
                        sfh, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n"
                    )
                    summary_out.writeheader()

                    n_genes = 0
                    fell_back: List[str] = []
                    for name, path, meta in bar(
                        todo, desc="calling genes", unit=" genome", total=len(todo)
                    ):
                        called = _call_and_write(name, path, meta, genes_dir)
                        genes_out.writerows(genecall.gene_rows(called))
                        summary_out.writerow(
                            {
                                "genome": called.name,
                                "n_contigs": called.n_contigs,
                                "total_bp": called.total_bp,
                                "n_genes": len(called.genes),
                                "n_partial_genes": called.n_partial,
                                "coding_density": "%.4f" % called.coding_density,
                                "calling_mode": called.mode,
                            }
                        )
                        n_genes += len(called.genes)
                        if called.mode == "meta" and not meta:
                            fell_back.append(called.name)

                comm("called %d genes in %d genome(s)" % (n_genes, len(todo)))
                if fell_back:
                    # Not a warning about a failure - a warning that these genomes were annotated
                    # by a different model than the rest, which matters when their gene counts are
                    # compared with their neighbours'.
                    warning(
                        "%d genome(s) had under %d bp and were called in metagenomic mode rather "
                        "than with their own coding model: %s"
                        % (
                            len(fell_back),
                            genecall.MIN_TRAIN_BP,
                            ", ".join(fell_back[:5]) + (" ..." if len(fell_back) > 5 else ""),
                        )
                    )
            ckpt.done("genecall")
        else:
            comm("skipping gene calling (already done; --resume)")

        if args.bakta:
            if ckpt.todo("bakta"):
                bakta_out = os.path.join(args.output, "bakta")
                ensure_dir(bakta_out)
                db = "--db %s" % bakta_db if bakta_db else ""
                for name, path, _meta in bar(todo, desc="bakta", unit=" genome", total=len(todo)):
                    rc = run(
                        BAKTA.format(
                            threads=args.threads,
                            outdir=os.path.join(bakta_out, name),
                            prefix=name,
                            db=db,
                            input=path,
                        ),
                        env=env,
                        tool="bakta",
                        check=False,
                    )
                    if rc != 0:
                        warning("Bakta failed on %s. Proceeding with the rest." % name)
                ckpt.done("bakta")
            else:
                comm("skipping bakta (already done; --resume)")

        if args.cazy:
            if ckpt.todo("cazy"):
                _run_cazy(args, genes_dir, todo)
                ckpt.done("cazy")
            else:
                comm("skipping CAZyme annotation (already done; --resume)")

        if args.bgc:
            if ckpt.todo("bgc"):
                _run_bgc(args, genes_dir, todo)
                ckpt.done("bgc")
            else:
                comm("skipping cluster prediction (already done; --resume)")

        if args.kofam:
            if ckpt.todo("kofam"):
                _run_kofam(args, genes_dir, todo)
                ckpt.done("kofam")
            else:
                comm("skipping KO assignment (already done; --resume)")

        if args.pfam:
            if ckpt.todo("pfam"):
                _run_pfam(args, genes_dir, todo, _ko_genes(args))
                ckpt.done("pfam")
            else:
                comm("skipping Pfam annotation (already done; --resume)")

        comm("per-gene table: %s" % os.path.join(args.output, GENES_TABLE))
        comm("per-genome summary: %s" % os.path.join(args.output, SUMMARY_TABLE))
        announcement("ANNOTATE PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
