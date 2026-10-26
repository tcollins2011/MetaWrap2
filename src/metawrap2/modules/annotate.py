"""MetaWrap2 annotate: call genes once over a bin set or a whole assembly, then annotate them.

This replaces ``annotate_bins``, which ran Prokka once per bin and sorted the output into three
folders. That design had two limits:

* **No joinable output.** Per-bin GFFs with product names do not give KEGG orthologs for pathway
  completeness, gene neighbourhoods for AMR, or stable gene ids for read mapping.
* **Bins only.** Most assembled sequence never reaches a MAG, and ``-b bins`` was the only input.

Genes are now called once, in :mod:`metawrap2.genecall`, from either a bin folder (``-b``) or an
assembly (``-a``). Orthologs, domains, pathways and resistance all read the gene table written here
instead of calling genes again, so every table refers to the same gene set.

The commands each external tool runs are in the COMMANDS block below, and the tunable numbers in the
CONSTANTS block.
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import Dict, List, Optional

from .. import bgc, cazy, eggnog, genecall, kofam, pfam
from ..checkm import tag_mimag
from ..cliformat import make_parser
from ..command import capture_in_env
from ..config import MODULE_ENVS, load_settings
from ..constants import FASTA_EXTENSIONS, STATS_SUFFIX, TSV_SUFFIX
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
    optional_env_for,
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
# Bakta re-predicts genes (--skip-cds with a --regions file is not stable across Bakta versions).
# It is used for product names and for non-coding features Prodigal does not find: tRNAs via
# tRNAscan-SE and rRNAs via Infernal/Rfam, which MIMAG's rRNA/tRNA criteria need.
BAKTA = (
    "bakta --threads {threads} --output {outdir} --prefix {prefix} --force "
    "--keep-contig-headers {db} {input}"
)
# CAZyme families over the proteins called here. --mode protein stops dbCAN calling its own genes
# with Prodigal, which would give CAZyme calls on genes absent from the other tables.
#
# All three methods run by default. dbCAN recommends trusting a family found by at least two, and on
# one bin a third of the calls were found by only one.
DBCAN = (
    "run_dbcan CAZyme_annotation --mode protein --input_raw_data {proteins}"
    " --output_dir {outdir} --db_dir {db} --threads {threads}"
)
# Biosynthetic gene clusters over the genes called here (--gff), so a BGC's member genes are the
# genes in genes.tsv. Without --gff GECCO calls its own genes.
#
# --force-tsv: otherwise a genome with no clusters (the common case) writes no table, which looks
# the same as a failed run.
GECCO = (
    "gecco run -g {genome} --gff {gff} -o {outdir} -j {threads} --force-tsv --no-color --no-markup"
)
# Orthology-based annotation over the proteins called here. `--itype proteins` is emapper's default
# but is passed explicitly; with `--itype genome` eggNOG would call its own genes with Prodigal.
#
# `--override` is required: emapper refuses to start if any of its outputs already exist, which is
# always the case on a re-run into the same -o.
#
# `--data_dir`, `--output_dir` and `--temp_dir` must exist (emapper's argparse uses
# `type=existing_dir`), so the output directory is created and the database checked beforehand.
EGGNOG = (
    "emapper.py -i {proteins} --itype proteins -m diamond --data_dir {db}"
    " --output_dir {outdir} --temp_dir {outdir} -o {prefix} --cpu {threads} --override"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_THREADS = 1
#: Where per-genome gene files go, one set per genome, named for the genome.
GENES_DIR = "genes"
#: The per-gene table every downstream module reads. Long, one row per gene.
GENES_TABLE = "genes" + TSV_SUFFIX
#: The per-gene KO table, written by --kofam. One row per accepted assignment, with the evidence
#: it qualified on (see metawrap2.kofam).
KOS_TABLE = "kos" + TSV_SUFFIX
#: One row per Pfam domain found, written by --pfam. Per domain because a protein can carry several,
#: including repeats, and the architecture is rebuilt from the rows.
PFAM_TABLE = "pfam_domains" + TSV_SUFFIX
#: One row per gene with its ordered domain architecture.
PFAM_ARCHITECTURE_TABLE = "pfam_architectures" + TSV_SUFFIX
#: The per-gene CAZyme table, written by --cazy. One row per gene with a CAZyme call, carrying which
#: of dbCAN's three methods found it.
CAZY_TABLE = "cazymes" + TSV_SUFFIX
#: Per-class CAZyme counts per genome.
CAZY_CLASS_TABLE = "cazyme_classes" + TSV_SUFFIX
#: Where dbCAN's own output is kept, one directory per genome.
CAZY_DIR = "dbcan"
#: The per-gene eggNOG table, written by --eggnog. One row per gene emapper placed, carrying the
#: COG category, GO terms and EC numbers that neither --kofam nor --pfam produces.
EGGNOG_TABLE = "eggnog" + TSV_SUFFIX
#: One row per (gene, namespace, term). Long format, since a gene can have any number of GO terms.
EGGNOG_TERM_TABLE = "eggnog_terms" + TSV_SUFFIX
#: Where emapper's own output is kept, one directory per genome: the .annotations table parsed here,
#: plus the .seed_orthologs and .hits files, which a re-annotation at a different --tax_scope needs
#: and which are expensive to regenerate.
EGGNOG_DIR = "eggnog"
#: Files `download_eggnog_data.py` writes and `-m diamond` needs: the annotation SQLite and the
#: DIAMOND database of eggNOG proteins. Checked by name to catch --eggnog-db pointing at the parent
#: directory; emapper's own error for that repeats per genome and names --data_dir.
EGGNOG_DB_FILES = ("eggnog.db", "eggnog_proteins.dmnd")
#: One row per predicted biosynthetic gene cluster, written by --bgc.
BGC_TABLE = "gene_clusters" + TSV_SUFFIX
#: Per-genome cluster counts, split by whether a cluster runs into unassembled sequence.
BGC_SUMMARY_TABLE = "gene_cluster_summary" + TSV_SUFFIX
#: Where GECCO's own output is kept, one directory per genome, including the per-cluster GenBank
#: files for an antiSMASH sideload or manual inspection.
BGC_DIR = "gecco"
#: Per-genome summary of gene calling. Check it when a genome's annotation looks wrong: the usual
#: causes, a coding density far from 0.85-0.90 or a bin that fell back to meta mode, show here.
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
        "profile's curated threshold, flagged as evidence=relative (default: off, curated "
        "thresholds only, as in KofamScan)",
    )
    p.add_argument(
        "--kofam-thresholdless-evalue",
        type=float,
        default=None,
        metavar="E",
        help="also accept hits to the 1866 profiles that have no KOfam threshold, when their "
        "E-value is at most E, flagged as evidence=no_threshold (default: off)",
    )
    p.add_argument(
        "--pfam",
        action="store_true",
        help="also assign Pfam domains to every called protein, at each family's gathering "
        "threshold. This annotates genes that get no KO (on real bins only 41%% of genes get a "
        "curated KO). Needs PFAM_DB (`metawrap2 install-db pfam`, 399 MB)",
    )
    p.add_argument(
        "--pfam-db",
        metavar="DIR",
        help="directory holding Pfam-A.hmm (default: PFAM_DB from metawrap2.toml)",
    )
    p.add_argument(
        "--eggnog",
        action="store_true",
        # No literal per-cent sign in this help: argparse treats it as a format template, and a
        # bare one once made `--help` fail for every command.
        help="also annotate the called proteins by orthology with eggNOG-mapper, which adds COG "
        "categories, GO terms and EC numbers (--kofam and --pfam produce none of these). Needs "
        "the eggNOG data (about 50 GB) via --eggnog-db",
    )
    p.add_argument(
        "--eggnog-db",
        metavar="DIR",
        help="directory holding eggnog.db and eggnog_proteins.dmnd, as "
        "download_eggnog_data.py writes them (default: EGGNOG_DB from metawrap2.toml)",
    )
    p.add_argument(
        "--cazy",
        action="store_true",
        help="also assign CAZyme families to the called proteins with dbCAN (needs DBCAN_DB). "
        "CAZy families give the glycosidic bond an enzyme acts on, which KOs do not",
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
        "(default: all three). The number of agreeing methods is reported per gene, so running "
        "fewer methods lowers every gene's count",
    )
    p.add_argument(
        "--bgc",
        action="store_true",
        help="also predict biosynthetic gene clusters with GECCO over the called genes (no "
        "database needed). The table gives each cluster's distance to a contig edge",
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

    One entry per bin, or a single entry for an assembly. An assembly is always called in
    metagenomic mode, whatever ``--meta`` says, because one coding model trained on many organisms
    fits none of them (see genecall).
    """
    if args.assembly:
        name = os.path.splitext(os.path.basename(args.assembly))[0]
        return [(name, args.assembly, True)]
    names = list_bins(args.bins)
    return [(os.path.splitext(name)[0], os.path.join(args.bins, name), args.meta) for name in names]


def _unit(args) -> str:
    """Name of the per-genome loop unit, so messages match the input flag given."""
    return "bins" if args.bins else "genomes"


def _require_some_output(tool: str, failed: int, total: int, unit: str) -> None:
    """Fail the step when *every* item in a per-item loop failed.

    One failing genome is warned about and skipped. When every item fails there is a common cause
    and no output, so the step errors out instead of exiting 0. Guards against a run that failed
    Bakta on all 26 bins (no database) and still reported success.
    """
    if total and failed >= total:
        error(
            "%s failed on %d of %d %s; no %s output was produced"
            % (tool, failed, total, unit, tool)
        )


def _call_and_write(name: str, path: str, meta: bool, genes_dir: str) -> genecall.GenomeGenes:
    called = genecall.call_genome(path, name=name, meta=meta)
    genecall.write_faa(called, os.path.join(genes_dir, name + ".faa"))
    genecall.write_ffn(called, os.path.join(genes_dir, name + ".ffn"))
    genecall.write_gff(called, os.path.join(genes_dir, name + ".gff"))
    return called


def _bakta_db(args, settings) -> Optional[str]:
    """The Bakta database to annotate against, checked before anything is launched.

    Bakta has no default database, so without BAKTA_DB every per-bin launch fails. Checking here
    reports that once, before the first launch. BAKTA_DB is a config key with a
    ``check_databases()`` entry.
    """
    db = args.bakta_db or settings.db("BAKTA_DB")
    if not args.bakta:
        return db
    if not db:
        error(
            "--bakta needs BAKTA_DB; set it in metawrap2.toml "
            "(bakta_db download --output <dir>), or pass --bakta-db"
        )
    if not os.path.exists(db):
        error("the Bakta database does not exist: %s (BAKTA_DB / --bakta-db)" % db)
    return db


def _kofam_profiles(args, settings) -> List[str]:
    """The HMM files to search, from a directory of them or a single concatenated file.

    KOfam ships one ``.hmm`` per KO (28,453 files) in a directory. A single concatenated file is
    also accepted; hmmsearch and PyHMMER both read it, and it opens faster.
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

    Read from disk so it also works on a `--resume` where the KOfam pass ran earlier. Empty when
    --kofam was not used; the Pfam report then omits the increment over KOs.
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


def _eggnog_db(args, settings) -> str:
    """The eggNOG data directory, checked before emapper is launched once per genome.

    ``EGGNOG_DB`` is not a config key yet, so ``settings.db`` returns "" for it and ``--eggnog-db``
    is currently the only way to pass it. Once the key is added the config default works here
    without changes.

    Checked here because emapper's ``--data_dir`` is ``type=existing_dir``: a missing directory
    would be an argparse error from emapper, repeated per genome, naming a flag the user did not
    type.
    """
    db = args.eggnog_db or settings.db("EGGNOG_DB")
    if not db:
        error(
            "--eggnog needs the eggNOG data (about 50 GB: eggnog.db and eggnog_proteins.dmnd): "
            "pass --eggnog-db, or set EGGNOG_DB in metawrap2.toml. Download it with "
            "eggNOG-mapper's `download_eggnog_data.py --data_dir <dir>`"
        )
    if not os.path.isdir(db):
        error("the eggNOG data directory does not exist: %s (EGGNOG_DB / --eggnog-db)" % db)
    missing = [name for name in EGGNOG_DB_FILES if not os.path.isfile(os.path.join(db, name))]
    if missing:
        error(
            "%s is not an eggNOG data directory: %s missing. Run `download_eggnog_data.py "
            "--data_dir %s`" % (db, ", ".join(missing), db)
        )
    return db


def _eggnog_version(env: Optional[str], db: str) -> tuple:
    """(what emapper reports about itself, the installed eggNOG DB version).

    ``emapper.py --version`` prints the tool version, the expected DB version and the installed DB
    version on one line. Recorded because eggNOG renumbers groups between releases, so a result is
    only reproducible against the same release.

    Returns "unknown" on failure; a failed version probe does not stop the annotation.
    """
    line = capture_in_env(env, ["emapper.py", "--data_dir", db, "--version"]).strip()
    line = line.splitlines()[0] if line else ""
    db_version = "unknown"
    marker = "Installed eggNOG DB version:"
    if marker in line:
        db_version = line.split(marker, 1)[1].split("/")[0].strip() or "unknown"
    return line or "unknown", db_version


def _kofam_kos(args) -> Dict[str, Dict[str, List[str]]]:
    """genome -> gene -> the KOs ``--kofam`` assigned it, from this run's own KO table.

    Read from disk, as in :func:`_ko_genes`. Keyed by gene with the K numbers kept, because the
    eggNOG comparison counts genes where both tools called the same KO, not just any KO.

    Empty when --kofam was not used; the eggNOG report then omits the overlap.
    """
    path = os.path.join(args.output, KOS_TABLE)
    if not os.path.isfile(path):
        return {}
    out: Dict[str, Dict[str, List[str]]] = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            genome = (row.get("genome") or "").strip()
            gene = (row.get("gene_id") or "").strip()
            ko = (row.get("ko") or "").strip()
            if genome and gene and ko:
                out.setdefault(genome, {}).setdefault(gene, []).append(ko)
    return out


def _run_pfam(args, genes_dir: str, todo: List[tuple], ko_genes_by_genome=None) -> None:
    """Search every genome's proteins against Pfam-A and write the domain tables.

    Runs in-process through PyHMMER, like the KOfam search. Per genome, as in the CAZyme pass,
    because tables here are keyed on (genome, gene) and genomes can share gene ids.
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
    failed = 0
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
                failed += 1
                warning("no proteins for %s; skipping Pfam annotation" % name)
                continue
            proteins = [(contig_id(header), sequence) for header, sequence in iter_fasta(faa)]
            if not proteins:
                failed += 1
                continue
            try:
                hits = pfam.search(proteins, hmm, threads=args.threads)
            except Exception as problem:  # noqa: BLE001 - one genome failing must not lose the rest
                failed += 1
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

    _require_some_output("Pfam", failed, len(todo), _unit(args))
    share = 100.0 * totals["annotated"] / totals["genes"] if totals["genes"] else 0.0
    comm(
        "%d domain(s) over %d of %d gene(s) (%.1f%%) -> %s"
        % (totals["domains"], totals["annotated"], totals["genes"], share, domain_path)
    )
    if ko_genes_by_genome:
        comm("%d gene(s) have a Pfam domain and no KEGG ortholog" % totals["rescued"])
        comm(pfam.INCREMENT_CAVEAT)
    comm("-> %s" % arch_path)
    comm(pfam.DOMAIN_CAVEAT)


def _run_cazy(args, genes_dir: str, todo: List[tuple]) -> None:
    """Run dbCAN over each genome's proteins and write the CAZyme tables.

    Per genome instead of over a concatenation: tables here are keyed on (genome, gene) but dbCAN
    keys on gene alone, so two genomes sharing a gene id would merge. The catalogue module prefixes
    its ids for the same reason.
    """
    settings = load_settings(args.config)
    # Looked up here instead of using main's `env`, which is None unless --bakta was given. dbCAN is
    # in the annotate env either way; None would pick up run_dbcan from the user's PATH.
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
    total = confident = failed = 0
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
                failed += 1
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
                failed += 1
                warning("dbCAN failed on %s; it contributes no CAZyme calls" % name)
                continue
            calls = cazy.parse_overview(overview, name)
            for call in calls:
                gene_writer.writerow(cazy.cazyme_row(call))
            for summary in cazy.summarise_classes(calls, name):
                class_writer.writerow(cazy.class_row(summary))
            total += len(calls)
            confident += sum(1 for call in calls if call.is_confident)

    _require_some_output("dbCAN", failed, len(todo), _unit(args))
    comm(
        "%d CAZyme gene(s), %d of them found by at least %d of dbCAN's three methods -> %s"
        % (total, confident, cazy.MIN_TOOLS, gene_path)
    )
    comm(
        "CBMs are binding modules, not enzymes, so %s lists them separately from the catalytic "
        "classes" % os.path.basename(class_path)
    )


def _run_bgc(args, genes_dir: str, todo: List[tuple]) -> None:
    """Predict gene clusters per genome and write both BGC tables.

    Per genome, as in the CAZyme pass: tables here are keyed on (genome, gene), and GECCO's output
    is keyed on contig and gene alone.
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
    total = truncated = failed = 0
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
                failed += 1
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
                # genome name (`bin.1.fa` has stem `bin.1`). Fall back to the table it wrote.
                found = sorted(
                    os.path.join(outdir, entry)
                    for entry in (os.listdir(outdir) if os.path.isdir(outdir) else [])
                    if entry.endswith(".clusters.tsv")
                )
                table = found[0] if found else ""
            if rc != 0 or not table:
                failed += 1
                warning("gecco failed on %s; it contributes no clusters" % name)
                continue
            lengths = {contig_id(header): len(seq) for header, seq in iter_fasta(path)}
            clusters = bgc.parse_clusters(table, name, lengths)
            for cluster in clusters:
                cluster_writer.writerow(bgc.cluster_row(cluster))
            summary_writer.writerow(bgc.summary_row(bgc.summarise(clusters, name)))
            total += len(clusters)
            truncated += sum(1 for cluster in clusters if cluster.is_truncated)

    _require_some_output("GECCO", failed, len(todo), _unit(args))
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
            "%s differs from the KOfam threshold list this release was tested with. KO calls "
            "use its cutoffs; each row records which list was used." % ko_list
        )
    profiles = kofam.load_ko_list(ko_list)
    ko_list_tag = kofam.ko_list_label(ko_list)
    out_path = os.path.join(args.output, KOS_TABLE)
    counts = {evidence: 0 for evidence in kofam.Evidence}
    failed = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=kofam.KO_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for name, _path, _meta in bar(todo, desc="assigning KOs", unit=" genome", total=len(todo)):
            faa = os.path.join(genes_dir, name + ".faa")
            if not os.path.isfile(faa):
                failed += 1
                warning("no proteins for %s; skipping KO assignment" % name)
                continue
            hits = kofam.search(faa, profile_files, threads=args.threads)
            for call in kofam.assign(hits, profiles, policy):
                writer.writerow(kofam.ko_row(name, call, policy, ko_list_tag))
                counts[call.evidence] += 1
    _require_some_output("KOfam", failed, len(todo), _unit(args))
    total = sum(counts.values())
    comm(
        "assigned %d KOs (%s)"
        % (
            total,
            ", ".join("%d %s" % (n, e.value) for e, n in counts.items() if n),
        )
    )
    comm("per-gene KO table: %s" % out_path)


def _run_eggnog(args, genes_dir: str, todo: List[tuple], kofam_by_genome=None) -> None:
    """Annotate every genome's proteins by orthology with eggNOG-mapper and write both tables.

    Per genome, as in the CAZyme and cluster passes: tables here are keyed on (genome, gene) and
    emapper keys on the query name alone, so genomes sharing a gene id would merge.

    The overlap with KOfam is reported but not merged into either table. eggNOG transfers a KO from
    an orthologous group; KOfam calls one from the profile's curated bit-score threshold. Keeping
    both preserves the comparison between them.
    """
    settings = load_settings(args.config)
    # A separate environment from annotate, looked up here because main's `env` is None unless
    # --bakta was given. eggNOG-mapper 2.1.x imports distutils, removed in Python 3.12, so it needs
    # a Python ceiling that bakta, dbCAN and GECCO should not share.
    env = optional_env_for("eggnog", settings, needed_by="annotate --eggnog")
    db = _eggnog_db(args, settings)
    if dry_run():
        comm("would annotate %d genome(s) by orthology with eggNOG-mapper" % len(todo))
        return

    version, db_version = _eggnog_version(env, db)
    comm("%s; annotations transferred from eggNOG DB %s" % (version, db_version))
    root = os.path.join(args.output, EGGNOG_DIR)
    ensure_dir(root)
    gene_path = os.path.join(args.output, EGGNOG_TABLE)
    term_path = os.path.join(args.output, EGGNOG_TERM_TABLE)
    totals = {"genes": 0, "annotated": 0, "orthogonal": 0, "terms": 0}
    overlaps = {"both": 0, "agreeing": 0, "eggnog_only": 0, "kofam_only": 0}
    failed = 0
    with (
        open(gene_path, "w", newline="") as gene_fh,
        open(term_path, "w", newline="") as term_fh,
    ):
        gene_writer = csv.DictWriter(
            gene_fh, fieldnames=eggnog.GENE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        gene_writer.writeheader()
        term_writer = csv.DictWriter(
            term_fh, fieldnames=eggnog.TERM_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        term_writer.writeheader()
        for name, _path, _meta in bar(todo, desc="eggNOG", unit=" genome", total=len(todo)):
            faa = os.path.join(genes_dir, name + ".faa")
            if not os.path.isfile(faa):
                failed += 1
                warning("no proteins for %s; skipping eggNOG annotation" % name)
                continue
            outdir = os.path.join(root, name)
            # emapper's --output_dir and --temp_dir are type=existing_dir, so create it first.
            ensure_dir(outdir)
            rc = run(
                EGGNOG.format(
                    proteins=faa, outdir=outdir, db=db, prefix=name, threads=args.threads
                ),
                env=env,
                tool="emapper.py",
                check=False,
            )
            annotations = os.path.join(outdir, "%s.emapper.annotations" % name)
            if rc != 0 or not os.path.isfile(annotations):
                failed += 1
                warning("eggNOG-mapper failed on %s; it contributes no annotations" % name)
                continue
            rows = eggnog.parse_annotations(annotations, name)
            for row in rows:
                gene_writer.writerow(eggnog.gene_row(row, db_version))
            terms = eggnog.term_rows(rows, db_version)
            for term in terms:
                term_writer.writerow(term)
            proteins = [contig_id(header) for header, _sequence in iter_fasta(faa)]
            found = eggnog.coverage(rows, proteins)
            totals["genes"] += found.genes
            totals["annotated"] += found.genes_annotated
            totals["orthogonal"] += found.genes_with_orthogonal_annotation
            totals["terms"] += len(terms)
            if kofam_by_genome:
                overlap = eggnog.compare_kos(rows, kofam_by_genome.get(name, {}))
                overlaps["both"] += overlap.both
                overlaps["agreeing"] += overlap.agreeing
                overlaps["eggnog_only"] += overlap.eggnog_only
                overlaps["kofam_only"] += overlap.kofam_only

    _require_some_output("eggNOG-mapper", failed, len(todo), _unit(args))
    share = 100.0 * totals["annotated"] / totals["genes"] if totals["genes"] else 0.0
    comm(
        "%d of %d gene(s) (%.1f%%) placed in an orthologous group -> %s"
        % (totals["annotated"], totals["genes"], share, gene_path)
    )
    comm(
        "%d of them carry a GO term, EC number or COG category, which --kofam and --pfam do "
        "not produce (%d term(s) -> %s)" % (totals["orthogonal"], totals["terms"], term_path)
    )
    comm(eggnog.ORTHOGONAL_CAVEAT)
    if kofam_by_genome:
        comm(
            eggnog.KoOverlap(
                both=overlaps["both"],
                eggnog_only=overlaps["eggnog_only"],
                kofam_only=overlaps["kofam_only"],
                agreeing=overlaps["agreeing"],
            ).describe()
        )
    comm(eggnog.TRANSFER_CAVEAT)


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
                    # These genomes used a different gene model than the rest, which affects
                    # comparisons of gene counts between genomes.
                    warning(
                        "%d genome(s) had under %d bp and were called in metagenomic mode "
                        "instead of with their own trained model: %s"
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
                failed = 0
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
                        failed += 1
                        warning("Bakta failed on %s. Proceeding with the rest." % name)
                _require_some_output("Bakta", failed, len(todo), _unit(args))
                # Only Bakta runs tRNAscan-SE, Aragorn and Infernal, so only its output can answer
                # MIMAG's rRNA and tRNA criteria. Tag the stats table beside the bins now. A genome
                # annotated without Bakta keeps `mimag_rrna=unknown` instead of being downgraded.
                for _name, genome_path, _meta in todo:
                    stats = os.path.dirname(genome_path.rstrip("/")) + STATS_SUFFIX
                    if os.path.isfile(stats):
                        tag_mimag(stats, evidence_dir=args.output)
                        break
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

        # Must run after KOfam: the overlap summary reads kos.tsv from disk. Run earlier, it finds
        # no KO table and reports that eggNOG called every KO alone. (The Pfam block once ran before
        # KOfam and its reported increment was always zero for this reason.)
        if args.eggnog:
            if ckpt.todo("eggnog"):
                _run_eggnog(args, genes_dir, todo, _kofam_kos(args))
                ckpt.done("eggnog")
            else:
                comm("skipping eggNOG annotation (already done; --resume)")

        comm("per-gene table: %s" % os.path.join(args.output, GENES_TABLE))
        comm("per-genome summary: %s" % os.path.join(args.output, SUMMARY_TABLE))
        announcement("ANNOTATE PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
