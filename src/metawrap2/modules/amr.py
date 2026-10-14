"""MetaWrap2 amr: resistance and virulence genes, where they sit, and what that licenses saying.

Runs AMRFinderPlus over bins or a whole assembly and turns its report into a table built around the
Cao et al. evidence ladder - see :mod:`metawrap2.amr` for the ladder, the contig-edge rule, and why
"no mobile element nearby" is a claim that needs a contig long enough to support it.

Two things this module does that a plain AMRFinderPlus wrapper does not:

* **It knows how long the contigs are**, so a gene's distance to the nearer contig end is on every
  row. AMRFinderPlus reports coordinates and not contig lengths, so a wrapper that passes its output
  through cannot tell a gene in the middle of a megabase from one 300 bp from the edge of what was
  assembled.
* **It runs before any length filter.** ARG-bearing contigs are short - one measurement puts
  MEGAHIT's at 284 bp on average - so a filter applied first removes the evidence.

``--plus`` is on by default, which also yields AMRFinderPlus' virulence and stress entries. They go
to their own table rather than into the resistance one: a virulence factor and a resistance gene are
not interchangeable units of risk, and nothing here should make it easy to add them up.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
from typing import Dict, List, Optional, Sequence, Tuple

from .. import amr as _amr
from ..cliformat import make_parser
from ..config import MODULE_ENVS, env_prefix, load_settings
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
    optional_env_for,
    require_nonempty_dir,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["amr"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# -n: nucleotide input, so the report carries contig coordinates. Protein input (-p) would lose the
# genomic context this module is built around.
#
# --plus adds AMRFinderPlus' stress-response and virulence references to the resistance ones. They
# are separated again on the way out, by Type.
#
# --database points at a specific dated release rather than letting the tool pick: an ARG call is
# only as datable as its database, and `amrfinder -u` moves "latest" under a user's feet.
AMRFINDER = (
    "amrfinder -n {input} --plus --database {db} --threads {threads} -o {out}" "{organism}{ident}"
)
# hAMRonization normalises a backend's own report into a 35-field schema shared across tools. The
# versions are passed explicitly rather than let it guess: it has no way to know which AMRFinderPlus
# database produced a report handed to it on disk, and "which database version" is one of Cao et
# al.'s L1 validity requirements rather than metadata.
#
# Worth running even with one backend. It is what makes a second one additive instead of producing a
# second incomparable table, and it maps every call onto a shared vocabulary - so a gene found by
# AMRFinderPlus and by RGI is recognisably one gene rather than two findings.
HAMRONIZE = (
    "hamronize amrfinderplus --analysis_software_version {tool_version}"
    " --reference_database_version {db_version} --input_file_name {genome}"
    " --format tsv --output {out} {report}"
)
# RGI, the opt-in second backend. `--local` keeps CARD in a directory this run owns rather than
# writing into the conda env, so two runs with two CARD versions do not overwrite each other's
# database - which matters because the version is what dates the calls.
#
# -t contig because the input is assembled sequence: RGI then calls genes with Prodigal itself and
# reports contig coordinates, which is what the context classification needs. -a DIAMOND rather than
# BLAST for the same reason the rest of the pipeline prefers it on protein search: same calls,
# minutes instead of hours. --clean removes the intermediate files, which are large and which
# nothing here reads.
RGI_LOAD = "rgi load --card_json {card_json} --local"
RGI_MAIN = "rgi main --local -i {genome} -o {prefix} -t contig -a DIAMOND -n {threads} --clean"
# hAMRonization again, for RGI's own schema. Passed the CARD version rather than letting it guess:
# it has no way to know which CARD produced a table handed to it on disk.
HAMRONIZE_RGI = (
    "hamronize rgi --analysis_software_version {tool_version}"
    " --reference_database_version {db_version} --input_file_name {genome}"
    " --format tsv --output {out} {report}"
)
# argNorm maps every normalised call onto an ARO term and onto the drugs that term confers
# resistance to. That is the second half of what makes a cross-tool table useful: hAMRonization makes
# the *columns* comparable, argNorm makes the *gene names* comparable, since one gene has a different
# symbol in every database. It also turns a gene list into something a reader can act on - "cepA"
# means little, "confers resistance to cephalosporin" is the finding.
#
# No --hamronization_skip_unsupported_tool: that flag downgrades an unsupported backend from an error
# to a warning, and a silently unmapped backend is exactly what this step exists to prevent. If a
# future backend is unsupported, failing here is the right outcome - the primary tables are already
# written by this point.
ARGNORM = "argnorm hamronization -i {input} -o {output}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_THREADS = 1
#: Resistance and stress calls.
AMR_TABLE = "amr" + TSV_SUFFIX
#: Virulence calls, deliberately a separate file.
VIRULENCE_TABLE = "virulence" + TSV_SUFFIX
#: The same calls in hAMRonization's cross-tool schema, for anyone combining this with another
#: pipeline's resistome output.
HARMONISED_TABLE = "amr_hamronized" + TSV_SUFFIX
#: The hAMRonization table with every gene mapped onto an ARO term and the drugs it confers
#: resistance to.
NORMALISED_TABLE = "amr_normalised" + TSV_SUFFIX
#: Co-located calls grouped into the unit that actually moves.
CLUSTERS_TABLE = "amr_clusters" + TSV_SUFFIX
#: Where RGI's own per-genome reports are kept, for the same reason AMRFinderPlus' are.
RGI_DIR = "rgi_reports"
#: Where the raw per-genome AMRFinderPlus reports are kept, because a derived table should never be
#: the only copy of a tool's own output.
REPORTS_DIR = "amrfinder_reports"
#: Columns of CLUSTERS_TABLE.
CLUSTER_COLUMNS = (
    "genome",
    "contig",
    "cluster_start",
    "cluster_end",
    "span_bp",
    "n_genes",
    "symbols",
    "classes",
    "distance_to_contig_end",
)
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 amr",
        usage="metawrap2 amr [options] -o output_dir (-b bin_folder | -a assembly.fasta)",
        description="Resistance and virulence genes with their genomic context and evidence level.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("-b", "--bins", help="folder of bins in fasta format")
    source.add_argument("-a", "--assembly", help="a whole assembly, screened without bins")
    threads_arg(p)
    p.add_argument(
        "--db",
        help="path to the AMRFinderPlus database (default: AMRFINDER_DB from metawrap2.toml). "
        "Point it at a dated release rather than 'latest' if you want the run to stay reproducible",
    )
    p.add_argument(
        "--organism",
        help="AMRFinderPlus organism for point-mutation screening and common-protein suppression "
        "(e.g. Escherichia). Only meaningful for a bin you have already classified",
    )
    p.add_argument(
        "--min-identity",
        type=float,
        default=None,
        metavar="FRACTION",
        help="minimum identity to the reference, 0-1. Default: AMRFinderPlus' own curated per-gene "
        "thresholds, which is the better choice - a single global cutoff replaces 7000 curated ones",
    )
    p.add_argument(
        "--rgi",
        action="store_true",
        help="also call resistance with RGI against CARD (needs `metawrap2 install-env rgi` and "
        "`install-db card`). The two backends disagree, and which genes they disagree about is the "
        "useful part; every row says which tool found it",
    )
    p.add_argument(
        "--card-db",
        metavar="DIR",
        help="CARD data directory holding card.json (default: CARD_DB from metawrap2.toml)",
    )
    p.add_argument(
        "--skip-amrfinder",
        action="store_true",
        help="run only the backends selected by flags. With --rgi that means RGI alone, which is a "
        "narrower screen: AMRFinderPlus is the only backend here that also reports virulence and "
        "stress genes",
    )
    p.add_argument(
        "--mobile-elements",
        metavar="DIR",
        help="a `metawrap2 viruses` output directory, or a geNomad run directory, to take mobile "
        "elements from. Without it every call's context is `context_unknown` and no call can reach "
        "evidence level L2 - which is accurate, since nothing looked",
    )
    p.add_argument(
        "--skip-hamronize",
        action="store_true",
        help="do not write the hAMRonization-normalised table. It is what makes this output "
        "comparable with another tool's resistome, so skipping it is a size decision rather than a "
        "quality one",
    )
    p.add_argument(
        "--mobility-window",
        type=int,
        default=_amr.MOBILITY_WINDOW,
        metavar="BP",
        help="how far from a gene to look for a mobile element, and how much flanking sequence a "
        "'nothing nearby' verdict requires (default: %(default)s). There is no field-wide standard "
        "for this, so it is recorded in every row",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _inputs(args) -> List[Tuple[str, str]]:
    if args.assembly:
        return [(os.path.splitext(os.path.basename(args.assembly))[0], args.assembly)]
    return [
        (os.path.splitext(name)[0], os.path.join(args.bins, name)) for name in list_bins(args.bins)
    ]


def _contig_lengths(path: str) -> Dict[str, int]:
    """Contig lengths, which are what make a context claim checkable."""
    return {contig_id(header): len(seq) for header, seq in iter_fasta(path)}


def _database(args, settings) -> str:
    db = args.db or settings.db("AMRFINDER_DB")
    if not db:
        error(
            "no AMRFinderPlus database: pass --db or set AMRFINDER_DB in metawrap2.toml "
            "(create one with `amrfinder_update --database DIR`)"
        )
    return db


def _database_version(db: str) -> str:
    """The dated release string, from the database directory itself.

    Read rather than asked of the tool, because the tool reports the version of whatever it resolves
    'latest' to at run time, and that is not necessarily the directory a rerun will see.
    """
    for candidate in (
        os.path.join(db, "version.txt"),
        os.path.join(os.path.dirname(db.rstrip("/")), "version.txt"),
    ):
        if os.path.isfile(candidate):
            with open(candidate) as fh:
                version = fh.readline().strip()
            if version:
                return version
    resolved = os.path.basename(os.path.realpath(db.rstrip("/")))
    return resolved or "unknown"


def env_for_rgi(args):
    """The RGI environment, for reading its version out of its own bin directory."""
    return optional_env_for("rgi", load_settings(args.config), needed_by="amr --rgi")


def _tool_version(env: Optional[str], executable: str) -> Optional[str]:
    """A tool's version string, read from inside *env* rather than from PATH.

    The same rule as :func:`_amrfinder_version`, for the same reason: a version read off PATH can
    belong to a different installation than the one that produced the calls, and it goes into the
    normalised table as provenance.
    """
    prefix = env_prefix(env) if env else None
    binary = os.path.join(prefix, "bin", executable) if prefix else executable
    if prefix and not os.path.isfile(binary):
        return None
    for argv in ([binary, "main", "--version"], [binary, "--version"]):
        try:
            out = subprocess.run(argv, capture_output=True, text=True, check=False, timeout=120)
        except (OSError, subprocess.SubprocessError):
            continue
        text = (out.stdout or out.stderr or "").strip().splitlines()
        if out.returncode == 0 and text:
            return text[-1].strip()
    return None


def _hamronize(args, env, version: str, todo: List[Tuple[str, str]]) -> List[_amr.NormalisedCall]:
    """Normalise every per-genome report into one hAMRonization table.

    hamronize takes one report at a time and each invocation writes a full header, so the per-genome
    outputs are concatenated with only the first header kept. Doing it here rather than asking
    hamronize for a combined run keeps the genome column - its ``input_file_name`` - meaningful: a
    combined run would name the file it was handed, which for a bin set is not the genome.
    """
    reports = os.path.join(args.output, REPORTS_DIR)
    pieces: List[str] = []
    work = os.path.join(args.output, "work_files")
    ensure_dir(work)
    for name, _path in todo:
        report = os.path.join(reports, name + TSV_SUFFIX)
        if not os.path.isfile(report):
            continue
        piece = os.path.join(work, "hamronized_%s%s" % (name, TSV_SUFFIX))
        rc = run(
            HAMRONIZE.format(
                tool_version=_amrfinder_version(env) or "unknown",
                db_version=version,
                genome=name,
                out=piece,
                report=report,
            ),
            env=env,
            tool="hamronize",
            check=False,
        )
        if rc == 0 and os.path.isfile(piece):
            pieces.append(piece)
        else:
            warning("hamronize failed on %s; it is left out of the normalised table" % name)

    if args.rgi:
        card = args.card_db or load_settings(args.config).db("CARD_DB")
        card_version = _amr.card_version(card) if card else "unknown"
        rgi_version = _tool_version(env_for_rgi(args), "rgi") or "unknown"
        for name, _path in todo:
            report = os.path.join(args.output, RGI_DIR, name + ".txt")
            if not os.path.isfile(report):
                continue
            piece = os.path.join(work, "hamronized_rgi_%s%s" % (name, TSV_SUFFIX))
            rc = run(
                HAMRONIZE_RGI.format(
                    tool_version=rgi_version,
                    db_version=card_version,
                    genome=name,
                    out=piece,
                    report=report,
                ),
                env=env,
                tool="hamronize rgi",
                check=False,
            )
            if rc == 0 and os.path.isfile(piece):
                pieces.append(piece)
            else:
                warning(
                    "hamronize failed on RGI's %s; it is left out of the normalised table" % name
                )

    if not pieces:
        warning("nothing to normalise; no hAMRonization table written")
        return []
    out_path = os.path.join(args.output, HARMONISED_TABLE)
    _concatenate(pieces, out_path)
    comm("hAMRonization schema -> %s" % out_path)

    # argNorm is given the AMRFinderPlus rows only. It refuses a hAMRonization table containing a
    # tool it cannot map - "argNorm can only map genes from the following tools: argsoap, deeparg,
    # resfinder, amrfinderplus, groot, abricate" - and raises on the whole file, so passing the
    # combined table produced no normalised output at all once --rgi was used. RGI loses nothing by
    # being left out: its calls are made against CARD, so they already carry their ARO accession,
    # which `amr.tsv` reports in the `aro` column.
    mappable = [piece for piece in pieces if "hamronized_rgi_" not in os.path.basename(piece)]
    if not mappable:
        comm(
            "only RGI calls to normalise, and those already carry their ARO terms from CARD; "
            "argNorm not run"
        )
        return []
    for_argnorm = out_path
    if len(mappable) != len(pieces):
        for_argnorm = os.path.join(work, "hamronized_mappable" + TSV_SUFFIX)
        _concatenate(mappable, for_argnorm)
    return _argnorm(args, env, for_argnorm)


def _concatenate(pieces: List[str], out_path: str) -> None:
    """Join per-genome tables keeping only the first header."""
    with open(out_path, "w") as out:
        for index, piece in enumerate(pieces):
            with open(piece) as fh:
                lines = fh.readlines()
            if not lines:
                continue
            out.writelines(lines if index == 0 else lines[1:])


def _argnorm(args, env, harmonised: str) -> List[_amr.NormalisedCall]:
    """Map every normalised call onto an ARO term, and report which drugs are implicated.

    A second pass over a table that already exists, so a failure here must not cost the tables
    already written - hence check=False and a warning rather than an error.
    """
    out_path = os.path.join(args.output, NORMALISED_TABLE)
    rc = run(
        ARGNORM.format(input=harmonised, output=out_path),
        env=env,
        tool="argnorm",
        check=False,
    )
    if rc != 0 or not os.path.isfile(out_path):
        warning(
            "argnorm failed, so the calls are not mapped onto ARO terms. %s is still written and "
            "is still comparable with another tool's hAMRonization output" % HARMONISED_TABLE
        )
        return []
    rows = _amr.parse_argnorm(out_path)
    mapped = [row for row in rows if row.aro]
    comm("%d of %d call(s) mapped to an ARO term -> %s" % (len(mapped), len(rows), out_path))
    drugs = _amr.drugs_implicated(rows)
    if drugs:
        # One line for the drugs and one per gene, rather than a line per gene-drug pair. A single
        # gene can confer resistance to nine aminoglycosides, and printing nine lines for it reads
        # as nine findings.
        comm(
            "%d drug(s) implicated: %s%s"
            % (
                len(drugs),
                ", ".join(drug for drug, _genes in drugs[:12]),
                " ..." if len(drugs) > 12 else "",
            )
        )
    for row in mapped[:10]:
        listed = ", ".join(row.drugs[:4]) + (
            " (+%d more)" % (len(row.drugs) - 4) if len(row.drugs) > 4 else ""
        )
        comm("  %s (%s, %s): %s" % (row.gene, row.genome, row.aro, listed or "no drug recorded"))
    if len(mapped) > 10:
        comm("  ... and %d more in the table" % (len(mapped) - 10))
    if len(mapped) < len(rows):
        comm(
            "an unmapped call is not a false one: argNorm maps what the ontology holds, and a "
            "recently described gene can be absent from it"
        )
    return rows


def _amrfinder_version(env: Optional[str]) -> Optional[str]:
    """AMRFinderPlus' own version, for the provenance hAMRonization records.

    Asked of the tool rather than hard-coded: the env is pinned to a range, so the version in a
    given install is whatever that range resolved to, and recording the wrong one makes the
    normalised table claim provenance it does not have.

    Run from inside the env's own bin directory, not from PATH. A bare
    ``subprocess.run(["amrfinder", ...])`` reads whatever is on PATH, which is the mistake `check`
    was making - it reported a tool present because the user had it installed system-wide while the
    run used a conda env that did not. Here it would be quieter and worse: the normalised table
    would record a version number from a different installation than the one that produced the
    calls.
    """
    executable = "amrfinder"
    if env:
        prefix = env_prefix(env)
        if not prefix:
            return None
        executable = os.path.join(prefix, "bin", "amrfinder")
        if not os.access(executable, os.X_OK):
            return None
    try:
        result = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
    except OSError:
        return None
    return result.stdout.strip() or None


def _run_rgi(args, settings, todo: List[Tuple[str, str]]) -> List[_amr.ArgCall]:
    """Screen every genome with RGI against CARD, and return its calls.

    CARD is loaded once into a localDB this run owns, not into the conda environment: two runs
    against two CARD versions would otherwise overwrite each other's database, and the version is
    what dates the calls. RGI is run from that directory because `--local` looks for localDB in the
    working directory.
    """
    card = args.card_db or settings.db("CARD_DB")
    if not card:
        error(
            "no CARD database: pass --card-db or set CARD_DB in metawrap2.toml "
            "(`metawrap2 install-db card`)"
        )
    card_json = card if card.endswith(".json") else os.path.join(card, "card.json")
    if not os.path.isfile(card_json):
        error("no card.json under %s; CARD_DB should name the directory holding it" % card)
    version = _amr.card_version(card_json)
    env = optional_env_for("rgi", settings, needed_by="amr --rgi")

    work = os.path.join(args.output, "work_files", "rgi")
    ensure_dir(work)
    reports = os.path.join(args.output, RGI_DIR)
    ensure_dir(reports)
    comm("RGI against CARD %s" % version)
    rc = run(
        RGI_LOAD.format(card_json=card_json),
        env=env,
        tool="rgi load",
        check=False,
        cwd=work,
    )
    if rc != 0:
        warning("`rgi load` failed, so RGI contributes no calls; AMRFinderPlus' stand")
        return []

    calls: List[_amr.ArgCall] = []
    for name, path in bar(todo, desc="rgi", unit=" genome", total=len(todo)):
        prefix = os.path.join(reports, name)
        rc = run(
            RGI_MAIN.format(genome=path, prefix=prefix, threads=args.threads),
            env=env,
            tool="rgi main",
            check=False,
            cwd=work,
        )
        table = prefix + ".txt"
        if rc != 0 or not os.path.isfile(table):
            warning("RGI failed on %s; that genome has AMRFinderPlus' calls only" % name)
            continue
        calls.extend(_amr.parse_rgi(table, name, _contig_lengths(path), version))
    comm("RGI made %d call(s) across %d genome(s)" % (len(calls), len(todo)))
    return calls


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
    validate_inputs("amr", checks)

    env = env_for("amr", settings)
    db = _database(args, settings)
    version = _database_version(db)
    rec = start_run("amr", args, env, settings, inputs=[source])

    ckpt = make_checkpoint(args.output)
    try:
        announcement("BEGIN AMR PIPELINE!")
        ensure_dir(args.output)
        reports = os.path.join(args.output, REPORTS_DIR)
        ensure_dir(reports)

        todo = _inputs(args)
        if not todo:
            error("No fasta input found (%s). Exiting." % "/".join(FASTA_EXTENSIONS))
        comm("%d genome(s) to screen; AMRFinderPlus database %s" % (len(todo), version))

        all_calls: List[_amr.ArgCall] = []
        if args.skip_amrfinder:
            comm(
                "skipping AMRFinderPlus (--skip-amrfinder). No virulence or stress calls will be "
                "made: it is the only backend here that reports them"
            )
        elif ckpt.todo("amrfinder"):
            organism = " --organism %s" % args.organism if args.organism else ""
            ident = " --ident_min %g" % args.min_identity if args.min_identity is not None else ""
            for name, path in bar(todo, desc="screening", unit=" genome", total=len(todo)):
                report = os.path.join(reports, name + TSV_SUFFIX)
                rc = run(
                    AMRFINDER.format(
                        input=path,
                        db=db,
                        threads=args.threads,
                        out=report,
                        organism=organism,
                        ident=ident,
                    ),
                    env=env,
                    tool="amrfinder",
                    check=False,
                )
                if rc != 0:
                    warning("AMRFinderPlus failed on %s. Proceeding with the rest." % name)
                    continue
                if dry_run():
                    continue
                all_calls.extend(_amr.parse_amrfinder(report, name, _contig_lengths(path), version))
            ckpt.done("amrfinder")
        elif not args.skip_amrfinder:
            comm("re-reading existing AMRFinderPlus reports (--resume)")
            for name, path in todo:
                report = os.path.join(reports, name + TSV_SUFFIX)
                if os.path.isfile(report):
                    all_calls.extend(
                        _amr.parse_amrfinder(report, name, _contig_lengths(path), version)
                    )

        if args.rgi and not dry_run():
            all_calls.extend(_run_rgi(args, settings, todo))

        if dry_run():
            if args.rgi:
                comm("would also screen every genome with RGI against CARD")
            announcement("AMR PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        # None rather than an empty list when no mobile-element caller has run: see
        # amr.classify_context on why those are different claims. With --mobile-elements the
        # question *was* asked, so an empty result is an empty list and can support
        # no_mge_in_window.
        elements: Optional[List[_amr.MobileElement]] = None
        whole_contig: Dict[str, str] = {}
        if args.mobile_elements:
            elements, whole_contig = _amr.mobile_elements_from_genomad(args.mobile_elements)
            comm(
                "mobile elements: %d provirus region(s) with coordinates, %d whole contig(s) "
                "called plasmid or viral, from %s"
                % (len(elements), len(whole_contig), args.mobile_elements)
            )
            if not elements and not whole_contig:
                warning(
                    "no mobile-element calls found under %s. Every call's context stays "
                    "`context_unknown`, which is the honest answer: nothing was found to compare "
                    "against, and that is not the same as nothing being there"
                    % args.mobile_elements
                )
                elements = None
        resistance = [c for c in all_calls if not c.is_virulence]
        virulence = [c for c in all_calls if c.is_virulence]

        for table, calls in ((AMR_TABLE, resistance), (VIRULENCE_TABLE, virulence)):
            path = os.path.join(args.output, table)
            with open(path, "w", newline="") as fh:
                writer = csv.DictWriter(
                    fh, fieldnames=_amr.ARG_COLUMNS, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                verdicts: List[Tuple[_amr.Context, Optional[int], Sequence[str]]] = []
                for call in sorted(calls, key=lambda c: (c.genome, c.contig, c.start)):
                    context, distance, kinds = _amr.classify_context(
                        call, elements, args.mobility_window, whole_contig=whole_contig or None
                    )
                    verdicts.append((context, distance, kinds))
                    writer.writerow(
                        _amr.arg_row(call, context, distance, kinds, args.mobility_window)
                    )
            comm("%d row(s) -> %s" % (len(calls), path))
            if calls and (elements is not None or whole_contig):
                counts = _amr.context_summary(verdicts)
                comm(
                    "  context: %d on or beside a mobile element, %d with none within %d bp, "
                    "%d undetermined%s"
                    % (
                        counts[_amr.Context.MGE_NEARBY.value],
                        counts[_amr.Context.NO_MGE_IN_WINDOW.value],
                        args.mobility_window,
                        counts[_amr.Context.UNKNOWN.value],
                        (
                            " ("
                            + ", ".join(
                                "%s %d" % (kind, counts[kind])
                                for kind in sorted(counts)
                                if kind not in {c.value for c in _amr.Context}
                            )
                            + ")"
                            if any(k not in {c.value for c in _amr.Context} for k in counts)
                            else ""
                        ),
                    )
                )

        clusters_path = os.path.join(args.output, CLUSTERS_TABLE)
        with open(clusters_path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=CLUSTER_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            n_clusters = 0
            for genome in sorted({c.genome for c in resistance}):
                genome_calls = [c for c in resistance if c.genome == genome]
                for group in _amr.cluster_calls(genome_calls, args.mobility_window):
                    if len(group) < 2:
                        continue
                    writer.writerow(
                        {
                            "genome": genome,
                            "contig": group[0].contig,
                            "cluster_start": group[0].start,
                            "cluster_end": group[-1].end,
                            "span_bp": group[-1].end - group[0].start + 1,
                            "n_genes": len(group),
                            "symbols": ",".join(c.symbol for c in group),
                            "classes": ",".join(sorted({c.drug_class for c in group})),
                            "distance_to_contig_end": min(c.distance_to_contig_end for c in group),
                        }
                    )
                    n_clusters += 1
        comm("%d multi-gene cluster(s) -> %s" % (n_clusters, clusters_path))

        normalised: List[_amr.NormalisedCall] = []
        if not args.skip_hamronize:
            normalised = _hamronize(args, env, version, todo)

        if len({call.tool for call in all_calls}) > 1:
            # Joined on ARO where it is known, which needs argNorm's output - hence after it. The
            # symbol join undercounts: AMRFinderPlus' `cepA` and RGI's `CepA-29` are the same
            # beta-lactamase at the same identity, two names and one ontology term.
            aro_for = {(row.genome, row.gene.lower()): row.aro for row in normalised if row.aro}
            agreement = _amr.backend_agreement(all_calls, aro_for=aro_for)
            comm(
                "backends: %s; %d gene(s) called by both%s"
                % (
                    ", ".join(
                        "%s %d" % (tool, agreement.per_tool[tool]) for tool in agreement.tools
                    ),
                    agreement.shared,
                    " (joined on ARO terms)" if aro_for else " (joined on gene symbols only)",
                )
            )
            for tool in agreement.tools:
                only = agreement.unique.get(tool, ())
                if only:
                    comm(
                        "  only %s: %s%s"
                        % (tool, ", ".join(only[:8]), " ..." if len(only) > 8 else "")
                    )
            comm(
                "  disagreement is the expected result rather than a fault: the backends model "
                "resistance differently and their thresholds are curated separately. RGI's extra "
                "calls here are mostly Strict-tier homologs well below 50% identity, which is a "
                "weaker claim than an exact match and is in the `method` column"
            )

        counts = _amr.summarise(all_calls)
        comm(
            "%d AMR, %d stress, %d virulence; %d called from a partial match, "
            "%d within %d bp of a contig end"
            % (
                counts["amr"],
                counts["stress"],
                counts["virulence"],
                counts["partial_method"],
                counts["within_window_of_contig_end"],
                _amr.MOBILITY_WINDOW,
            )
        )
        comm(
            "every row is evidence level L1 or L2 (Cao et al. 2026): detection, and detection with "
            "a genomic context. Nothing here supports 'confers resistance'."
        )
        announcement("AMR PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
