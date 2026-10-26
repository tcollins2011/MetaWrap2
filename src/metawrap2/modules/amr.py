"""MetaWrap2 amr: resistance and virulence genes, their genomic context, and evidence level.

Runs AMRFinderPlus over bins or a whole assembly and builds a table around the Cao et al. evidence
ladder. See :mod:`metawrap2.amr` for the ladder, the contig-edge rule, and the contig length needed
to report "no mobile element nearby".

Compared with a plain AMRFinderPlus wrapper:

* Every row has the gene's distance to the nearer contig end. AMRFinderPlus reports coordinates
  but not contig lengths, so its output alone cannot separate a gene in a long contig from one
  300 bp from a contig end.
* It runs before any length filter, because ARG-bearing contigs are short (one measurement puts
  MEGAHIT's at 284 bp on average).

``--plus`` is on by default and adds AMRFinderPlus' virulence and stress entries. Virulence calls go
to a separate table so they are not summed with resistance genes.
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
# -n: nucleotide input, so the report has contig coordinates (protein input, -p, would not).
#
# --plus adds AMRFinderPlus' stress-response and virulence references. Output is split again by
# Type.
#
# --database names a dated release, because `amrfinder -u` changes what "latest" points to and the
# database version dates each call.
AMRFINDER = (
    "amrfinder -n {input} --plus --database {db} --threads {threads} -o {out}" "{organism}{ident}"
)
# hAMRonization converts a backend's report into a 35-field schema shared across tools. Versions
# are passed explicitly because it cannot tell which database produced a report on disk, and the
# database version is one of Cao et al.'s L1 validity requirements.
#
# Run even with one backend, so a second backend's calls share the schema and a gene found by both
# AMRFinderPlus and RGI can be matched.
HAMRONIZE = (
    "hamronize amrfinderplus --analysis_software_version {tool_version}"
    " --reference_database_version {db_version} --input_file_name {genome}"
    " --format tsv --output {out} {report}"
)
# RGI, the optional second backend. `--local` keeps CARD in a directory owned by this run instead
# of the conda env, so runs with different CARD versions do not overwrite each other.
#
# -t contig: the input is assembled sequence, so RGI calls genes with Prodigal and reports contig
# coordinates for the context classification. -a DIAMOND: same calls as BLAST in minutes instead of
# hours. --clean removes large intermediate files that are not used.
RGI_LOAD = "rgi load --card_json {card_json} --local"
RGI_MAIN = "rgi main --local -i {genome} -o {prefix} -t contig -a DIAMOND -n {threads} --clean"
# hAMRonization for RGI's schema, with the CARD version passed explicitly.
HAMRONIZE_RGI = (
    "hamronize rgi --analysis_software_version {tool_version}"
    " --reference_database_version {db_version} --input_file_name {genome}"
    " --format tsv --output {out} {report}"
)
# argNorm maps each normalised call to an ARO term and the drugs it confers resistance to.
# hAMRonization makes the columns comparable across tools; argNorm makes the gene names comparable,
# since each database uses its own symbols, and adds the drug classes.
#
# No --hamronization_skip_unsupported_tool: it turns an unsupported backend into a warning and
# leaves its calls unmapped. An unsupported backend should fail here; the primary tables are already
# written by then.
ARGNORM = "argnorm hamronization -i {input} -o {output}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_THREADS = 1
#: Resistance and stress calls.
AMR_TABLE = "amr" + TSV_SUFFIX
#: Virulence calls, kept in a separate file.
VIRULENCE_TABLE = "virulence" + TSV_SUFFIX
#: The same calls in hAMRonization's cross-tool schema, for combining with other resistome output.
HARMONISED_TABLE = "amr_hamronized" + TSV_SUFFIX
#: The hAMRonization table with every gene mapped onto an ARO term and the drugs it confers
#: resistance to.
NORMALISED_TABLE = "amr_normalised" + TSV_SUFFIX
#: Co-located calls grouped into clusters on one contig.
CLUSTERS_TABLE = "amr_clusters" + TSV_SUFFIX
#: Directory for RGI's raw per-genome reports.
RGI_DIR = "rgi_reports"
#: Directory for the raw per-genome AMRFinderPlus reports, kept alongside the derived tables.
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
        "Use a dated release instead of 'latest' for a reproducible run",
    )
    p.add_argument(
        "--organism",
        help="AMRFinderPlus organism for point-mutation screening and common-protein suppression "
        "(e.g. Escherichia). Use it only for a bin that is already classified",
    )
    p.add_argument(
        "--min-identity",
        type=float,
        default=None,
        metavar="FRACTION",
        help="minimum identity to the reference, 0-1. Default: AMRFinderPlus' curated per-gene "
        "thresholds. A single global cutoff replaces about 7000 curated ones",
    )
    p.add_argument(
        "--rgi",
        action="store_true",
        help="also call resistance with RGI against CARD (needs `metawrap2 install-env rgi` and "
        "`install-db card`). The two backends call different genes; every row says which tool "
        "found it",
    )
    p.add_argument(
        "--card-db",
        metavar="DIR",
        help="CARD data directory holding card.json (default: CARD_DB from metawrap2.toml)",
    )
    p.add_argument(
        "--skip-amrfinder",
        action="store_true",
        help="run only the backends selected by flags. With --rgi that means RGI alone, which "
        "reports no virulence or stress genes (only AMRFinderPlus reports them)",
    )
    p.add_argument(
        "--mobile-elements",
        metavar="DIR",
        help="a `metawrap2 viruses` output directory, or a geNomad run directory, to take mobile "
        "elements from. Without it every call's context is `context_unknown` and no call can reach "
        "evidence level L2",
    )
    p.add_argument(
        "--skip-hamronize",
        action="store_true",
        help="do not write the hAMRonization-normalised table, which makes this output "
        "comparable with other tools' resistome tables. Skipping it saves space and does not "
        "change the calls",
    )
    p.add_argument(
        "--mobility-window",
        type=int,
        default=_amr.MOBILITY_WINDOW,
        metavar="BP",
        help="how far from a gene to look for a mobile element, and how much flanking sequence a "
        "'nothing nearby' call requires (default: %(default)s). There is no standard value, so it "
        "is recorded in every row",
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
    """Contig lengths by contig id, used for distance to contig end."""
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

    Read from the directory, because the tool reports whatever 'latest' resolves to at run time,
    which a rerun may not see.
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

    As in :func:`_amrfinder_version`: a version found on PATH may come from a different
    installation than the one that made the calls, and it is recorded as provenance.
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

    hamronize takes one report per call and writes a header each time, so the outputs are
    concatenated keeping the first header. Running per genome keeps ``input_file_name`` set to the
    genome name.
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

    # argNorm gets the AMRFinderPlus rows only. It rejects the whole file if it contains a tool it
    # cannot map ("argNorm can only map genes from the following tools: argsoap, deeparg, resfinder,
    # amrfinderplus, groot, abricate"), so with --rgi it produced no output. RGI calls are made
    # against CARD and already have an ARO accession (the `aro` column of `amr.tsv`).
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

    Runs with check=False and warns on failure, so the tables already written are kept.
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
            "can be compared with other tools' hAMRonization output" % HARMONISED_TABLE
        )
        return []
    rows = _amr.parse_argnorm(out_path)
    mapped = [row for row in rows if row.aro]
    comm("%d of %d call(s) mapped to an ARO term -> %s" % (len(mapped), len(rows), out_path))
    drugs = _amr.drugs_implicated(rows)
    if drugs:
        # One line for all drugs, then one per gene, not one per gene-drug pair: a single gene can
        # cover nine aminoglycosides.
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
            "unmapped calls are not false calls: argNorm only maps genes in the ontology, and "
            "recently described genes may be missing from it"
        )
    return rows


def _amrfinder_version(env: Optional[str]) -> Optional[str]:
    """AMRFinderPlus' own version, for the provenance hAMRonization records.

    Queried from the tool, because the env pins a version range and the installed version varies.

    Runs the binary in the env's bin directory, not from PATH, which may hold a different
    installation (as happened with `check`, which found a system-wide copy while the run used an
    env without it). A PATH lookup would record the wrong version in the normalised table.
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

    CARD is loaded once into a localDB owned by this run, not into the conda environment, so runs
    with different CARD versions do not overwrite each other. RGI runs from that directory because
    `--local` looks for localDB in the working directory.
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
        warning("`rgi load` failed, so RGI adds no calls; the AMRFinderPlus calls are kept")
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
                "made, since only AMRFinderPlus reports them"
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

        # None when no mobile-element caller has run, and an empty list when one ran and found
        # nothing; only the latter can support no_mge_in_window (see amr.classify_context).
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
                    "`context_unknown`: with no mobile elements to compare against, the run "
                    "cannot tell whether any are nearby" % args.mobile_elements
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
            # Joined on ARO where known (from argNorm, so this runs after it). A symbol join
            # undercounts: AMRFinderPlus' `cepA` and RGI's `CepA-29` are one ARO term.
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
                "  the backends are expected to disagree: they model resistance differently and "
                "curate their thresholds separately. RGI's extra calls are mostly Strict-tier "
                "homologs below 50% identity, a weaker match than an exact hit; the `method` "
                "column records it"
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
