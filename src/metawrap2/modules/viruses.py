"""MetaWrap2 viruses: recover viral contigs, estimate their completeness, and group them into vOTUs.

Binning keeps prokaryotes and discards phage, which is often the most abundant thing in a sample.
geNomad finds the viral contigs, and CheckV estimates how complete each one is and trims host
sequence off proviruses: the viral counterpart of ``binning`` plus ``bin_refinement``.

* geNomad runs on the assembly (``-a``), because phage contigs are usually in no prokaryotic bin.
  Bins (``-b``) are accepted for finding prophages integrated in binned genomes.
* geNomad's plasmid calls from the same run are written out too, for the ``amr`` module's mobility
  analysis.

See :mod:`metawrap2.viruses` for the quality tiers, why an unreliable completeness estimate is
reported as ``undetermined``, and CRISPR host linking.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from .. import viruses as _viruses
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
    require_nonempty_dir,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["viruses"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# end-to-end runs the whole geNomad chain: gene calling, marker annotation, the classifier, and the
# summary. --cleanup removes its large intermediates, which nothing here reads.
#
# --min-score is not set by default, so geNomad's FDR-calibrated cutoff applies. --min-virus-score
# sets it.
GENOMAD = "genomad end-to-end --cleanup --threads {threads} {score} {input} {outdir} {db}"
# skani over the viral contigs, for vOTU clustering. Same tool and table format as `dereplicate`,
# with viral thresholds: 95% ANI over 85% of the shorter sequence instead of 95% over 10%. -E writes
# one row per pair instead of a full matrix.
SKANI_VIRAL = (
    "skani triangle -E -l {genome_list} -o {output} -t {threads} -s {screen} --min-af {min_af}"
)
# cctyper over the prokaryotic bins, to find CRISPR arrays. --prodigal meta because a MAG has many
# contigs, while cctyper's single-genome gene calling assumes one closed chromosome.
#
# --db is required: bioconda's cctyper 1.9 activation script exports
# CCTYPER_DB=$CONDA_PREFIX/cct_data, a 1.8-era directory without the profile database 1.9 needs,
# so cctyper exits with "Could not find HMM profile database" unless pointed at the data inside
# the Python package.
CCTYPER = "cctyper {input} {outdir} --prodigal meta -t {threads} --db {db} --no_plot"
# CheckV's end_to_end runs contamination, completeness, repeats and the quality summary in order.
# -d names the database directory, so metawrap2.toml decides it instead of a CHECKVDB variable in
# the environment.
CHECKV = "checkv end_to_end {input} {outdir} -t {threads} -d {db}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_THREADS = 1
#: One row per viral sequence, with its tier and the flags that qualify it.
VIRAL_TABLE = "viral_contigs" + TSV_SUFFIX
#: geNomad's plasmid calls, kept because the amr module's mobility work consumes them.
PLASMID_TABLE = "plasmid_contigs" + TSV_SUFFIX
#: The viral sequences themselves, after CheckV's provirus trimming.
VIRAL_FASTA = "viruses.fasta"
#: One row per viral population (vOTU). Counting contigs instead would count one population once
#: per contig.
VOTU_TABLE = "votus" + TSV_SUFFIX
#: One row per piece of evidence that a host was infected: a CRISPR spacer from a host genome found
#: in a viral contig.
HOST_LINK_TABLE = "host_links" + TSV_SUFFIX
#: Where cctyper's own output is kept, one directory per host genome.
CRISPR_DIR = "crispr"
#: Where each tool's own output is kept, alongside the derived tables.
GENOMAD_DIR = "genomad"
CHECKV_DIR = "checkv"
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 viruses",
        usage="metawrap2 viruses [options] -o output_dir (-a assembly.fasta | -b bin_folder)",
        description="Recover viral contigs, estimate their completeness, and group them into vOTUs.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "-a",
        "--assembly",
        help="an assembly. Use this to find free phage, which is most of the viral signal: "
        "phage contigs are usually left out of prokaryotic bins",
    )
    source.add_argument(
        "-b",
        "--bins",
        help="a bin folder, for finding proviruses integrated in binned genomes. Bins alone "
        "miss free phage",
    )
    threads_arg(p)
    p.add_argument(
        "--genomad-db",
        help="path to the geNomad database (default: GENOMAD_DB from metawrap2.toml)",
    )
    p.add_argument(
        "--checkv-db",
        help="path to the CheckV database (default: CHECKV_DB from metawrap2.toml)",
    )
    p.add_argument(
        "--min-length",
        type=int,
        default=_viruses.MIN_VIRAL_LENGTH,
        metavar="BP",
        help="ignore contigs shorter than this (default: %(default)s). Shorter virus calls are "
        "fragments CheckV cannot place",
    )
    p.add_argument(
        "--min-virus-score",
        type=float,
        default=None,
        metavar="SCORE",
        help="geNomad minimum virus score, 0-1 (default: geNomad's own FDR-calibrated "
        "cutoff, which is recommended)",
    )
    p.add_argument(
        "--skip-votus",
        action="store_true",
        help="do not cluster the viral contigs into populations. Counts are then of contigs, "
        "which overestimates the number of distinct viruses",
    )
    p.add_argument(
        "--crispr-hosts",
        metavar="BIN_FOLDER",
        help="predict hosts by matching CRISPR spacers from these prokaryotic genomes against the "
        "viral contigs. A spacer records a past infection, so a match is direct evidence of the "
        "host, unlike a prediction from sequence composition",
    )
    p.add_argument(
        "--untrusted-arrays",
        action="store_true",
        help="also use spacers from CRISPR arrays cctyper does not trust. A spurious array "
        "contributes host sequence as 'spacers', which can give false host matches",
    )
    p.add_argument(
        "--skip-checkv",
        action="store_true",
        help="run geNomad only. The viral set then has no completeness estimates, so every "
        "contig is tier `undetermined`",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _inputs(args) -> List[Tuple[str, str]]:
    if args.assembly:
        return [(os.path.splitext(os.path.basename(args.assembly))[0], args.assembly)]
    return [
        (os.path.splitext(name)[0], os.path.join(args.bins, name)) for name in list_bins(args.bins)
    ]


def _database(path: Optional[str], settings, key: str, what: str) -> str:
    db = path or settings.db(key)
    if not db:
        error("no %s database: pass the flag or set %s in metawrap2.toml" % (what, key))
    return db


def _genomad_summary(outdir: str, name: str, kind: str) -> Optional[str]:
    """geNomad's summary path. It nests under ``<name>_summary/``, which is easy to get wrong."""
    candidate = os.path.join(outdir, "%s_summary" % name, "%s_%s_summary.tsv" % (name, kind))
    return candidate if os.path.isfile(candidate) else None


def _write_plasmids(path: str, rows: Dict[str, Dict[str, str]], genome: str, append: bool) -> int:
    # No taxonomy column: geNomad's plasmid summary does not have one, unlike its virus summary.
    columns = (
        "genome",
        "seq_name",
        "length",
        "topology",
        "plasmid_score",
        "n_genes",
        "n_hallmarks",
    )
    mode = "a" if append else "w"
    with open(path, mode, newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        if not append:
            writer.writeheader()
        for seq_name, row in rows.items():
            writer.writerow(
                {
                    "genome": genome,
                    "seq_name": seq_name,
                    "length": row.get("length", ""),
                    "topology": row.get("topology", ""),
                    "plasmid_score": row.get("plasmid_score", ""),
                    "n_genes": row.get("n_genes", ""),
                    "n_hallmarks": row.get("n_hallmarks", ""),
                }
            )
    return len(rows)


def _cluster_votus(args, env, viral_fasta: str, contigs: List[_viruses.ViralContig]) -> None:
    """Cluster the viral sequences into populations and write the vOTU table.

    skani compares files, so each sequence is written to its own fasta first (161 files on one
    assembly). This reuses `dereplicate`'s skani parser.
    """
    from ..modules.dereplicate import parse_skani

    work = os.path.join(args.output, "work_files", "votu")
    ensure_dir(work)
    lengths: Dict[str, int] = {}
    paths: List[str] = []
    for header, sequence in iter_fasta(viral_fasta):
        name = contig_id(header)
        lengths[name] = len(sequence)
        # Named by index, since contig ids can contain characters awkward in a filename. Names
        # are mapped back from the files after skani runs.
        path = os.path.join(work, "seq_%05d.fna" % len(paths))
        with open(path, "w") as out:
            out.write(">%s\n%s\n" % (name, sequence))
        paths.append(path)

    if len(paths) < 2:
        comm("only %d viral sequence(s); nothing to cluster" % len(paths))
        return

    listing = os.path.join(work, "sequences.txt")
    with open(listing, "w") as out:
        out.writelines(path + "\n" for path in paths)
    matrix = os.path.join(work, "ani.tsv")
    rc = run(
        SKANI_VIRAL.format(
            genome_list=listing,
            output=matrix,
            threads=args.threads,
            screen=_viruses.VOTU_ANI - 10,
            min_af=_viruses.VOTU_COVERAGE,
        ),
        env=env_for("dereplicate", load_settings(args.config)),
        tool="skani triangle",
        check=False,
    )
    if rc != 0 or not os.path.isfile(matrix):
        warning(
            "skani failed, so the viral contigs are not clustered. Counts in %s are of contigs, "
            "which overestimates the number of distinct viruses" % VIRAL_TABLE
        )
        return

    # skani reports the file paths it was given; map them back to the sequence names.
    by_path = {os.path.realpath(p): contig_id(next(iter_fasta(p))[0]) for p in paths}

    def rename(value: str) -> str:
        return by_path.get(os.path.realpath(value), value)

    edges = [
        edge._replace(ref=rename(edge.ref), query=rename(edge.query))
        for edge in parse_skani(matrix)
    ]
    votus = _viruses.build_votus(sorted(lengths), edges, lengths)
    table = os.path.join(args.output, VOTU_TABLE)
    with open(table, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=_viruses.VOTU_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for votu in votus:
            writer.writerow(_viruses.votu_row(votu))
    multi = sum(1 for v in votus if v.size > 1)
    comm(
        "%d viral sequence(s) -> %d vOTU(s) at %.0f%% ANI over %.0f%% of the shorter; %d population(s) "
        "have more than one contig"
        % (len(lengths), len(votus), _viruses.VOTU_ANI, _viruses.VOTU_COVERAGE, multi)
    )
    comm("-> %s" % table)


def _cctyper_db(env: Optional[str]) -> Optional[str]:
    """Where cctyper's packaged data lives inside *env*, asked of the env's own Python.

    Needed because bioconda's activation script points CCTYPER_DB at a stale data directory (see
    the CCTYPER command above). Found by importing cctyper in the environment, since the path
    includes the Python version.
    """
    prefix = env_prefix(env) if env else None
    python = os.path.join(prefix, "bin", "python") if prefix else sys.executable
    if not os.path.isfile(python):
        return None
    try:
        out = subprocess.run(
            [
                python,
                "-c",
                "import cctyper, os; print(os.path.join(os.path.dirname(cctyper.__file__), 'data'))",
            ],
            capture_output=True,
            text=True,
            check=False,
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    path = out.stdout.strip()
    return path if path and os.path.isdir(path) else None


def _crispr_host_links(args, env, viral_fasta: str) -> None:
    """Run cctyper over the host genomes and write every spacer match against the viral set.

    Covers only some of the viruses: a link needs a host with a CRISPR system that has recorded
    this virus. A missing link is not evidence of a missing host, so the table has one row per
    match, not per virus.
    """
    bins = list_bins(args.crispr_hosts)
    if not bins:
        warning("no genomes in %s, so no host linking" % args.crispr_hosts)
        return
    viral_sequences = {contig_id(header): seq for header, seq in iter_fasta(viral_fasta)}
    if not viral_sequences:
        comm("no viral sequences, so nothing to link a host to")
        return

    db = _cctyper_db(env)
    if db is None:
        warning(
            "could not locate cctyper's packaged database inside the environment, so host linking "
            "is skipped. (bioconda's cctyper build sets CCTYPER_DB to a data directory without "
            "the profiles version 1.9 needs, so cctyper would fail here anyway)"
        )
        return

    root = os.path.join(args.output, CRISPR_DIR)
    ensure_dir(root)
    spacers_by_host: Dict[str, Mapping[str, Sequence[str]]] = {}
    subtypes: Dict[str, str] = {}
    arrays_found = 0
    for name in bar(bins, desc="cctyper", unit=" genome", total=len(bins)):
        host = os.path.splitext(name)[0]
        outdir = os.path.join(root, host)
        rc = run(
            CCTYPER.format(
                input=os.path.join(args.crispr_hosts, name),
                outdir=outdir,
                threads=args.threads,
                db=db,
            ),
            env=env,
            tool="cctyper",
            check=False,
        )
        if rc != 0:
            warning("cctyper failed on %s; that genome contributes no spacers" % host)
            continue
        arrays = _viruses.read_crispr_arrays(outdir, host=host)
        arrays_found += sum(1 for a in arrays if a.trusted or args.untrusted_arrays)
        subtypes.update({a.array: a.subtype for a in arrays})
        spacers = _viruses.read_spacers(outdir, trusted_only=not args.untrusted_arrays)
        if spacers:
            spacers_by_host[host] = spacers

    total_spacers = sum(len(s) for arrays in spacers_by_host.values() for s in arrays.values())
    comm(
        "%d CRISPR array(s) with %d spacer(s) across %d of %d genome(s)"
        % (arrays_found, total_spacers, len(spacers_by_host), len(bins))
    )
    if not spacers_by_host:
        return

    links = _viruses.host_links(spacers_by_host, viral_sequences, subtypes=subtypes)
    table = os.path.join(args.output, HOST_LINK_TABLE)
    with open(table, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=_viruses.HOST_LINK_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for link in sorted(links, key=lambda link: (link.host, link.viral_contig)):
            writer.writerow(_viruses.host_link_row(link))
    linked_viruses = {link.viral_contig for link in links}
    comm(
        "%d spacer match(es) linking %d host(s) to %d of %d viral sequence(s) -> %s"
        % (
            len(links),
            len({link.host for link in links}),
            len(linked_viruses),
            len(viral_sequences),
            table,
        )
    )
    if not links:
        comm(
            "no matches. This is common for a single sample: a link needs the host to carry a "
            "CRISPR system with a spacer from this virus"
        )


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
    validate_inputs("viruses", checks)

    env = env_for("viruses", settings)
    genomad_db = _database(args.genomad_db, settings, "GENOMAD_DB", "geNomad")
    checkv_db = (
        None if args.skip_checkv else _database(args.checkv_db, settings, "CHECKV_DB", "CheckV")
    )
    rec = start_run("viruses", args, env, settings, inputs=[source])

    ckpt = make_checkpoint(args.output)
    try:
        announcement("BEGIN VIRUSES PIPELINE!")
        ensure_dir(args.output)
        genomad_root = os.path.join(args.output, GENOMAD_DIR)
        ensure_dir(genomad_root)

        todo = _inputs(args)
        if not todo:
            error("No fasta input found (%s). Exiting." % "/".join(FASTA_EXTENSIONS))
        comm("%d input(s) to screen" % len(todo))

        score = "--min-score %g" % args.min_virus_score if args.min_virus_score is not None else ""
        if ckpt.todo("genomad"):
            for name, path in bar(todo, desc="genomad", unit=" input", total=len(todo)):
                rc = run(
                    GENOMAD.format(
                        threads=args.threads,
                        score=score,
                        input=path,
                        outdir=os.path.join(genomad_root, name),
                        db=genomad_db,
                    ),
                    env=env,
                    tool="genomad",
                    check=False,
                )
                if rc != 0:
                    warning("geNomad failed on %s. Proceeding with the rest." % name)
            ckpt.done("genomad")
        else:
            comm("skipping geNomad (already done; --resume)")

        if dry_run():
            announcement("VIRUSES PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        viral_fasta = os.path.join(args.output, VIRAL_FASTA)
        plasmid_table = os.path.join(args.output, PLASMID_TABLE)
        per_input: Dict[str, Dict[str, Dict[str, str]]] = {}
        n_plasmids = 0
        with open(viral_fasta, "w") as out:
            for index, (name, _path) in enumerate(todo):
                outdir = os.path.join(genomad_root, name)
                summary = _genomad_summary(outdir, name, "virus")
                if not summary:
                    warning("no geNomad virus summary for %s; skipping it" % name)
                    continue
                rows = _viruses.parse_genomad(summary)
                per_input[name] = {
                    seq: row
                    for seq, row in rows.items()
                    if int(float(row.get("length", 0) or 0)) >= args.min_length
                }
                fasta = os.path.join(outdir, "%s_summary" % name, "%s_virus.fna" % name)
                if os.path.isfile(fasta):
                    for header, sequence in iter_fasta(fasta):
                        if contig_id(header) in per_input[name]:
                            out.write(">%s\n%s\n" % (contig_id(header), sequence))
                plasmids = _genomad_summary(outdir, name, "plasmid")
                if plasmids:
                    n_plasmids += _write_plasmids(
                        plasmid_table, _viruses.parse_genomad(plasmids), name, index > 0
                    )

        total_viral = sum(len(rows) for rows in per_input.values())
        comm(
            "%d viral sequence(s) at least %d bp -> %s"
            % (total_viral, args.min_length, viral_fasta)
        )
        if n_plasmids:
            comm(
                "%d plasmid call(s) -> %s (the amr module's mobility input)"
                % (n_plasmids, plasmid_table)
            )
        if not total_viral:
            comm("nothing viral to score; finishing here")
            announcement("VIRUSES PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        checkv_rows: Dict[str, Dict[str, str]] = {}
        if not args.skip_checkv:
            checkv_out = os.path.join(args.output, CHECKV_DIR)
            if ckpt.todo("checkv"):
                rc = run(
                    CHECKV.format(
                        input=viral_fasta,
                        outdir=checkv_out,
                        threads=args.threads,
                        db=checkv_db,
                    ),
                    env=env,
                    tool="checkv",
                    check=False,
                )
                if rc != 0:
                    warning(
                        "CheckV failed. The viral set is still written, but every contig will be "
                        "tier `undetermined` because nothing estimated its completeness."
                    )
                ckpt.done("checkv")
            else:
                comm("skipping CheckV (already done; --resume)")
            quality = os.path.join(checkv_out, "quality_summary.tsv")
            if os.path.isfile(quality):
                checkv_rows = _viruses.parse_checkv(quality)

        table = os.path.join(args.output, VIRAL_TABLE)
        contigs: List[_viruses.ViralContig] = []
        with open(table, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=_viruses.VIRAL_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for name, rows in per_input.items():
                merged = _viruses.merge(name, rows, checkv_rows)
                contigs.extend(merged)
                for contig in sorted(merged, key=lambda c: -c.length):
                    writer.writerow(_viruses.viral_row(contig))
        comm("%d row(s) -> %s" % (len(contigs), table))

        if not args.skip_votus and contigs:
            _cluster_votus(args, env, viral_fasta, contigs)

        if args.crispr_hosts:
            _crispr_host_links(args, env, viral_fasta)

        counts = _viruses.summarise(contigs)
        comm(
            "tiers: %d complete, %d high, %d medium, %d low, %d undetermined"
            % (
                counts["complete"],
                counts["high_quality"],
                counts["medium_quality"],
                counts["low_quality"],
                counts["undetermined"],
            )
        )
        comm(
            "%d provirus(es) trimmed of host sequence; %d contig(s) have a completeness estimate "
            "CheckV does not trust, which usually means a virus with no close relative in the "
            "database" % (counts["proviruses"], counts["unreliable_completeness"])
        )
        announcement("VIRUSES PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
