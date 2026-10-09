"""MetaWrap2 coverage module: one cross-sample depth table for a whole cohort, in N alignments.

Multi-sample binning wants every contig's depth in every sample, and the obvious way to get that
costs N^2 alignments: map all N samples' reads against each of the N per-sample assemblies. For a
twenty-sample cohort that is four hundred alignments to produce twenty depth tables.

The multi-split design gets the same tables for N. Concatenate the per-sample assemblies into one
contig catalogue, map each sample's reads once against that catalogue, build one depth table over
it, then split the table back into one slice per source assembly. Each sample's binning then runs
on its own contigs with every sample's coverage, via ``metawrap2 binning --depth``.

What makes the split safe is the naming. Contig names collide across assemblies - every assembler
numbers from one - so each contig is prefixed with its assembly's label before concatenation, and
the prefix is what the split demultiplexes on. That is only unambiguous if no original name
contains the separator, so this module refuses to build a catalogue where one does rather than
guessing later; see SEPARATOR.

The slices carry the original contig names, not the prefixed ones, because they are handed to a
binning run over the original assembly. A table of prefixed names against an unprefixed assembly is
exactly the well-formed-but-wrong input that `binning`'s own check exists to catch.

How the depth itself is measured is a separate choice from all of that: ``--coverage align`` maps
the reads (the default, and what every binner's published accuracy was measured with), while
``--coverage fairy`` estimates depth from k-mer sketches and never aligns anything. The table, the
split and everything downstream are identical either way, because fairy writes the same
contig x sample layout metaBAT2 reads; see COVERAGE_ROUTES for what the
cheaper route gives up.

This file is meant to be read and edited: the commands live in the COMMANDS block below and the
tunable values in the CONSTANTS block.
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple

from .. import depth as _depth
from ..aligners import ALIGNERS, DEFAULT_ALIGNER, INDEX_SUFFIXES
from ..cliformat import make_parser, percent
from ..config import MODULE_ENVS, load_settings
from ..io.seqio import contig_id, is_binnable, iter_fasta
from ._common import (
    absolutize_paths,
    align_read_sets,
    announcement,
    check_fasta,
    check_fastq,
    collect_read_pairs,
    collect_reads,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    expect_produced,
    finish_run,
    fs_remove,
    listdir,
    make_checkpoint,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)
from .binning import DEPTH_HEADER, depth_means

#: Shared with `binning`: the aligner and samtools this needs are
#: exactly the ones that module installs, and a second environment holding the same four tools would
#: be a second solve to keep current for no gain. `install-env coverage` therefore builds (or finds)
#: metawrap2-bins - see _env_yaml, which resolves the yaml from the env name for this reason.
CONDA_ENV = MODULE_ENVS["coverage"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
#: The k-mer route (--coverage fairy): sketch each sample's reads once, then estimate every contig's
#: depth from the sketches. ``-S`` names the sketches after the samples, which is what puts sample
#: names in the table's columns where the alignment route puts BAM paths. Defaults left alone on
#: purpose: ``-m 95`` is fairy's own warning not to lower it, and ``-c``/``-k`` are what its published
#: accuracy was measured at.
FAIRY_SKETCH_PAIRED = "fairy sketch -1 {r1} -2 {r2} -S {names} -d {sketches} -t {threads}"
FAIRY_SKETCH_SINGLE = "fairy sketch -r {reads} -S {names} -d {sketches} -t {threads}"
FAIRY_COVERAGE = "fairy coverage -t {threads} {sketches} {catalogue} -o {depth}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: How coverage is measured. ``align`` maps the reads and reads the alignments, which is what every
#: binner's published accuracy was measured with. ``fairy`` estimates depth from k-mer sketches
#: instead and skips alignment entirely - minutes rather than hours for a cohort - at the cost of
#: resolving nothing below ~95% ANI, so a strain that differs from the contig it came from is
#: credited to it or to nothing.
#:
#: ``align`` is the default on measurement, not on caution, and the measurement is now the one that
#: matters. Fed to metaBAT2 on one CAMI2 GI sample, single-sample fairy scored AMBER F1 0.212 against
#: minimap2's 0.452 - half the bins - and raising its ANI threshold made agreement worse, not better.
#: What that could not settle was fairy's actual claim, a multi-sample one, which needs the
#: cross-sample table this module exists to build.
#:
#: Settled by benchmark B5 over three CAMI2 samples, every other token held identical: fairy builds the
#: cohort table for 0.04 cpu-hours against minimap2's 2.26, 57x cheaper, and the catalogue that comes
#: out holds 27 high-quality MAGs against 37, with 76% of its bins passing GUNC against 87%. In the
#: cohort setting it was supposed to win in, it loses 10 genomes out of 37.
#:
#: So the route stays, opt-in, and is closed as a default candidate rather than merely unmeasured. Note
#: what the loss is *not*: the variance ablation showed metaBAT2's variance column carries no signal at
#: this sample count, so what fairy gets wrong is its mean depths (Pearson 0.237 against bwa), and an
#: alignment-free estimator with accurate means would still be worth having.
COVERAGE_ROUTES = ("align", "fairy")
DEFAULT_ROUTE = "align"
#: Where the fairy route keeps its per-sample sketches, inside the work directory.
SKETCH_DIR = "sketches"
#: Separates an assembly's label from the original contig name in the catalogue. ``~`` because it is
#: legal in a fasta header, has no meaning to any aligner or to samtools, and essentially never
#: appears in an assembler's contig names - SPAdes, MEGAHIT and Flye all use letters, digits, ``_``
#: and ``.``. It is a flag rather than a constant only so that a cohort whose names do contain it has
#: a way out; the names are checked either way, because a separator inside a name would make the
#: split silently cut in the wrong place.
SEPARATOR = "~"
#: What the catalogue and the tables are called inside the output directory.
CATALOGUE = "catalogue.fa"
CATALOGUE_DEPTH = "catalogue_depth.txt"
DEPTH_DIR = "depth"
WORK_DIR = "work_files"
#: Written beside the slices so a later run - or a human - can map a prefixed catalogue contig back
#: to its assembly without re-deriving the labels from the filenames.
ASSEMBLIES_TABLE = "assemblies.tsv"
# ───────────────────────────────────────────────────────────────────────────────────────────


class Assembly(NamedTuple):
    """One sample's assembly: the label its contigs are prefixed with, and where it is."""

    label: str
    path: str


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 coverage",
        usage="metawrap2 coverage [options] -o output_dir -a asm_1 [-a asm_2 ...] reads...",
        description="Build one cross-sample depth table over a cohort's assemblies, and split it "
        "back into one table per assembly for `metawrap2 binning --depth`.",
    )
    p.add_argument(
        "-a",
        "--assembly",
        action="append",
        default=[],
        required=True,
        metavar="[LABEL=]FASTA",
        help="an assembly to include; repeat once per sample. Give LABEL=path to name it, "
        "otherwise the filename is the label",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    threads_arg(p)
    p.add_argument(
        "--aligner",
        choices=sorted(ALIGNERS),
        default=DEFAULT_ALIGNER,
        help="aligner used to measure coverage (default %s; see ALIGNERS for why)"
        % DEFAULT_ALIGNER,
    )
    p.add_argument(
        "--coverage",
        choices=COVERAGE_ROUTES,
        default=DEFAULT_ROUTE,
        help="how depth is measured: align the reads (default) or estimate it from k-mer "
        "sketches with fairy, which is far faster and resolves nothing below ~95%% ANI",
    )
    p.add_argument(
        "--separator",
        default=SEPARATOR,
        help="character joining an assembly's label to its contig names in the catalogue "
        "(default %s)" % SEPARATOR,
    )
    p.add_argument(
        "--keep-bams",
        action="store_true",
        help="keep the per-sample BAMs in work_files/ after the depth table is built (about 5 GB "
        "per sample; discarded by default)",
    )
    p.add_argument(
        "--keep-catalogue",
        action="store_true",
        help="keep the concatenated catalogue and its aligner index (the size of all the "
        "assemblies together; discarded by default once the tables are written)",
    )
    p.add_argument(
        "--min-identity",
        type=percent(),
        default=_depth.PERCENT_IDENTITY,
        metavar="PERCENT",
        help="minimum alignment identity for a read to count toward a contig's depth. 90 beats jgi's "
        "97 on every measure taken: 88 bins against 74, more pure, and 47 clean high-quality MAGs "
        "against 35 - because on strain-rich data a tight floor removes signal rather than separating "
        "strains. Pass 97 to reproduce a jgi depth table, or 0 for no filter",
    )
    p.add_argument("--single-end", action="store_true", help="reads are single-end (*.fastq)")
    p.add_argument("--interleaved", action="store_true", help="reads are interleaved paired-end")
    p.add_argument("--config", help="path to metawrap2.toml")
    p.add_argument(
        "reads",
        nargs="*",
        help="every sample's read files (paired *_1.fastq/*_2.fastq, or single/interleaved)",
    )
    return p.parse_args(argv)


# --- the parts that are pure, and therefore testable without an aligner -------------------


def parse_assembly(spec: str) -> Assembly:
    """Read one ``-a`` value into a label and a path.

    ``label=path`` is explicit; a bare path is labelled by its filename with the extension taken
    off, which for a cohort laid out as ``sample_0/final_assembly.fasta`` would label every sample
    the same - hence the explicit form, and hence the duplicate check in :func:`parse_assemblies`
    rather than silently appending a number. A label that is wrong is worth an error, because it
    ends up in the contig names of every table downstream.
    """
    if "=" in spec:
        label, _sep, path = spec.partition("=")
        label, path = label.strip(), path.strip()
        if not label or not path:
            raise ValueError(
                "%r is not a LABEL=path assembly: both sides of the '=' have to be given." % spec
            )
        return Assembly(label, path)
    base = os.path.basename(os.path.normpath(spec))
    for ext in (".gz", ".bz2", ".zst"):
        base = base.removesuffix(ext)
    label = os.path.splitext(base)[0]
    if not label:
        raise ValueError("cannot label the assembly %r; give it one as LABEL=%s" % (spec, spec))
    return Assembly(label, spec)


def parse_assemblies(specs: Sequence[str], separator: str = SEPARATOR) -> List[Assembly]:
    """Every ``-a`` value, with the labels checked for the things that would break the split.

    Two assemblies sharing a label cannot be told apart afterwards, and a label containing the
    separator makes the prefix ambiguous in the same way a contig name containing it does. Both are
    errors here rather than surprises in the split, which would have no way to report them.
    """
    assemblies = [parse_assembly(spec) for spec in specs]
    seen: Dict[str, str] = {}
    for assembly in assemblies:
        if separator in assembly.label:
            raise ValueError(
                "the label %r contains %r, which separates a label from a contig name in the "
                "catalogue.\nGive the assembly another label as LABEL=%s, or pick a different "
                "--separator." % (assembly.label, separator, assembly.path)
            )
        if assembly.label in seen:
            raise ValueError(
                "two assemblies are both labelled %r (%s and %s), so their contigs could not be "
                "told apart.\nName them with -a LABEL=path."
                % (assembly.label, seen[assembly.label], assembly.path)
            )
        seen[assembly.label] = assembly.path
    return assemblies


def catalogue_name(label: str, contig: str, separator: str = SEPARATOR) -> str:
    """The name a contig takes in the catalogue."""
    return "%s%s%s" % (label, separator, contig)


def split_catalogue_name(name: str, separator: str = SEPARATOR) -> Tuple[str, str]:
    """Undo :func:`catalogue_name`. Returns ``(label, original contig)``.

    Split on the *first* separator, which is the one this module wrote: a name is only allowed into
    the catalogue if it does not contain the separator, so everything after the first one belongs to
    the contig. A name with no separator at all belongs to no label and comes back as ``("", name)``,
    which the caller reports rather than guesses about.
    """
    label, sep, contig = name.partition(separator)
    if not sep:
        return "", name
    return label, contig


def write_catalogue(
    assemblies: Sequence[Assembly],
    out_path: str,
    *,
    separator: str = SEPARATOR,
) -> Dict[str, int]:
    """Concatenate *assemblies* into one fasta, prefixing every contig with its label.

    Returns each label's contig count. Prefixing unconditionally, including for a single assembly:
    the split then has one code path instead of two, and the one-assembly case is the one a test
    would otherwise never exercise.

    A contig name containing the separator is refused, because the split demultiplexes on the first
    one - a name like ``k99~3`` under label ``s1`` would come back as belonging to label ``s1`` and
    contig ``k99~3``, which is right, but a name containing it in the *label* position could not be
    distinguished at all. Checking the names is cheap and the alternative is a depth table quietly
    attributed to the wrong sample.

    Records with no sequence, and records with no A, C, G, T or U anywhere in them, are left out -
    the same two classes ``binning.copy_assembly`` drops, for the same reasons, and the agreement is
    the point. A depth table carrying rows for contigs the binner's own copy of the assembly does not
    have is a table whose rows do not match its fasta, and every binner reacts to that differently.
    """
    counts: Dict[str, int] = {}
    with open(out_path, "w") as out:
        for assembly in assemblies:
            count = 0
            for header, seq in iter_fasta(assembly.path):
                name = contig_id(header)
                if separator in name:
                    raise ValueError(
                        "the contig %r in %s contains %r, which this module uses to join an "
                        "assembly's label to its contig names.\nPass another one with "
                        "--separator (any character the names do not use)."
                        % (name, assembly.path, separator)
                    )
                if not is_binnable(seq):
                    continue
                out.write(">%s\n%s\n" % (catalogue_name(assembly.label, name, separator), seq))
                count += 1
            if not count:
                raise ValueError("%s contains no contigs." % assembly.path)
            counts[assembly.label] = count
    return counts


def split_depth_table(
    table: str,
    out_dir: str,
    labels: Iterable[str],
    *,
    separator: str = SEPARATOR,
) -> Dict[str, str]:
    """Write one depth table per label from the catalogue-wide *table*.

    Each slice keeps the header - the sample columns are the whole point, and every slice carries
    all of them - and restores the original contig names, so the slice describes the assembly it
    will be binned against rather than the catalogue.

    Rows are streamed and dispatched by prefix in one pass rather than read per label, because the
    catalogue table has a row for every contig of every sample and reading it N times is the cost
    this module exists to avoid. Labels with no rows are reported by the caller: an empty slice
    means the catalogue did not contain that assembly, which is a real failure and not an empty file
    to hand to a binner.
    """
    ensure_dir(out_dir)
    handles = {}
    paths: Dict[str, str] = {}
    written: Dict[str, int] = {label: 0 for label in labels}
    try:
        with open(table) as fh:
            header = fh.readline()
            if tuple(header.rstrip("\n").split("\t")[:3]) != DEPTH_HEADER:
                raise ValueError(
                    "%s is not a contig depth table: its first three columns "
                    "are %s, not %s."
                    % (
                        table,
                        ", ".join(header.rstrip("\n").split("\t")[:3]) or "(empty)",
                        ", ".join(DEPTH_HEADER),
                    )
                )
            for label in written:
                path = os.path.join(out_dir, "%s.txt" % label)
                paths[label] = path
                # Held open for the whole pass: the point is one read of the catalogue table.
                handles[label] = open(path, "w")  # noqa: SIM115
                handles[label].write(header)
            unknown = set()
            for line in fh:
                if not line.strip():
                    continue
                name, tab, rest = line.partition("\t")
                label, contig = split_catalogue_name(name, separator)
                if label not in handles:
                    unknown.add(label)
                    continue
                handles[label].write(contig + tab + rest)
                written[label] += 1
    finally:
        for handle in handles.values():
            handle.close()
    if unknown:
        # Not fatal: a catalogue can legitimately hold more than this run's assemblies if somebody
        # reused one. Worth saying, because the usual cause is a label typo and the symptom of that
        # is a slice that is short rather than missing.
        warning(
            "%d contig(s) in %s carry a label this run does not know (%s) and were not written to "
            "any slice" % (sum(1 for _ in unknown), table, ", ".join(sorted(unknown)[:5]) or "none")
        )
    empty = [label for label, rows in written.items() if not rows]
    if empty:
        raise ValueError(
            "no contig in %s belongs to %s, so there is no depth table to bin %s with."
            % (table, ", ".join(empty), "them" if len(empty) > 1 else "it")
        )
    for label, rows in sorted(written.items()):
        comm("%s: %d contig(s)" % (label, rows))
    return paths


# --- the orchestration --------------------------------------------------------------------


def _write_assemblies_table(
    assemblies: Sequence[Assembly], path: str, depth: Dict[str, str]
) -> None:
    with open(path, "w") as fh:
        fh.write("label\tassembly\tdepth_table\n")
        fh.writelines(
            "%s\t%s\t%s\n" % (assembly.label, assembly.path, depth.get(assembly.label, ""))
            for assembly in assemblies
        )


def sketch_reads(
    reads: Sequence[str],
    sketch_dir: str,
    *,
    threads: int,
    env: Optional[str] = None,
    single_end: bool = False,
    interleaved: bool = False,
) -> List[str]:
    """Sketch every sample's reads for the fairy route, returning the sketch paths in sample order.

    One invocation for every sample rather than one each: fairy sketches samples concurrently, and
    the names it writes come from ``-S``, which has to list them all in the same order as the reads.

    An interleaved file is sketched as single-end, which is not a compromise: the sketch is a
    subsample of the read k-mers, and which mate a k-mer came from is not part of it. The alignment
    route needs the pairing because an aligner does; this does not.
    """
    ensure_dir(sketch_dir)
    sets = (
        collect_reads(list(reads))
        if (single_end or interleaved)
        else collect_read_pairs(list(reads))
    )
    names = [entry[0] for entry in sets]
    sketches = [os.path.join(sketch_dir, "%s.paired.bcsp" % name) for name in names]
    if single_end or interleaved:
        # fairy names a single-end sketch differently from a paired one, and the coverage call takes
        # the paths rather than the directory, so the two cases cannot share one suffix.
        sketches = [os.path.join(sketch_dir, "%s.single.bcsp" % name) for name in names]
        command = FAIRY_SKETCH_SINGLE.format(
            reads=" ".join(entry[1] for entry in sets),
            names=" ".join(names),
            sketches=sketch_dir,
            threads=threads,
        )
    else:
        command = FAIRY_SKETCH_PAIRED.format(
            r1=" ".join(entry[1] for entry in sets),
            r2=" ".join(entry[2] for entry in sets),
            names=" ".join(names),
            sketches=sketch_dir,
            threads=threads,
        )
    comm("sketching %d sample(s) once each; nothing is aligned on this route" % len(sets))
    run(command, env=env, tool="fairy sketch")
    return sketches


def _bams(work: str) -> List[str]:
    return [os.path.join(work, f) for f in sorted(listdir(work)) if f.endswith(".bam")]


def _sketches(sketch_dir: str) -> List[str]:
    """The sketches a previous run left behind, for a --resume that skips the sketching stage."""
    return [os.path.join(sketch_dir, f) for f in sorted(listdir(sketch_dir)) if f.endswith(".bcsp")]


def _sample_count(table: str) -> int:
    with open(table) as fh:
        return len(depth_means(fh.readline().rstrip("\n").split("\t")))


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    try:
        assemblies = parse_assemblies(args.assembly, args.separator)
    except ValueError as exc:
        error(str(exc))
    if not args.reads:
        error(
            "No read files given. Coverage is measured from reads; pass every sample's reads as "
            "positional arguments."
        )

    settings = load_settings(args.config)
    resolve_threads(args, settings)
    absolutize_paths(args)
    assemblies = [Assembly(a.label, os.path.abspath(a.path)) for a in assemblies]
    validate_inputs(
        "coverage",
        [(check_fasta, a.path, "assembly (-a)") for a in assemblies]
        + [(check_fastq, path, "read file") for path in args.reads],
    )
    env = env_for("coverage", settings)
    # fairy ships in this module's env. It used to have one of its own, which meant the route needed
    # a second install command to try; it is one static binary, so that bought nothing.
    fairy_env = env if args.coverage == "fairy" else None
    rec = start_run(
        "coverage",
        args,
        env,
        settings,
        inputs=[a.path for a in assemblies] + list(args.reads),
    )

    ckpt = make_checkpoint(args.output)
    try:
        announcement("BUILDING ONE CONTIG CATALOGUE FROM %d ASSEMBL(IES)" % len(assemblies))
        ensure_dir(args.output)
        work = os.path.join(args.output, WORK_DIR)
        ensure_dir(work)
        catalogue = os.path.join(args.output, CATALOGUE)
        if ckpt.todo("catalogue"):
            if dry_run():
                comm(
                    "(dry run) would concatenate %d assembl(ies) into %s"
                    % (len(assemblies), catalogue)
                )
            else:
                try:
                    counts = write_catalogue(assemblies, catalogue, separator=args.separator)
                except (OSError, ValueError) as exc:
                    error(str(exc))
                comm(
                    "%d contig(s) from %d assembl(ies), prefixed with their label and %r"
                    % (sum(counts.values()), len(counts), args.separator)
                )
            ckpt.done("catalogue")
        else:
            comm("skipping the catalogue (already built; --resume)")

        depth_table = os.path.join(args.output, CATALOGUE_DEPTH)
        if args.coverage == "fairy":
            announcement("SKETCHING EVERY SAMPLE AND ESTIMATING DEPTH FROM THE SKETCHES")
            comm("no alignment on this route: nothing below ~95%% ANI is resolved")
            sketch_dir = os.path.join(work, SKETCH_DIR)
            if ckpt.todo("fairy_sketch"):
                sketches = sketch_reads(
                    args.reads,
                    sketch_dir,
                    threads=args.threads,
                    env=fairy_env,
                    single_end=args.single_end,
                    interleaved=args.interleaved,
                )
                ckpt.done("fairy_sketch")
            else:
                comm("skipping the sketches (already built; --resume)")
                sketches = _sketches(sketch_dir)
            # The stage is named for the route rather than shared with the alignment one, so that a
            # run resumed with the other --coverage rebuilds the table instead of inheriting a table
            # the new route did not produce.
            if ckpt.todo("fairy_depth_table"):
                run(
                    FAIRY_COVERAGE.format(
                        threads=args.threads,
                        sketches=" ".join(sketches),
                        catalogue=catalogue,
                        depth=depth_table,
                    ),
                    env=fairy_env,
                    tool="fairy coverage",
                )
                ckpt.done("fairy_depth_table")
            else:
                comm("skipping the depth table (already built; --resume)")
        else:
            announcement("MAPPING EVERY SAMPLE AGAINST THE CATALOGUE")
            comm(
                "one alignment per sample against one catalogue, rather than one per sample per "
                "assembly"
            )
            if ckpt.todo("align_reads"):
                align_read_sets(
                    catalogue,
                    args.reads,
                    work,
                    aligner=args.aligner,
                    threads=args.threads,
                    env=env,
                    single_end=args.single_end,
                    interleaved=args.interleaved,
                )
                ckpt.done("align_reads")
            else:
                comm("skipping alignment (already done; --resume)")

            if ckpt.todo("depth_table"):
                announcement("MEASURING EVERY CONTIG'S DEPTH IN EVERY SAMPLE")
                # One pass over every BAM, producing the one table that is then sliced. No minimum
                # contig length is applied: filtering belongs to the binning run that reads the slice,
                # and a contig dropped here could not be recovered by a later run with a lower
                # threshold. See metawrap2.depth for the algorithm and for why it is not a shell-out.
                bams = sorted(_bams(work))
                _depth.write_depth_table(
                    bams,
                    depth_table,
                    samples=[os.path.splitext(os.path.basename(b))[0] for b in bams],
                    env=env,
                    threads=args.threads,
                    percent_identity=args.min_identity,
                )
                ckpt.done("depth_table")
            else:
                comm("skipping the depth table (already built; --resume)")

        if dry_run():
            comm(
                "(dry run) would split %s into one table per assembly under %s"
                % (depth_table, os.path.join(args.output, DEPTH_DIR))
            )
            announcement("COVERAGE PIPELINE FINISHED SUCCESSFULLY!!!")
            return 0

        expect_produced(
            depth_table,
            "the catalogue depth table",
            "fairy coverage" if args.coverage == "fairy" else "samtools depth",
        )
        announcement("SPLITTING THE TABLE BACK INTO ONE TABLE PER ASSEMBLY")
        comm(
            "%d sample column(s) in the table; every slice carries all of them"
            % _sample_count(depth_table)
        )
        depth_dir = os.path.join(args.output, DEPTH_DIR)
        try:
            paths = split_depth_table(
                depth_table,
                depth_dir,
                [a.label for a in assemblies],
                separator=args.separator,
            )
        except (OSError, ValueError) as exc:
            error(str(exc))
        _write_assemblies_table(assemblies, os.path.join(args.output, ASSEMBLIES_TABLE), paths)

        if args.coverage == "fairy":
            # The sketches are kept: a few MB per sample, and they are not specific to this
            # catalogue, so a cohort that gains a sample re-sketches one sample rather than all of
            # them. The BAMs are swept for the opposite reason - gigabytes, and specific to the
            # catalogue that produced them.
            comm("kept the read sketches in %s; they are reusable for another catalogue" % work)
        elif not args.keep_bams:
            for bam in _bams(work):
                fs_remove(bam)
            comm("removed the per-sample BAMs (pass --keep-bams to keep them)")
        if not args.keep_catalogue:
            # By suffix rather than by prefix: the index names are known (INDEX_SUFFIXES), and a
            # sweep of everything starting with the catalogue's name would eventually take something
            # a later version of this module put beside it.
            for path in [catalogue] + [catalogue + suffix for suffix in INDEX_SUFFIXES]:
                fs_remove(path)
            comm("removed the catalogue and its index (pass --keep-catalogue to keep them)")

        announcement("COVERAGE PIPELINE FINISHED SUCCESSFULLY!!!")
        comm("bin each sample with its own slice:")
        for assembly in assemblies:
            comm(
                "  metawrap2 binning -a %s --depth %s -o <out> --metabat2 --maxbin2"
                % (assembly.path, paths[assembly.label])
            )
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":  # pragma: no cover - module entry point
    import sys

    raise SystemExit(main(sys.argv[1:]))
