"""MetaWrap2 coverage module: one cross-sample depth table for a cohort, in N alignments.

Multi-sample binning needs every contig's depth in every sample. Mapping all N samples against
each of the N per-sample assemblies costs N^2 alignments (400 for a 20-sample cohort).

This module does it in N: concatenate the per-sample assemblies into one contig catalogue, map
each sample's reads once against it, build one depth table, then split the table into one slice
per source assembly. Each sample is then binned on its own contigs with every sample's coverage,
via ``metawrap2 binning --depth``.

Contig names collide across assemblies, so each contig is prefixed with its assembly's label
before concatenation and the split demultiplexes on that prefix. A catalogue is refused if any
original name contains the separator; see SEPARATOR. The slices carry the original contig names,
because they are used with the original assemblies.

``--coverage align`` (default) maps the reads; ``--coverage fairy`` estimates depth from k-mer
sketches without alignment. fairy writes the same contig x sample layout metaBAT2 reads, so the
split and everything downstream are the same; see COVERAGE_ROUTES for the accuracy cost.

Commands are in the COMMANDS block below and tunable values in the CONSTANTS block.
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple

from .. import depth as _depth
from ..aligners import ALIGNERS, DEFAULT_ALIGNER, INDEX_SUFFIXES
from ..cliformat import make_parser, percent
from ..command import capture_in_env
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
    fs_move,
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
from .binning import DEPTH_HEADER, SAMTOOLS_INDEX, depth_means

#: Shared with `binning`, which installs the same aligner and samtools. `install-env coverage`
#: builds or finds metawrap2-bins (see _env_yaml).
CONDA_ENV = MODULE_ENVS["coverage"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
#: The k-mer route (--coverage fairy): sketch each sample's reads once, then estimate every contig's
#: depth from the sketches. ``-S`` names the sketches after the samples, so the table's columns
#: carry sample names where the alignment route has BAM paths. fairy's defaults are kept: its
#: documentation warns against lowering ``-m 95``, and ``-c``/``-k`` are the values its published
#: accuracy was measured at.
FAIRY_SKETCH_PAIRED = "fairy sketch -1 {r1} -2 {r2} -S {names} -d {sketches} -t {threads}"
FAIRY_SKETCH_SINGLE = "fairy sketch -r {reads} -S {names} -d {sketches} -t {threads}"
FAIRY_COVERAGE = "fairy coverage -t {threads} {sketches} {catalogue} -o {depth}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: How coverage is measured. ``align`` maps the reads, which is how every binner's published
#: accuracy was measured. ``fairy`` estimates depth from k-mer sketches without alignment (minutes
#: instead of hours for a cohort) but resolves nothing below ~95% ANI, so a strain that differs
#: from its contig is credited to that contig or to nothing.
#:
#: ``align`` is the default on measured results. With metaBAT2 on one CAMI2 GI sample,
#: single-sample fairy scored AMBER F1 0.212 against minimap2's 0.452, and raising its ANI
#: threshold made agreement worse. In benchmark B5 (three CAMI2 samples, everything else
#: identical), fairy built the cohort table for 0.04 CPU-hours against minimap2's 2.26 (57x
#: cheaper), but the catalogue held 27 high-quality MAGs against 37, with 76% of bins passing GUNC
#: against 87%.
#:
#: The route remains opt-in. The variance ablation showed metaBAT2's variance column carries no
#: signal at this sample count, so fairy's loss comes from its mean depths (Pearson 0.237 against
#: bwa).
COVERAGE_ROUTES = ("align", "fairy")
DEFAULT_ROUTE = "align"
#: Where the fairy route keeps its per-sample sketches, inside the work directory.
SKETCH_DIR = "sketches"
#: Separates an assembly's label from the original contig name in the catalogue. ``~`` is legal in
#: a fasta header, means nothing to aligners or samtools, and does not appear in SPAdes, MEGAHIT or
#: Flye contig names. It is configurable for cohorts whose names contain it; names are checked
#: either way, because a separator inside a name would make the split cut in the wrong place.
SEPARATOR = "~"
#: What the catalogue and the tables are called inside the output directory.
CATALOGUE = "catalogue.fa"
CATALOGUE_DEPTH = "catalogue_depth.txt"
DEPTH_DIR = "depth"
#: Where --split-bams writes one directory per assembly label.
BAMS_DIR = "bams"
WORK_DIR = "work_files"
#: Written beside the slices; maps each assembly label to its assembly and depth table.
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
        "sketches with fairy (much faster, but does not resolve anything below ~95%% ANI)",
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
        "--split-bams",
        action="store_true",
        help="also write one BAM per assembly per sample under <out>/" + BAMS_DIR + "/<label>/, "
        "with contigs named as in that assembly. SemiBin2, VAMB, COMEBin and CONCOCT read BAMs "
        "and cannot use a depth table; pass these to them with `binning --bams`. Uses about as "
        "much disk as the catalogue BAMs; extraction reads each record once, through the index",
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
        help="minimum alignment identity for a read to count toward a contig's depth (default "
        "90). Against jgi's 97, 90 gave 88 bins instead of 74, purer bins, and 47 clean "
        "high-quality MAGs instead of 35; on strain-rich data a high cutoff removes signal "
        "without separating strains. Pass 97 to reproduce a jgi depth table, or 0 for no filter",
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

    ``label=path`` is explicit. A bare path is labelled by its filename without extension, so a
    cohort laid out as ``sample_0/final_assembly.fasta`` gives every sample the same label; the
    duplicate check in :func:`parse_assemblies` reports that. Labels end up in contig names
    downstream, so they are not renamed automatically.
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

    Duplicate labels, or a label containing the separator, would make the split ambiguous, so
    both raise ValueError.
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

    Splits on the first separator. Names containing the separator are not allowed into the
    catalogue, so everything after it is the contig name. A name with no separator returns
    ``("", name)``.
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

    Returns each label's contig count. Contigs are prefixed even for a single assembly, so the
    split has one code path.

    A contig name containing the separator raises ValueError, because the split demultiplexes on
    the separator.

    Records with no sequence, or with no A, C, G, T or U, are skipped. ``binning.copy_assembly``
    drops the same records, so the depth table rows match the binner's copy of the assembly.
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

    Each slice keeps the header (all sample columns) and restores the original contig names.

    Rows are dispatched by prefix in a single pass, since the catalogue table has a row for every
    contig of every sample. A label with no rows raises ValueError: it means the catalogue did not
    contain that assembly.
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
                # Held open for the whole pass, so the catalogue table is read once.
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
        # Not fatal: a reused catalogue can hold more than this run's assemblies. The usual cause
        # is a label typo, which shows up as a short slice.
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

    One fairy invocation covers all samples: fairy sketches them concurrently, and ``-S`` lists the
    names in read order.

    Interleaved files are sketched as single-end. The sketch is a subsample of read k-mers and does
    not record which mate a k-mer came from.
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
        # fairy names single-end and paired sketches differently, and the coverage call takes the
        # sketch paths, so each case needs its own suffix.
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


def plain_header(header: str, label: str, separator: str) -> Tuple[str, int]:
    """*header* with ``label<sep>`` stripped from every ``@SQ SN:``. Returns it and the count.

    Only SN fields with that prefix are changed. Read names and @PG command lines can contain the
    same text and are left alone.
    """
    prefix = label + separator
    out: List[str] = []
    renamed = 0
    for line in header.splitlines():
        if line.startswith("@SQ\t"):
            fields = line.split("\t")
            for i, field in enumerate(fields):
                if field.startswith("SN:" + prefix):
                    fields[i] = "SN:" + field[len("SN:") + len(prefix) :]
                    renamed += 1
            line = "\t".join(fields)
        out.append(line)
    return "\n".join(out) + "\n", renamed


def _index_kept_bams(bams: Sequence[str], threads: int, env: Optional[str]) -> None:
    """Index the BAMs that are being kept.

    The BAMs are already sorted, and BAM readers (the BAM binners, region queries with
    `samtools view`) need the index. Indexing is cheap compared with alignment.
    """
    for bam in bams:
        if os.path.isfile(bam + ".bai"):
            continue
        run(SAMTOOLS_INDEX.format(threads=threads, bam=bam), env=env, tool="samtools index")


def bam_references(bam: str, env: Optional[str]) -> List[Tuple[str, int]]:
    """``(name, length)`` for every reference in *bam*'s header, in header order."""
    out: List[Tuple[str, int]] = []
    for line in capture_in_env(env, ["samtools", "view", "-H", bam]).splitlines():
        if not line.startswith("@SQ\t"):
            continue
        name = length = None
        for field in line.split("\t"):
            if field.startswith("SN:"):
                name = field[3:]
            elif field.startswith("LN:"):
                length = int(field[3:]) if field[3:].isdigit() else None
        if name and length:
            out.append((name, length))
    return out


def write_region_bed(references: Sequence[Tuple[str, int]], prefix: str, path: str) -> int:
    """Write a BED of every reference starting with *prefix*. Returns how many it held.

    A BED file is used because a label can have hundreds of thousands of contigs, more than fits
    in ARG_MAX on a command line. `samtools view -M -L` uses the index with it, so each record is
    read once across all labels.
    """
    kept = 0
    with open(path, "w") as handle:
        for name, length in references:
            if name.startswith(prefix):
                handle.write("%s\t0\t%d\n" % (name, length))
                kept += 1
    return kept


def split_bams(
    assemblies: Sequence[Assembly],
    work: str,
    out_dir: str,
    separator: str,
    env: Optional[str],
    threads: int,
) -> Dict[str, List[str]]:
    """Write one BAM per assembly per sample, named as that assembly names its contigs.

    This lets the BAM-reading binners (SemiBin2, VAMB, COMEBin, CONCOCT), which do not accept
    `--depth`, use the cohort's single alignment pass. `binning --bams <dir>` reads the output.

    * Extraction uses a region file and the index, so the whole split costs about one pass over the
      catalogue BAMs.
    * Only this label's `@SQ` names are de-prefixed. Other labels' `@SQ` lines must stay: reference
      IDs are positional, so dropping lines would renumber every record, and de-prefixing them
      would merge identically named contigs (`k141_1`) from different assemblies. They carry no
      records in the extracted file, so binners ignore them.
    """
    ensure_dir(out_dir)
    bams = _bams(work)
    written: Dict[str, List[str]] = {}
    if not bams:
        return written
    # `samtools view -M -L` needs the index to read only one label's records; without it, samtools
    # refuses the region file.
    _index_kept_bams(bams, threads, env)
    references = bam_references(bams[0], env)
    for assembly in assemblies:
        label_dir = os.path.join(out_dir, assembly.label)
        ensure_dir(label_dir)
        bed = os.path.join(label_dir, "regions.bed")
        count = write_region_bed(references, assembly.label + separator, bed)
        if not count:
            warning(
                "no contig in the catalogue BAMs is labelled %r; %s gets no per-assembly BAM"
                % (assembly.label, assembly.label)
            )
            continue
        comm("extracting %d contig(s) of %s from %d BAM(s)" % (count, assembly.label, len(bams)))
        produced = []
        for bam in bams:
            sample = os.path.splitext(os.path.basename(bam))[0]
            target = os.path.join(label_dir, sample + ".bam")
            run(
                "samtools view -b -M -L %s -@ %d -o %s %s"
                % (bed, max(1, threads - 1), target, bam),
                env=env,
                tool="samtools view",
            )
            header = capture_in_env(env, ["samtools", "view", "-H", target])
            plain, renamed = plain_header(header, assembly.label, separator)
            if renamed:
                header_path = target + ".plain.sam"
                with open(header_path, "w") as handle:
                    handle.write(plain)
                reheadered = target + ".plain.bam"
                run(
                    "samtools reheader %s %s" % (header_path, target),
                    env=env,
                    tool="samtools reheader",
                    stdout_path=reheadered,
                )
                fs_move(reheadered, target)
                fs_remove(header_path)
            produced.append(target)
        _index_kept_bams(produced, threads, env)
        written[assembly.label] = produced
        comm("wrote %s" % label_dir)
    return written


def _deprefix_kept_bams(
    assemblies: Sequence[Assembly],
    work: str,
    separator: str,
    env: Optional[str],
    threads: int = 1,
) -> None:
    """Make kept BAMs name contigs the way the assembly does, or say why they cannot.

    The alignment is against the catalogue, whose contigs are ``label<sep>contig``, while the split
    depth tables and the assemblies use plain names. BAM-reading binners (SemiBin2, COMEBin, VAMB,
    CONCOCT) then match zero contigs and produce no bins without an error; metaBAT2 reads the table
    and is unaffected.

    With one assembly the headers are rewritten. `samtools reheader` copies each BAM, so this runs
    only with `--keep-bams`.

    With several assemblies the headers are left in catalogue space, because contig names repeat
    between assemblies and stripping the prefix would merge different sequences. A warning gives
    the command to extract one assembly's BAM.
    """
    bams = _bams(work)
    if not bams:
        return
    if len(assemblies) != 1:
        warning(
            "the kept BAMs are aligned against the catalogue, so their contigs are named "
            f"label{separator}contig, while the depth tables and assemblies use plain names. "
            "metaBAT2 and MaxBin2 read the table and are unaffected; SemiBin2, VAMB, COMEBin and "
            "CONCOCT read the BAM and will match no contigs in a plain assembly.\n"
            f"The headers are not rewritten because {len(assemblies)} assemblies share this "
            "catalogue and contig names repeat between them, so stripping the prefix would merge "
            "different sequences. To get one assembly's BAM, select its contigs first:\n"
            "  samtools view -b -o <sample>.<label>.bam <sample>.bam $(grep '^>' <assembly> | "
            f"sed 's/^>/<label>{separator}/;s/ .*//')\n"
            "then reheader that file to drop the prefix."
        )
        _index_kept_bams(bams, threads, env)
        return
    label = assemblies[0].label
    comm(
        "rewriting the kept BAM headers to drop the '%s%s' catalogue prefix, so BAM-reading "
        "binners see the contig names of the assembly and depth table" % (label, separator)
    )
    for bam in bams:
        header = capture_in_env(env, ["samtools", "view", "-H", bam])
        rewritten, renamed = plain_header(header, label, separator)
        if not renamed:
            continue
        header_path = bam + ".plain.sam"
        with open(header_path, "w") as handle:
            handle.write(rewritten)
        plain = bam + ".plain.bam"
        run(
            "samtools reheader %s %s" % (header_path, bam),
            env=env,
            tool="samtools reheader",
            stdout_path=plain,
        )
        # In-place reheader is CRAM-only, so write a new file and move it over the original once it
        # exists; an interrupted rewrite leaves the BAM intact.
        fs_move(plain, bam)
        fs_remove(header_path)
    comm("%d BAM(s) now name contigs as the assembly does" % len(bams))
    # Index after the rewrite, since reheadering writes a new file.
    _index_kept_bams(bams, threads, env)


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
    # fairy is installed in this module's env.
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
            # Checkpoint names differ by route, so resuming with the other --coverage rebuilds the
            # table.
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
                # One pass over every BAM into one table. No minimum contig length is applied here;
                # the binning run that reads the slice filters, so a later run can use a lower
                # threshold. See metawrap2.depth for the algorithm.
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
            # Sketches are kept: a few MB per sample and not specific to this catalogue, so adding
            # a sample means sketching only that sample. BAMs are large and catalogue-specific, so
            # they are removed.
            comm("kept the read sketches in %s; they are reusable for another catalogue" % work)
        else:
            # Split before the catalogue BAMs are removed, since the per-assembly BAMs are extracted
            # from them.
            if args.split_bams:
                announcement("SPLITTING THE ALIGNMENTS INTO ONE BAM PER ASSEMBLY")
                split = split_bams(
                    assemblies,
                    work,
                    os.path.join(args.output, BAMS_DIR),
                    args.separator,
                    env,
                    args.threads,
                )
                if split:
                    comm(
                        "bin each assembly's BAMs with `binning -a <assembly> --bams %s/<label>`"
                        % os.path.join(args.output, BAMS_DIR)
                    )
            if not args.keep_bams:
                for bam in _bams(work):
                    fs_remove(bam)
                    fs_remove(bam + ".bai")
                comm("removed the catalogue BAMs (pass --keep-bams to keep them)")
            else:
                _deprefix_kept_bams(assemblies, work, args.separator, env, args.threads)
        if not args.keep_catalogue:
            # Remove by known index suffix (INDEX_SUFFIXES), not by name prefix, so unrelated files
            # beside the catalogue are left alone.
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
