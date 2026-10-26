"""MetaWrap2 strains: strain sharing between samples and within-sample population diversity.

A MAG is a population consensus, and two samples whose MAGs are 99.9% identical can carry
different strains, so consensus ANI cannot support persistence, transmission or engraftment
claims. inStrain works from the alignments and per-position allele frequencies, so e.g. two
variants at 30/70 are visible. The module maps each sample's reads against one genome set,
profiles each alignment, and compares every pair of samples for every genome present in both.

Rules, both in :mod:`metawrap2.strains`:

* popANI is reported with its breadth. It is computed only over positions both samples covered
  deeply enough, so two samples sharing 2% of a genome at 100% identity give popANI 100%. Below
  50% of the genome compared, the verdict is `undetermined` instead of `same_strain`.
* Nucleotide diversity is left blank (with the breadth beside it) for a genome with low breadth,
  e.g. 10%.

Every sample is mapped to the same genome set; popANI between samples mapped to different
references is not meaningful.
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import Dict, List, Optional, Sequence, Tuple

from .. import aligners as _aligners
from .. import strains as _strains
from ..cliformat import make_parser
from ..config import MODULE_ENVS, load_settings
from ..constants import TSV_SUFFIX
from ..io.seqio import contig_id, iter_fasta
from ..progress import bar
from ._common import (
    ALIGN_SORT,
    SORT_THREADS,
    absolutize_paths,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    list_bins,
    make_checkpoint,
    optional_env_for,
    require_file,
    require_nonempty_dir,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["strains"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# Mapping uses the same aligner registry as `binning` and `coverage`, piped into sort through
# _common.ALIGN_SORT (pipefail, sort temp directory), so the SAM is never written (tens of GB for
# a 12 M-read sample against a 26-genome set).
INDEX_BAM = "samtools index -@ {threads} {bam}"
# --database_mode because the reference is a collection of genomes: it skips the unused mm-level
# profiling, sets the read-ANI floor to 0.92 and only profiles scaffolds with coverage (hours
# instead of days on a bin set).
#
# -s is the scaffold-to-bin file; without it inStrain reports per scaffold only. This module writes
# it from the bin folder, so it matches the mapped fasta.
INSTRAIN_PROFILE = (
    "inStrain profile {bam} {fasta} -o {outdir} -p {threads} -s {stb}"
    " --min_read_ani {min_read_ani} --skip_plot_generation"
)
INSTRAIN_COMPARE = (
    "inStrain compare -i {profiles} -o {outdir} -p {threads} -s {stb} --skip_plot_generation"
)
# TRACS, the second backend. It screens samples against a sourmash signature database of the
# reference genomes before aligning, so build-db runs first, once per run.
TRACS_BUILD_DB = "tracs build-db -i {genomes} -o {db} -t {threads}"
# `pipe` is align + distance + cluster in one call. -i is a file listing the read paths, which
# the module writes.
TRACS_PIPE = "tracs pipe -i {inputs} --database {db} --refseqs {refs} -o {outdir} -t {threads}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
#: One row per genome per sample: coverage, breadth, diversity.
PROFILE_TABLE = "populations" + TSV_SUFFIX
#: One row per genome per pair of samples: popANI, breadth, verdict.
COMPARE_TABLE = "strain_comparisons" + TSV_SUFFIX
#: TRACS' distances, kept separate from inStrain's table because they measure different
#: quantities on different scales.
TRACS_TABLE = "strain_distances_tracs" + TSV_SUFFIX
#: TRACS' output directory.
TRACS_DIR = "tracs"
#: Directory for the concatenated genome set used for mapping, and for inStrain's output.
WORK_DIR = "work_files"
PROFILE_DIR = "instrain"
#: Scaffold-to-bin mapping, written by this module so it matches the mapped fasta.
STB_FILE = "scaffold_to_bin" + TSV_SUFFIX
#: inStrain's default read-ANI floor for a single genome. `--database_mode` sets 0.92; the value
#: is passed explicitly so it is recorded in provenance.
DEFAULT_MIN_READ_ANI = 0.95
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 strains",
        usage="metawrap2 strains [options] -o out -b bin_folder -r s1=r1.fq,r2.fq -r s2=...",
        description="Within-population variation and strain sharing between samples, by inStrain.",
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    p.add_argument(
        "-b",
        "--bins",
        required=True,
        help="the genome set every sample is compared against: a bin folder or a dereplicated "
        "catalogue. All samples use the same set, because popANI is only defined between "
        "samples mapped to the same references",
    )
    p.add_argument(
        "-r",
        "--reads",
        nargs="+",
        metavar="SAMPLE=R1[,R2]",
        help="one argument per sample: `gut1=gut1_1.fq,gut1_2.fq`, or `gut1=gut1.fq` for "
        "single-end. QC'd reads",
    )
    p.add_argument(
        "--bam",
        nargs="+",
        metavar="SAMPLE=FILE",
        help="already-sorted BAMs of each sample against this genome set, instead of --reads. "
        "They must be mapped to the same genomes; the run cannot check this",
    )
    p.add_argument(
        "--aligner",
        choices=sorted(_aligners.ALIGNERS),
        default=_aligners.DEFAULT_ALIGNER,
        help="read aligner (default: %(default)s), from the same registry binning uses",
    )
    p.add_argument(
        "--min-read-ani",
        type=float,
        default=DEFAULT_MIN_READ_ANI,
        metavar="FRACTION",
        help="minimum identity of a read to its reference for inStrain to use it (default: "
        "%(default)s). Lower it for genomes more divergent than the strain level; raise it to "
        "narrow what counts as the same population",
    )
    p.add_argument(
        "--database-mode",
        action="store_true",
        help="inStrain's --database_mode, for a reference made of many genomes: skips profiling "
        "steps whose output this module does not use, and is much faster. It also lowers the "
        "read-ANI floor to 0.92, so it is off by default",
    )
    p.add_argument(
        "--backend",
        choices=("instrain", "tracs", "both"),
        default="instrain",
        help="which strain caller to use (default: %(default)s). TRACS was the only tool in Shaw & "
        "Yu's evaluation that did not overestimate SNP distances. It reports a transmission "
        "distance, while inStrain reports a population identity; `both` runs each and reports "
        "where they agree, without averaging them",
    )
    p.add_argument(
        "--skip-compare",
        action="store_true",
        help="profile each sample but do not compare them. The per-sample diversity is still "
        "written; strain sharing between samples is not",
    )
    threads_arg(p)
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def _parse_pairs(specs: Sequence[str], what: str) -> List[Tuple[str, List[str]]]:
    """``sample=path[,path]`` -> (sample, [paths]), validated before anything runs."""
    out: List[Tuple[str, List[str]]] = []
    seen = set()
    for spec in specs:
        sample, _, paths = spec.partition("=")
        sample = sample.strip()
        files = [path.strip() for path in paths.split(",") if path.strip()]
        if not sample or not files:
            error("%s takes sample=path, e.g. gut1=gut1_1.fq,gut1_2.fq (got %r)" % (what, spec))
        if sample in seen:
            error(
                "sample %r is given twice to %s; sample names must be unique, since each is a "
                "row of the output" % (sample, what)
            )
        seen.add(sample)
        files = [os.path.abspath(path) for path in files]
        for path in files:
            require_file(path, "%s for %s" % (what, sample))
        out.append((sample, files))
    return out


def _concatenate_genomes(bins_dir: str, fasta: str, stb: str) -> Tuple[int, int]:
    """Write the genome set as one fasta and its scaffold-to-bin table. Returns (genomes, contigs).

    Both are written in one pass over the same files, so the table matches the fasta (inStrain
    assigns every scaffold to a genome through it).

    Contig names are kept unchanged. Bins from this pipeline already have unique contig names (one
    assembly), and a duplicate is an error, since it would make inStrain's genome assignment
    arbitrary.
    """
    names = list_bins(bins_dir)
    if not names:
        error("no genomes in %s" % bins_dir)
    seen: Dict[str, str] = {}
    contigs = 0
    with open(fasta, "w") as out, open(stb, "w", newline="") as table:
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        for name in names:
            genome = os.path.splitext(name)[0]
            for header, sequence in iter_fasta(os.path.join(bins_dir, name)):
                contig = contig_id(header)
                if contig in seen:
                    error(
                        "contig %r is in both %s and %s. inStrain assigns every scaffold to one "
                        "genome through the scaffold-to-bin table, so contig names must be "
                        "unique across genomes" % (contig, seen[contig], genome)
                    )
                seen[contig] = genome
                out.write(">%s\n%s\n" % (contig, sequence))
                writer.writerow([contig, genome])
                contigs += 1
    return len(names), contigs


def _map_sample(args, env, fasta: str, sample: str, reads: List[str], work: str) -> Optional[str]:
    """Map one sample against the genome set, returning its sorted BAM."""
    bam = os.path.join(work, "%s.bam" % sample)
    if os.path.isfile(bam) and os.path.getsize(bam):
        comm("%s: alignment already present" % sample)
        return bam
    spec = _aligners.ALIGNERS[args.aligner]
    index = _aligners.index_path(fasta, args.aligner)
    if len(reads) == 2:
        align = spec["paired"].format(threads=args.threads, index=index, r1=reads[0], r2=reads[1])
    else:
        align = spec["single"].format(threads=args.threads, index=index, reads=reads[0])
    rc = run(
        [
            "bash",
            "-c",
            ALIGN_SORT.format(
                align=align,
                tmp=os.path.join(work, "tmp-samtools"),
                threads=min(args.threads, SORT_THREADS),
                bam=bam,
            ),
        ],
        env=env,
        tool="%s | samtools sort" % args.aligner,
        check=False,
    )
    if rc != 0 or not os.path.isfile(bam):
        warning("mapping failed for %s; that sample is left out" % sample)
        return None
    run(
        INDEX_BAM.format(threads=args.threads, bam=bam),
        env=env,
        tool="samtools index",
        check=False,
    )
    return bam


def _run_tracs(args, settings, samples, work: str) -> List[_strains.TracsDistance]:
    """Screen, align and compare with TRACS, returning its distances.

    Writes what TRACS needs: a sourmash signature database of the reference genomes (used to
    screen which genomes a sample carries), a file listing the read paths, and an output
    directory. All are per run, since the database must match this genome set.
    """
    env = optional_env_for("tracs", settings, needed_by="strains --backend tracs")
    root = os.path.join(args.output, TRACS_DIR)
    ensure_dir(root)
    genomes = [os.path.join(args.bins, name) for name in list_bins(args.bins)]
    db = os.path.join(work, "tracs_db")
    rc = run(
        TRACS_BUILD_DB.format(genomes=" ".join(genomes), db=db, threads=args.threads),
        env=env,
        tool="tracs build-db",
        check=False,
    )
    database = db + ".zip"
    if rc != 0 or not os.path.isfile(database):
        warning("`tracs build-db` failed, so TRACS contributes no distances")
        return []

    listing = os.path.join(work, "tracs_inputs" + TSV_SUFFIX)
    with open(listing, "w", newline="") as fh:
        inputs_writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        # Header row: `tracs pipe` skips the first line with `next(infile)` (undocumented), so
        # without it the first sample is dropped (an empty distance table from two samples).
        inputs_writer.writerow(["sample", "read_1", "read_2"])
        for sample, reads in samples:
            inputs_writer.writerow([sample, *reads])
    rc = run(
        TRACS_PIPE.format(
            inputs=listing,
            db=database,
            refs=args.bins,
            outdir=os.path.join(root, "pipe"),
            threads=args.threads,
        ),
        env=env,
        tool="tracs pipe",
        check=False,
    )
    table = os.path.join(root, "pipe", "transmission_distances.csv")
    if rc != 0 or not os.path.isfile(table):
        warning("`tracs pipe` produced no distance table")
        return []
    distances = _strains.parse_tracs(table, samples=[sample for sample, _ in samples])
    if not distances:
        comm(
            "TRACS compared nothing: it skips any reference with less than %.0f%% of its genome "
            "covered, so the samples lack the depth for a comparison"
            % (_strains.TRACS_MIN_BREADTH * 100)
        )
        return []

    path = os.path.join(args.output, TRACS_TABLE)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=_strains.TRACS_TABLE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for distance in sorted(distances, key=lambda d: (d.genome, d.sample_a, d.sample_b)):
            writer.writerow(_strains.tracs_row(distance))
    same = [d for d in distances if d.verdict is _strains.StrainVerdict.SAME]
    comm(
        "%d TRACS distance(s) -> %s; %d within %d SNPs"
        % (len(distances), path, len(same), _strains.TRACS_SNP_THRESHOLD)
    )
    for distance in same[:10]:
        comm(
            "  %s: %s vs %s, %.0f SNPs over %s sites (transmission distance %.3f)"
            % (
                distance.genome,
                distance.sample_a,
                distance.sample_b,
                (
                    distance.filtered_snp_distance
                    if distance.filtered_snp_distance is not None
                    else (distance.snp_distance or 0.0)
                ),
                distance.sites_considered,
                distance.transmission_distance or 0.0,
            )
        )
    comm(_strains.TRACS_CAVEAT)
    return distances


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    resolve_threads(args, settings)
    # Parse sample specs before absolutize_paths, which treats `ERR011347=reads_1.fq` as a path
    # and would make the sample name `/abs/path/ERR011347`.
    specs = list(args.reads or args.bam or ())
    absolutize_paths(args)
    if not args.reads and not args.bam:
        error("nothing to profile: pass --reads sample=... or --bam sample=...")
    if args.reads and args.bam:
        error("pass --reads or --bam, not both (--bam skips the read mapping)")
    validate_inputs("strains", [(require_nonempty_dir, args.bins, "genome set (-b)")])
    samples = _parse_pairs(specs, "--reads" if args.reads else "--bam")
    if len(samples) < 2 and not args.skip_compare:
        comm(
            "only one sample: its populations are profiled and the comparison between samples "
            "is skipped"
        )
        args.skip_compare = True

    env = env_for("strains", settings)
    rec = start_run("strains", args, env, settings, inputs=[args.bins])
    ckpt = make_checkpoint(args.output)
    try:
        announcement("BEGIN STRAINS PIPELINE!")
        ensure_dir(args.output)
        work = os.path.join(args.output, WORK_DIR)
        ensure_dir(work)
        fasta = os.path.join(work, "genomes.fasta")
        stb = os.path.join(args.output, STB_FILE)

        if dry_run():
            comm(
                "would profile %d sample(s) against the genomes in %s and compare every pair"
                % (len(samples), args.bins)
            )
            announcement("STRAINS PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        tracs_distances: List[_strains.TracsDistance] = []
        if args.backend in ("tracs", "both"):
            tracs_distances = _run_tracs(args, settings, samples, work)
            if args.backend == "tracs":
                announcement("STRAINS PIPELINE SUCCESSFULLY FINISHED!!!")
                return 0

        genomes, contigs = _concatenate_genomes(args.bins, fasta, stb)
        comm("%d genome(s), %d contig(s) -> %s" % (genomes, contigs, os.path.basename(fasta)))

        bams: List[Tuple[str, str]] = []
        if args.bam:
            bams = [(sample, files[0]) for sample, files in samples]
        else:
            if ckpt.todo("index"):
                # minimap2's template uses {index}; bwa's only {assembly}, since its index files
                # sit beside it. Both are supplied so no branch is needed.
                rc = run(
                    _aligners.ALIGNERS[args.aligner]["index_cmd"].format(
                        assembly=fasta, index=_aligners.index_path(fasta, args.aligner)
                    ),
                    env=env,
                    tool="%s index" % args.aligner,
                    check=False,
                )
                if rc != 0:
                    error("could not index the genome set with %s" % args.aligner)
                ckpt.done("index")
            for sample, reads in bar(samples, desc="mapping", unit=" sample", total=len(samples)):
                bam = _map_sample(args, env, fasta, sample, reads, work)
                if bam:
                    bams.append((sample, bam))
        if not bams:
            error("no sample could be aligned, so there is nothing to profile")

        profiles: List[_strains.GenomeProfile] = []
        profile_dirs: List[Tuple[str, str]] = []
        for sample, bam in bar(bams, desc="profiling", unit=" sample", total=len(bams)):
            outdir = os.path.join(args.output, PROFILE_DIR, sample)
            rc = run(
                INSTRAIN_PROFILE.format(
                    bam=bam,
                    fasta=fasta,
                    outdir=outdir,
                    threads=args.threads,
                    stb=stb,
                    min_read_ani=args.min_read_ani,
                )
                + (" --database_mode" if args.database_mode else ""),
                env=env,
                tool="inStrain profile",
                check=False,
            )
            table = os.path.join(outdir, "output", "%s_genome_info.tsv" % os.path.basename(outdir))
            if rc != 0 or not os.path.isfile(table):
                warning("inStrain profile failed on %s; that sample has no populations" % sample)
                continue
            profiles.extend(_strains.parse_profile(table, sample))
            profile_dirs.append((sample, outdir))

        if not profiles:
            error("no sample produced a population profile")

        path = os.path.join(args.output, PROFILE_TABLE)
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=_strains.PROFILE_TABLE_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for profile in sorted(profiles, key=lambda p: (p.sample, p.genome)):
                writer.writerow(_strains.profile_row(profile))
        reportable = [p for p in profiles if p.diversity_is_reportable]
        comm(
            "%d genome-sample pair(s) -> %s; %d covered well enough for a diversity figure"
            % (len(profiles), path, len(reportable))
        )
        comm(_strains.DIVERSITY_CAVEAT)

        if args.skip_compare or len(profile_dirs) < 2:
            announcement("STRAINS PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        comparable = _strains.comparable_genomes(profiles, [sample for sample, _ in profile_dirs])
        for (first, second), count in sorted(comparable.items()):
            comm(
                "%s vs %s: %d of %d genome(s) covered deeply enough in both to be comparable"
                % (first, second, count, genomes)
            )
        compare_out = os.path.join(args.output, PROFILE_DIR, "compare")
        rc = run(
            INSTRAIN_COMPARE.format(
                profiles=" ".join(outdir for _sample, outdir in profile_dirs),
                outdir=compare_out,
                threads=args.threads,
                stb=stb,
            ),
            env=env,
            tool="inStrain compare",
            check=False,
        )
        table = os.path.join(compare_out, "output", "compare_genomeWide_compare.tsv")
        if rc != 0 or not os.path.isfile(table):
            warning(
                "inStrain compare failed, so no strain-sharing table was written. The per-sample "
                "populations above are unaffected"
            )
            announcement("STRAINS PIPELINE SUCCESSFULLY FINISHED!!!")
            return 0

        # inStrain labels rows with the profile directory name; relabel with the pipeline's sample
        # names so the table joins with the rest of the output.
        labels = {os.path.basename(outdir): sample for sample, outdir in profile_dirs}
        comparisons = _strains.relabel(_strains.parse_comparison(table), labels)
        path = os.path.join(args.output, COMPARE_TABLE)
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=_strains.COMPARE_TABLE_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for comparison in sorted(comparisons, key=lambda c: (c.genome, c.sample_a, c.sample_b)):
                writer.writerow(_strains.comparison_row(comparison))
        counts = _strains.summarise(comparisons)
        comm(
            "%d comparison(s) -> %s: %d same strain, %d different, %d undetermined, %d not in both"
            % (
                len(comparisons),
                path,
                counts[_strains.StrainVerdict.SAME.value],
                counts[_strains.StrainVerdict.DIFFERENT.value],
                counts[_strains.StrainVerdict.UNDETERMINED.value],
                counts[_strains.StrainVerdict.ABSENT.value],
            )
        )
        for shared in _strains.shared_strains(comparisons)[:10]:
            comm(
                "  %s shared between %s and %s (popANI %.5f%% over %.1f%% of the genome)"
                % (
                    shared.genome,
                    shared.sample_a,
                    shared.sample_b,
                    shared.pop_ani or 0.0,
                    shared.percent_compared or 0.0,
                )
            )
        comm(_strains.BREADTH_CAVEAT)
        if tracs_distances:
            agreement = _strains.compare_backends(comparisons, tracs_distances)
            comm(
                "backends: %s of %s comparable genome-pairs agree on same versus different"
                % (agreement["agree"], agreement["comparable"])
            )
            for line in agreement["disagree"][:10]:  # type: ignore[index]
                comm("  %s" % line)
            comm(
                "  the tools can disagree without either being wrong: inStrain measures a "
                "population identity and TRACS a consensus SNP distance, so both are reported "
                "side by side and never averaged"
            )
        announcement("STRAINS PIPELINE SUCCESSFULLY FINISHED!!!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
