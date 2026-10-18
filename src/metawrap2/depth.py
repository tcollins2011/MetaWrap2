"""Per-contig, per-sample coverage depth and intra-contig depth variance, computed here.

This replaces ``jgi_summarize_bam_contig_depths``, which MetaWrap2 used to shell out to. The reason
is not that the tool is wrong - it is that the depth table is the single most load-bearing number in
the whole pipeline, and it was the one thing the pipeline could not see inside. When MaxBin2 was
found to be binning on invented depths, the defect was in how that table was being *converted*, and
diagnosing it meant reasoning about a binary with no test coverage here and no way to check an
intermediate. The table every binner depends on should be code this project owns and tests.

The algorithm is the same one, deliberately, so numbers stay comparable with everything measured
before the switch:

* **Reads are filtered** by mapping quality (``min_map_qual``) and by alignment identity
  (``percent_identity``), and unmapped / secondary / QC-failed / duplicate records are excluded. The
  identity filter is the reason a depth table is not just ``samtools depth``. Its floor is **90%**,
  which is the one deliberate departure from jgi's defaults in this module - see
  :data:`PERCENT_IDENTITY` for the measurement behind that, and for why jgi's 97 costs 26% of the
  clean bins on strain-rich data.
* **Contig ends are trimmed** before averaging - ``edge_bases`` (75) at each end. Coverage falls off
  at the end of a contig because reads cannot hang over the edge, so including the ends biases every
  contig's mean downward, short contigs much more than long ones. A contig too short to trim
  (``<= 2 * edge_bases``) is measured whole rather than dropped.
* **Variance is the variance of per-base depth** inside that trimmed window - not variance across
  samples. metaBAT2 reads it as a second signal per contig, and a table whose variance column is
  constant is one metaBAT2 rejects outright.

  Worth knowing before building a cheaper contig-abundance estimator: **that column carries no usable
  signal once there are several samples.** Measured by handing metaBAT2 four depth tables with identical
  means and identical contigs, differing only in the variance column - real, zeroed, shuffled between
  contigs, and resampled at random - and scoring the bins with AMBER. At three samples the spread across
  all four is 1.5% of F1, and the nominally best arm is the one where every contig wears another contig's
  variance. At one sample it is worth about 10% (F1 0.457 real against 0.416 zeroed).

  So the variance is a weak substitute for cross-sample coverage rather than an independent signal: given
  several samples, the covariance between them says what the within-contig variance was standing in for.
  The consequence for a future estimator is that it does **not** have to produce a real variance to be
  usable on a cohort - emitting a plausible constant-free column costs nothing measurable, and what such
  an estimator has to get right is the means. That is where fairy loses (see modules/coverage.py), not
  here. One dataset, metaBAT2 only; MaxBin2 and CONCOCT read no variance column at all.

Implementation note: the per-base depths come from ``samtools depth``, which does the pileup in C and
emits one row per position with one column per BAM - the same shape as the output. Everything after
that is numpy over chunks, so the cost is dominated by the pileup rather than by Python. Nothing is
held per contig except two running sums per sample, so a 200k-contig assembly costs no more memory
than a small one.
"""

from __future__ import annotations

import os
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .command import runner
from .modules._common import comm

#: Bases trimmed from each end of a contig before averaging. 75 is jgi's ``--maxEdgeBases`` default
#: and is about a read length: the window where coverage is depressed purely because an alignment
#: cannot start before the contig does.
EDGE_BASES = 75
#: Minimum alignment identity, in percent, for a read to count toward a contig's depth. Below this a
#: read is treated as coming from something else.
#:
#: **90, not jgi's 97.** This is the one place MetaWrap2 departs from
#: ``jgi_summarize_bam_contig_depths``' defaults, and it is measured rather than preferred. 97 was
#: inherited, never evaluated, and turned out to be a 2x lever: on CAMI2 GI it discards half the mapped
#: bases *against the true reference genomes*, since the gold assembly is those genomes and a read's
#: divergence from it is sequencing error plus real strain variation.
#:
#: One set of minimap2 alignments, six depth tables from it, metaBAT2 on each, AMBER and GUNC on the
#: results:
#:
#:     floor   mean depth   bins   AMBER F1   purity   HQ   clean HQ
#:     none        15.71      90     0.6607   0.8374     48      45
#:     90          15.70      88     0.6624   0.8375     48      47
#:     95          13.46      86     0.6574   0.8264     43      40
#:     97           8.35      74     0.6371   0.8175     37      35
#:     99           1.48      21     0.2452   0.7904      -       -
#:     100            -        7     0.0525   0.9998      -       -
#:
#: "clean HQ" is bins that are both >=70/<=10 and non-chimeric by GUNC, counted per bin. 97 loses 14 bins
#: of 88 and gains nothing: at 90 the bins are more complete, *more pure*, and less often chimeric, for
#: **47 clean high-quality MAGs against 35** - a 34% gain. That last column is what made this safe to
#: change. The objection to a looser floor is that a relative's reads inflate a contig's depth,
#: two strains come to look like one organism and the binner merges them - so it should buy completeness
#: at purity's expense. Both AMBER purity and GUNC say it does not. On this data the filter was not
#: separating strains, it was removing signal the binner was using.
#:
#: 90 rather than 0 because a floor still has a job - a read 20% divergent from a contig is evidence
#: about some other organism - and 0 and 90 are within 0.002 AMBER F1 of each other, so the cheap
#: insurance costs nothing. Measured on one dataset with one binner, which is why this is a module
#: constant and an argument rather than a hardcoded number.
PERCENT_IDENTITY = 90.0
#: jgi's own default, kept named because it is what metaBAT2's published behaviour was validated
#: against, and anyone reproducing a metaWRAP 1 result or a jgi depth table needs it.
JGI_PERCENT_IDENTITY = 97.0
#: Minimum MAPQ. 0, as jgi's ``--minMapQual`` default: multi-mapping reads are not thrown away here,
#: because a repeat shared between two genomes is real coverage of both and the binner is what has to
#: decide. Raise it if you would rather lose repeats than share them.
MIN_MAP_QUAL = 0
#: SAM flags excluded from the pileup: unmapped, secondary alignment, QC fail, PCR/optical duplicate.
#: Secondary alignments are excluded and multi-mapping is not: the *primary* record of a
#: multi-mapping read still counts, so a repeat contributes once per read rather than once per place
#: it could have gone.
EXCLUDED_FLAGS = 0x4 | 0x100 | 0x200 | 0x400
#: Rows parsed per chunk. Large enough that the parse is amortised, small enough that the frame is a
#: few hundred MB at ten samples.
CHUNK_ROWS = 4_000_000

#: How many samtools threads one BAM's pileup gets, and so how many BAMs are measured at once.
#:
#: `samtools depth` does not scale: a user measured ~2 minutes per BAM using every thread of a large
#: machine, with the machine at ~25% utilisation for the hours it took to walk 50 BAMs one at a time.
#: The work is embarrassingly parallel across BAMs - each one fills its own column of the running
#: sums and nothing is shared - so the threads are spent on several BAMs at once instead. Four per
#: pileup rather than one: `samtools view` does the identity filtering and decompression, which
#: does thread, and the pandas C parser and numpy accumulation both release the GIL.
DEPTH_THREADS_PER_BAM = 4

#: A ceiling on how many pileups run at once, whatever the thread budget allows. Each worker holds
#: a chunk of up to CHUNK_ROWS rows while it folds it in, so this bounds memory as well as the
#: number of samtools processes.
MAX_DEPTH_WORKERS = 16

#: The first three columns, then two per sample. Written out here because three separate places used
#: to spell this layout and one of them had it wrong.
LEADING_COLUMNS = ("contigName", "contigLen", "totalAvgDepth")


def _samtools(env: Optional[str]) -> List[str]:
    """``samtools`` as it must be invoked - inside *env* when conda envs are in use."""
    if not env:
        return ["samtools"]
    from .command import _run_flags

    return [runner.env_manager, "run"] + _run_flags(runner.env_manager) + ["-n", env, "samtools"]


def contig_lengths(bam: str, env: Optional[str] = None) -> Dict[str, int]:
    """Contig names and lengths from *bam*'s header, in header order.

    From the BAM rather than from the assembly fasta on purpose: the header records what the reads
    were actually aligned against, so a table built from it cannot disagree with the alignment it
    came from. If the two differ, the alignment is the fact.
    """
    out = subprocess.run(
        _samtools(env) + ["view", "-H", bam],
        capture_output=True,
        text=True,
        check=True,
    )
    lengths: Dict[str, int] = {}
    for line in out.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        name = length = None
        for field in line.split("\t")[1:]:
            if field.startswith("SN:"):
                name = field[3:]
            elif field.startswith("LN:"):
                length = int(field[3:])
        if name is not None and length is not None:
            lengths[name] = length
    if not lengths:
        raise ValueError("%s has no @SQ lines in its header, so it names no contigs." % bam)
    return lengths


def _identity_expression(percent_identity: float) -> List[str]:
    """A ``samtools view`` filter expression keeping reads at or above *percent_identity*.

    Identity is ``1 - NM/qlen``: NM counts mismatches and indel bases. Written as a multiplication
    rather than a division so that a read with no aligned length cannot divide by zero, and
    expressed over ``[NM]`` so a BAM without NM tags fails loudly at the filter instead of silently
    counting every read.

    Known difference from jgi, not yet settled. ``qlen`` in a samtools filter expression is the
    query length implied by the CIGAR, which includes soft-clipped and inserted bases; jgi's
    ``percentIdentity`` is over the *aligned* length alone. For a 150 bp read with 90 bp
    soft-clipped, this permits ``NM <= 15`` over 60 aligned bases - about 75% identity passing a
    nominal 90% floor - and minimap2 ``-ax sr``, the default aligner, soft-clips routinely. So the
    effective floor is looser than the number asked for, and looser than the sweep in
    PERCENT_IDENTITY's note was measuring. The aggregate effect on the means is evidently small
    (Pearson 0.9992 against jgi over all 210,632 CAMI2 contigs), which is why this is recorded here
    rather than changed: tightening it to the aligned length would move every depth table and so
    every binning result, and that is a change to re-measure deliberately, not to slip in.
    """
    if percent_identity <= 0:
        return []
    keep = (100.0 - percent_identity) / 100.0
    return ["-e", "[NM] <= (qlen * %g)" % keep]


def write_depth_table(
    bams: Sequence[str],
    out_path: str,
    *,
    samples: Optional[Sequence[str]] = None,
    edge_bases: int = EDGE_BASES,
    percent_identity: float = PERCENT_IDENTITY,
    min_map_qual: int = MIN_MAP_QUAL,
    min_contig_len: int = 1,
    env: Optional[str] = None,
    threads: int = 1,
) -> str:
    """Write a contig x sample depth table for *bams* to *out_path*. Returns *out_path*.

    *samples* names the columns; it defaults to each BAM's filename without its extension, which is
    how every caller in this project names them anyway. The column order follows *bams*, so two
    tables built from the same list are comparable column by column.
    """
    if not bams:
        raise ValueError("no BAMs given, so there is nothing to measure depth from.")
    names = list(samples) if samples else [_sample_name(b) for b in bams]
    if len(names) != len(bams):
        raise ValueError(
            "%d sample name(s) for %d BAM(s); they have to correspond one to one."
            % (len(names), len(bams))
        )

    lengths = contig_lengths(bams[0], env=env)
    order = [c for c, n in lengths.items() if n >= min_contig_len]
    index = {contig: i for i, contig in enumerate(order)}

    # Two running sums per contig per sample. Mean and variance come out of these at the end, so no
    # contig's per-base depths are ever all in memory at once.
    total = np.zeros((len(order), len(bams)), dtype=np.float64)
    total_sq = np.zeros((len(order), len(bams)), dtype=np.float64)
    # The divisor is the contig's trimmed length, from the header - not the number of positions the
    # pileup emitted. ``samtools depth -a`` reports all positions of every reference that has at least
    # one read, and says nothing at all about a reference with none, so counting emitted positions
    # makes the divisor depend on which samples happened to cover the contig. Taking it from the
    # length instead makes it identical across samples, which is both what makes the columns
    # comparable and what jgi does. (Found by comparing against jgi: a contig with no reads in the
    # first BAM got a divisor of 1 in all the others, and its depth came out as the sum over the whole
    # contig - a 4 Mbp contig reporting a depth of 1.6 million.)
    widths = trimmed_widths(np.array([lengths[c] for c in order], dtype=np.int64), edge_bases)

    workers = depth_workers(len(bams), threads)
    per_bam = max(1, threads // workers) if workers > 1 else threads

    def measure(job: Tuple[int, str]) -> None:
        column, bam = job
        comm("measuring per-base depth in %s" % os.path.basename(bam))
        _accumulate_bam(
            bam,
            column,
            lengths,
            index,
            edge_bases,
            percent_identity,
            min_map_qual,
            per_bam,
            env,
            total,
            total_sq,
        )

    jobs = list(enumerate(bams))
    if workers == 1:
        for job in jobs:
            measure(job)
    else:
        comm(
            "measuring %d BAM(s) %d at a time, %d samtools thread(s) each: one pileup does not "
            "use a whole machine" % (len(bams), workers, per_bam)
        )
        # Each BAM writes its own column of `total`/`total_sq` and reads nothing another writes,
        # so no lock is needed. Exceptions surface on result(), and the first one stops the rest
        # rather than letting a half-measured table be written.
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for future in [pool.submit(measure, job) for job in jobs]:
                future.result()

    _write(out_path, order, lengths, names, total, total_sq, widths)
    return out_path


def depth_workers(bam_count: int, threads: int) -> int:
    """How many pileups to run at once, given *threads* for the whole stage.

    Never more than there are BAMs, never more than MAX_DEPTH_WORKERS, and never so many that a
    pileup gets less than one thread. One worker when the budget is small, which is also what keeps
    a single-BAM call exactly as it was.
    """
    if bam_count <= 1 or threads <= DEPTH_THREADS_PER_BAM:
        return 1
    return max(1, min(bam_count, MAX_DEPTH_WORKERS, threads // DEPTH_THREADS_PER_BAM))


def trimmed_widths(contig_len: np.ndarray, edge_bases: int) -> np.ndarray:
    """How many bases of each contig are averaged over, after trimming both ends.

    A contig too short to have a middle left (``<= 2 * edge_bases``) is measured whole: dropping it
    would lose a contig the binner was asked to bin, and measuring 0 bases of it is not an option.
    """
    trimmable = contig_len > 2 * edge_bases
    return np.where(trimmable, contig_len - 2 * edge_bases, contig_len).astype(np.float64)


def _accumulate_bam(
    bam: str,
    column: int,
    lengths: Dict[str, int],
    index: Dict[str, int],
    edge_bases: int,
    percent_identity: float,
    min_map_qual: int,
    threads: int,
    env: Optional[str],
    total: np.ndarray,
    total_sq: np.ndarray,
) -> None:
    """Pile up one BAM and fold it into *column* of the running sums.

    One BAM at a time rather than all of them in one ``samtools depth`` call, because the identity
    filter has to happen first and it is a ``samtools view`` expression - and view reads one file.
    Per-BAM is the better shape anyway: the row sequence no longer has to line up across samples, and
    a BAM that is corrupt names itself in the error.
    """
    identity = _identity_expression(percent_identity)
    extra_threads = ["-@", str(max(0, threads - 1))] if threads > 1 else []
    view = subprocess.Popen(
        _samtools(env) + ["view", "-u"] + extra_threads + identity + [bam],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    pileup = None
    try:
        pileup = subprocess.Popen(
            _samtools(env)
            + ["depth", "-a", "-Q", str(min_map_qual), "-G", str(EXCLUDED_FLAGS)]
            + extra_threads
            + ["-"],
            stdin=view.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    finally:
        # Closed on both paths: if the second Popen raised, view is left writing into a pipe
        # nothing will ever read, so it has to see EPIPE. On the normal path this is what lets
        # view notice that the pileup died first.
        if view.stdout is not None:
            view.stdout.close()
        if pileup is None:
            view.kill()
            if view.stderr is not None:
                view.stderr.close()
            view.wait()
    assert pileup.stdout is not None
    # Both stderrs are drained by their own thread, started before anything reads stdout. These
    # are per-read tools over hundreds of millions of reads: one warning per record fills the 64 KiB
    # pipe buffer, at which point samtools blocks on the write, stops producing rows, and
    # pd.read_csv below blocks forever on a stream that will never end. Reading the pipes after the
    # depth stream is consumed is too late - the deadlock happens while it is being consumed.
    errors: Dict[str, str] = {}
    drains = [
        threading.Thread(target=_drain, args=(pipe, errors, name), daemon=True)
        for name, pipe in (("view", view.stderr), ("pileup", pileup.stderr))
    ]
    for thread in drains:
        thread.start()
    try:
        frames = pd.read_csv(
            pileup.stdout,
            sep="\t",
            header=None,
            names=["contig", "pos", "depth"],
            chunksize=CHUNK_ROWS,
            dtype={"contig": str, "pos": np.int64, "depth": np.float64},
            engine="c",
        )
        for chunk in frames:
            _accumulate(chunk, column, lengths, index, edge_bases, total, total_sq)
    finally:
        pileup.stdout.close()
        pileup_rc, view_rc = pileup.wait(), view.wait()
        for thread in drains:
            thread.join()
        pileup_err = errors.get("pileup", "")
        view_err = errors.get("view", "")
    if view_rc != 0:
        raise RuntimeError(
            "samtools view failed on %s (exit %d).\n%s" % (bam, view_rc, view_err.strip())
        )
    if pileup_rc != 0:
        raise RuntimeError(
            "samtools depth failed on %s (exit %d).\n%s" % (bam, pileup_rc, pileup_err.strip())
        )


def _drain(pipe, into: Dict[str, str], name: str) -> None:
    """Read *pipe* to EOF into ``into[name]``, so the writer never blocks on a full buffer."""
    if pipe is None:
        return
    try:
        into[name] = pipe.read().decode(errors="replace")
    except (OSError, ValueError):
        into[name] = ""
    finally:
        pipe.close()


def _accumulate(
    chunk: pd.DataFrame,
    column: int,
    lengths: Dict[str, int],
    index: Dict[str, int],
    edge_bases: int,
    total: np.ndarray,
    total_sq: np.ndarray,
) -> None:
    """Fold one chunk of ``samtools depth`` rows into *column* of the running sums.

    A chunk can straddle contigs and the trim window has to be applied per contig, so the positions
    are masked by each row's own contig length rather than by one global rule.
    """
    rows = chunk["contig"].map(index)
    keep = rows.notna()
    if not keep.any():
        return
    chunk = chunk[keep.to_numpy()]
    rows = rows[keep].to_numpy(dtype=np.int64)

    contig_len = chunk["contig"].map(lengths).to_numpy(dtype=np.int64)
    pos = chunk["pos"].to_numpy(dtype=np.int64)
    # Trim both ends, except on a contig too short to have a middle left: there, measuring the whole
    # thing is less wrong than measuring none of it.
    trimmable = contig_len > 2 * edge_bases
    low = np.where(trimmable, edge_bases, 0)
    high = np.where(trimmable, contig_len - edge_bases, contig_len)
    inside = (pos > low) & (pos <= high)
    if not inside.any():
        return
    rows = rows[inside]
    depths = chunk["depth"].to_numpy(dtype=np.float64)[inside]

    np.add.at(total[:, column], rows, depths)
    np.add.at(total_sq[:, column], rows, depths * depths)


def _write(
    out_path: str,
    order: Sequence[str],
    lengths: Dict[str, int],
    names: Sequence[str],
    total: np.ndarray,
    total_sq: np.ndarray,
    widths: np.ndarray,
) -> None:
    """Write the table. Six significant figures, which is what the tool this replaces emitted."""
    positions = np.maximum(widths, 1.0)[:, None]
    mean = total / positions
    # E[x^2] - E[x]^2, floored at zero: a contig whose depth never varies gives a tiny negative
    # number here from rounding, and a negative variance in the table is a value metaBAT2 would read.
    variance = np.maximum(total_sq / positions - mean * mean, 0.0)
    totals = mean.sum(axis=1)

    header = list(LEADING_COLUMNS)
    for name in names:
        header += [name, "%s-var" % name]
    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for i, contig in enumerate(order):
            cells = [contig, str(lengths[contig]), "%g" % totals[i]]
            for j in range(len(names)):
                cells += ["%g" % mean[i, j], "%g" % variance[i, j]]
            out.write("\t".join(cells) + "\n")


def _sample_name(bam: str) -> str:
    return os.path.splitext(os.path.basename(bam))[0]


def read_depth_table(path: str) -> Tuple[List[str], pd.DataFrame]:
    """Read a depth table, returning its sample names and the frame. For tests and for scoring."""
    frame = pd.read_csv(path, sep="\t")
    samples = [c for c in frame.columns[3:] if not c.endswith("-var")]
    return samples, frame
