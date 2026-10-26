"""Per-contig, per-sample coverage depth and intra-contig depth variance.

Replaces ``jgi_summarize_bam_contig_depths``, which MetaWrap2 previously called. Every binner
depends on the depth table, and keeping it in tested project code makes intermediate values
checkable. (A MaxBin2 defect, binning on wrong depths, was in how the jgi table was converted, and
was hard to diagnose through the external binary.)

The algorithm matches jgi's so numbers stay comparable with earlier measurements:

* Reads are filtered by mapping quality (``min_map_qual``) and alignment identity
  (``percent_identity``); unmapped, secondary, QC-failed and duplicate records are excluded. The
  identity floor is 90%, the only change from jgi's defaults (97%); see :data:`PERCENT_IDENTITY`
  for the measurement, including 26% fewer clean bins at 97% on strain-rich data.
* ``edge_bases`` (75) are trimmed from each contig end before averaging, because coverage drops
  where reads cannot overhang the edge, which biases short contigs most. A contig of
  ``<= 2 * edge_bases`` is measured whole.
* Variance is the variance of per-base depth inside the trimmed window, not across samples.
  metaBAT2 uses it as a second signal and rejects a table whose variance column is constant.

  With several samples the variance column carries little signal. metaBAT2 was given four depth
  tables with the same means and contigs and a variance column that was real, zeroed, shuffled
  between contigs, or randomly resampled, and the bins were scored with AMBER. At three samples
  all four were within 1.5% F1 (the shuffled one scored highest). At one sample real variance
  was worth about 10% (F1 0.457 against 0.416 zeroed). Cross-sample covariance carries the same
  information, so a cheaper estimator for cohorts mainly has to get the means right (where fairy
  falls short, see modules/coverage.py). One dataset, metaBAT2 only; MaxBin2 and CONCOCT do not
  read the variance column.

Per-base depths come from ``samtools depth``, which piles up in C and writes one row per position.
The rest is numpy over chunks, so the pileup dominates the cost. Only two running sums per contig
per sample are kept, so memory does not grow with per-base data.
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

#: Bases trimmed from each end of a contig before averaging. jgi's ``--maxEdgeBases`` default, about
#: one read length: the region where coverage drops because alignments cannot start before the
#: contig does.
EDGE_BASES = 75
#: Minimum alignment identity, in percent, for a read to count toward a contig's depth.
#:
#: 90 instead of jgi's 97, the only change from ``jgi_summarize_bam_contig_depths``' defaults.
#: 97 was inherited without evaluation. On CAMI2 GI it discards half the mapped bases against the
#: true reference genomes (the gold assembly), where divergence is sequencing error plus strain
#: variation.
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
#: "clean HQ" counts bins that are >=70/<=10 and non-chimeric by GUNC. Compared with 90, the 97
#: floor gives 14 fewer bins and lower completeness, purity and clean HQ count (35 against 47, a 34%
#: gain at 90). A looser floor could let a relative's reads merge two strains and lower purity, but
#: AMBER purity and GUNC both show no such loss on this data.
#:
#: A floor is kept (90 instead of 0) because a read 20% divergent from a contig likely comes from
#: another organism; 0 and 90 are within 0.002 AMBER F1. Measured on one dataset with one binner,
#: so it is a module constant and an argument.
PERCENT_IDENTITY = 90.0
#: jgi's default, used by metaBAT2's published validation; needed to reproduce metaWRAP 1 results
#: or jgi depth tables.
JGI_PERCENT_IDENTITY = 97.0
#: Minimum MAPQ. 0, as jgi's ``--minMapQual`` default: multi-mapping reads are kept, since a repeat
#: shared between two genomes is coverage of both. Raise it to drop repeats instead.
MIN_MAP_QUAL = 0
#: SAM flags excluded from the pileup: unmapped, secondary, QC fail, PCR/optical duplicate.
#: The primary record of a multi-mapping read still counts, so each read contributes once.
EXCLUDED_FLAGS = 0x4 | 0x100 | 0x200 | 0x400
#: Rows parsed per chunk; keeps a chunk to a few hundred MB at ten samples.
CHUNK_ROWS = 4_000_000

#: How many samtools threads one BAM's pileup gets, and so how many BAMs are measured at once.
#:
#: `samtools depth` does not scale with threads: a user measured ~2 minutes per BAM on a large
#: machine at ~25% utilisation, over hours for 50 BAMs in sequence. Each BAM fills its own column
#: of the running sums, so several BAMs are measured at once. Four threads per pileup because
#: `samtools view` (identity filter, decompression) does use threads, and the pandas C parser and
#: numpy accumulation release the GIL.
DEPTH_THREADS_PER_BAM = 4

#: Maximum pileups at once. Each worker holds up to CHUNK_ROWS rows, so this bounds memory as well
#: as the number of samtools processes.
MAX_DEPTH_WORKERS = 16

#: The first three columns, then two per sample. Defined once; three places used to spell it out
#: and one had it wrong.
LEADING_COLUMNS = ("contigName", "contigLen", "totalAvgDepth")


def _samtools(env: Optional[str]) -> List[str]:
    """``samtools`` as it must be invoked - inside *env* when conda envs are in use."""
    if not env:
        return ["samtools"]
    from .command import _run_flags

    return [runner.env_manager, "run"] + _run_flags(runner.env_manager) + ["-n", env, "samtools"]


def contig_lengths(bam: str, env: Optional[str] = None) -> Dict[str, int]:
    """Contig names and lengths from *bam*'s header, in header order.

    Taken from the BAM instead of the assembly FASTA, so the table always matches the reference the
    reads were aligned to.
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

    Identity is ``1 - NM/qlen``, where NM counts mismatches and indel bases. A multiplication avoids
    division by zero, and ``[NM]`` makes a BAM without NM tags fail at the filter instead of
    passing every read.

    Open difference from jgi: ``qlen`` here is the CIGAR query length, including soft-clipped and
    inserted bases, while jgi's ``percentIdentity`` uses the aligned length. A 150 bp read with
    90 bp soft-clipped passes with ``NM <= 15`` over 60 aligned bases (about 75% identity), and
    minimap2 ``-ax sr`` soft-clips often, so the effective floor is looser than requested. The
    effect on means is small (Pearson 0.9992 against jgi over 210,632 CAMI2 contigs). Changing it
    would shift every depth table and binning result, so it needs a separate re-measurement.
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

    *samples* names the columns and defaults to each BAM's filename without extension. Column order
    follows *bams*.
    """
    if not bams:
        raise ValueError("no BAMs given; there is nothing to measure depth from.")
    names = list(samples) if samples else [_sample_name(b) for b in bams]
    if len(names) != len(bams):
        raise ValueError(
            "%d sample name(s) for %d BAM(s); they have to correspond one to one."
            % (len(names), len(bams))
        )

    lengths = contig_lengths(bams[0], env=env)
    order = [c for c, n in lengths.items() if n >= min_contig_len]
    index = {contig: i for i, contig in enumerate(order)}

    # Two running sums per contig per sample, from which mean and variance are computed at the end.
    total = np.zeros((len(order), len(bams)), dtype=np.float64)
    total_sq = np.zeros((len(order), len(bams)), dtype=np.float64)
    # The divisor is the trimmed contig length from the header, as in jgi, not the number of
    # positions emitted. ``samtools depth -a`` emits nothing for a reference with no reads, so
    # counting positions made the divisor vary by sample: a contig with no reads in the first BAM
    # got a divisor of 1 in the others, and a 4 Mbp contig reported a depth of 1.6 million.
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
            "measuring %d BAM(s), %d at a time, %d samtools thread(s) "
            "per BAM" % (len(bams), workers, per_bam)
        )
        # Each BAM writes only its own column of `total`/`total_sq`, so no lock is needed. The
        # first exception raised by result() stops the stage before a partial table is written.
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for future in [pool.submit(measure, job) for job in jobs]:
                future.result()

    _write(out_path, order, lengths, names, total, total_sq, widths)
    return out_path


def depth_workers(bam_count: int, threads: int) -> int:
    """How many pileups to run at once, given *threads* for the whole stage.

    At most the number of BAMs and MAX_DEPTH_WORKERS, and each worker gets DEPTH_THREADS_PER_BAM
    threads. One worker for a single BAM or a small thread budget.
    """
    if bam_count <= 1 or threads <= DEPTH_THREADS_PER_BAM:
        return 1
    return max(1, min(bam_count, MAX_DEPTH_WORKERS, threads // DEPTH_THREADS_PER_BAM))


def trimmed_widths(contig_len: np.ndarray, edge_bases: int) -> np.ndarray:
    """How many bases of each contig are averaged over, after trimming both ends.

    A contig of ``<= 2 * edge_bases`` is measured whole so it is not lost to the binner.
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

    One BAM per call because the identity filter is a ``samtools view`` expression, and view reads
    one file. Rows then need not line up across samples, and errors name the failing BAM.
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
        # Closed on both paths so view gets EPIPE if the pileup failed to start or exits first.
        if view.stdout is not None:
            view.stdout.close()
        if pileup is None:
            view.kill()
            if view.stderr is not None:
                view.stderr.close()
            view.wait()
    assert pileup.stdout is not None
    # Drain both stderrs in threads started before stdout is read. One warning per read can fill
    # the 64 KiB pipe buffer, which blocks samtools and then pd.read_csv below (a deadlock).
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

    A chunk can span contigs, so the trim window is applied using each row's contig length.
    """
    rows = chunk["contig"].map(index)
    keep = rows.notna()
    if not keep.any():
        return
    chunk = chunk[keep.to_numpy()]
    rows = rows[keep].to_numpy(dtype=np.int64)

    contig_len = chunk["contig"].map(lengths).to_numpy(dtype=np.int64)
    pos = chunk["pos"].to_numpy(dtype=np.int64)
    # Trim both ends, except on contigs too short to trim, which are measured whole.
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
    """Write the table with six significant figures, as jgi's tool did."""
    positions = np.maximum(widths, 1.0)[:, None]
    mean = total / positions
    # E[x^2] - E[x]^2, floored at zero because rounding gives small negatives for constant depth.
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
