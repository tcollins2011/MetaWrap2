"""Core bin-refinement algorithms: consolidation and contig dereplication.

This is the scientific heart of the bin_refinement module, kept faithful to the
behaviour it has always had. Two changes over the original scripts:

* Contig identity is the full FASTA header string. The old code parsed a length out
  of the header with ``line.split('_')[3]`` (a SPAdes-only assumption) in one code
  path; here lengths always come from the actual sequence, so any assembler works.
* Where a source bin overlapped several candidate replacements, the old code's choice
  depended on Python's dict iteration order (which changed between Python 2 and 3).
  We deterministically pick the single best-scoring qualifying replacement, matching
  the documented intent ("decide which version is best").

Bin quality score (unchanged): completeness - 5 * contamination.
"""

from __future__ import annotations

import os
import shutil
import sys
from typing import Dict, List, Optional, Tuple

from .constants import BIN_EXTENSION, OVERLAP_THRESHOLD, STATS_SUFFIX, bin_filename, bin_stem
from .io.seqio import contig_id, iter_fasta
from .utils import bin_name_from_filename, quality_score


def _read_stats(path: str) -> Tuple[str, Dict[str, str], Dict[str, Tuple[float, ...]]]:
    """Parse a CheckM ``.stats`` file.

    Returns ``(header_line, {bin.fa: raw_line}, {bin.fa: (completeness, contamination, extra)})``
    where ``extra`` is the N50 column when present, used only as a dereplication tiebreaker
    (weighted 1e-10, so it separates otherwise-identical bins without shifting ranking).
    """
    header = ""
    raw: Dict[str, str] = {}
    stats: Dict[str, Tuple[float, ...]] = {}
    with open(path) as fh:
        for line in fh:
            if "completeness" in line:
                header = line
                continue
            cut = line.rstrip("\n").split("\t")
            if len(cut) < 3:
                continue
            name = cut[0] + BIN_EXTENSION
            raw[name] = line
            extra = float(cut[5]) if len(cut) > 5 else 0.0
            stats[name] = (float(cut[1]), float(cut[2]), extra)
    return header, raw, stats


def _bin_contig_lengths(bin_path: str) -> Dict[str, int]:
    """Map each contig id -> length for one bin FASTA (compression-transparent).

    Keyed on :func:`contig_id`, not the raw header, so a binner's extra header annotations
    cannot stop two copies of the same contig from matching.
    """
    return {contig_id(name): len(seq) for name, seq in iter_fasta(bin_path)}


def _overlap_percent(a: Dict[str, int], b: Dict[str, int]) -> float:
    """Return the max directional length-overlap percentage between two bins."""
    match_a = sum(b[c] for c in a if c in b)
    mismatch_a = sum(a[c] for c in a if c not in b)
    match_b = sum(a[c] for c in b if c in a)
    mismatch_b = sum(b[c] for c in b if c not in a)
    ratio_a = 100.0 * match_a / (match_a + mismatch_a) if (match_a + mismatch_a) else 0.0
    ratio_b = 100.0 * match_b / (match_b + mismatch_b) if (match_b + mismatch_b) else 0.0
    return max(ratio_a, ratio_b)


def consolidate(
    folder_1: str,
    folder_2: str,
    stats_1: str,
    stats_2: str,
    out_folder: str,
    min_completion: float,
    max_contamination: float,
) -> int:
    """Merge two bin sets, keeping the best version of each bin. Returns the bin count.

    Mirrors ``consolidate_two_sets_of_bins.py``: for each good bin in set 1, find bins in
    set 2 with >=80% overlap; if the best such match scores higher, take it, else keep the
    original. Finally add good bins from set 2 that never matched anything in set 1.
    """

    def good(stats):
        return {
            name: v
            for name, v in stats.items()
            if v[0] > min_completion and v[1] < max_contamination
        }

    header1, raw1, stats1 = _read_stats(stats_1)
    _header2, raw2, stats2 = _read_stats(stats_2)
    good1 = good(stats1)
    good2 = good(stats2)

    lengths1 = {b: _bin_contig_lengths(os.path.join(folder_1, b)) for b in good1}
    lengths2 = {b: _bin_contig_lengths(os.path.join(folder_2, b)) for b in good2}

    os.makedirs(out_folder, exist_ok=True)
    new_summary = [header1] if header1 else []
    matched_in_2: Dict[str, None] = {}
    bin_ct = 1

    for bin_1 in sorted(good1):
        base_score = quality_score(*stats1[bin_1][:2])
        # find the single best qualifying replacement in set 2 (deterministic)
        best_bin_2: Optional[str] = None
        best_bin_2_score = base_score
        for bin_2 in sorted(good2):
            if _overlap_percent(lengths1[bin_1], lengths2[bin_2]) < OVERLAP_THRESHOLD:
                continue
            matched_in_2[bin_2] = None
            score_2 = quality_score(*stats2[bin_2][:2])
            if score_2 > best_bin_2_score:
                best_bin_2 = bin_2
                best_bin_2_score = score_2

        if best_bin_2 is not None:
            src, raw = os.path.join(folder_2, best_bin_2), raw2[best_bin_2]
        else:
            src, raw = os.path.join(folder_1, bin_1), raw1[bin_1]
        shutil.copy(src, os.path.join(out_folder, bin_filename(bin_ct)))
        new_summary.append("%s\t%s" % (bin_stem(bin_ct), "\t".join(raw.split("\t")[1:])))
        bin_ct += 1

    for bin_2 in sorted(good2):
        if bin_2 in matched_in_2:
            continue
        shutil.copy(os.path.join(folder_2, bin_2), os.path.join(out_folder, bin_filename(bin_ct)))
        new_summary.append("%s\t%s" % (bin_stem(bin_ct), "\t".join(raw2[bin_2].split("\t")[1:])))
        bin_ct += 1

    with open(out_folder + STATS_SUFFIX, "w") as fh:
        fh.write("".join(new_summary))
    return bin_ct


def dereplicate(stats_file: str, bins_folder: str, out_folder: str, mode: str = "best") -> None:
    """Remove contigs shared between bins, keeping each in its best bin (or removing it).

    Mirrors ``dereplicate_contigs_in_bins.py``. ``mode="remove"`` drops any contig that
    appears in more than one bin; otherwise a shared contig is kept only in the
    highest-scoring bin. Score: completeness - 5*contamination + 1e-10*extra (tiebreak).
    """
    bin_scores: Dict[str, float] = {}
    with open(stats_file) as fh:
        for line in fh:
            if "completeness" in line:
                continue
            cut = line.rstrip("\n").split("\t")
            if len(cut) < 3:
                continue
            extra = float(cut[5]) if len(cut) > 5 else 0.0
            bin_scores[cut[0]] = quality_score(float(cut[1]), float(cut[2]), extra)

    def score(bin_name: str) -> float:
        # A bin file with no row in the .stats file (the two can drift apart if a run was
        # interrupted between writing bins and re-scoring them) used to raise KeyError here
        # and abort the whole refinement. Treat it as the worst possible bin instead, so it
        # never wins a contig but also never crashes the run.
        if bin_name not in bin_scores:
            missing_stats.add(bin_name)
            return float("-inf")
        return bin_scores[bin_name]

    missing_stats: set = set()
    contig_owner: Dict[str, Optional[str]] = {}
    bin_files = sorted(os.listdir(bins_folder))
    for bin_file in bin_files:
        name = bin_name_from_filename(bin_file)
        for header, _seq in iter_fasta(os.path.join(bins_folder, bin_file)):
            contig = contig_id(header)
            if contig not in contig_owner:
                contig_owner[contig] = name
            elif mode == "remove":
                contig_owner[contig] = None
            else:
                owner = contig_owner[contig]
                if owner is not None and score(name) > score(owner):
                    contig_owner[contig] = name
    if missing_stats:
        sys.stderr.write(
            "WARNING: %d bin(s) in %s have no row in %s, so they lost every contested "
            "contig: %s\n"
            % (len(missing_stats), bins_folder, stats_file, ", ".join(sorted(missing_stats)[:10]))
        )

    os.makedirs(out_folder, exist_ok=True)
    for bin_file in bin_files:
        name = bin_name_from_filename(bin_file)
        kept: List[Tuple[str, str]] = [
            (contig_id(h), s)
            for h, s in iter_fasta(os.path.join(bins_folder, bin_file))
            if contig_owner.get(contig_id(h)) == name
        ]
        if not kept:
            continue
        with open(os.path.join(out_folder, bin_file), "w") as out:
            for contig, seq in kept:
                out.write(">%s\n%s\n" % (contig, seq))
