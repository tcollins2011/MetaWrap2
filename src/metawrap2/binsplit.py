"""Split a bin set back into per-sample bins, for multi-split binning.

Default cohort binning bins each assembly separately with a cross-sample depth table:
``coverage`` concatenates the assemblies, maps every sample against the catalogue, then slices the
depth table back per assembly, so each binning run sees N samples' coverage of its own contigs.

Multi-split binning bins the catalogue itself, so the model is fitted on every sample's contigs at
once, then splits the resulting clusters by the sample each contig came from. The split happens
after the model instead of before it. Söding's analysis reports that splitting the embedding space
by sample gives 19% more bins and 79% more high-quality bins than splitting the final clusters.
Fitted on the catalogue, the model sees an organism's contigs from every sample when it clusters.

VAMB (``-o SEP``) and SemiBin2 (``multi_easy_bin -s SEP``) can do the split themselves. The current
code passes VAMB ``-o ''``, which is correct for a single assembly but not for a catalogue. For the
other binners this module does the split.
"""

from __future__ import annotations

import os
import shutil
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

from .io.seqio import contig_id, iter_fasta

#: The separator ``coverage`` writes between a sample label and a contig name. Defined once so a
#: catalogue built by one module can be split by another.
SEPARATOR = "~"


@dataclass(frozen=True)
class SplitBin:
    """One per-sample piece of one catalogue bin."""

    #: The bin as the binner produced it, over the whole catalogue.
    source_bin: str
    #: The sample label its contigs came from.
    label: str
    #: Contig names with the label prefix removed, so the piece describes that sample's assembly.
    contigs: Tuple[str, ...]
    total_bp: int

    @property
    def name(self) -> str:
        """``<source bin>.<label>``, so a split piece is traceable to the cluster it came from."""
        return "%s.%s" % (self.source_bin, self.label)


def split_contig_names(
    names: Iterable[str], separator: str = SEPARATOR
) -> Tuple[Dict[str, List[str]], List[str]]:
    """Group prefixed contig names by sample label.

    Returns ``(label -> original contig names, names with no prefix)``. Unprefixed names are
    returned, not dropped or assigned: they mean the catalogue was not built by ``coverage``, and
    assigning them could attribute sequence to the wrong sample.
    """
    grouped: Dict[str, List[str]] = defaultdict(list)
    unprefixed: List[str] = []
    for name in names:
        label, found, contig = name.partition(separator)
        if not found:
            unprefixed.append(name)
            continue
        grouped[label].append(contig)
    return dict(grouped), unprefixed


def read_split(
    path: str, separator: str = SEPARATOR
) -> Tuple[Dict[str, List[Tuple[str, str]]], List[str]]:
    """One pass over a bin's fasta: ``(label -> [(contig, sequence)], unprefixed names)``."""
    per_label: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    unprefixed: List[str] = []
    for header, sequence in iter_fasta(path):
        name = contig_id(header)
        label, found, contig = name.partition(separator)
        if not found:
            unprefixed.append(name)
            continue
        per_label[label].append((contig, sequence))
    return dict(per_label), unprefixed


def _pieces_of(bin_name: str, per_label: Dict[str, List[Tuple[str, str]]]) -> List[SplitBin]:
    return [
        SplitBin(
            source_bin=bin_name,
            label=label,
            contigs=tuple(contig for contig, _seq in records),
            total_bp=sum(len(seq) for _contig, seq in records),
        )
        for label, records in sorted(per_label.items())
    ]


def split_bin(
    path: str, bin_name: str, separator: str = SEPARATOR
) -> Tuple[List[SplitBin], List[str]]:
    """Split one catalogue bin's fasta into per-sample pieces, without writing anything.

    A bin whose contigs all come from one sample yields a single piece; this is not special-cased.
    """
    per_label, unprefixed = read_split(path, separator)
    return _pieces_of(bin_name, per_label), unprefixed


def write_split_bins(
    path: str,
    bin_name: str,
    out_dir: str,
    separator: str = SEPARATOR,
    min_bp: int = 0,
) -> Tuple[List[SplitBin], List[str]]:
    """Split one bin and write each piece into ``out_dir/<label>/``.

    The original contig names are restored, because downstream tools pair a piece with that
    sample's own (unprefixed) assembly.

    *min_bp* drops pieces below a size, such as the few kilobases a multi-sample cluster often
    takes from a sample that barely contains the organism. The default 0 keeps everything and
    leaves size filtering to whatever scores the bins.
    """
    # One read of the fasta, used for both the piece summaries and the sequences written out.
    by_label, unprefixed = read_split(path, separator)
    pieces = _pieces_of(bin_name, by_label)
    kept: List[SplitBin] = []

    for piece in pieces:
        if piece.total_bp < min_bp:
            continue
        label_dir = os.path.join(out_dir, piece.label)
        os.makedirs(label_dir, exist_ok=True)
        with open(os.path.join(label_dir, piece.name + ".fasta"), "w") as out:
            out.writelines(
                ">%s\n%s\n" % (contig, sequence) for contig, sequence in by_label[piece.label]
            )
        kept.append(piece)
    return kept, unprefixed


def bin_files(bins_dir: str, extensions: Sequence[str] = (".fa", ".fasta", ".fna")) -> List[str]:
    """Bin fastas in *bins_dir*, sorted, so a split is reproducible."""
    if not os.path.isdir(bins_dir):
        return []
    return sorted(
        os.path.join(bins_dir, name)
        for name in os.listdir(bins_dir)
        if name.endswith(tuple(extensions))
    )


@dataclass
class SplitReport:
    """What a split did, for the log and for the caller to check."""

    bins_in: int
    pieces_out: int
    labels: Tuple[str, ...]
    #: Bins whose contigs came from more than one sample. If zero, the mode found nothing that
    #: per-sample binning would not.
    multi_sample_bins: int
    unprefixed_contigs: int

    @property
    def worthwhile(self) -> bool:
        return self.multi_sample_bins > 0


def flatten_split_bins(split_dir: str, flat_dir: str) -> int:
    """Link every per-sample piece under *split_dir* into one directory, keeping its own name.

    The nested layout (`<out>/<sample>/bin_001.fasta`) suits per-sample refinement. Tools that treat
    a bin set as one folder of FASTA files need the flat copy; the benchmark scorer collects
    top-level `*_bins` directories and would otherwise score the unsplit catalogue bins.

    Uses hard links where the filesystem allows, and copies otherwise.
    """
    if not os.path.isdir(split_dir):
        return 0
    os.makedirs(flat_dir, exist_ok=True)
    linked = 0
    for label in sorted(os.listdir(split_dir)):
        inner = os.path.join(split_dir, label)
        if not os.path.isdir(inner):
            continue
        for name in sorted(os.listdir(inner)):
            if not name.endswith((".fa", ".fasta", ".fna")):
                continue
            source = os.path.join(inner, name)
            # The filename already carries the sample (`bin_001.sample_0.fasta`), so it is unique.
            dest = os.path.join(flat_dir, name)
            if os.path.exists(dest):
                continue
            try:
                os.link(source, dest)
            except OSError:
                shutil.copyfile(source, dest)
            linked += 1
    return linked


def split_bin_set(
    bins_dir: str,
    out_dir: str,
    separator: str = SEPARATOR,
    min_bp: int = 0,
) -> SplitReport:
    """Split every bin in *bins_dir* into per-sample pieces under *out_dir*."""
    paths = bin_files(bins_dir)
    labels: set = set()
    pieces_out = 0
    multi = 0
    unprefixed = 0
    for path in paths:
        bin_name = os.path.splitext(os.path.basename(path))[0]
        pieces, missing = write_split_bins(path, bin_name, out_dir, separator, min_bp)
        pieces_out += len(pieces)
        unprefixed += len(missing)
        labels.update(piece.label for piece in pieces)
        if len(pieces) > 1:
            multi += 1
    return SplitReport(
        bins_in=len(paths),
        pieces_out=pieces_out,
        labels=tuple(sorted(labels)),
        multi_sample_bins=multi,
        unprefixed_contigs=unprefixed,
    )
