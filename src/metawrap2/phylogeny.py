"""Trees from alignments the pipeline already produces, and what a support value does not mean.

Two alignments exist by the time anyone wants a tree, and neither needs new computation:

* **GTDB-Tk's** concatenated marker alignment (``*.user_msa.fasta``), 120 bacterial or 53 archaeal
  single-copy markers. This is the across-species tree - where a catalogue's genomes sit relative to
  each other and to GTDB's reference set.
* **PPanGGOLiN's** core-gene alignment for one species group, written with ``--phylo``. This is the
  within-species tree - which strain in which sample is closest to which.

They answer different questions and the module does not pretend otherwise: a tree from 120 universal
markers cannot resolve strains of one species, and a tree from one species' core genes says nothing
about where that species sits.

IQ-TREE rather than FastTree by default. ModelFinder picks a substitution model from the alignment
instead of assuming one, and the ultrafast bootstrap gives support values calibrated to be
interpreted at 95% rather than at FastTree's 0.8-ish SH-like values. FastTree stays available
because at a few hundred genomes IQ-TREE's model selection becomes the dominant cost.

**What a support value means, and what it does not.** Ultrafast bootstrap support is the frequency
with which a branch appears in resampled alignments: it measures whether the *alignment* supports
the branch, not whether the branch is true. On a concatenation of marker genes from MAGs, three
things it cannot see: genes missing because the MAG is incomplete (so the branch rests on fewer
columns than the alignment's width suggests), horizontal transfer, and contamination placing another
organism's marker in this genome's row. So :func:`support_summary` reports the occupancy of the
alignment beside the support, and :data:`SUPPORT_CAVEAT` travels with every tree.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

#: Ultrafast bootstrap support at or above which a branch is conventionally called supported.
#: IQ-TREE's own guidance for UFBoot, which is calibrated differently from a standard bootstrap -
#: 95 here is roughly as strict as 70 there, and reading UFBoot at 70 accepts branches the data does
#: not support.
SUPPORTED_BRANCH = 95.0

#: Replicates for the ultrafast bootstrap. IQ-TREE's minimum is 1000 and it refuses fewer.
BOOTSTRAP_REPLICATES = 1000

#: Taxa below which a tree is not a tree: three is the fewest that has a topology at all, and with
#: three there is exactly one unrooted topology, so the model and the support are decoration.
MIN_TAXA = 4

#: Fraction of an alignment column that must be a residue rather than a gap for the column to carry
#: information about every taxon. Reported rather than enforced - trimming is a separate decision -
#: because a concatenation of markers from incomplete MAGs is mostly gaps in places, and a tree
#: built on it rests on fewer columns than its width suggests.
OCCUPANCY_FLOOR = 0.5


@dataclass(frozen=True)
class Alignment:
    """An alignment's shape, which is what decides whether a tree from it means anything."""

    path: str
    taxa: Tuple[str, ...]
    #: Alignment width in columns.
    columns: int
    #: Per-taxon fraction of non-gap characters.
    occupancy: Dict[str, float]
    #: ``"nucleotide"`` or ``"protein"``, measured from the residues. Which one decides the
    #: substitution model, and getting it wrong is not a subtle error: FastTree run with ``-nt``
    #: over a protein alignment reads every residue that is not A, C, G or T as ambiguous and
    #: returns a tree built from almost no information.
    kind: str = "protein"

    @property
    def size(self) -> int:
        return len(self.taxa)

    @property
    def mean_occupancy(self) -> Optional[float]:
        if not self.occupancy:
            return None
        return sum(self.occupancy.values()) / len(self.occupancy)

    @property
    def sparse_taxa(self) -> Tuple[str, ...]:
        """Taxa whose row is more gap than residue - the ones a branch cannot really be about."""
        return tuple(
            sorted(name for name, filled in self.occupancy.items() if filled < OCCUPANCY_FLOOR)
        )

    @property
    def is_a_tree(self) -> bool:
        return self.size >= MIN_TAXA


#: Residues counted as nucleotide when deciding an alignment's type: the four bases, U for RNA,
#: N for an unknown base, and the IUPAC ambiguity codes. T and N are also amino acids, which is why
#: the test is a *fraction* rather than a membership check.
NUCLEOTIDE_CHARS = frozenset("ACGTUNRYSWKMBDHV")

#: How much of an alignment has to be nucleotide-like to call it nucleotide. A protein alignment
#: reaches roughly 20-30% of this set by chance (A, C, G, T, N, W, Y, K, M, H, R, S, V, D are all
#: amino-acid codes too), and a real nucleotide alignment is above 99%, so anything in between is
#: unlike both and treated as protein - the safer default, since a protein alignment is what every
#: aligner upstream of this module writes.
NUCLEOTIDE_FRACTION = 0.9


def alignment_kind(counts: Dict[str, int]) -> str:
    """``"nucleotide"`` or ``"protein"``, from a residue histogram."""
    total = sum(counts.values())
    if not total:
        return "protein"
    nucleotide = sum(n for char, n in counts.items() if char.upper() in NUCLEOTIDE_CHARS)
    return "nucleotide" if nucleotide / total >= NUCLEOTIDE_FRACTION else "protein"


def read_alignment(path: str) -> Alignment:
    """Measure a FASTA alignment: its taxa, its width, and how much of each row is sequence.

    Gzip is handled because GTDB-Tk writes its MSA compressed, and a module that silently failed on
    the file its own documentation points at would be useless.
    """
    opener: object = open
    if path.endswith(".gz"):
        import gzip

        opener = gzip.open
    names: List[str] = []
    histogram: Dict[str, int] = {}
    residues: Dict[str, int] = {}
    lengths: Dict[str, int] = {}
    current: Optional[str] = None
    with opener(path, "rt") as handle:  # type: ignore[operator]
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                names.append(current)
                residues.setdefault(current, 0)
                lengths.setdefault(current, 0)
                continue
            if current is None:
                raise ValueError("%s does not start with a FASTA header" % path)
            lengths[current] += len(line)
            residues[current] += sum(1 for char in line if char not in "-.~?")
            for char in line:
                if char not in "-.~?":
                    histogram[char] = histogram.get(char, 0) + 1
    if not names:
        raise ValueError("%s holds no sequences" % path)
    width = max(lengths.values())
    return Alignment(
        path=path,
        taxa=tuple(names),
        columns=width,
        occupancy={name: (residues[name] / width if width else 0.0) for name in names},
        kind=alignment_kind(histogram),
    )


#: IQ-TREE writes the model it chose into its log and its .iqtree report, in this form.
_MODEL = re.compile(r"Best-fit model(?: according to BIC)?:\s*(\S+)")
#: Newick support values after a closing bracket, in the two forms IQ-TREE writes. With --alrt it
#: writes *both* tests as `)97.5/100:` - SH-aLRT first, ultrafast bootstrap second - and with only
#: -B it writes one value as `)100:`. A pattern that matches a single number reads the combined form
#: as zero branches, which is what the first real tree reported: 4 taxa, one internal branch
#: labelled 100/100, counted as none.
_SUPPORT = re.compile(r"\)(\d+(?:\.\d+)?)(?:/(\d+(?:\.\d+)?))?[:,;)]")


def model_from_report(path: str) -> Optional[str]:
    """The substitution model IQ-TREE's ModelFinder chose, from its report file.

    Reported because it is a decision the run made on the user's behalf, and two trees built under
    different models are not comparable even from the same alignment.
    """
    if not os.path.isfile(path):
        return None
    with open(path, errors="replace") as handle:
        for line in handle:
            match = _MODEL.search(line)
            if match:
                return match.group(1)
    return None


@dataclass(frozen=True)
class SupportSummary:
    """How much of a tree's branching the alignment actually supports."""

    branches: int
    supported: int
    threshold: float = SUPPORTED_BRANCH

    @property
    def fraction_supported(self) -> Optional[float]:
        if not self.branches:
            return None
        return self.supported / self.branches


def support_summary(newick_path: str, threshold: float = SUPPORTED_BRANCH) -> SupportSummary:
    """Count internal branches and how many clear *threshold*, from the Newick string.

    Parsed out of the tree rather than taken from IQ-TREE's log because the log reports the
    bootstrap settings, not the resulting distribution - and "how many branches does this tree
    actually resolve" is the first question to ask of a MAG phylogeny, where the usual answer is
    "fewer than it has".
    """
    if not os.path.isfile(newick_path):
        return SupportSummary(branches=0, supported=0, threshold=threshold)
    with open(newick_path) as handle:
        newick = handle.read()
    # The ultrafast bootstrap is the second value where both are present: SUPPORTED_BRANCH is
    # calibrated for UFBoot, and reading an SH-aLRT value at 95 is a different and stricter test.
    values = [float(second or first) for first, second in _SUPPORT.findall(newick)]
    return SupportSummary(
        branches=len(values),
        supported=sum(1 for value in values if value >= threshold),
        threshold=threshold,
    )


#: Columns of the per-tree summary.
SUMMARY_COLUMNS: Tuple[str, ...] = (
    "tree",
    "alignment",
    "taxa",
    "alignment_columns",
    "mean_occupancy",
    "taxa_below_occupancy_floor",
    "model",
    "internal_branches",
    "branches_supported",
    "percent_branches_supported",
    "support_threshold",
    "bootstrap_replicates",
)


def summary_row(
    name: str,
    alignment: Alignment,
    model: Optional[str],
    support: SupportSummary,
    replicates: int = BOOTSTRAP_REPLICATES,
) -> Dict[str, str]:
    fraction = support.fraction_supported
    return {
        "tree": name,
        "alignment": os.path.basename(alignment.path),
        "taxa": str(alignment.size),
        "alignment_columns": str(alignment.columns),
        "mean_occupancy": (
            "" if alignment.mean_occupancy is None else "%.3f" % alignment.mean_occupancy
        ),
        # Named rather than counted: with a handful of genomes the identity of the sparse one is
        # what a reader needs, and with many the count is in the length of this field.
        "taxa_below_occupancy_floor": ";".join(alignment.sparse_taxa),
        "model": model or "",
        "internal_branches": str(support.branches),
        "branches_supported": str(support.supported),
        "percent_branches_supported": "" if fraction is None else "%.1f" % (100.0 * fraction),
        "support_threshold": "%g" % support.threshold,
        "bootstrap_replicates": str(replicates),
    }


#: Printed with every tree.
SUPPORT_CAVEAT = (
    "ultrafast bootstrap support says whether the alignment supports a branch, not whether the "
    "branch is true. Three things it cannot see on a MAG phylogeny: markers absent because a "
    "genome is incomplete, so a branch rests on fewer columns than the alignment is wide; "
    "horizontal transfer; and contamination putting another organism's marker in a genome's row. "
    "Read it at %g, which for UFBoot is about as strict as 70 for a standard bootstrap"
    % SUPPORTED_BRANCH
)

#: Printed when the alignment is mostly gaps for some taxa.
OCCUPANCY_CAVEAT = (
    "some genomes fill less than %.0f%% of the alignment, which on a marker concatenation means "
    "their markers were not recovered rather than that they differ. Their placement is the least "
    "reliable part of the tree" % (OCCUPANCY_FLOOR * 100)
)

#: Printed when the tree is of marker genes rather than of one species' core.
SCOPE_CAVEAT = (
    "a tree from universal marker genes places species relative to each other; it cannot resolve "
    "strains of one species. For that, build the tree from one dereplication cluster's core genes "
    "- `metawrap2 pangenome --msa` writes that alignment"
)


def too_few_taxa(alignment: Alignment) -> Optional[str]:
    """Why this alignment cannot give a tree, or None if it can."""
    if alignment.size < MIN_TAXA:
        return (
            "%d taxa in %s: a tree needs at least %d. With three there is exactly one unrooted "
            "topology, so a model and a support value would be decoration"
            % (alignment.size, os.path.basename(alignment.path), MIN_TAXA)
        )
    return None


def find_alignments(root: str) -> List[str]:
    """Alignments under *root* that this module knows how to use.

    Both shapes it accepts, searched for rather than demanded by exact path: GTDB-Tk nests its MSA
    under ``align/`` with a name that carries the domain, and `pangenome` writes one per group under
    ``groups/<cluster>/msa/``.
    """
    found: List[str] = []
    for base, _dirs, files in os.walk(root):
        for name in sorted(files):
            if (
                name.endswith(("user_msa.fasta", "user_msa.fasta.gz"))
                or name.startswith("core_genome_alignment")
                and name.endswith((".aln", ".fasta"))
            ):
                found.append(os.path.join(base, name))
    return found
