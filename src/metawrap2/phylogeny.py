"""Trees from alignments the pipeline already produces, with support and occupancy summaries.

Two alignments are used:

* GTDB-Tk's concatenated marker alignment (``*.user_msa.fasta``; 120 bacterial or 53 archaeal
  single-copy markers), for placing a catalogue's genomes relative to each other and to GTDB's
  reference set.
* PPanGGOLiN's core-gene alignment for one species group (written with ``--phylo``), for
  relationships between strains of one species.

A marker tree cannot resolve strains, and a core-gene tree does not place the species.

IQ-TREE is the default. ModelFinder selects the substitution model, and ultrafast bootstrap support
is read at 95% (FastTree gives SH-like values around 0.8). FastTree remains available because at a
few hundred genomes IQ-TREE's model selection dominates the run time.

Ultrafast bootstrap support measures whether the alignment supports a branch. On a MAG marker
concatenation it does not account for markers missing from incomplete MAGs, horizontal transfer, or
contamination. :func:`support_summary` is therefore reported with alignment occupancy, and
:data:`SUPPORT_CAVEAT` is printed with every tree.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

#: Ultrafast bootstrap support at or above which a branch is called supported, per IQ-TREE's
#: guidance. UFBoot 95 is about as strict as 70 for a standard bootstrap.
SUPPORTED_BRANCH = 95.0

#: Replicates for the ultrafast bootstrap. IQ-TREE's minimum is 1000 and it refuses fewer.
BOOTSTRAP_REPLICATES = 1000

#: Minimum number of taxa. Three taxa have only one unrooted topology.
MIN_TAXA = 4

#: Per-taxon fraction of non-gap characters below which a taxon is reported as sparse. Reported,
#: not enforced; trimming is a separate step. Marker concatenations from incomplete MAGs can be
#: mostly gaps for some taxa.
OCCUPANCY_FLOOR = 0.5


@dataclass(frozen=True)
class Alignment:
    """Taxa, width, occupancy and residue type of an alignment."""

    path: str
    taxa: Tuple[str, ...]
    #: Alignment width in columns.
    columns: int
    #: Per-taxon fraction of non-gap characters.
    occupancy: Dict[str, float]
    #: ``"nucleotide"`` or ``"protein"``, measured from the residues; it selects the substitution
    #: model. FastTree with ``-nt`` on a protein alignment treats most residues as ambiguous.
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
        """Taxa whose occupancy is below OCCUPANCY_FLOOR."""
        return tuple(
            sorted(name for name, filled in self.occupancy.items() if filled < OCCUPANCY_FLOOR)
        )

    @property
    def is_a_tree(self) -> bool:
        return self.size >= MIN_TAXA


#: Residues counted as nucleotide: the four bases, U, N and the IUPAC ambiguity codes. Several are
#: also amino-acid codes, so the test uses a fraction.
NUCLEOTIDE_CHARS = frozenset("ACGTUNRYSWKMBDHV")

#: Fraction of nucleotide-like residues needed to call an alignment nucleotide. Protein
#: alignments reach about 20-30% (A, C, G, T, N, W, Y, K, M, H, R, S, V, D are amino-acid codes
#: too) and nucleotide alignments exceed 99%. Values in between are treated as protein, which is
#: what the upstream aligners write.
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

    Reads gzip, since GTDB-Tk writes its MSA compressed.
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
#: Newick support values after a closing bracket, in both forms IQ-TREE writes: `)97.5/100:` with
#: --alrt (SH-aLRT, then ultrafast bootstrap) and `)100:` with only -B. A single-number pattern
#: would count no branches in the combined form.
_SUPPORT = re.compile(r"\)(\d+(?:\.\d+)?)(?:/(\d+(?:\.\d+)?))?[:,;)]")


def model_from_report(path: str) -> Optional[str]:
    """The substitution model IQ-TREE's ModelFinder chose, from its report file.

    Trees built under different models are not comparable, even from the same alignment.
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
    """Number of internal branches and how many reach the support threshold."""

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

    Parsed from the tree because IQ-TREE's log reports the bootstrap settings, not the resulting
    values.
    """
    if not os.path.isfile(newick_path):
        return SupportSummary(branches=0, supported=0, threshold=threshold)
    with open(newick_path) as handle:
        newick = handle.read()
    # Where both values are present, use the second (UFBoot); SUPPORTED_BRANCH is calibrated for
    # UFBoot, not SH-aLRT.
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
        # Listed by name so the sparse genomes can be identified.
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
    "ultrafast bootstrap support measures how well the alignment supports a branch, not whether "
    "the branch is correct. On a MAG phylogeny it does not account for markers missing from "
    "incomplete genomes (fewer columns behind a branch), horizontal transfer, or contamination "
    "that puts another organism's marker in a genome's row. Use a cutoff of %g; for UFBoot this "
    "is about as strict as 70 for a standard bootstrap" % SUPPORTED_BRANCH
)

#: Printed when the alignment is mostly gaps for some taxa.
OCCUPANCY_CAVEAT = (
    "some genomes fill less than %.0f%% of the alignment. In a marker concatenation this means "
    "their markers were not recovered, so their placement in the tree is the least "
    "reliable" % (OCCUPANCY_FLOOR * 100)
)

#: Printed when the tree is of marker genes rather than of one species' core.
SCOPE_CAVEAT = (
    "a tree from universal marker genes places species relative to each other but cannot resolve "
    "strains of one species. For strains, build the tree from one dereplication cluster's core "
    "genes: `metawrap2 pangenome --msa` writes that alignment"
)


def too_few_taxa(alignment: Alignment) -> Optional[str]:
    """Why this alignment cannot give a tree, or None if it can."""
    if alignment.size < MIN_TAXA:
        return (
            "%d taxa in %s: a tree needs at least %d (three taxa have only one possible "
            "unrooted topology)" % (alignment.size, os.path.basename(alignment.path), MIN_TAXA)
        )
    return None


def find_alignments(root: str) -> List[str]:
    """Alignments under *root* that this module can use.

    Found by name anywhere under *root*: GTDB-Tk writes its MSA under ``align/`` with the domain in
    the name, and `pangenome` writes one per group under ``groups/<cluster>/msa/``.
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
