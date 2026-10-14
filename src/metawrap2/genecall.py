"""Gene calling, once, over whatever sequence the pipeline has - a bin set or a whole assembly.

Every functional question downstream of this file starts from the same list of genes: what a genome
carries, which KEGG modules its KOs complete, which resistance genes it has and what sits next to
them, and how much of each is in the sample. Calling genes separately for each of those is how
``annotate_bins`` worked, and it costs more than the duplicated compute: two callers that disagree by
one gene produce two answers that cannot be joined.

So this is a library rather than a module. It calls genes with **Pyrodigal** - Cython bindings to
Prodigal 2.6.3, so the calls are Prodigal's calls - and returns them as data, with the two fields a
fragmented assembly makes necessary:

* ``partial``, Prodigal's own flag for a gene whose start or stop codon was not found because the
  contig ended first. A partial gene is a real observation with an unknown boundary, not a bad call.
* ``distance_to_contig_end``, which is what makes a *negative* claim auditable. "No mobile element
  within 10 kb of this resistance gene" means nothing if the gene sits 2 kb from the end of its
  contig; the absence is of sequence, not of mobile elements. Every consumer of these genes that
  reports a neighbourhood has to carry this number, so it is computed here once rather than
  rediscovered per module.

**Single-genome mode or metagenomic mode is chosen per input, not globally.** Prodigal learns its
coding model from the sequence it is given, and that model is sharper when the sequence is one
organism. A bin is one organism by construction, so each bin trains its own model; a whole assembly
is not, so it gets Prodigal's pre-trained metagenomic bins. Prodigal's own floor for training is
20 kb (:data:`pyrodigal.MIN_SINGLE_GENOME`), below which the model is worse than the pre-trained
ones, so a bin under that falls back to meta mode and says so in the record rather than silently
producing a different kind of call than its neighbours.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import pyrodigal

from .io.seqio import contig_id, iter_fasta, smart_open

#: Translation table used to **train** a single-genome coding model: 11, bacterial/archaeal/plant
#: plastid, Prodigal's own default and the right one for everything this pipeline bins. Stated rather
#: than left implicit because a model trained under a different table finds different genes.
#:
#: Proteins are *not* translated with this table unconditionally. Prodigal's metagenomic models carry
#: their own table - some of the pre-trained bins use table 4, where TGA codes for tryptophan rather
#: than stopping - and translating a gene called under table 4 with table 11 truncates it at the
#: first TGA. So each gene is translated with the table it was called under, and that table is
#: recorded per gene. pyrodigal warns about the mismatch, which is how this was found.
TRAINING_TABLE = 11

#: Minimum sequence a bin must have before its own coding model is trained rather than borrowing
#: Prodigal's metagenomic ones. Prodigal's own floor; below it the trained model is the worse option.
MIN_TRAIN_BP = pyrodigal.MIN_SINGLE_GENOME

#: Written into the GFF header so a file says which caller and which mode produced it.
SOURCE = "Pyrodigal_%s" % pyrodigal.PRODIGAL_VERSION


@dataclass(frozen=True)
class Gene:
    """One predicted gene, with the provenance a fragmented assembly needs.

    Coordinates are 1-based and inclusive, matching GFF and Prodigal's own output rather than
    Python's half-open convention, so that a number here can be compared with a number in a GFF
    without a silent off-by-one.
    """

    gene_id: str
    contig: str
    start: int
    end: int
    strand: int
    #: Prodigal's two-character flag: "00" complete, "10" missing a start, "01" missing a stop,
    #: "11" running off both ends of its contig.
    partial: str
    #: Bases between this gene and the nearer end of its contig. 0 for a gene that reaches the end.
    distance_to_contig_end: int
    contig_length: int
    #: The genetic code this gene was called and translated under. Usually 11; Prodigal's
    #: metagenomic models include table-4 bins, and a protein is only comparable to another
    #: translated the same way.
    translation_table: int
    protein: str
    nucleotide: str

    @property
    def is_partial(self) -> bool:
        return self.partial != "00"

    @property
    def length_bp(self) -> int:
        return self.end - self.start + 1


def _partial_flag(gene: pyrodigal.Gene) -> str:
    """Prodigal's partial flag, in Prodigal's own order: begin then end, in genome orientation."""
    return "%d%d" % (int(gene.partial_begin), int(gene.partial_end))


def _genes_of_contig(finder: pyrodigal.GeneFinder, name: str, sequence: str) -> Iterator[Gene]:
    contig_length = len(sequence)
    for index, found in enumerate(finder.find_genes(sequence), start=1):
        # min() over both ends: a gene is only as well contextualised as its nearer edge, and a
        # consumer asking "is there room for a mobile element beside this" cares about the short
        # side. Both ends are recoverable from start/end and contig_length if a caller wants them.
        distance = min(found.begin - 1, contig_length - found.end)
        yield Gene(
            gene_id="%s_%d" % (name, index),
            contig=name,
            start=found.begin,
            end=found.end,
            strand=found.strand,
            partial=_partial_flag(found),
            distance_to_contig_end=distance,
            contig_length=contig_length,
            translation_table=found.translation_table,
            # No translation_table= argument: translate under the code the gene was called with.
            protein=found.translate(),
            nucleotide=found.sequence(),
        )


def _trained_finder(sequences: Sequence[str]) -> Tuple[pyrodigal.GeneFinder, str]:
    """A finder for one organism, trained on it, or the metagenomic one if there is too little.

    Returns the finder and the mode it is in, because which one ran is part of the result: a bin
    annotated in meta mode has systematically different gene boundaries from its neighbours, and a
    reader comparing gene counts across bins needs to know that happened.
    """
    total = sum(len(s) for s in sequences)
    if total < MIN_TRAIN_BP:
        return pyrodigal.GeneFinder(meta=True), "meta"
    finder = pyrodigal.GeneFinder(meta=False)
    # The literal is what pyrodigal's overloads accept; TRAINING_TABLE is the documented name for it
    # and is asserted equal to 11 in the tests, so the two cannot drift apart silently.
    finder.train(*sequences, translation_table=11)
    return finder, "single"


def _covered_bases(genes: Sequence[Gene]) -> int:
    """Bases inside at least one gene, counting each base once.

    Per contig, because gene coordinates restart on every contig and merging across them would
    join intervals that are nowhere near each other.
    """
    by_contig: Dict[str, List[Tuple[int, int]]] = {}
    for gene in genes:
        by_contig.setdefault(gene.contig, []).append((gene.start, gene.end))
    covered = 0
    for spans in by_contig.values():
        spans.sort()
        current_start, current_end = spans[0]
        for start, end in spans[1:]:
            if start <= current_end + 1:
                current_end = max(current_end, end)
            else:
                covered += current_end - current_start + 1
                current_start, current_end = start, end
        covered += current_end - current_start + 1
    return covered


@dataclass(frozen=True)
class GenomeGenes:
    """The genes of one genome (or one assembly), with how they were called."""

    name: str
    genes: List[Gene]
    #: "single" if a coding model was trained on this sequence, "meta" if Prodigal's pre-trained
    #: metagenomic bins were used instead.
    mode: str
    n_contigs: int
    total_bp: int

    @property
    def n_partial(self) -> int:
        return sum(1 for g in self.genes if g.is_partial)

    @property
    def coding_density(self) -> float:
        """Fraction of the genome inside a predicted gene, as a sanity signal.

        A bacterial genome runs about 0.85-0.90. A number far below that usually means the sequence
        is not what the caller assumed - a viral contig with overlapping reading frames, or a bin
        made mostly of very short fragments.

        Computed over the **union** of gene intervals rather than the sum of their lengths. Genes
        overlap - Prodigal allows up to ``max_overlap`` bases between neighbours - so summing
        lengths double-counts the overlap and can exceed the genome's length. A density above 1.0 is
        how that was found here.
        """
        if not self.total_bp:
            return 0.0
        return _covered_bases(self.genes) / float(self.total_bp)


def call_genome(path: str, name: Optional[str] = None, meta: bool = False) -> GenomeGenes:
    """Call genes in one FASTA.

    *meta* forces metagenomic mode; otherwise a coding model is trained on the file's own sequence
    when there is enough of it. Pass ``meta=True`` for a whole assembly, where "the file's own
    sequence" is many organisms and training on it trains on an average of them.
    """
    records = [(contig_id(header), seq) for header, seq in iter_fasta(path)]
    sequences = [seq for _name, seq in records]
    if meta:
        finder, mode = pyrodigal.GeneFinder(meta=True), "meta"
    else:
        finder, mode = _trained_finder(sequences)

    genes: List[Gene] = []
    for contig, sequence in records:
        genes.extend(_genes_of_contig(finder, contig, sequence))
    return GenomeGenes(
        name=name if name is not None else os.path.splitext(os.path.basename(path))[0],
        genes=genes,
        mode=mode,
        n_contigs=len(records),
        total_bp=sum(len(s) for s in sequences),
    )


# ── writing ──────────────────────────────────────────────────────────────────────────────────


def write_faa(called: GenomeGenes, path: str) -> int:
    """Translated genes. The header carries the fields a downstream join needs, not just an id.

    ``>gene_id # start # end # strand # partial=00;dist_to_end=1234`` keeps Prodigal's own
    ``#``-delimited shape - which eggNOG-mapper, PyHMMER and every other consumer treat as an
    opaque description - while adding the two fields this pipeline added to Gene. Anything parsing
    only the first whitespace-delimited token is unaffected.
    """
    with smart_open(path, "wt") as out:
        for gene in called.genes:
            out.write(
                ">%s # %d # %d # %d # partial=%s;dist_to_end=%d\n"
                % (
                    gene.gene_id,
                    gene.start,
                    gene.end,
                    gene.strand,
                    gene.partial,
                    gene.distance_to_contig_end,
                )
            )
            out.write(gene.protein + "\n")
    return len(called.genes)


def write_ffn(called: GenomeGenes, path: str) -> int:
    """Untranslated genes, with the same headers as the .faa so the two files join on gene_id."""
    with smart_open(path, "wt") as out:
        for gene in called.genes:
            out.write(
                ">%s # %d # %d # %d # partial=%s;dist_to_end=%d\n"
                % (
                    gene.gene_id,
                    gene.start,
                    gene.end,
                    gene.strand,
                    gene.partial,
                    gene.distance_to_contig_end,
                )
            )
            out.write(gene.nucleotide + "\n")
    return len(called.genes)


def write_gff(called: GenomeGenes, path: str) -> int:
    """GFF3 of the calls, with partial and distance as attributes rather than lost in a comment."""
    with smart_open(path, "wt") as out:
        out.write("##gff-version 3\n")
        out.write("# source: %s (%s mode)\n" % (SOURCE, called.mode))
        for gene in called.genes:
            out.write(
                "\t".join(
                    (
                        gene.contig,
                        SOURCE,
                        "CDS",
                        str(gene.start),
                        str(gene.end),
                        ".",
                        "+" if gene.strand > 0 else "-",
                        "0",
                        "ID=%s;partial=%s;dist_to_contig_end=%d"
                        % (gene.gene_id, gene.partial, gene.distance_to_contig_end),
                    )
                )
                + "\n"
            )
    return len(called.genes)


#: Columns of the per-gene table, in order. One row per gene per genome, long rather than wide,
#: because the things joined onto it later (KOs, Pfams, ARG calls) are many-to-one per gene.
GENE_COLUMNS: Tuple[str, ...] = (
    "genome",
    "gene_id",
    "contig",
    "start",
    "end",
    "strand",
    "length_bp",
    "partial",
    "distance_to_contig_end",
    "contig_length",
    "translation_table",
    "calling_mode",
)


def gene_rows(called: GenomeGenes) -> Iterator[Dict[str, object]]:
    for gene in called.genes:
        yield {
            "genome": called.name,
            "gene_id": gene.gene_id,
            "contig": gene.contig,
            "start": gene.start,
            "end": gene.end,
            "strand": "+" if gene.strand > 0 else "-",
            "length_bp": gene.length_bp,
            "partial": gene.partial,
            "distance_to_contig_end": gene.distance_to_contig_end,
            "contig_length": gene.contig_length,
            "translation_table": gene.translation_table,
            "calling_mode": called.mode,
        }
