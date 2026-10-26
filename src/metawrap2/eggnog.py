"""eggNOG-mapper output parsing, for the COG categories, GO terms and EC numbers that KOfam and
Pfam do not produce.

In 2.1.2, of 3,891 genes in one bin, 1,578 had a Pfam domain and no KO, and 13 had a KO and no
Pfam. Pfam covers nearly all of KEGG's genes; KEGG adds reaction-level calls. Neither gives GO
terms, COG functional categories or EC numbers.

eggNOG annotates by orthology: a gene is placed in an orthologous group, and GO, EC, COG and KEGG
annotations are transferred from characterised members of the group. This is a different inference
from a profile match, with different failure modes, and is reported beside KOfam and Pfam.

eggNOG's KEGG orthologs are kept separate from KOfam's. eggNOG transfers a KO from an orthologous
group; KOfam calls one from a curated per-profile bit-score threshold. :func:`compare_kos` reports
genes called by both, by eggNOG only and by KOfam only, and neither source is preferred.

Nothing is downloaded here. The eggNOG data (~50 GB: ``eggnog.db`` and ``eggnog_proteins.dmnd``)
is given on the command line or in ``metawrap2.toml``.

Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J. eggNOG-mapper v2:
functional annotation, orthology assignments, and domain prediction at the metagenomic scale.
*Mol Biol Evol* 38(12):5825-5829 (2021). doi:10.1093/molbev/msab293
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

#: eggNOG-mapper writes ``-`` in any field with no value, in any column. Unhandled, it would be
#: counted as a GO term, EC number or COG category. All values are read through :func:`_split` and
#: :func:`_scalar`, which drop it. (dbCAN output had the same problem with ``-|-``.)
EMPTY = "-"

#: The columns of ``emapper.annotations``, in order, taken from eggNOG-mapper 2.1.13's writer
#: (``eggnogmapper/annotation/output.py``: ``HIT_HEADER + ANNOTATIONS_HEADER``), since the
#: documentation lags the code.
ANNOTATION_COLUMNS: Tuple[str, ...] = (
    "query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
)

#: Extra column appended by ``--md5``; accepted but not required.
MD5_COLUMN = "md5"

#: emapper prefixes the header line with ``#`` (``#query<TAB>...``), and writes ``## `` provenance
#: lines above it and ``## N queries scanned`` footer lines after the last row. Comment lines are
#: skipped at both ends so footers are not read as data rows.
COMMENT = "#"

#: Comma-joined fields (emapper writes them sorted). They are split and deduplicated, keeping
#: emapper's order so re-runs give byte-identical tables.
MULTI_VALUED: Tuple[str, ...] = (
    "eggNOG_OGs",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "BRITE",
    "PFAMs",
)

#: ``KEGG_ko`` values carry a namespace (``ko:K00059``). The package's KO tables
#: (``annotate --kofam``'s ``kos.tsv``, ``pathways``' BRITE hierarchy) use bare K numbers, so the
#: prefix is stripped; otherwise the join with KOfam returns no rows.
KO_PREFIX = "ko:"

#: Not stripped, unlike ``KO_PREFIX``: ``GO:`` is part of the accession (``GO:0003674``).
#: ``KEGG_Pathway`` values such as ``ko00010`` have no colon and are left as is.
GO_PREFIX = "GO:"


def _scalar(value: Optional[str]) -> str:
    """One field as text, with eggNOG's absent-value sentinel resolved to empty."""
    text = (value or "").strip()
    return "" if text == EMPTY else text


def _split(value: Optional[str]) -> Tuple[str, ...]:
    """One comma-joined field as its members, sentinels and empties dropped, order preserved.

    emapper's order is kept so the output depends only on the input.
    """
    text = _scalar(value)
    if not text:
        return ()
    members = [piece.strip() for piece in text.split(",")]
    return tuple(dict.fromkeys(piece for piece in members if piece and piece != EMPTY))


def _number(value: Optional[str]) -> Optional[float]:
    """An e-value or bit score, or None.

    None for an unparseable field: emapper's ``-m no_search`` and cache modes can leave the hit
    columns as ``-``, and 0.0 would read as a perfect hit.
    """
    text = _scalar(value)
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


#: COG categories are concatenated, not comma-separated: a gene in three categories has ``EGP``.
#: The column is therefore not in MULTI_VALUED and is split per character.
COG_CATEGORY_WIDTH = 1


@dataclass(frozen=True)
class EggnogRow:
    """One gene's eggNOG annotation, as one line of ``emapper.annotations`` describes it."""

    genome: str
    gene: str
    #: The best-matching eggNOG protein, ``taxid.protein``; the row's annotations were transferred
    #: from its orthologous groups.
    seed_ortholog: str
    evalue: Optional[float]
    score: Optional[float]
    #: Orthologous groups, finest first, each ``OG@taxid|taxon``.
    og_names: Tuple[str, ...]
    #: The taxonomic level the annotations were transferred at; a genus-level transfer is
    #: stronger evidence than a phylum-level one.
    max_annot_lvl: str
    #: Concatenated one-letter categories, as written; see COG_CATEGORY_WIDTH.
    cog_category: str
    description: str
    preferred_name: str
    #: GO terms, with the ``GO:`` prefix kept as part of the accession.
    gos: Tuple[str, ...]
    ecs: Tuple[str, ...]
    #: Bare K numbers with ``ko:`` stripped; see KO_PREFIX.
    kegg_kos: Tuple[str, ...]
    kegg_pathways: Tuple[str, ...]
    brites: Tuple[str, ...]
    pfams: Tuple[str, ...]

    @property
    def cog_categories(self) -> Tuple[str, ...]:
        """``"EGP"`` -> ``("E", "G", "P")``; ``""`` -> ``()``.

        A gene with no COG category yields an empty tuple.
        """
        return tuple(
            self.cog_category[i : i + COG_CATEGORY_WIDTH]
            for i in range(0, len(self.cog_category), COG_CATEGORY_WIDTH)
        )

    @property
    def has_orthogonal_annotation(self) -> bool:
        """Whether this row has a GO term, EC number or COG category.

        These are the annotations KOfam and Pfam do not produce. A row with only a seed ortholog
        and description does not count.
        """
        return bool(self.gos or self.ecs or self.cog_category)


class AnnotationFormatError(ValueError):
    """``emapper.annotations`` did not have the expected columns.

    A separate class so callers can tell a schema change from an IO error. Raised instead of
    returning no rows, which would look like a genome with no annotations.
    """


def _check_header(line: str, path: str) -> Tuple[str, ...]:
    """Validate emapper's header line and return the column names it declares.

    On a mismatch (an eggNOG-mapper version whose columns differ from 2.1.13's), the error lists
    the missing and unexpected columns.
    """
    found = tuple(line.lstrip(COMMENT).rstrip("\n").split("\t"))
    expected = ANNOTATION_COLUMNS
    if found == expected or found == expected + (MD5_COLUMN,):
        return found
    missing = [name for name in expected if name not in found]
    extra = [name for name in found if name not in expected + (MD5_COLUMN,)]
    raise AnnotationFormatError(
        "%s does not look like eggNOG-mapper %d-column output: found %d column(s)%s%s"
        % (
            path,
            len(expected),
            len(found),
            "; missing %s" % ", ".join(missing) if missing else "",
            "; unexpected %s" % ", ".join(extra) if extra else "",
        )
    )


def parse_annotations(path: str, genome: str) -> List[EggnogRow]:
    """Read one ``emapper.annotations`` file.

    Only genes emapper placed appear in it, so the annotated share of a genome must be computed
    against the gene table, as in the Pfam and CAZyme passes.
    """
    rows: List[EggnogRow] = []
    columns: Optional[Tuple[str, ...]] = None
    with open(path) as handle:
        for line in handle:
            if line.startswith(COMMENT):
                # The header starts with "#"; other "#" lines are provenance or footer counts.
                if columns is None and line.startswith(COMMENT + ANNOTATION_COLUMNS[0]):
                    columns = _check_header(line, path)
                continue
            if columns is None:
                raise AnnotationFormatError(
                    "%s has data rows before any %s%s header line"
                    % (path, COMMENT, ANNOTATION_COLUMNS[0])
                )
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(ANNOTATION_COLUMNS):
                continue
            cell = dict(zip(columns, fields))
            rows.append(
                EggnogRow(
                    genome=genome,
                    gene=_scalar(cell.get("query")),
                    seed_ortholog=_scalar(cell.get("seed_ortholog")),
                    evalue=_number(cell.get("evalue")),
                    score=_number(cell.get("score")),
                    og_names=_split(cell.get("eggNOG_OGs")),
                    max_annot_lvl=_scalar(cell.get("max_annot_lvl")),
                    cog_category=_scalar(cell.get("COG_category")),
                    description=_scalar(cell.get("Description")),
                    preferred_name=_scalar(cell.get("Preferred_name")),
                    gos=_split(cell.get("GOs")),
                    ecs=_split(cell.get("EC")),
                    kegg_kos=strip_ko_prefixes(_split(cell.get("KEGG_ko"))),
                    kegg_pathways=_split(cell.get("KEGG_Pathway")),
                    brites=_split(cell.get("BRITE")),
                    pfams=_split(cell.get("PFAMs")),
                )
            )
    if columns is None:
        raise AnnotationFormatError(
            "%s has no %s%s header line; it is not an eggNOG-mapper annotations file"
            % (path, COMMENT, ANNOTATION_COLUMNS[0])
        )
    return rows


def strip_ko_prefixes(kos: Iterable[str]) -> Tuple[str, ...]:
    """``("ko:K00059",)`` -> ``("K00059",)``, leaving anything already bare alone.

    Public because joins with ``kos.tsv`` need bare K numbers; with the prefix, every comparison
    would come back empty.
    """
    out: List[str] = []
    for ko in kos:
        bare = ko.removeprefix(KO_PREFIX)
        if bare and bare != EMPTY and bare not in out:
            out.append(bare)
    return tuple(out)


# ── this pipeline's own tables ────────────────────────────────────────────────────────────────

#: Columns of the per-gene table. COG, GO and EC come first, before the hit statistics.
GENE_COLUMNS: Tuple[str, ...] = (
    "genome",
    "gene",
    "cog_category",
    "go_terms",
    "ec_numbers",
    "kegg_ko",
    "description",
    "preferred_name",
    "orthologous_groups",
    "max_annot_lvl",
    "seed_ortholog",
    "evalue",
    "score",
    "kegg_pathway",
    "brite",
    "pfams",
    "eggnog_db_version",
)


def gene_row(row: EggnogRow, db_version: str = "unknown") -> Dict[str, str]:
    """One row of the per-gene table.

    Lists are ``;``-separated, as in the CAZyme and Pfam tables.
    """
    return {
        "genome": row.genome,
        "gene": row.gene,
        "cog_category": row.cog_category,
        "go_terms": ";".join(row.gos),
        "ec_numbers": ";".join(row.ecs),
        # Bare K numbers, so this column joins kos.tsv directly; see KO_PREFIX.
        "kegg_ko": ";".join(row.kegg_kos),
        "description": row.description,
        "preferred_name": row.preferred_name,
        "orthologous_groups": ";".join(row.og_names),
        "max_annot_lvl": row.max_annot_lvl,
        "seed_ortholog": row.seed_ortholog,
        "evalue": "" if row.evalue is None else "%.2g" % row.evalue,
        "score": "" if row.score is None else "%.1f" % row.score,
        "kegg_pathway": ";".join(row.kegg_pathways),
        "brite": ";".join(row.brites),
        "pfams": ";".join(row.pfams),
        # Recorded per row (as kofam.py records its threshold list) because eggNOG renumbers
        # groups between releases.
        "eggnog_db_version": db_version,
    }


#: The namespaces in the long table, in row order.
TERM_NAMESPACES: Tuple[str, ...] = ("COG", "GO", "EC", "KEGG_ko")

#: Columns of the long per-term table: one row per (genome, gene, namespace, term), for counting
#: with a group-by. GENE_COLUMNS is the wide table for reading.
TERM_COLUMNS: Tuple[str, ...] = (
    "genome",
    "gene",
    "namespace",
    "term",
    "eggnog_db_version",
)


def _terms(row: EggnogRow, namespace: str) -> Tuple[str, ...]:
    if namespace == "COG":
        return row.cog_categories
    if namespace == "GO":
        return row.gos
    if namespace == "EC":
        return row.ecs
    return row.kegg_kos


def term_rows(rows: Sequence[EggnogRow], db_version: str = "unknown") -> List[Dict[str, str]]:
    """The long table: one row per term a gene carries, no row for a gene that carries none.

    Rows follow input order, then TERM_NAMESPACES order, so output is byte-identical across runs.
    """
    out: List[Dict[str, str]] = []
    for row in rows:
        for namespace in TERM_NAMESPACES:
            for term in _terms(row, namespace):
                out.append(
                    {
                        "genome": row.genome,
                        "gene": row.gene,
                        "namespace": namespace,
                        "term": term,
                        "eggnog_db_version": db_version,
                    }
                )
    return out


# ── what eggNOG added, and where it disagrees ─────────────────────────────────────────────────


@dataclass(frozen=True)
class Coverage:
    """How much of a gene set eggNOG annotated, and how much with GO, EC or COG terms.

    ``genes_with_orthogonal_annotation`` counts genes with GO, EC or COG terms, which KOfam and
    Pfam do not provide. Total coverage adds little, since Pfam already covers nearly every gene
    KEGG covers (1,578 Pfam-only genes against 13 KEGG-only in the bin measured in 2.1.2).
    """

    genes: int
    genes_annotated: int
    genes_with_orthogonal_annotation: int
    go_terms: int
    ec_numbers: int
    cog_categories: int

    @property
    def fraction_annotated(self) -> Optional[float]:
        return self.genes_annotated / self.genes if self.genes else None


def coverage(rows: Sequence[EggnogRow], genes: Sequence[str]) -> Coverage:
    """Count what eggNOG placed over a known gene set.

    *genes* is every called gene, from the gene table, since emapper omits genes it could not
    place.
    """
    gos: Set[str] = set()
    ecs: Set[str] = set()
    cogs: Set[str] = set()
    for row in rows:
        gos.update(row.gos)
        ecs.update(row.ecs)
        cogs.update(row.cog_categories)
    return Coverage(
        genes=len(set(genes)),
        genes_annotated=len({row.gene for row in rows}),
        genes_with_orthogonal_annotation=len(
            {row.gene for row in rows if row.has_orthogonal_annotation}
        ),
        go_terms=len(gos),
        ec_numbers=len(ecs),
        cog_categories=len(cogs),
    )


@dataclass(frozen=True)
class KoOverlap:
    """Per-gene agreement between eggNOG's and KOfam's KEGG orthologs.

    The two use different inferences (orthologous-group transfer versus a curated bit-score
    threshold), so they are counted, not merged. Both tables are kept on disk.
    """

    #: Genes both annotators gave at least one KO to, whether or not it is the same KO.
    both: int
    eggnog_only: int
    kofam_only: int
    #: Of `both`, the genes where the two KO sets intersect. Genes in `both` but not here have
    #: different KOs from the two sources and are worth checking by hand.
    agreeing: int

    @property
    def total(self) -> int:
        return self.both + self.eggnog_only + self.kofam_only

    def describe(self) -> str:
        """One console line with the counts, preferring neither source."""
        return (
            "KEGG orthologs: %d gene(s) called by both eggNOG and KOfam (%d on the same KO), "
            "%d by eggNOG only, %d by KOfam only; both tables are kept and neither is preferred"
            % (self.both, self.agreeing, self.eggnog_only, self.kofam_only)
        )


def compare_kos(rows: Sequence[EggnogRow], kofam_kos: Mapping[str, Iterable[str]]) -> KoOverlap:
    """Compare this genome's eggNOG KOs against KOfam's, per gene.

    *kofam_kos* is gene -> the K numbers ``--kofam`` assigned it. Both sides are reduced to bare K
    numbers with :func:`strip_ko_prefixes`, since eggNOG writes ``ko:K00059`` and ``kos.tsv`` holds
    ``K00059``.
    """
    mine: Dict[str, Set[str]] = {}
    for row in rows:
        if row.kegg_kos:
            mine.setdefault(row.gene, set()).update(row.kegg_kos)
    theirs = {
        gene: set(strip_ko_prefixes(kos))
        for gene, kos in kofam_kos.items()
        if set(strip_ko_prefixes(kos))
    }
    shared = set(mine) & set(theirs)
    return KoOverlap(
        both=len(shared),
        eggnog_only=len(set(mine) - set(theirs)),
        kofam_only=len(set(theirs) - set(mine)),
        agreeing=sum(1 for gene in shared if mine[gene] & theirs[gene]),
    )


#: Printed with every eggNOG result.
TRANSFER_CAVEAT = (
    "eggNOG annotations are transferred from an orthologous group: a gene gets the GO, EC and COG "
    "terms of the characterised members of its group. max_annot_lvl gives the taxonomic level the "
    "transfer came from; a phylum-level transfer is weaker evidence than a genus-level one"
)

#: Printed when eggNOG is reported beside KOfam and Pfam.
ORTHOGONAL_CAVEAT = (
    "the useful count is the genes eggNOG gave a GO term, EC number or COG category, since KOfam "
    "and Pfam do not produce these. Pfam already covers nearly every gene KEGG covers: on one "
    "test bin, 1578 genes had a Pfam domain and no KO, and 13 had a KO and no Pfam domain"
)
