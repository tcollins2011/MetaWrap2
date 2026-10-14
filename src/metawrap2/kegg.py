"""KEGG module definitions, and how complete a genome's KOs make each one.

"Complete" is not one definition. Five published tools implement five mutually incompatible
algorithms over the same boolean expressions, and they disagree structurally rather than at the
margin - on the same KO set, anvi'o reports 11-12% more modules at >=80% completeness than
MicrobeAnnotator or KofamScan do. A user comparing a number from here against a number in a paper is
comparing algorithms unless both are named.

So this file computes completeness rather than wrapping a tool that does, for the same reason
:mod:`metawrap2.depth` computes the depth table rather than shelling out to
``jgi_summarize_bam_contig_depths``: it is a load-bearing number, and a wrapper cannot be looked
inside. Two scores come out of every module, always:

* **pathwise** - unroll every alternative path through the definition, score each as the fraction of
  its steps whose KOs are present, and take the best. This is anvi'o's primary score. It is
  optimistic by construction, because the best path is chosen after seeing the data.
* **stepwise** - take the definition's top-level steps, call each one complete if any of its
  alternatives is complete, and report the fraction of steps that are. This is the stricter reading
  and the one that answers "is the pathway there", where pathwise answers "how much of the best
  route is there".

Reporting one of them hides a disagreement that is as large as the biological signal, so both are
emitted and the caller is not given a way to ask for "the" completeness.

Definition grammar, loosest to tightest binding:

====================  ==========================================================================
``K00844 K01810``     a space is AND: sequential steps, both required
``K00844,K12407``     a comma is OR: alternatives for one step
``K00844+K12407``     a plus is AND *within* a step: subunits of a complex, all essential
``K00844-K12407``     a minus marks a non-essential component, which does not count toward the step
``(...)``             grouping
``--``                a reaction with no KO assigned to it
====================  ==========================================================================

``--`` is the one that needs a decision rather than a parse. A step with no KO cannot be satisfied by
any genome, so a module containing one can never reach 100% under pathwise scoring. anvi'o scores it
0 and lets the module cap below 1; EBI's tool substitutes it at weight 0, so it does not cap the
module at all. Neither is wrong - they answer different questions - and this file takes anvi'o's
reading, because a module reported as complete when one of its reactions has no known gene is a
claim the data cannot support. Which reading ran is recorded in the output, and
:data:`ModuleScore.has_no_ko_step` marks every module where the choice could matter.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

#: The vendored snapshot. See data/kegg/SNAPSHOT.md for its provenance and how to move it.
MODULES_TABLE = os.path.join(os.path.dirname(__file__), "data", "kegg", "modules_table.tsv")

#: The date of that snapshot, emitted on every scored row. Hard-coded rather than read from the
#: file's mtime, which a checkout resets.
DEFINITION_SNAPSHOT = "2026-10-12"

#: Stand-in for a reaction with no KO (``--`` in a definition). Not a real K number, so no genome
#: can carry it, which is exactly the intended behaviour - see the module docstring.
NO_KO = "--"

#: Path enumeration is capped. Only a handful of modules have many alternative routes (M00001 has
#: the most at ~22k), so this is not normally reached; it exists so that a future definition with a
#: pathological expression degrades to a reported partial answer instead of hanging. When it is hit,
#: ``paths_truncated`` says so rather than the score quietly meaning something else.
MAX_PATHS = 200_000

_TOKEN = re.compile(r"K\d{5}|M\d{5}|--|[(),+\- ]")


# ── the expression tree ──────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Ko:
    """One KO, or the ``--`` sentinel, or a reference to another module."""

    name: str

    @property
    def is_module_reference(self) -> bool:
        return self.name.startswith("M")


@dataclass(frozen=True)
class Seq:
    """Sequential steps, space-separated. Each child is a step of the pathway in its own right."""

    children: Tuple[Node, ...]


@dataclass(frozen=True)
class All:
    """Subunits of one complex, ``+``-joined. All essential, but together they are **one** step.

    Distinct from :class:`Seq` only so that :func:`top_level_steps` does not split a complex into its
    subunits. It matters: M00142 is NADH dehydrogenase, seven subunits joined by ``+``. Read as seven
    sequential steps, a genome with three of them scores 3/7 complete; read correctly as one step
    requiring all seven, it scores 0 - which is the biologically meaningful answer, because a
    partially assembled complex is not a working one. 31 of the 573 modules are a single complex.
    """

    children: Tuple[Node, ...]


@dataclass(frozen=True)
class Or:
    """Any one child suffices: comma-separated alternatives."""

    children: Tuple[Node, ...]


@dataclass(frozen=True)
class Optional_:
    """A non-essential component (``-K12345``). Never required, never counted."""

    child: Node


Node = object  # Ko | Seq | All | Or | Optional_, named for readability in signatures


class DefinitionError(ValueError):
    """A definition that does not parse. Raised rather than guessed at."""


def _tokenize(definition: str) -> List[str]:
    out: List[str] = []
    position = 0
    for found in _TOKEN.finditer(definition):
        if found.start() != position:
            raise DefinitionError(
                "unparsable text %r at %d in %r"
                % (definition[position : found.start()], position, definition)
            )
        position = found.end()
        token = found.group(0)
        if token != " ":
            out.append(token)
        else:
            out.append(" ")
    if position != len(definition):
        raise DefinitionError("trailing %r in %r" % (definition[position:], definition))
    return out


class _Parser:
    """Recursive descent over the grammar in the module docstring.

    Precedence, loosest first: space (AND), comma (OR), plus/minus (complex). Space being looser
    than comma is what makes ``K1 K2,K3`` mean "K1, then either K2 or K3" rather than "either K1 and
    K2, or K3" - the reading every published parser uses, and the one the definitions are written
    for.
    """

    def __init__(self, tokens: Sequence[str]) -> None:
        self.tokens = list(tokens)
        self.i = 0

    def peek(self) -> Optional[str]:
        return self.tokens[self.i] if self.i < len(self.tokens) else None

    def take(self) -> str:
        token = self.tokens[self.i]
        self.i += 1
        return token

    def skip_spaces(self) -> None:
        while self.peek() == " ":
            self.take()

    def parse(self) -> Node:
        node = self.sequence()
        self.skip_spaces()
        if self.peek() is not None:
            raise DefinitionError("unexpected %r" % self.peek())
        return node

    def sequence(self) -> Node:
        """Space-separated steps: the loosest binding."""
        parts = [self.alternatives()]
        while True:
            # A space only separates steps if something follows it; trailing spaces and the space
            # before a closing paren are not operators.
            save = self.i
            if self.peek() != " ":
                break
            self.skip_spaces()
            if self.peek() in (None, ")", ","):
                self.i = save
                break
            parts.append(self.alternatives())
        return parts[0] if len(parts) == 1 else Seq(tuple(parts))

    def alternatives(self) -> Node:
        parts = [self.complex_()]
        while self.peek() == ",":
            self.take()
            parts.append(self.complex_())
        return parts[0] if len(parts) == 1 else Or(tuple(parts))

    def complex_(self) -> Node:
        """``+`` joins essential subunits; ``-`` marks a non-essential one.

        A ``-`` also appears at the *start* of a step rather than attached to a preceding KO -
        ``K11023 -K11024 K11025`` - where it marks the whole step non-essential rather than one
        subunit of a complex. Seven of the 573 modules are written that way, and reading the leading
        minus as a binary operator rejects all of them.
        """
        if self.peek() == "-":
            self.take()
            return Optional_(self.atom())
        parts: List[Node] = [self.atom()]
        while self.peek() in ("+", "-"):
            operator = self.take()
            child = self.atom()
            parts.append(child if operator == "+" else Optional_(child))
        return parts[0] if len(parts) == 1 else All(tuple(parts))

    def atom(self) -> Node:
        token = self.peek()
        if token == "(":
            self.take()
            self.skip_spaces()
            inner = self.sequence()
            self.skip_spaces()
            if self.peek() != ")":
                raise DefinitionError("expected ) but found %r" % self.peek())
            self.take()
            return inner
        if token is None or token in (")", ",", "+", "-", " "):
            raise DefinitionError("expected a KO but found %r" % token)
        return Ko(self.take())


def parse_definition(definition: str) -> Node:
    """Parse a KEGG module definition into an expression tree."""
    return _Parser(_tokenize(definition.strip())).parse()


# ── scoring ──────────────────────────────────────────────────────────────────────────────────


def top_level_steps(node: Node) -> List[Node]:
    """The definition's top-level steps.

    This is what stepwise completeness is a fraction of, and the count is reported per module
    because it decides how much a completeness figure can mean: 8% of modules have exactly one
    top-level step, which makes their completeness binary and a single false-positive KO call enough
    to report the whole module present.
    """
    if isinstance(node, Seq):
        return list(node.children)
    return [node]


def _required_sets(node: Node) -> List[FrozenSet[str]]:
    """Every alternative route through *node*, each as the set of KOs it requires.

    An AND multiplies its children's routes together; an OR concatenates them; an Optional
    contributes nothing, which is what "non-essential" means. Returned as sets rather than lists so
    that a KO appearing twice in one route is required once.
    """
    if isinstance(node, Ko):
        return [frozenset((node.name,))]
    if isinstance(node, Optional_):
        return [frozenset()]
    if isinstance(node, Or):
        out: List[FrozenSet[str]] = []
        for child in node.children:
            out.extend(_required_sets(child))
        return out[:MAX_PATHS]
    if isinstance(node, (Seq, All)):
        combined: List[FrozenSet[str]] = [frozenset()]
        for child in node.children:
            child_sets = _required_sets(child)
            grown: List[FrozenSet[str]] = []
            for left in combined:
                for right in child_sets:
                    grown.append(left | right)
                    if len(grown) >= MAX_PATHS:
                        return grown
            combined = grown
        return combined
    raise TypeError("not a definition node: %r" % (node,))


def all_kos(node: Node) -> Set[str]:
    """Every KO mentioned anywhere in *node*, essential or not."""
    if isinstance(node, Ko):
        return {node.name}
    if isinstance(node, Optional_):
        return all_kos(node.child)
    found: Set[str] = set()
    for child in node.children:  # type: ignore[attr-defined]
        found |= all_kos(child)
    return found


@dataclass(frozen=True)
class ModuleScore:
    """One module scored against one genome's KOs, with the flags that qualify the score.

    Each flag corresponds to a documented way completeness misleads, which is why they are part of
    the result rather than something a reader is expected to look up per module.
    """

    module: str
    name: str
    module_class: str
    #: Best alternative route's fraction of present steps. anvi'o's primary score.
    completeness_pathwise: float
    #: Fraction of top-level steps that can be fully satisfied. The stricter reading.
    completeness_stepwise: float
    #: KOs of the best pathwise route that are present, and the route's size.
    matched_kos: Tuple[str, ...]
    route_size: int
    n_top_level_steps: int
    #: True when a reaction in this module has no KO, so pathwise completeness cannot reach 1.0.
    has_no_ko_step: bool
    #: True when the definition refers to another module, so this score is of the reference rather
    #: than of what the referenced module expands to.
    references_modules: bool
    #: Of the KOs this module needs, the fraction that appear in no other module. anvi'o's
    #: mitigation for promiscuous KOs: a module made of KOs shared with many others is easy to score
    #: as present for reasons that have nothing to do with it. 45 KOs appear in five or more modules.
    proportion_unique_enzymes: float
    #: True if path enumeration hit MAX_PATHS, in which case pathwise is a lower bound.
    paths_truncated: bool
    definition_snapshot: str = DEFINITION_SNAPSHOT


@dataclass
class Module:
    """A module definition, parsed once and scored many times."""

    module: str
    name: str
    module_class: str
    definition: str
    tree: Node = field(repr=False)

    @property
    def is_signature(self) -> bool:
        """Signature modules are markers of a phenotype, not metabolic capability.

        KEGG's own class string says so. Reporting "Carbapenem resistance" beside "Glycolysis" as
        though both were pathways a genome can perform is the mistake this flag exists to prevent -
        and that module in particular is a flat OR of 21 KOs, so any single hit calls it complete.
        """
        return self.module_class.startswith("Signature")


def load_modules(path: str = MODULES_TABLE) -> Dict[str, Module]:
    """Parse the vendored module table. Keyed by module id, in file order."""
    modules: Dict[str, Module] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            definition = (row.get("definition") or "").strip()
            if not row.get("module") or not definition:
                continue
            modules[row["module"]] = Module(
                module=row["module"],
                name=(row.get("name") or "").strip(),
                module_class=(row.get("class") or "").strip(),
                definition=definition,
                tree=parse_definition(definition),
            )
    return modules


def ko_module_counts(modules: Dict[str, Module]) -> Dict[str, int]:
    """How many modules each KO appears in, for ``proportion_unique_enzymes``."""
    counts: Dict[str, int] = {}
    for module in modules.values():
        for ko in all_kos(module.tree):
            counts[ko] = counts.get(ko, 0) + 1
    return counts


def _step_is_complete(step: Node, present: Set[str]) -> bool:
    """True when some alternative route through *step* has every KO it requires."""
    for route in _required_sets(step):
        if route and route <= present:
            return True
        if not route:
            # An all-optional step requires nothing, so there is nothing missing from it. Rare, but
            # treating it as incomplete would penalise a module for a component KEGG called
            # non-essential.
            return True
    return False


def score_module(
    module: Module,
    present: Set[str],
    ko_counts: Optional[Dict[str, int]] = None,
) -> ModuleScore:
    """Score one module against the KOs *present* in one genome."""
    routes = _required_sets(module.tree)
    truncated = len(routes) >= MAX_PATHS

    best_fraction = 0.0
    best_matched: Tuple[str, ...] = ()
    best_size = 0
    for route in routes:
        if not route:
            continue
        matched = route & present
        fraction = len(matched) / float(len(route))
        if fraction > best_fraction or (fraction == best_fraction and len(route) > best_size):
            best_fraction = fraction
            best_matched = tuple(sorted(matched))
            best_size = len(route)

    steps = top_level_steps(module.tree)
    complete_steps = sum(1 for step in steps if _step_is_complete(step, present))

    needed = all_kos(module.tree) - {NO_KO}
    if ko_counts and needed:
        unique = sum(1 for ko in needed if ko_counts.get(ko, 0) <= 1)
        proportion_unique = unique / float(len(needed))
    else:
        proportion_unique = 0.0

    return ModuleScore(
        module=module.module,
        name=module.name,
        module_class=module.module_class,
        completeness_pathwise=best_fraction,
        completeness_stepwise=complete_steps / float(len(steps)) if steps else 0.0,
        matched_kos=best_matched,
        route_size=best_size,
        n_top_level_steps=len(steps),
        has_no_ko_step=NO_KO in all_kos(module.tree),
        references_modules=any(ko.startswith("M") for ko in all_kos(module.tree) if ko != NO_KO),
        proportion_unique_enzymes=proportion_unique,
        paths_truncated=truncated,
    )


def score_all(
    present: Iterable[str],
    modules: Optional[Dict[str, Module]] = None,
    ko_counts: Optional[Dict[str, int]] = None,
) -> Iterator[ModuleScore]:
    """Score every module against one genome's KOs."""
    mods = modules if modules is not None else load_modules()
    counts = ko_counts if ko_counts is not None else ko_module_counts(mods)
    have = set(present)
    for module in mods.values():
        yield score_module(module, have, counts)


#: Columns of the per-module table, in order.
MODULE_COLUMNS: Tuple[str, ...] = (
    "genome",
    "module",
    "name",
    "class",
    "completeness_pathwise",
    "completeness_stepwise",
    "n_top_level_steps",
    "route_size",
    "n_matched_kos",
    "matched_kos",
    "is_signature_module",
    "has_no_ko_step",
    "references_modules",
    "proportion_unique_enzymes",
    "paths_truncated",
    "definition_snapshot",
)


def module_row(genome: str, score: ModuleScore, is_signature: bool) -> Dict[str, object]:
    return {
        "genome": genome,
        "module": score.module,
        "name": score.name,
        "class": score.module_class,
        "completeness_pathwise": "%.4f" % score.completeness_pathwise,
        "completeness_stepwise": "%.4f" % score.completeness_stepwise,
        "n_top_level_steps": score.n_top_level_steps,
        "route_size": score.route_size,
        "n_matched_kos": len(score.matched_kos),
        "matched_kos": ",".join(score.matched_kos),
        "is_signature_module": int(is_signature),
        "has_no_ko_step": int(score.has_no_ko_step),
        "references_modules": int(score.references_modules),
        "proportion_unique_enzymes": "%.4f" % score.proportion_unique_enzymes,
        "paths_truncated": int(score.paths_truncated),
        "definition_snapshot": score.definition_snapshot,
    }
