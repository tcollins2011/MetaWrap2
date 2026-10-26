"""KEGG module definitions, and how complete a genome's KOs make each one.

Published tools use different completeness algorithms over the same definitions. On the same KO
set, anvi'o reports 11-12% more modules at >=80% completeness than MicrobeAnnotator or KofamScan,
so completeness values are only comparable when the algorithm is named.

Completeness is computed here instead of by an external tool (as :mod:`metawrap2.depth` does for
the depth table), and every module gets two scores:

* **pathwise**: enumerate every alternative path through the definition, score each as the
  fraction of its steps whose KOs are present, and take the best. This is anvi'o's primary score;
  it is optimistic because the best path is chosen after seeing the data.
* **stepwise**: the fraction of top-level steps that are complete, where a step is complete if
  any of its alternatives is. This is stricter and measures whether the pathway is present.

Both are always reported; there is no single "completeness".

Definition grammar, loosest to tightest binding:

====================  ==========================================================================
``K00844 K01810``     a space is AND: sequential steps, both required
``K00844,K12407``     a comma is OR: alternatives for one step
``K00844+K12407``     a plus is AND *within* a step: subunits of a complex, all essential
``K00844-K12407``     a minus marks a non-essential component, which does not count toward the step
``(...)``             grouping
``--``                a reaction with no KO assigned to it
====================  ==========================================================================

A ``--`` step cannot be satisfied by any genome. anvi'o scores it 0, so the module stays below
100% pathwise; EBI's tool gives it weight 0, so it does not limit the module. This file follows
anvi'o, so a module with a reaction of no known gene is never reported complete.
:data:`ModuleScore.has_no_ko_step` marks the modules affected.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

#: The vendored snapshot. See data/kegg/SNAPSHOT.md for its provenance and how to move it.
MODULES_TABLE = os.path.join(os.path.dirname(__file__), "data", "kegg", "modules_table.tsv")

#: The snapshot date, written on every scored row. Hard-coded because a checkout resets mtimes.
DEFINITION_SNAPSHOT = "2026-10-12"

#: Stand-in for a reaction with no KO (``--`` in a definition). No genome can carry it; see the
#: module docstring.
NO_KO = "--"

#: Cap on path enumeration. M00001 has the most routes (~22k), so the cap is not normally reached;
#: it prevents a pathological definition from hanging. ``paths_truncated`` records when it is hit.
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
    """Subunits of one complex, ``+``-joined. All are essential and together form one step.

    Separate from :class:`Seq` so that :func:`top_level_steps` does not split a complex. M00142
    (NADH dehydrogenase) is seven ``+``-joined subunits: a genome with three of them scores 0, not
    3/7. 31 of the 573 modules are a single complex.
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
    """Raised for a definition that does not parse."""


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

    Precedence, loosest first: space (AND), comma (OR), plus/minus (complex). So ``K1 K2,K3``
    means "K1, then K2 or K3", as in published parsers.
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

        A leading ``-`` (``K11023 -K11024 K11025``) marks the whole step non-essential. Seven of
        the 573 modules use this form.
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

    Stepwise completeness is a fraction of these. The count is reported per module: 8% of modules
    have one top-level step, so their completeness is binary and one false-positive KO call marks
    the module present.
    """
    if isinstance(node, Seq):
        return list(node.children)
    return [node]


def _required_sets(node: Node) -> List[FrozenSet[str]]:
    """Every alternative route through *node*, each as the set of KOs it requires.

    AND takes the product of its children's routes, OR concatenates them, and an Optional
    contributes nothing. Sets, so a KO repeated in one route is required once.
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

    Each flag marks a known way completeness can mislead.
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
    #: True when the definition refers to another module; the reference is not expanded.
    references_modules: bool
    #: Fraction of this module's KOs found in no other module (anvi'o's measure for shared KOs,
    #: which can make a module look present). 45 KOs appear in five or more modules.
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
        """Signature modules mark a phenotype, not a metabolic capability (per KEGG's class).

        For example, "Carbapenem resistance" is a flat OR of 21 KOs, so any single hit makes it
        complete.
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
            # An all-optional step requires nothing, so it counts as complete.
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
