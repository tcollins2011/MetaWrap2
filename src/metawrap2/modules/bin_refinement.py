"""MetaWrap2 bin_refinement module: consolidate several bin sets into one refined set.

Given any number of folders of bins produced from the same assembly (e.g. metaBAT2, MaxBin2,
CONCOCT outputs), this hybridizes them with Binning_refiner, scores every candidate set with CheckM,
consolidates the best version of each bin, and dereplicates shared contigs. The scientific
core lives in :mod:`metawrap2.refinement` (consolidate/dereplicate) and
:mod:`metawrap2.checkm` (CheckM); this module orchestrates them and preserves the original
output layout.

``--consolidator`` replaces that middle part with an external tool (currently DAS_Tool), which
decides from single-copy genes in one pass instead of from CheckM on every combination. Both
routes leave one bin set in ``binsO``, and everything after that - re-scoring, the -c/-x cut, the
plots, the final layout, GUNC - is the same code either way, which is what makes the two
comparable on the NC MAGs each recovers, with each route's cost reported beside that.

This file is meant to be read and edited: the commands each external tool runs live in the
COMMANDS block below and the tunable numbers in the CONSTANTS block.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import os
import shlex
import shutil
from typing import Dict, List, Optional, Sequence, Tuple

from .. import checkm as _checkm
from .. import gunc as _gunc
from .. import mimag as _mimag
from .. import recovery as _recovery
from .. import refinement as _refinement
from ..cliformat import bounded_float, make_parser, percent
from ..config import MODULE_ENVS, load_settings
from ..constants import (
    BIN_EXTENSION,
    CONTIGS_SUFFIX,
    FASTA_EXTENSIONS,
    MAX_BIN_SIZE,
    MIN_BIN_SIZE,
    STATS_SUFFIX,
    TSV_SUFFIX,
    bin_filename,
    is_refined_bins_dirname,
    is_unbinned_name,
    refined_bins_dirname,
)
from ..io.seqio import contig_id, fasta_files, iter_fasta
from ..progress import bar
from ..pyrun import PY
from ..utils import strip_trailing_slash
from ..validate import check_fasta
from ._common import (
    absolutize_paths,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    force,
    fs_replace_dir,
    listdir,
    make_checkpoint,
    memory_arg,
    optional_env_for,
    require_nonempty_dir,
    resolve_memory_gb,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = MODULE_ENVS["bin_refinement"]

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# These are MetaWrap2's own helpers, so they run in the host interpreter (env=None), not in
# the module's conda env - that env holds CheckM, not python+metawrap2. See metawrap2.pyrun.
BINNING_REFINER = PY + " -m metawrap2.vendor.binning_refiner -1 {b1} -2 {b2} -o {out}"
BINNING_REFINER_3 = PY + " -m metawrap2.vendor.binning_refiner -1 {b1} -2 {b2} -3 {b3} -o {out}"
BINNING_REFINER_N = PY + " -m metawrap2.vendor.binning_refiner {bins} -o {out}"
PLOT_BINNING = PY + " -m metawrap2.scripts.plot_binning_results {comp} {cont} {stats}"
# --search_engine diamond, not the DAS_Tool default: the default is usearch, which is not free
# software and is not in the env, so the default command would simply fail to find it.
DAS_TOOL = (
    "DAS_Tool -i {tables} -l {labels} -c {assembly} -o {prefix} -t {threads}"
    " --search_engine diamond --score_threshold {score} --write_bins"
)
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_MIN_COMPLETION = 70  # -c; bins must exceed this % completion to be "good"
DEFAULT_MAX_CONTAM = 10  # -x; bins must stay below this % contamination
# bin_refinement ignores bins outside this size range (bytes); see constants.py
MIN_SIZE = MIN_BIN_SIZE  # 50 kb
MAX_SIZE = MAX_BIN_SIZE  # 20 Mb
#: Upper limit on input bin sets. Refinement cost grows with the number of combinations, and
#: CheckM is run on every one of them, so this is a guard against an accidental 20-set run
#: rather than an algorithmic limit.
MAX_BIN_SETS = 6
#: How many bin-set combinations to build. "pairs-and-all" reproduces what the 2- and 3-set
#: code always did (every pair, plus the all-way combination); see refinement_plan().
COMBINATION_MODES = ("pairs-and-all", "pairs", "all")
DEFAULT_COMBINATIONS = "pairs-and-all"
#: How the bin sets are turned into one set. "metawrap" is this module's own algorithm - hybridise
#: every combination with Binning_refiner, score each with CheckM, then pick the best version of
#: each bin. "dastool" scores the input bins directly on single-copy genes and aggregates in one
#: pass, which is a different accuracy/cost point rather than a different implementation of the same
#: idea. The default is "metawrap" because it is the incumbent, not because it has been shown to
#: win: that is decided on NC MAGs recovered once the comparison has run, with cost beside it. On
#: per-sample bins so far, metawrap recovers 13 to DAS Tool's 10, for 4.52 cpu-hours against 0.62.
CONSOLIDATORS = ("metawrap", "dastool")
DEFAULT_CONSOLIDATOR = "metawrap"
#: DAS_Tool's own default score threshold. Named here because it is a quality floor on the bins it
#: keeps, so it belongs next to -c/-x rather than buried in a command string.
DEFAULT_DASTOOL_SCORE = 0.5
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = make_parser(
        prog="metawrap2 bin_refinement",
        usage="metawrap2 bin_refinement [options] -o output_dir -A bin_folderA "
        "[-B bin_folderB ... -F bin_folderF]  |  --bins DIR [--bins DIR ...]",
        description="Consolidate up to %d bin sets (from the same assembly) into a refined "
        "set. Contig names must be consistent across the bin folders." % MAX_BIN_SETS,
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    threads_arg(p)
    memory_arg(p)
    p.add_argument(
        "-c",
        "--completeness",
        type=percent(),
        default=DEFAULT_MIN_COMPLETION,
        help="minimum %% completion of bins; below 50 is not a genome by any MIMAG tier",
    )
    p.add_argument(
        "-x",
        "--contamination",
        type=percent(),
        default=DEFAULT_MAX_CONTAM,
        help="maximum acceptable %% contamination of bins",
    )
    # -A..-F name the bin sets positionally, as they always have. --bins may be repeated
    # instead, which is easier to generate from a script or a loop.
    p.add_argument("-A", "--bins-a", help="folder with metagenomic bins (.fa/.fasta)")
    for letter in "BCDEF":
        p.add_argument(
            "-" + letter, "--bins-" + letter.lower(), help="another folder with metagenomic bins"
        )
    p.add_argument(
        "--bins",
        action="append",
        default=[],
        dest="bins_extra",
        metavar="DIR",
        help="a folder with metagenomic bins; repeat once per bin set "
        "(alternative to -A/-B/-C...)",
    )
    p.add_argument(
        "--refine-combinations",
        choices=COMBINATION_MODES,
        default=DEFAULT_COMBINATIONS,
        help="which bin-set combinations to build with Binning_refiner: "
        "'pairs-and-all' (default) every pair plus the all-way combination; "
        "'pairs' every pair only (cheaper); 'all' every subset of size >=2 "
        "(thorough, exponential). 2- and 3-set runs are unaffected.",
    )
    p.add_argument(
        "--skip-refinement",
        dest="refine",
        action="store_false",
        help="don't use Binning_refiner to make refined bins from binner combinations",
    )
    p.add_argument(
        "--score-everything",
        action="store_true",
        help="score every bin with CheckM, including ones whose size or fragmentation already rules "
        "them out. Slower, and the only reason to pass it is to check the filter is not wrong",
    )
    p.add_argument(
        "--skip-checkm",
        dest="run_checkm",
        action="store_false",
        help="don't run CheckM to assess bins",
    )
    p.add_argument(
        "--skip-gunc",
        dest="run_gunc",
        action="store_false",
        help="don't run GUNC on the final bins; skips the chimerism check CheckM cannot make "
        "(and its reference database)",
    )
    p.add_argument(
        "--skip-consolidation",
        dest="consolidate",
        action="store_false",
        help="don't choose the best version of each bin; just pick the best folder",
    )
    amb = p.add_mutually_exclusive_group()
    amb.add_argument(
        "--keep-ambiguous",
        dest="dereplicate",
        action="store_const",
        const="keep",
        help="for contigs in more than one bin, keep them in all bins",
    )
    amb.add_argument(
        "--remove-ambiguous",
        dest="dereplicate",
        action="store_const",
        const="remove",
        help="for contigs in more than one bin, remove them from all bins",
    )
    p.set_defaults(dereplicate="partial")  # default: keep only in the best bin
    p.add_argument(
        "--quick",
        action="store_true",
        help="add --reduced_tree to CheckM, reducing runtime (helps with low memory)",
    )
    p.add_argument(
        "--checkm2",
        action="store_true",
        help="score bins with CheckM2 instead of CheckM1 (fast, no pplacer memory "
        "wall; needs the metawrap2-checkm2 env)",
    )
    p.add_argument(
        "--consolidator",
        choices=CONSOLIDATORS,
        default=DEFAULT_CONSOLIDATOR,
        help="how to turn the bin sets into one set: 'metawrap' (default) hybridises every "
        "combination and scores each with CheckM; 'dastool' scores the input bins on "
        "single-copy genes and aggregates in one pass, which needs -a and the "
        "metawrap2-binqc env, which `metawrap2 install-env bin_refinement` creates",
    )
    p.add_argument(
        "-a",
        "--assembly",
        help="assembly the bins came from (required by --consolidator dastool, which reads the "
        "contigs itself; MetaWrap2's own consolidation does not need it)",
    )
    p.add_argument(
        "--dastool-score",
        type=bounded_float(0, 1),
        default=DEFAULT_DASTOOL_SCORE,
        help="DAS_Tool --score_threshold: minimum bin score to keep",
    )
    p.add_argument("--config", help="path to metawrap2.toml")
    return p.parse_args(argv)


def collect_bin_sets(args) -> List[str]:
    """The input bin-set directories, in order: -A..-F first, then any repeated --bins.

    Order decides each set's label (binsA, binsB, ...) and therefore the refinement plan and
    the bin numbering inside each refined set, so it is part of the interface rather than an
    implementation detail.
    """
    ordered: List[str] = []
    for letter in "abcdef":
        value = getattr(args, "bins_" + letter, None)
        if value:
            ordered.append(value)
    for value in getattr(args, "bins_extra", []) or []:
        if value not in ordered:
            ordered.append(value)
    if not ordered:
        error(
            "No bin sets given. Pass at least one with -A (and -B, -C, ... for more), or "
            "repeat --bins DIR."
        )
    if len(ordered) > MAX_BIN_SETS:
        error(
            "%d bin sets given, but at most %d are supported.\n"
            "Every combination of them is refined and scored with CheckM, so the cost grows "
            "quickly. Refine in groups instead: consolidate the first few, then consolidate "
            "that result with the rest." % (len(ordered), MAX_BIN_SETS)
        )
    return ordered


def drop_empty_bin_sets(bin_sets: Sequence[str]) -> List[str]:
    """The bin sets that hold at least one bin, warning about each one that holds none.

    A binner that finds nothing to cluster writes an empty folder and exits 0 (metaBAT2 on a small
    assembly does), and that is a result, not a failure: refining the other sets is still worth
    doing. Only a folder that exists and has no FASTA in it is dropped. A missing path is left for
    the input checks to report, since that is a typo rather than an empty result. When every set is
    empty the result is an empty list, and the caller writes an empty final set rather than failing.
    """
    kept: List[str] = []
    for index, path in enumerate(bin_sets):
        if os.path.isdir(path) and not any(
            name.endswith(tuple(FASTA_EXTENSIONS)) for name in os.listdir(path)
        ):
            warning("bin set %d (%s) has no bins in it; refining without it." % (index + 1, path))
            continue
        kept.append(path)
    if not kept:
        warning(
            "Every bin set given is empty (%s), so there is nothing to refine; writing an empty "
            "final bin set.\nCheck the binning logs: a binner that makes no bins usually had too "
            "few contigs above its length cut-off, or no coverage for them." % ", ".join(bin_sets)
        )
    elif len(kept) < len(bin_sets):
        comm(
            "Refining %d of %d bin sets; they are relabelled %s in order."
            % (len(kept), len(bin_sets), ", ".join(bin_set_label(i) for i in range(len(kept))))
        )
    return kept


def write_empty_result(out: str, completeness: float, contamination: float) -> str:
    """Write the final bin set for a run with no bins: an empty folder and a header-only table."""
    final = os.path.join(out, refined_bins_dirname(completeness, contamination))
    os.makedirs(final, exist_ok=True)
    with open(final + STATS_SUFFIX, "w") as fh:
        fh.write("\t".join(_checkm.STATS_HEADER) + "\n")
    return final


def bin_set_label(index: int) -> str:
    """The internal name of the *index*-th input bin set: binsA, binsB, ... binsF."""
    return "bins" + chr(ord("A") + index)


def combined_label(labels: Sequence[str]) -> str:
    """The name of the set refined from *labels*, e.g. (binsA, binsC) -> binsAC."""
    return "bins" + "".join(label[len("bins") :] for label in labels)


#: The exact refinement plan the 2- and 3-input code always used, kept verbatim.
#:
#: Two details of it are load-bearing and would not fall out of a generic scheme: the order
#: the combinations are built in, and the *orientation* of each pair. Binning_refiner numbers
#: its output bins from the order the input folders are given, so building binsBC as (C, B) -
#: which is what this module has always done - produces different bin numbering from (B, C).
#: Any change here changes results for every existing 2- or 3-binner analysis.
_LEGACY_PLANS: Dict[int, List[Tuple[Tuple[str, ...], str]]] = {
    2: [(("binsA", "binsB"), "binsAB")],
    3: [
        (("binsA", "binsB"), "binsAB"),
        (("binsC", "binsB"), "binsBC"),
        (("binsA", "binsC"), "binsAC"),
        (("binsA", "binsB", "binsC"), "binsABC"),
    ],
}


def refinement_plan(
    labels: Sequence[str], combinations: str = DEFAULT_COMBINATIONS
) -> List[Tuple[Tuple[str, ...], str]]:
    """The bin-set combinations to refine, as [(input labels, output label), ...].

    For two or three inputs this returns :data:`_LEGACY_PLANS` unchanged, so those runs
    reproduce exactly what they did before this function existed. For more inputs:

    ``pairs-and-all`` (default)
        every pair, then the combination of all of them - the natural extension of what the
        3-input case did (3 pairs + the triple).
    ``pairs``
        every pair only. Cheaper: the number of CheckM runs grows quadratically, not
        exponentially.
    ``all``
        every subset of size >= 2. Thorough and expensive: 2**n - n - 1 combinations.
    """
    if combinations not in COMBINATION_MODES:
        raise ValueError(
            "combinations must be one of %s, got %r" % (", ".join(COMBINATION_MODES), combinations)
        )
    labels = list(labels)
    if len(labels) < 2:
        return []
    if combinations == DEFAULT_COMBINATIONS and len(labels) in _LEGACY_PLANS:
        # Only valid for the canonical A/B/C labels, which is how main() names them.
        legacy = _LEGACY_PLANS[len(labels)]
        if labels == [bin_set_label(i) for i in range(len(labels))]:
            return list(legacy)

    plan: List[Tuple[Tuple[str, ...], str]] = []
    if combinations == "all":
        sizes = range(2, len(labels) + 1)
    else:
        sizes = range(2, 3)
    for size in sizes:
        for combo in itertools.combinations(labels, size):
            plan.append((combo, combined_label(combo)))
    if combinations == DEFAULT_COMBINATIONS and len(labels) > 2:
        all_way = tuple(labels)
        plan.append((all_way, combined_label(all_way)))
    return plan


def _fix_contig_naming(path: str) -> bool:
    """Normalise a bin fasta's headers in place. Returns True if the file was rewritten.

    Two things happen here: the header is reduced to the contig id (:func:`contig_id`),
    dropping description text that binners append - current metaBAT2 writes
    ``total_depth=.. sample_depths=..`` after the name, which made its bins match nothing -
    and any remaining '=' becomes '_' (ports fix_config_naming.py, which existed because
    '=' breaks some downstream tools).

    Files whose headers are already clean are left untouched. Rewriting every bin on every
    invocation meant re-reading and re-writing the entire bin set each time the module ran,
    which on a large study is a lot of pointless I/O - and it needlessly changed mtimes that
    a user (or a resumed run) might be relying on.
    """
    records = list(iter_fasta(path))
    fixed = [(contig_id(h).replace("=", "_"), seq) for h, seq in records]
    if all(new == old for (new, _), (old, _) in zip(fixed, records)):
        return False
    with open(path, "w") as out:
        out.writelines(">%s\n%s\n" % (header, seq) for header, seq in fixed)
    return True


def _copy_bins(src: str, dest: str) -> int:
    """Copy size-filtered bins from *src* into *dest* as bin-name.fa. Returns count kept.

    The unbinned catch-all is not a bin and is skipped by name. A bin set is allowed to carry one:
    `binning` writes CONCOCT's leftovers to ``concoct_bins/unbinned.fasta`` on purpose, so that a
    user can see what no binner claimed, and the AMBER scorer skips it for the same reason - scoring
    the leftovers as a genome would be scoring a bin nobody proposed. Here it was being copied in as
    a candidate genome, which puts it through CheckM, into every pairwise Binning_refiner
    combination, and into the dereplication that decides which bin each contig belongs in. On the
    CAMI2 GI run it happened not to matter, because CONCOCT's leftovers came to 69.8 MB and the size
    filter below drops anything over 20 MB - the right outcome reached for the wrong reason, and only
    while the catch-all stays large. A tighter CONCOCT run leaves a smaller one, and a smaller one
    gets refined as if it were a genome.

    What the size filter drops is summarised at the end, separately for the two ends, because the
    two mean different things and the per-bin lines above do not add up to either. Over the ceiling
    is a catch-all bin - too much sequence to be one genome, so it could never have become a MAG and
    losing it costs nothing but the recall figure. Under the floor is a fragment, which might have
    grown into something had a different binner been in the set. On a CAMI2 CONCOCT set this filter
    silently took 10 of 99 bins holding 41% of the set's sequence - 7 over the ceiling, the largest
    131.8 Mb - which is most of the recall that refinement appears to lose; a run that does not say
    so leaves that looking like refinement's doing. Measured on file size rather than base count,
    which is what the filter itself compares.
    """
    ensure_dir(dest)
    n = 0
    dropped: dict[str, list[int]] = {"over": [], "under": []}
    total = 0
    for name in sorted(os.listdir(src)):
        full = os.path.join(src, name)
        if not os.path.isfile(full):
            continue
        if is_unbinned_name(name):
            comm("Skipping %s because it is the unbinned leftovers, not a bin" % full)
            continue
        size = os.path.getsize(full)  # portable (no `stat -c`)
        total += size
        if MIN_SIZE < size < MAX_SIZE:
            base = os.path.splitext(name)[0]
            shutil.copy(full, os.path.join(dest, base + BIN_EXTENSION))
            n += 1
        else:
            dropped["over" if size >= MAX_SIZE else "under"].append(size)
            comm("Skipping %s because the bin size is not between 50kb and 20Mb" % full)
    for end, label in (("over", "above the 20Mb ceiling"), ("under", "below the 50kb floor")):
        sizes = dropped[end]
        if sizes:
            comm(
                "  %d bin(s) %s, holding %.1f%% of the set's sequence (largest %.1fMb)"
                % (len(sizes), label, 100.0 * sum(sizes) / total, max(sizes) / 1e6)
            )
    return n


def _count_bins(path: str) -> int:
    return (
        len([f for f in os.listdir(path) if f.endswith(BIN_EXTENSION)])
        if os.path.isdir(path)
        else 0
    )


def write_contig_table(bins_dir: str, path: str) -> int:
    """Write *bins_dir*'s contig-to-bin membership to *path* as TSV. Returns the row count.

    Two consumers want the same table for different reasons. Anvio and similar tools read the one
    laid down beside each final bin set, and DAS_Tool takes one per input bin set as its only view
    of what each binner decided. Writing it in one place keeps those identical, which matters
    because a consolidator reading a differently-shaped table silently consolidates nothing.
    """
    rows = 0
    with open(path, "w") as out_fh:
        for name in sorted(os.listdir(bins_dir)):
            if not name.endswith(BIN_EXTENSION):
                continue
            bin_name = os.path.splitext(name)[0]
            for header, _seq in iter_fasta(os.path.join(bins_dir, name)):
                out_fh.write("%s\t%s\n" % (contig_id(header), bin_name))
                rows += 1
    return rows


def assembly_for_bins(assembly: str, work: str) -> str:
    """The assembly to hand a consolidator, with headers that match the copied bins.

    The bins in *out* have been through :func:`_fix_contig_naming`, which reduces a header to its
    contig id and turns any remaining '=' into '_'. A consolidator that reads the assembly itself
    matches its contig names against the names in the bin tables, so if normalisation changed
    anything the two no longer agree - and the failure mode is not an error but an empty
    consolidation, because every contig simply fails to be found.

    Records with no sequence are dropped here too, for a reason that has nothing to do with naming.
    DAS_Tool reads the assembly with R's ``fread`` and a header with an empty sequence yields a line
    with one field where it expects two, which stops the read part way and then reports every contig
    after that point as "not found in assembly" - 36,854 of them, in a message that names none of the
    real problem. CAMI2's gold-standard assembly ships two such records, so this is not hypothetical:
    `binning` already filters them out of its own working copy, and a consolidator reading the
    original was getting what binning had been careful not to pass on.

    A copy is written only when something actually needed changing, which for metaSPAdes and MEGAHIT
    names with no empty records it does not. Assemblies are gigabytes, so a copy nothing needed would
    be the most expensive no-op in the module.
    """
    renamed = [(contig_id(h).replace("=", "_"), h, seq) for h, seq in iter_fasta(assembly)]
    empty = [new for new, _old, seq in renamed if not seq]
    if not empty and all(new == contig_id(old) for new, old, _seq in renamed):
        return assembly
    ensure_dir(work)
    normalized = os.path.join(work, "assembly_normalized.fa")
    reasons = []
    if any(new != contig_id(old) for new, old, _seq in renamed):
        reasons.append("headers need the same normalisation the bins got")
    if empty:
        reasons.append(
            "%d record(s) have no sequence, which DAS_Tool's assembly reader stops on" % len(empty)
        )
    comm("writing %s: %s" % (os.path.basename(normalized), "; ".join(reasons)))
    with open(normalized, "w") as fh:
        fh.writelines(">%s\n%s\n" % (new, seq) for new, _old, seq in renamed if seq)
    return normalized


def _run_dastool(args, env: Optional[str], out: str, labels: Sequence[str]) -> None:
    """Consolidate the bin sets in *out* with DAS_Tool, leaving the result as ``binsO``.

    binsO is where this module's own consolidation leaves its answer, so landing here means the
    whole downstream half - re-scoring with CheckM, dropping bins under -c/-x, the ranking plots,
    the final layout and GUNC - runs identically whichever consolidator produced it. That is what
    makes the arms comparable: the only thing that differs is how the one set was arrived at.
    """
    work = os.path.join(out, "work_files")
    ensure_dir(work)
    tables = []
    for label in labels:
        table = os.path.join(work, label + ".dastool.tsv")
        if not dry_run():
            rows = write_contig_table(os.path.join(out, label), table)
            comm("%s contributes %d contigs to DAS_Tool" % (label, rows))
        tables.append(table)

    assembly = args.assembly if dry_run() else assembly_for_bins(args.assembly, work)
    prefix = os.path.join(work, "dastool", "DASTool")
    ensure_dir(os.path.dirname(prefix))
    run(
        DAS_TOOL.format(
            tables=",".join(tables),
            labels=",".join(labels),
            assembly=assembly,
            prefix=prefix,
            threads=args.threads,
            score=args.dastool_score,
        ),
        env=env,
        tool="DAS_Tool",
    )
    if dry_run():
        return
    produced = prefix + "_DASTool_bins"
    if not os.path.isdir(produced):
        error(
            "DAS_Tool did not produce %s.\n"
            "Its own log is in %s. A DAS_Tool run that finds no single-copy genes ends this way, "
            "which usually means the contig names in the bins do not match the assembly given "
            "with -a." % (os.path.basename(produced), os.path.dirname(prefix))
        )
    binsO = os.path.join(out, "binsO")
    ensure_dir(binsO)
    # Renamed to bin_001.fasta, not copied under DAS_Tool's own names. Every bin set MetaWrap2
    # produces is numbered this way, the .stats.tsv rows are keyed on these stems, and a
    # consolidator whose output is named differently would make the two arms distinguishable by
    # something other than their contents. DAS_Tool's names say which input set each bin came
    # from, which is worth keeping, so the mapping is written out beside the bins rather than
    # discarded.
    kept = 0
    provenance = []
    for name in sorted(os.listdir(produced)):
        if os.path.splitext(name)[1] in FASTA_EXTENSIONS:
            kept += 1
            renamed = bin_filename(kept)
            shutil.copy(os.path.join(produced, name), os.path.join(binsO, renamed))
            provenance.append((os.path.splitext(renamed)[0], name))
    with open(os.path.join(work, "dastool_bin_names.tsv"), "w") as fh:
        fh.write("bin\tdastool_name\n")
        fh.writelines("%s\t%s\n" % row for row in provenance)
    comm("DAS_Tool kept %d bins (score threshold %s)" % (kept, args.dastool_score))
    if kept == 0:
        error(
            "DAS_Tool kept no bins at --dastool-score %s. Its per-bin scores are in "
            "%s_DASTool_summary.tsv; lower the threshold if the best score there is acceptable."
            % (args.dastool_score, prefix)
        )


def _refine(out: str, env: Optional[str], inputs: Sequence[str], dest: str) -> None:
    """Run Binning_refiner on any number of bin folders (relative to *out*) into *dest*.

    Two and three folders are passed as -1/-2/-3, exactly as before, so the command line for
    those cases is byte-identical to what it always was. More than three are passed with the
    repeatable -i flag.
    """
    refined = "Refined_" + dest[len("bins") :]  # binsAB -> Refined_AB
    if len(inputs) == 2:
        cmd = BINNING_REFINER.format(b1=inputs[0], b2=inputs[1], out=refined)
    elif len(inputs) == 3:
        cmd = BINNING_REFINER_3.format(b1=inputs[0], b2=inputs[1], b3=inputs[2], out=refined)
    else:
        cmd = BINNING_REFINER_N.format(
            bins=" ".join("-i %s" % shlex.quote(b) for b in inputs), out=refined
        )
    run(cmd, env=None, tool="binning_refiner", cwd=out)
    produced = os.path.join(out, refined, "Refined")
    if dry_run():
        return
    if not os.path.isdir(produced):
        error("Bin_refiner did not finish correctly for %s. Exiting..." % dest)
    # Destination-safe: this is reached on a --resume redoing the refine step after a crash, and what
    # is already there is that crashed run's partial version of the set about to be written.
    fs_replace_dir(produced, os.path.join(out, dest))
    comm("there are %d refined bins in %s" % (_count_bins(os.path.join(out, dest)), dest))
    shutil.rmtree(os.path.join(out, refined))


def _normalize_bin_extensions(bin_set_dir: str) -> None:
    """Give every bin in a bin-set directory the canonical BIN_EXTENSION.

    Input bin sets come from whatever produced them, so the same directory can hold ``.fa``,
    ``.fna`` and ``.fasta``. CheckM is invoked with a single ``-x`` extension and would silently
    score only the files matching it, so they are normalised first rather than half-processed.
    """
    for f in listdir(bin_set_dir):
        stem, ext = os.path.splitext(f)
        if ext and ext in FASTA_EXTENSIONS and ext != BIN_EXTENSION:
            os.replace(
                os.path.join(bin_set_dir, f), os.path.join(bin_set_dir, stem + BIN_EXTENSION)
            )


def _good_bins(stats_path: str, comp: float, cont: float) -> List[str]:
    """Return names of bins in a .stats file meeting the completion/contamination thresholds."""
    good: List[str] = []
    if not os.path.isfile(stats_path):
        return good
    with open(stats_path) as fh:
        for line in fh:
            if "completeness" in line or "compl" in line:
                continue
            cut = line.rstrip("\n").split("\t")
            if len(cut) < 3:
                continue
            if comp <= float(cut[1]) <= 100 and 0 <= float(cut[2]) <= cont:
                good.append(cut[0])
    return good


def _run_checkm_set(
    bin_set_dir: str,
    args,
    env: Optional[str],
    checkm2_env: Optional[str] = None,
    cache: Optional[_checkm.BinScoreCache] = None,
) -> str:
    """Run CheckM on one bin set and report the number of 'good' bins. Returns the .stats path.

    Through `run_checkm_cached`, because this is called once per hybrid bin set and the sets overlap: on
    a three-binner refinement, 32% of the bins across the seven sets are the same bin scored again. The
    reported numbers are identical either way - the same contigs have the same completeness.
    """
    comm(
        "Running %s on %s bins"
        % ("CheckM2" if args.checkm2 else "CheckM", os.path.basename(bin_set_dir))
    )
    stats = _checkm.run_checkm_cached(
        bin_set_dir,
        cache=cache,
        prefilter=not args.score_everything,
        say=comm,
        threads=args.threads,
        mem_gb=args.memory,
        env=env,
        quick=args.quick,
        checkm2=args.checkm2,
        checkm2_env=checkm2_env,
    )
    n = len(_good_bins(stats, args.completeness, args.contamination))
    comm(
        "There are %d 'good' bins found in %s! (>%s%% completion and <%s%% contamination)"
        % (n, os.path.basename(bin_set_dir), args.completeness, args.contamination)
    )
    return stats


def _bin_set_dirs(out: str) -> List[str]:
    """Bin-set directories in *out* (binsA, binsB, ..., binsAB, ...), sorted, excluding binsM/binsO."""
    dirs = []
    for name in sorted(listdir(out)):
        if (
            name.startswith("bins")
            and name not in ("binsM", "binsO")
            and os.path.isdir(os.path.join(out, name))
        ):
            dirs.append(name)
    return dirs


def main(argv: List[str]) -> int:
    args = _parse_args(argv)
    settings = load_settings(args.config)
    resolve_threads(args, settings)
    resolve_memory_gb(args, settings)
    absolutize_paths(args)
    env = env_for("bin_refinement", settings)
    # CheckM2 is the one tool here that cannot share: it has no common solution with CheckM1, so it
    # keeps an env of its own and this is the only env lookup in the module that can fail separately.
    checkm2_env = optional_env_for("checkm2", settings, "--checkm2") if args.checkm2 else None
    # GUNC and DAS Tool are in this module's own env. They used to have one each, which made a step
    # that runs by default - GUNC, unless --skip-gunc - depend on a second install command.
    gunc_env = env
    # An alternative consolidator replaces the middle of the pipeline rather than adding to it, so
    # the precondition it still has is checked before any bin is copied: the assembly it cannot work
    # without. Finding out after prodigal has run on every bin would be the whole cost of the arm for
    # nothing. Its env is no longer a separate precondition - DAS Tool ships with CheckM1.
    external = args.consolidator != "metawrap"
    external_env = env if external else None
    if external and not args.assembly:
        error(
            "--consolidator %s reads the assembly contigs itself, so it needs -a/--assembly.\n"
            "Pass the same assembly the bins were produced from. MetaWrap2's own consolidation "
            "(--consolidator metawrap, the default) works from the bins alone and does not."
            % args.consolidator
        )
    if external and not args.consolidate:
        # Named as a contradiction rather than resolved by precedence. Either order of precedence
        # gives the user something they did not ask for, and both flags were typed deliberately.
        error(
            "--consolidator %s and --skip-consolidation ask for opposite things: one names which "
            "consolidator to use, the other says not to consolidate at all.\n"
            "Drop --skip-consolidation to run %s, or drop --consolidator to pick the best single "
            "bin set without mixing bins." % (args.consolidator, args.consolidator)
        )
    given = collect_bin_sets(args)
    bin_sets = drop_empty_bin_sets(given)
    if not bin_sets:
        rec = start_run("bin_refinement", args, env, settings, inputs=list(given))
        try:
            if not dry_run():
                final = write_empty_result(args.output, args.completeness, args.contamination)
                comm("0 bins -> %s" % final)
            announcement("BIN_REFINEMENT FINISHED WITH NO BINS")
        finally:
            finish_run(rec, args.output)
        return 0
    validate_inputs(
        "bin_refinement",
        [
            (
                lambda path, what: require_nonempty_dir(path, what, FASTA_EXTENSIONS),
                path,
                "bin set %s" % bin_set_label(i),
            )
            for i, path in enumerate(bin_sets)
        ]
        + ([(check_fasta, args.assembly, "assembly (-a)")] if external else []),
    )
    rec = start_run("bin_refinement", args, env, settings, inputs=list(bin_sets))

    ckpt = make_checkpoint(args.output)
    try:
        # Checked here rather than at the point of use: GUNC runs last, and finding out that its
        # database is absent after CheckM has spent hours on every bin-set combination is the worst
        # possible moment to find out.
        gunc_db = settings.db("GUNC_DB") if args.run_gunc else ""
        if args.run_gunc and not (gunc_db and os.path.exists(gunc_db)):
            # A dry run is asked for precisely when the tools and databases are not in place yet, so
            # it reports the problem and carries on printing the plan. A real run stops.
            (warning if dry_run() else error)(
                "GUNC runs by default and its database is not configured (GUNC_DB=%s). %s"
                % (gunc_db or "<not set>", _gunc.DB_HINT)
            )

        comp, cont = args.completeness, args.contamination
        folders = [(bin_set_label(i), src) for i, src in enumerate(bin_sets)]
        for _label, src in folders:
            if not os.path.isdir(src):
                error("%s is not a valid directory. Exiting." % src)

        announcement("BEGIN PIPELINE!")
        out = args.output
        comm("setting up output folder and copying over bins...")
        # Only under --force. --resume arrives here with the same already-existing output folder and
        # means the opposite: it is about to skip the steps whose checkpoints are set, so clearing
        # their output leaves it with nothing to skip to. Unconditional, this deleted out/bins* on
        # every resume - and since the bin sets only move into work_files at the finalize step, that
        # was every resume between "refine" and the end.
        if os.path.isdir(out) and force():
            warning("Warning: %s already exists. Attempting to clean." % out)
            # Remove every bins* working directory, whatever the combination names are for
            # this number of inputs, plus the consolidated (M) and dereplicated (O) sets.
            #
            # Under work_files as well as in the output root, because that is where a finished run
            # leaves them: cleaning only the root left the previous run's sets sitting in work_files,
            # and a rerun with fewer inputs then kept a binsC nobody in this run produced. The
            # staleness is the reason to do it - every set this run does produce is overwritten at
            # the end anyway.
            # Sidecars as well as directories: a binsC.stats.tsv left behind by a 3-set run is as
            # misleading as the folder was, and it is the file the stats tables and the plots are
            # built from.
            for parent in (out, os.path.join(out, "work_files")):
                for name in listdir(parent):
                    path = os.path.join(parent, name)
                    if not name.startswith("bins"):
                        continue
                    if os.path.isdir(path):
                        shutil.rmtree(path, ignore_errors=True)
                    else:
                        os.remove(path)
            # The previous run's final set goes too, with its sidecars. Not by name: the directory is
            # named for the thresholds it was cut at, so a rerun at different -c/-x would leave a
            # second final set beside the new one with nothing to say which is this run's. Keeping it
            # was also the next thing to fail - the layout step copytrees onto this path, and
            # copytree will not write into a directory that exists.
            for name in listdir(out):
                if is_refined_bins_dirname(name) and os.path.isdir(os.path.join(out, name)):
                    shutil.rmtree(os.path.join(out, name), ignore_errors=True)
                    for suffix in (STATS_SUFFIX, CONTIGS_SUFFIX):
                        sidecar = os.path.join(out, name + suffix)
                        if os.path.isfile(sidecar):
                            os.remove(sidecar)
            # And the copies of the input sets that the last run left under their input folder names.
            # Those are only ever written when the destination does not already exist - so once a
            # folder has been refined, a rerun whose inputs have changed keeps presenting the old
            # ones. That is how a refined/ directory came to advertise a 188-bin MaxBin2 set next to
            # the results of refining the corrected 129-bin one: the superseded set, under the name a
            # reader would take for this run's input, with a .stats.tsv to match.
            for _label, src in folders:
                base = os.path.basename(strip_trailing_slash(src))
                copy = os.path.join(out, base)
                # Not if it *is* the input: nothing stops someone pointing -A at a folder inside the
                # output directory, and deleting a user's bins would be a far worse bug than the one
                # being fixed.
                if os.path.realpath(copy) == os.path.realpath(src):
                    continue
                if os.path.isdir(copy):
                    shutil.rmtree(copy, ignore_errors=True)
                for suffix in (STATS_SUFFIX, CONTIGS_SUFFIX):
                    sidecar = copy + suffix
                    if os.path.isfile(sidecar):
                        os.remove(sidecar)
        ensure_dir(out)

        n_binnings = 0
        for dest, src in folders:
            if not os.path.isdir(src):
                error("%s is not a valid directory. Exiting." % src)
            kept = _copy_bins(src, os.path.join(out, dest))
            comm("there are %d bins in %s" % (kept, dest))
            if kept == 0 and not dry_run():
                error("Please provide valid input. Exiting...")
            n_binnings += 1
        comm("There are %d bin sets!" % n_binnings)

        # One CheckM score cache for the whole refinement. The hybrid sets this builds overlap heavily -
        # measured at 32% identical bins across the seven sets of a three-binner run - and a bin's score
        # depends only on its contigs, so scoring one twice is pure waste. Kept under work_files so a
        # --resume does not re-pay either.
        score_cache = _checkm.BinScoreCache(os.path.join(out, "work_files", "checkm_scores.json"))

        comm("Fix contig naming by removing special characters...")
        rewritten = 0
        for dest, _src in folders:
            for f in listdir(os.path.join(out, dest)):
                if _fix_contig_naming(os.path.join(out, dest, f)):
                    rewritten += 1
        comm("normalised contig headers in %d bin file(s); the rest were already clean" % rewritten)

        if external:
            comm(
                "Skipping bin refinement: %s consolidates the input bin sets directly, so there "
                "are no combinations to build." % args.consolidator
            )
        elif args.refine and n_binnings >= 2:
            if ckpt.todo("refine"):
                announcement("BEGIN BIN REFINEMENT")
                plan = refinement_plan(
                    [label for label, _src in folders], combinations=args.refine_combinations
                )
                comm(
                    "There are %d bin folders, so there are %d way(s) to refine them (%s). "
                    "Trying all of them!"
                    % (
                        n_binnings,
                        len(plan),
                        ", ".join(
                            "+".join(lbl[len("bins") :] for lbl in inputs) for inputs, _dest in plan
                        ),
                    )
                )
                for inputs, dest in bar(
                    plan, desc="refining bin-set combinations", unit=" combo", total=len(plan)
                ):
                    _refine(out, env, inputs, dest)
                comm("Bin refinement finished successfully!")
                ckpt.done("refine")
            else:
                comm("skipping bin refinement (already done; --resume)")
        elif not args.refine:
            comm("Skipping bin refinement. Will proceed with the %d bins specified." % n_binnings)
        else:
            comm("There is only one bin folder, so no refinement of bins possible. Moving on...")

        comm("fixing bin naming to .fa convention for consistency...")
        for name in _bin_set_dirs(out):
            _normalize_bin_extensions(os.path.join(out, name))

        comm("making sure every refined bin set contains bins...")
        for name in [] if dry_run() else _bin_set_dirs(out):
            if _count_bins(os.path.join(out, name)) == 0:
                comm("Removing bin set %s because it yielded 0 refined bins ..." % name)
                shutil.rmtree(os.path.join(out, name))

        # ── CheckM on all bin sets ──────────────────────────────────────────────────────────
        if external:
            # Deliberately not run. CheckM on every input set exists to tell this module's own
            # consolidation which version of each bin to keep; an external consolidator makes that
            # decision from its own scoring. Running it anyway would add the single largest cost in
            # the module to the arm that was supposed to be the cheap one, and compare nothing.
            comm(
                "Skipping CheckM on the input bin sets: %s scores the bins itself. The "
                "consolidated set is still scored with CheckM below, so the arms stay comparable."
                % args.consolidator
            )
        elif args.run_checkm:
            if ckpt.todo("checkm"):
                announcement("RUNNING CHECKM ON ALL SETS OF BINS")
                # CheckM on each of up to 7 bin sets; each one is many minutes.
                sets = _bin_set_dirs(out) or ([d for d, _ in folders] if dry_run() else [])
                for name in bar(sets, desc="CheckM on bin sets", unit=" set", total=len(sets)):
                    _run_checkm_set(
                        os.path.join(out, name), args, env, checkm2_env, cache=score_cache
                    )
                ckpt.done("checkm")
            else:
                comm("skipping CheckM on bin sets (already done; --resume)")
        else:
            comm("Skipping CheckM. Warning: bin consolidation will not be possible.")

        # ── consolidate ─────────────────────────────────────────────────────────────────────
        best_bin_set: Optional[str] = None
        if external:
            if ckpt.todo("consolidate"):
                announcement("CONSOLIDATING ALL BIN SETS WITH " + args.consolidator.upper())
                _run_dastool(args, external_env, out, [label for label, _src in folders])
                ckpt.done("consolidate")
            else:
                comm("skipping %s consolidation (already done; --resume)" % args.consolidator)
            best_bin_set = "binsO"
        elif args.consolidate and args.run_checkm:
            announcement("CONSOLIDATING ALL BIN SETS BY CHOOSING THE BEST VERSION OF EACH BIN")
            if n_binnings == 1:
                comm(
                    "There is only one original bin folder, so no consolidation possible. Moving on..."
                )
                best_bin_set = "binsA"
            else:
                if ckpt.todo("consolidate") and not dry_run():
                    comm("There are %d original bin folders, plus the refined bins." % n_binnings)
                    binsM = os.path.join(out, "binsM")
                    # Both of these are destination-safe for the same reason: a --resume whose
                    # consolidate checkpoint is unset got here because a previous attempt crashed
                    # partway through consolidating, and left a binsM (or, below, a binsO) behind.
                    shutil.rmtree(binsM, ignore_errors=True)
                    shutil.copytree(os.path.join(out, "binsA"), binsM)
                    shutil.copy(os.path.join(out, "binsA" + STATS_SUFFIX), binsM + STATS_SUFFIX)
                    for name in _bin_set_dirs(out):
                        if name == "binsA":
                            continue
                        stats = os.path.join(out, name + STATS_SUFFIX)
                        if not os.path.isfile(stats):
                            continue
                        comm("merging %s and binsM" % name)
                        tmp = os.path.join(out, "binsM1")
                        _refinement.consolidate(
                            binsM,
                            os.path.join(out, name),
                            binsM + STATS_SUFFIX,
                            stats,
                            tmp,
                            comp,
                            cont,
                        )
                        shutil.rmtree(binsM)
                        os.remove(binsM + STATS_SUFFIX)
                        os.replace(tmp, binsM)
                        os.replace(tmp + STATS_SUFFIX, binsM + STATS_SUFFIX)

                    binsO = os.path.join(out, "binsO")
                    # A crashed attempt at this step leaves a partial binsO behind, and dereplicate
                    # writes into the folder rather than replacing it - so its leftover bins would
                    # survive into the final set. The move in the "keep" branch would have refused
                    # the directory outright. Cleared once, for all three dereplication modes.
                    shutil.rmtree(binsO, ignore_errors=True)
                    if args.dereplicate == "keep":
                        comm("Skipping dereplication of contigs between bins...")
                        os.replace(binsM, binsO)
                        os.replace(binsM + STATS_SUFFIX, binsO + STATS_SUFFIX)
                    elif args.dereplicate == "partial":
                        comm(
                            "Scanning to find duplicate contigs between bins and only keep them in "
                            "the best bin..."
                        )
                        _refinement.dereplicate(binsM + STATS_SUFFIX, binsM, binsO, mode="best")
                        shutil.copy(binsM + STATS_SUFFIX, binsO + STATS_SUFFIX)
                    elif args.dereplicate == "remove":
                        comm(
                            "Scanning to find duplicate contigs between bins and deleting them in all "
                            "bins..."
                        )
                        _refinement.dereplicate(binsM + STATS_SUFFIX, binsM, binsO, mode="remove")
                        shutil.copy(binsM + STATS_SUFFIX, binsO + STATS_SUFFIX)
                    ckpt.done("consolidate")
                elif dry_run():
                    comm(
                        "(dry run) would consolidate the bin sets into binsM, then "
                        "dereplicate contigs into binsO"
                    )
                else:
                    comm("skipping bin consolidation (already done; --resume)")
                best_bin_set = "binsO"
        elif not args.consolidate:
            comm(
                "Skipping bin consolidation. Will pick the best binning folder without mixing bins."
            )
            if not args.run_checkm:
                comm("cannot decide on best bin set because CheckM was not run. Assuming binsA.")
                best_bin_set = "binsA"
            else:
                best, max_good = None, -1
                for name in _bin_set_dirs(out):
                    n = len(_good_bins(os.path.join(out, name + STATS_SUFFIX), comp, cont))
                    comm("There are %d 'good' bins found in %s!" % (n, name))
                    if n > max_good:
                        max_good, best = n, name
                best_bin_set = best
                comm("looks like the best bin set is %s" % best_bin_set)
        else:
            best_bin_set = "binsA"
        comm("You will find the best non-reassembled versions of the bins in %s" % best_bin_set)

        # ── finalize: re-QC the consolidated set and lay out the final outputs ───────────────
        if ckpt.todo("finalize") and not dry_run():
            # --keep-ambiguous is what skips the re-score for this module's own consolidation,
            # because that is the one mode where binsO is binsM unchanged and already scored. An
            # external consolidator's output has never been through CheckM at all, and its bins are
            # non-overlapping by construction, so the flag has nothing to say about it.
            rescore = args.run_checkm and best_bin_set == "binsO"
            if rescore and not external and args.dereplicate == "keep":
                rescore = False
            if rescore:
                binsO = os.path.join(out, "binsO")
                if _checkm.count_bins(binsO) == 0:
                    error(
                        "Consolidation produced no bins: not one bin from any bin set met "
                        "your quality thresholds of >%s%% completeness and <%s%% "
                        "contamination.\n"
                        "The per-bin-set scores are in %s/bins*.stats - look at the best "
                        "completeness there and re-run with a lower -c (and/or a higher "
                        "-x) if that is acceptable for your analysis. A very low ceiling "
                        "usually means the assembly or the coverage is too thin to recover "
                        "genomes, rather than anything wrong with the thresholds."
                        % (comp, cont, out)
                    )
                comm("Re-running CheckM on binsO bins")
                stats = _checkm.run_checkm_cached(
                    binsO,
                    cache=score_cache,
                    prefilter=not args.score_everything,
                    say=comm,
                    threads=args.threads,
                    mem_gb=args.memory,
                    env=env,
                    quick=args.quick,
                    checkm2=args.checkm2,
                    checkm2_env=checkm2_env,
                )
                comm("Removing bins that are inadequate quality...")
                good = set(_good_bins(stats, comp, cont))
                for f in list(os.listdir(binsO)):
                    if f.endswith(BIN_EXTENSION) and os.path.splitext(f)[0] not in good:
                        print(
                            "%s will be removed because it fell below the quality threshold after "
                            "de-replication of contigs..." % os.path.splitext(f)[0]
                        )
                        os.remove(os.path.join(binsO, f))
                # rewrite stats to keep only the good bins
                with open(stats) as fh:
                    lines = fh.readlines()
                with open(stats, "w") as fh:
                    fh.write(lines[0])
                    for line in lines[1:]:
                        cut = line.rstrip("\n").split("\t")
                        if (
                            len(cut) >= 3
                            and comp <= float(cut[1]) <= 100
                            and 0 <= float(cut[2]) <= cont
                        ):
                            fh.write(line)
                comm(
                    "Re-evaluating bin quality after contig de-replication is complete! There are still "
                    "%d high quality bins." % len(good)
                )

            if args.run_checkm:
                comm("making completion and contamination ranking plots of final outputs")
                # This runs with cwd=out (the plot script writes its png into the working
                # directory), so the arguments must be absolute. Passing paths built from a
                # relative -o meant "cd out && ... out/binsA.stats", i.e. out/out/binsA.stats,
                # and the step failed for every user who passed a relative output directory.
                stats_files = " ".join(
                    sorted(
                        os.path.abspath(os.path.join(out, f))
                        for f in os.listdir(out)
                        if f.endswith(STATS_SUFFIX)
                    )
                )
                ensure_dir(os.path.join(out, "figures"))
                run(
                    PLOT_BINNING.format(comp=int(comp), cont=int(cont), stats=stats_files),
                    env=None,
                    tool="plot_binning_results",
                    cwd=out,
                )
                png = os.path.join(out, "binning_results.png")
                if os.path.isfile(png):
                    os.replace(png, os.path.join(out, "figures", "binning_results.png"))

            # ── move intermediate files aside and lay out the final outputs ──────────────────
            announcement("MOVING OVER TEMPORARY FILES")
            if n_binnings != 1:
                work = os.path.join(out, "work_files")
                ensure_dir(work)
                # Every working bin set and its sidecar files, whatever the combination is
                # named. This was a hardcoded ("binsA", "binsB", "binsC", "binsM", "binsO")
                # prefix list, so with more than three inputs binsD/binsE and their .stats.tsv
                # were left behind in the output root instead of being tidied away.
                for f in list(os.listdir(out)):
                    if f.startswith("bins"):
                        # Destination-safe: the destination already exists whenever this output
                        # folder has been refined before, so this tidy-up used to raise
                        # FileExistsError at the very end of an otherwise complete run, with every
                        # final bin already written, and leave the remaining sets half-moved for the
                        # next --resume to trip over.
                        fs_replace_dir(os.path.join(out, f), os.path.join(work, f))

            final = os.path.join(out, refined_bins_dirname(comp, cont))
            if best_bin_set is None:
                # Every branch above sets this; assert it for the type checker and so that a
                # future branch that forgets fails loudly rather than writing to "None".
                error("internal error: no best bin set was chosen")
            best_dir = (
                os.path.join(out, "work_files", best_bin_set)
                if n_binnings != 1
                else os.path.join(out, best_bin_set)
            )
            if os.path.isdir(best_dir):
                # With two or more bin sets, consolidation has already dropped everything below the
                # thresholds. With exactly one there was nothing to consolidate, so the thresholds had
                # never been applied to anything and this folder - named for them - was getting every
                # input bin, 0%-complete ones included.
                single_set_stats = best_dir + STATS_SUFFIX
                if n_binnings == 1 and args.run_checkm and os.path.isfile(single_set_stats):
                    kept = _refinement.select_passing(best_dir, single_set_stats, final, comp, cont)
                    comm(
                        "%d of the %d bins pass %d%% completeness and %d%% contamination"
                        % (kept, len(fasta_files(best_dir)), int(comp), int(cont))
                    )
                else:
                    shutil.copytree(best_dir, final)
                    if os.path.isfile(single_set_stats):
                        shutil.copy(single_set_stats, final + STATS_SUFFIX)

            # restore the original input bin sets under their input folder names
            for dest, src in folders:
                base = os.path.basename(strip_trailing_slash(src))
                src_dir = (
                    os.path.join(out, "work_files", dest)
                    if n_binnings != 1
                    else os.path.join(out, dest)
                )
                if os.path.isdir(src_dir) and not os.path.exists(os.path.join(out, base)):
                    shutil.copytree(src_dir, os.path.join(out, base))
                    if os.path.isfile(src_dir + STATS_SUFFIX):
                        shutil.copy(src_dir + STATS_SUFFIX, os.path.join(out, base + STATS_SUFFIX))

            if args.run_checkm:
                comm("making contig membership files (for Anvio and other applications)")
                for name in os.listdir(out):
                    d = os.path.join(out, name)
                    if name.endswith("_bins") and os.path.isdir(d):
                        write_contig_table(d, d + CONTIGS_SUFFIX)
            ckpt.done("finalize")
        elif dry_run():
            if args.run_checkm:
                ensure_dir(os.path.join(out, "figures"))
                run(
                    PLOT_BINNING.format(
                        comp=int(comp), cont=int(cont), stats=os.path.join(out, "<bin-set>.stats")
                    ),
                    env=None,
                    tool="plot_binning_results",
                    cwd=out,
                )
            comm("(dry run) would lay out the final bins in %s" % refined_bins_dirname(comp, cont))
        else:
            comm("skipping final output layout (already done; --resume)")

        # ── GUNC: the chimerism question CheckM's marker genes cannot answer ─────────────────
        # Its own checkpoint step rather than part of "finalize", so a resume of a run that
        # finalized before --skip-gunc was dropped still gets the annotation without redoing the
        # layout. It annotates and never deletes: a chimeric bin is reported, and what to do about
        # it is the reader's call.
        final_bins = os.path.join(out, refined_bins_dirname(comp, cont))
        if args.run_gunc and ckpt.todo("gunc"):
            announcement("ASSESSING THE FINAL BINS FOR CHIMERISM WITH GUNC")
            if dry_run():
                # A placeholder rather than the empty string, so the printed command reads as
                # "this is where your database goes" instead of as a malformed --db_file.
                _gunc.run_gunc(
                    final_bins,
                    threads=args.threads,
                    db=gunc_db or "<GUNC_DB>",
                    env=gunc_env,
                    out_dir=os.path.join(out, "work_files", "gunc"),
                )
            elif _checkm.count_bins(final_bins) == 0:
                comm("no final bins to assess; skipping GUNC")
                ckpt.done("gunc")
            else:
                report = _gunc.run_gunc(
                    final_bins,
                    threads=args.threads,
                    db=gunc_db,
                    env=gunc_env,
                    out_dir=os.path.join(out, "work_files", "gunc"),
                )
                verdicts = _gunc.parse_report(str(report))
                chimeric = sorted(name for name, v in verdicts.items() if not v.passed)
                n = _gunc.annotate_stats(final_bins + STATS_SUFFIX, verdicts)
                comm(
                    "GUNC scored %d bins; %d annotated in %s"
                    % (len(verdicts), n, os.path.basename(final_bins) + STATS_SUFFIX)
                )
                if chimeric:
                    warning(
                        "GUNC calls %d of %d final bins chimeric (clade separation score above "
                        "%s): %s\nThey are kept, with gunc_pass=no in the stats table - a bin two "
                        "genomes wide can look complete and clean to CheckM."
                        % (len(chimeric), len(verdicts), _gunc.CSS_THRESHOLD, ", ".join(chimeric))
                    )
                else:
                    comm("GUNC found no chimeric bins among the final %d" % len(verdicts))
                ckpt.done("gunc")
        elif args.run_gunc:
            comm("skipping GUNC (already done; --resume)")

        # ── MIMAG tiers for the final bin set ───────────────────────────────────────────────
        # Here rather than in checkm.run_checkm because this is the module whose output is the
        # answer: the hybrid sets scored on the way in are candidates, and tagging each of them
        # would put a quality claim on bins that were discarded.
        #
        # No rRNA evidence exists at this point and that is not a defect - refinement runs before
        # annotate, so nothing has searched these bins for a 16S yet. Bins that clear 90/5 are
        # therefore tagged high_quality_unverified_rrna, which is the distinction this step is for:
        # it says "the completeness half passed and the rRNA half was not measured", where a plain
        # "medium" would report a measurement that was never made. Re-running after
        # `annotate --bakta` resolves them.
        if not dry_run() and os.path.isdir(final_bins):
            stats_path = final_bins + STATS_SUFFIX
            assessments = _checkm.tag_mimag(stats_path)
            if assessments:
                comm(_mimag.summarise(assessments))
                counts = _mimag.tally(assessments)
                if counts[_mimag.MimagTag.HIGH_UNVERIFIED]:
                    comm(_mimag.NO_EVIDENCE_HINT)

        # What the bin set does not say about itself. Both cheap, and both change what a reader
        # concludes rather than changing a bin - see metawrap2.recovery.
        if not dry_run() and args.assembly and os.path.isdir(final_bins):
            recovery = _recovery.measure(args.assembly, final_bins)
            table = os.path.join(args.output, "recovery" + TSV_SUFFIX)
            with open(table, "w", newline="") as fh:
                writer = csv.DictWriter(
                    fh,
                    fieldnames=_recovery.RECOVERY_COLUMNS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow(_recovery.recovery_row(recovery))
            comm(recovery.describe())
            comm(_recovery.UNBINNED_CAVEAT)
            comm("-> %s" % table)
        comm(_recovery.CHIMERISM_CAVEAT)

        announcement("BIN_REFINEMENT PIPELINE FINISHED SUCCESSFULLY!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
