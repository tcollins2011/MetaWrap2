"""MetaWrap2 bin_refinement module: consolidate several bin sets into one refined set.

Given any number of folders of bins produced from the same assembly (e.g. metaBAT2, MaxBin2,
CONCOCT outputs), this hybridizes them with Binning_refiner, scores every candidate set with CheckM,
consolidates the best version of each bin, and dereplicates shared contigs. The scientific
core lives in :mod:`metawrap2.refinement` (consolidate/dereplicate) and
:mod:`metawrap2.checkm` (CheckM); this module orchestrates them and preserves the original
output layout.

This file is meant to be read and edited: the commands each external tool runs live in the
COMMANDS block below and the tunable numbers in the CONSTANTS block.
"""

from __future__ import annotations

import argparse
import itertools
import os
import shlex
import shutil
from typing import Dict, List, Optional, Sequence, Tuple

from .. import checkm as _checkm
from .. import refinement as _refinement
from ..config import OPTIONAL_ENVS, load_settings
from ..constants import (
    BIN_EXTENSION,
    CONTIGS_SUFFIX,
    FASTA_EXTENSIONS,
    MAX_BIN_SIZE,
    MIN_BIN_SIZE,
    STATS_SUFFIX,
)
from ..io.seqio import contig_id, iter_fasta
from ..progress import bar
from ..pyrun import PY
from ..utils import strip_trailing_slash
from ._common import (
    absolutize_paths,
    announcement,
    comm,
    dry_run,
    ensure_dir,
    env_for,
    error,
    finish_run,
    listdir,
    make_checkpoint,
    require_nonempty_dir,
    resolve_threads,
    run,
    start_run,
    threads_arg,
    validate_inputs,
    warning,
)

CONDA_ENV = "metawrap2-bin_refinement"

# ─── COMMANDS (edit flags here) ──────────────────────────────────────────────────────────
# These are MetaWrap2's own helpers, so they run in the host interpreter (env=None), not in
# the module's conda env - that env holds CheckM, not python+metawrap2. See metawrap2.pyrun.
BINNING_REFINER = PY + " -m metawrap2.vendor.binning_refiner -1 {b1} -2 {b2} -o {out}"
BINNING_REFINER_3 = PY + " -m metawrap2.vendor.binning_refiner -1 {b1} -2 {b2} -3 {b3} -o {out}"
BINNING_REFINER_N = PY + " -m metawrap2.vendor.binning_refiner {bins} -o {out}"
PLOT_BINNING = PY + " -m metawrap2.scripts.plot_binning_results {comp} {cont} {stats}"
# ─── CONSTANTS ───────────────────────────────────────────────────────────────────────────
DEFAULT_MEM_GB = 40
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
# ───────────────────────────────────────────────────────────────────────────────────────────


def _parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="metawrap2 bin_refinement",
        usage="metawrap2 bin_refinement [options] -o output_dir -A bin_folderA "
        "[-B bin_folderB ... -F bin_folderF]  |  --bins DIR [--bins DIR ...]",
        description="Consolidate up to %d bin sets (from the same assembly) into a refined "
        "set. Contig names must be consistent across the bin folders." % MAX_BIN_SETS,
    )
    p.add_argument("-o", "--output", required=True, help="output directory")
    threads_arg(p)
    p.add_argument(
        "-m", "--memory", type=int, default=DEFAULT_MEM_GB, help="RAM in GB (default 40)"
    )
    p.add_argument(
        "-c",
        "--completeness",
        type=float,
        default=DEFAULT_MIN_COMPLETION,
        help="minimum %% completion of bins [should be >50] (default 70)",
    )
    p.add_argument(
        "-x",
        "--contamination",
        type=float,
        default=DEFAULT_MAX_CONTAM,
        help="maximum acceptable %% contamination of bins (default 10)",
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
        "--skip-checkm",
        dest="run_checkm",
        action="store_false",
        help="don't run CheckM to assess bins",
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
        "wall; needs the metawrap2-bin_refinement-checkm2 env)",
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
    """Copy size-filtered bins from *src* into *dest* as bin-name.fa. Returns count kept."""
    ensure_dir(dest)
    n = 0
    for name in sorted(os.listdir(src)):
        full = os.path.join(src, name)
        if not os.path.isfile(full):
            continue
        size = os.path.getsize(full)  # portable (no `stat -c`)
        if MIN_SIZE < size < MAX_SIZE:
            base = os.path.splitext(name)[0]
            shutil.copy(full, os.path.join(dest, base + BIN_EXTENSION))
            n += 1
        else:
            comm("Skipping %s because the bin size is not between 50kb and 20Mb" % full)
    return n


def _count_bins(path: str) -> int:
    return (
        len([f for f in os.listdir(path) if f.endswith(BIN_EXTENSION)])
        if os.path.isdir(path)
        else 0
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
    os.replace(produced, os.path.join(out, dest))
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
    bin_set_dir: str, args, env: Optional[str], checkm2_env: Optional[str] = None
) -> str:
    """Run CheckM on one bin set and report the number of 'good' bins. Returns the .stats path."""
    comm(
        "Running %s on %s bins"
        % ("CheckM2" if args.checkm2 else "CheckM", os.path.basename(bin_set_dir))
    )
    stats = _checkm.run_checkm(
        bin_set_dir,
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
    absolutize_paths(args)
    env = env_for("bin_refinement", settings)
    # CheckM2 lives in its own opt-in env; everything else in this module still uses the
    # default bin_refinement env.
    checkm2_env = OPTIONAL_ENVS["checkm2"] if (args.checkm2 and settings.use_conda_envs) else None
    bin_sets = collect_bin_sets(args)
    validate_inputs(
        "bin_refinement",
        [
            (
                lambda path, what: require_nonempty_dir(path, what, FASTA_EXTENSIONS),
                path,
                "bin set %s" % bin_set_label(i),
            )
            for i, path in enumerate(bin_sets)
        ],
    )
    rec = start_run("bin_refinement", args, env, settings, inputs=list(bin_sets))

    ckpt = make_checkpoint(args.output)
    try:

        comp, cont = args.completeness, args.contamination
        folders = [(bin_set_label(i), src) for i, src in enumerate(bin_sets)]
        for _label, src in folders:
            if not os.path.isdir(src):
                error("%s is not a valid directory. Exiting." % src)

        announcement("BEGIN PIPELINE!")
        out = args.output
        comm("setting up output folder and copying over bins...")
        if os.path.isdir(out):
            warning("Warning: %s already exists. Attempting to clean." % out)
            # Remove every bins* working directory, whatever the combination names are for
            # this number of inputs, plus the consolidated (M) and dereplicated (O) sets.
            for name in listdir(out):
                if name.startswith("bins") and os.path.isdir(os.path.join(out, name)):
                    shutil.rmtree(os.path.join(out, name), ignore_errors=True)
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

        comm("Fix contig naming by removing special characters...")
        rewritten = 0
        for dest, _src in folders:
            for f in listdir(os.path.join(out, dest)):
                if _fix_contig_naming(os.path.join(out, dest, f)):
                    rewritten += 1
        comm("normalised contig headers in %d bin file(s); the rest were already clean" % rewritten)

        if args.refine and n_binnings >= 2:
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
        for name in ([] if dry_run() else _bin_set_dirs(out)):
            if _count_bins(os.path.join(out, name)) == 0:
                comm("Removing bin set %s because it yielded 0 refined bins ..." % name)
                shutil.rmtree(os.path.join(out, name))

        # ── CheckM on all bin sets ──────────────────────────────────────────────────────────
        if args.run_checkm:
            if ckpt.todo("checkm"):
                announcement("RUNNING CHECKM ON ALL SETS OF BINS")
                # CheckM on each of up to 7 bin sets; each one is many minutes.
                sets = _bin_set_dirs(out) or ([d for d, _ in folders] if dry_run() else [])
                for name in bar(sets, desc="CheckM on bin sets", unit=" set", total=len(sets)):
                    _run_checkm_set(os.path.join(out, name), args, env, checkm2_env)
                ckpt.done("checkm")
            else:
                comm("skipping CheckM on bin sets (already done; --resume)")
        else:
            comm("Skipping CheckM. Warning: bin consolidation will not be possible.")

        # ── consolidate ─────────────────────────────────────────────────────────────────────
        best_bin_set: Optional[str] = None
        if args.consolidate and args.run_checkm:
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
            if args.run_checkm and best_bin_set == "binsO" and args.dereplicate != "keep":
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
                stats = _checkm.run_checkm(
                    binsO,
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
                        os.replace(os.path.join(out, f), os.path.join(work, f))

            final = os.path.join(out, "metawrap_%d_%d_bins" % (int(comp), int(cont)))
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
                shutil.copytree(best_dir, final)
                if os.path.isfile(best_dir + STATS_SUFFIX):
                    shutil.copy(best_dir + STATS_SUFFIX, final + STATS_SUFFIX)

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
                        with open(d + CONTIGS_SUFFIX, "w") as out_fh:
                            for f in sorted(os.listdir(d)):
                                if f.endswith(BIN_EXTENSION):
                                    bin_name = os.path.splitext(f)[0]
                                    out_fh.writelines(
                                        "%s\t%s\n" % (contig_id(header), bin_name)
                                        for header, _seq in iter_fasta(os.path.join(d, f))
                                    )
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
            comm(
                "(dry run) would lay out the final bins in metawrap_%d_%d_bins"
                % (int(comp), int(cont))
            )
        else:
            comm("skipping final output layout (already done; --resume)")

        announcement("BIN_REFINEMENT PIPELINE FINISHED SUCCESSFULLY!")
    finally:
        finish_run(rec, args.output)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
