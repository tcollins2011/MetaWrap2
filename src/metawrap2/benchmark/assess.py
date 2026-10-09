"""Turning a variant's bins into numbers, at score time rather than inside the measured run.

A binning comparison leaves behind directories of FASTA files. Nothing about those files says how
good they are, so something has to assess them, and where that assessment happens decides what the
cost column means.

It happens here - after the run, outside the measurement - for one reason: the number MetaWrap2
ranks defaults on is the cost of the thing being chosen. ``binner_set`` asks which binners should
run by default "given each one costs a full pass over the reads", so the cost of a variant is the
cost of binning. Folding CheckM into the measured command would add a charge that the real pipeline
does not levy at that step (``metawrap2 binning`` does not run CheckM unless asked; refinement
does), and would charge it unevenly, because CheckM's cost scales with how many bins a variant
happened to emit. A variant would then be penalised for being productive.

Two assessors, because they answer different questions and the plan asks for both:

- **AMBER** (Meyer et al., *Gigascience* 2018) compares bins against the dataset's known
  contig-to-genome truth. It is marker-gene independent, which matters because several modern
  binners use single-copy markers internally - scoring those with a marker-based tool asks the tool
  to grade its own homework.
- **CheckM** gives the completeness/contamination tiers papers report (NC, HQ, MQ), so a MetaWrap2
  number can be put next to a published one. It does not know the truth and cannot notice a bin
  assembled from halves of two genomes, which is exactly why AMBER is not redundant with it.

Neither is run unless the comparison's declared metrics ask for it, so scoring a trimmer comparison
does not start a CheckM run.
"""

from __future__ import annotations

import os
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .. import checkm as _checkm
from .. import config as _config
from .. import gunc as _gunc
from ..command import run as _run
from ..constants import BIN_EXTENSION, is_refined_bins_dirname
from ..io.seqio import fasta_files

__all__ = [
    "AMBER_ENV_KEY",
    "AMBER_ENV_YAML",
    "AMBER_RESULTS",
    "AssessError",
    "amber_env",
    "assess_bins",
    "bin_dirs",
    "bins_to_bioboxes",
    "bioboxes_sample_ids",
    "checkm_envs",
    "gold_sample_id",
    "gold_to_bioboxes",
    "gunc_settings",
    "needs_amber",
    "needs_checkm",
    "needs_gunc",
    "run_amber",
    "run_bin_quality",
    "run_chimerism",
]

#: The key under which the AMBER env is registered in :data:`metawrap2.config.OPTIONAL_ENVS`. Its
#: own env, not a module's: AMBER is a benchmark-only tool and has no business in the environment of
#: a module a user runs.
AMBER_ENV_KEY = "amber"

#: What to type after ``metawrap2 install-env``, which is the yaml's stem and not the key above. The
#: two differ, and printing the key in an error message sends people to a command that does not
#: exist ("No env file for module 'amber'").
AMBER_ENV_YAML = "benchmark-amber"

#: AMBER's summary file, and the name :func:`~metawrap2.commands.benchmark._collect_scores` already
#: knows how to read. Written under a per-variant subdirectory so several bin sets can be scored
#: without overwriting each other.
AMBER_RESULTS = "results.tsv"

#: Sentinel for "the caller said nothing about the env", so that passing None explicitly still means
#: "run whatever is on PATH" rather than being indistinguishable from not passing anything.
_UNSET = "\0unset"

#: Bioboxes binning-format header. AMBER reads this format and nothing else, so both the gold
#: standard and every query set have to be converted into it.
_BIOBOXES_HEADER = "@Version:0.9.1\n@SampleID:%s\n\n@@SEQUENCEID\tBINID\n"

#: The gold standard's header, which carries a third column. AMBER's bp-weighted metrics - the ones
#: worth ranking on, since a 5 Mb genome recovered is not the same result as a 50 kb one - need a
#: length per sequence, and it takes them from the gold standard only (a query is read as id and bin
#: id alone). Without the column it dies with ``KeyError: "['LENGTH'] not in index"``.
_BIOBOXES_GOLD_HEADER = "@Version:0.9.1\n@SampleID:%s\n\n@@SEQUENCEID\tBINID\t_LENGTH\n"

#: Metrics that need AMBER, and metrics that need CheckM. Used to decide whether to run an assessor
#: at all: a comparison declares what it is ranked on, and anything it does not declare is not
#: worth an hour of somebody's machine.
_AMBER_METRICS = ("amber_f1_bp", "amber_ari_bp", "amber_completeness_bp", "amber_purity_bp")
_CHECKM_METRICS = ("bins_near_complete", "bins_high_quality", "bins_medium_quality", "bins_total")
_GUNC_METRICS = ("gunc_pass_fraction",)


def amber_env() -> Optional[str]:
    """The conda env to run AMBER in, or None when the config says not to use conda envs.

    Resolved at call time rather than at import, because ``use_conda_envs`` is a setting and a
    setting read at import is a setting that ignores the config file the user just edited.
    """
    settings = _config.load_settings()
    if not settings.use_conda_envs:
        return None
    return _config.OPTIONAL_ENVS[AMBER_ENV_KEY]


#: The module whose env carries CheckM1. Scoring is not that module, but it runs that module's tool,
#: and a benchmark env holding a second copy of CheckM would be a second thing to keep in step.
CHECKM_MODULE = "bin_refinement"


def checkm_envs() -> Tuple[Optional[str], Optional[str]]:
    """The envs to run CheckM1 and CheckM2 in, or ``(None, None)`` when conda envs are off.

    Resolved here rather than left to the caller because leaving it unset means running whatever
    ``checkm`` is on PATH, and on a machine that installed MetaWrap2 the way the docs say, that is
    nothing at all: scoring then reported ``No such file or directory: 'checkm'`` for a tool it had
    installed itself.
    """
    settings = _config.load_settings()
    if not settings.use_conda_envs:
        return None, None
    return _config.MODULE_ENVS.get(CHECKM_MODULE), _config.OPTIONAL_ENVS.get("checkm2")


class AssessError(RuntimeError):
    """Raised when a variant's bins cannot be assessed, with the reason and what to do about it."""


def needs_amber(metrics: Iterable[str]) -> bool:
    """Do any of *metrics* come from AMBER?"""
    return any(m in _AMBER_METRICS for m in metrics)


def needs_checkm(metrics: Iterable[str]) -> bool:
    """Do any of *metrics* come from CheckM?"""
    return any(m in _CHECKM_METRICS for m in metrics)


def needs_gunc(metrics: Iterable[str]) -> bool:
    """Do any of *metrics* come from GUNC?

    Worth asking separately from CheckM even though both score bins, because GUNC answers a question
    neither CheckM nor AMBER can. CheckM counts single-copy markers, so a bin that is 95% one genome
    and 5% another reads as barely contaminated; AMBER weights by base pairs, so the same bin costs
    almost nothing in purity. GUNC asks whether the genes agree on a lineage, and fails that bin
    outright. On this project's own measurements a 1.3-point AMBER purity gap between two bin sets
    came with a 2.8-fold difference in chimeric bins, which is the gap this metric exists to see.
    """
    return any(m in _GUNC_METRICS for m in metrics)


def bin_dirs(root: str) -> List[str]:
    """Every bin directory under *root*, in a stable order.

    A binning run names its output ``<binner>_bins``, one directory per binner, so a variant that
    ran three binners leaves three. Refinement's output is ``metawrap_<c>_<x>_bins``, which matches
    the same suffix, so this finds it too without a second rule. Which of them is the variant's own
    answer is a separate question - see :func:`metawrap2.constants.is_refined_bins_dirname` and the
    caller that uses it.

    A consolidated set one level down is found too, and that is not a convenience. A binning
    comparison's variants are refined as a second step, which writes into a subdirectory of the
    variant - so with a flat scan every multi-binner arm reported "left several raw bin sets and no
    consolidated one" and went unranked, when its consolidated set was sitting right there. Six arms of
    the binner comparison, including every ensemble worth having, were in that state.

    Only refined names are taken from a subdirectory. Refinement restores its input sets beside its
    output under their original names, so accepting any ``*_bins`` from one level down would pull in a
    second copy of the raw sets and double-count every genome in them.
    """
    if not os.path.isdir(root):
        return []
    found = []
    for name in sorted(os.listdir(root)):
        path = os.path.join(root, name)
        if not os.path.isdir(path):
            continue
        if name.endswith("_bins") and _checkm.count_bins(path):
            found.append(path)
            continue
        for inner in sorted(os.listdir(path)):
            nested = os.path.join(path, inner)
            if (
                inner.endswith("_bins")
                and is_refined_bins_dirname(inner)
                and os.path.isdir(nested)
                and _checkm.count_bins(nested)
            ):
                found.append(nested)
    return found


def _fasta_ids(path: str) -> List[str]:
    """Sequence ids in a FASTA file, cut at the first whitespace.

    Cut because a bin's FASTA header keeps whatever description the assembler wrote, while the gold
    standard names the contig alone. Comparing the two without trimming silently scores every
    contig as unknown - AMBER reports that as a purity of zero rather than as an error, which is
    the sort of wrong number that gets believed.
    """
    ids = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split(None, 1)[0].strip())
    return ids


def fasta_lengths(path: str) -> Dict[str, int]:
    """Sequence id to sequence length for every record in the FASTA at *path*.

    Read from the assembly the binners were given, because that is the sequence set the gold standard
    describes - lengths taken from the bins instead would only cover the contigs somebody binned, and
    a completeness denominator built from those is the fraction of what was recovered that was
    recovered, which is 1 by construction.
    """
    lengths: Dict[str, int] = {}
    current = None
    total = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = total
                current = line[1:].split(None, 1)[0].strip()
                total = 0
            elif current is not None:
                total += len(line.strip())
    if current is not None:
        lengths[current] = total
    return lengths


def bins_to_bioboxes(bins_dir: str, out_path: str, *, sample_id: str = "benchmark") -> str:
    """Write *bins_dir* as a bioboxes binning file at *out_path*. Returns *out_path*.

    Each bin's FASTA name (without extension) becomes the bin id, which keeps AMBER's output
    readable against the files on disk: ``bin_007`` in the table is ``bin_007.fasta`` in the
    directory. Contigs a binner left unbinned are simply absent, which is what the format means by
    unbinned - writing them under a bin id named "unbinned" would score the leftovers as a genome.
    """
    # fasta_files is what drops the unbinned catch-all, and it is shared with the modules so that the
    # scorer and the pipeline agree on what a bin folder contains. The dotfile check stays here: an
    # editor's ``.bin_001.fasta.swp`` is not something the pipeline writes, but it is something a
    # directory being looked at by a person has in it.
    bins = [
        os.path.basename(p)
        for p in fasta_files(bins_dir, (BIN_EXTENSION,))
        if not os.path.basename(p).startswith(".")
    ]
    if not bins:
        raise AssessError("no %s files in %s" % (BIN_EXTENSION, bins_dir))
    with open(out_path, "w", encoding="utf-8") as out:
        out.write(_BIOBOXES_HEADER % sample_id)
        for name in bins:
            bin_id = name[: -len(BIN_EXTENSION)]
            out.writelines(
                "%s\t%s\n" % (seq_id, bin_id) for seq_id in _fasta_ids(os.path.join(bins_dir, name))
            )
    return out_path


def _declares_lengths(head: str) -> bool:
    """Does a bioboxes header already carry a length column, under either of its two spellings?"""
    for line in head.splitlines():
        if line.startswith("@@"):
            columns = [c.strip().upper() for c in line[2:].split("\t")]
            return "LENGTH" in columns or "_LENGTH" in columns
    return False


#: Columns of a CAMI2 contig mapping, past the two every mapping has.
_CAMI_CONTIG_ID, _CAMI_START, _CAMI_END = 3, 5, 6


def _cami_source_contig(fields: List[str]) -> Optional[str]:
    """The name a CAMI2 gold-standard assembly gives the contig described by *fields*, or None.

    ``gsa.fa`` names a contig after the region it was cut from -
    ``CP013019.1_from_3_to_4351861_total_4351859`` - while the mapping keeps those three facts in
    separate columns and leads with an anonymous id instead. The total is inclusive of both ends,
    which is why it is ``end - start + 1`` and not the difference.

    None when the row is not that shape, so a mapping from another dataset falls through to being
    reported as unmatched rather than turned into a plausible-looking wrong name.
    """
    if len(fields) <= _CAMI_END:
        return None
    try:
        start = int(float(fields[_CAMI_START]))
        end = int(float(fields[_CAMI_END]))
    except ValueError:
        return None
    contig_id = fields[_CAMI_CONTIG_ID].strip()
    if not contig_id:
        return None
    return "%s_from_%d_to_%d_total_%d" % (contig_id, start, end, end - start + 1)


def contig_to_genome(path: str, lengths: Dict[str, int]) -> Tuple[Dict[str, str], int]:
    """Parse a gold-standard mapping into contig-to-genome, plus how many contigs went unmatched.

    Shared by the two things that need this truth: :func:`gold_to_bioboxes`, which writes it out for
    AMBER, and :mod:`metawrap2.benchmark.abundance`, which needs to know which genome a bin mostly
    is. Both have to resolve contig names the same way or they would disagree about the same dataset.

    Only contigs *lengths* knows are kept - that is the assembly the binners were given, and a contig
    no binner could have seen does not belong in anybody's denominator. Insertion order is the file's
    order, and the first of a repeated contig wins: CAMI2 publishes both a raw mapping that repeats a
    contig once per source region and a de-duplicated one, and the two should give the same answer.

    A CAMI2 mapping's first column is an *anonymous* contig id (``PC0``) that the gold-standard
    assembly does not use: ``gsa.fa`` names the same contig by where it came from,
    ``CP013019.1_from_3_to_4351861_total_4351859``, which the mapping spells out across its
    ``contig_id``, ``start_position`` and ``end_position`` columns. On the gastrointestinal
    challenge, not one of the 210,618 anonymous ids appears in the 210,633-contig assembly, so keying
    on the first column alone matched nothing at all. Reconstructing the name matches all 210,618.
    Both spellings are tried per row, first column first, so a dataset that does use the anonymous
    ids is unaffected.
    """
    mapping: Dict[str, str] = {}
    missing = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line[0] in "#@":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            contig, genome = fields[0].strip(), fields[1].strip()
            if not contig or not genome:
                continue
            if not lengths.get(contig):
                # The anonymous id is not what the assembly calls this contig; try the name the
                # remaining columns spell out. See this function's docstring.
                source = _cami_source_contig(fields)
                if source is None or not lengths.get(source):
                    missing += 1
                    continue
                contig = source
            mapping.setdefault(contig, genome)
    return mapping, missing


def gold_to_bioboxes(
    path: str,
    out_path: str,
    *,
    sample_id: str = "benchmark",
    lengths: Optional[Dict[str, int]] = None,
) -> str:
    """Convert a gold-standard mapping at *path* into bioboxes format at *out_path*.

    Handles the two shapes the datasets ship. CAMI2's per-body-site mapping is a ``#``-commented
    table whose first two columns are the anonymous contig id and the genome it came from; the
    remaining columns are taxonomy and coordinates. A file that is already bioboxes *and* already
    carries a length column is returned as-is rather than copied, so a dataset that ships a complete
    gold standard is used exactly as published.

    *lengths* is contig id to length, from the assembly the binners were given. AMBER's bp-weighted
    metrics need a length per sequence and read it from the gold standard alone, so a bioboxes gold
    standard without that column is rewritten with it rather than passed through - otherwise AMBER
    exits on a ``KeyError`` deep inside pandas, which says nothing about what is missing. A contig
    the assembly does not mention is dropped: it is a contig no binner could have seen, and keeping
    it with a length of zero would leave it in the denominator of every recall.

    The parse itself - including CAMI2's anonymous contig ids, which its own gold assembly does not
    use, and the repeated rows its raw mapping carries - is :func:`contig_to_genome`, shared with the
    abundance scoring so the two cannot end up disagreeing about the same dataset.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        head = fh.read(4096)
    already_bioboxes = "@@SEQUENCEID" in head
    if already_bioboxes and _declares_lengths(head):
        return path
    if lengths is None:
        raise AssessError(
            "%s carries no sequence lengths and none were supplied, so AMBER cannot compute the "
            "bp-weighted metrics it is being asked for. Lengths come from the assembly the binners "
            "were given; pass it with --assembly." % path
        )
    if already_bioboxes:
        # Keep the sample id the dataset published rather than imposing the default: it is the name
        # the rest of that dataset's files use, and it is what shows up in AMBER's own report.
        declared = bioboxes_sample_ids(path)
        sample_id = declared[0] if declared else sample_id

    mapping, missing = contig_to_genome(path, lengths)
    with open(out_path, "w", encoding="utf-8") as out:
        out.write(_BIOBOXES_GOLD_HEADER % sample_id)
        out.writelines(
            "%s\t%s\t%d\n" % (contig, genome, lengths[contig]) for contig, genome in mapping.items()
        )
    if len(mapping) == 0:
        raise AssessError(
            "%s yielded no contig-to-genome pairs the assembly also contains, so there is nothing "
            "to score against (%d contig(s) were named by the gold standard but absent from the "
            "assembly). Expected either a bioboxes binning file or a tab-separated table whose "
            "first two columns are the contig id and its genome, over the same contigs as the "
            "assembly the binners were given." % (path, missing)
        )
    return out_path


def bioboxes_sample_ids(path: str) -> List[str]:
    """Every ``@SampleID`` a bioboxes file declares, in the order it declares them."""
    ids = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("@@"):
                continue
            if line.lower().startswith("@sampleid:"):
                value = line.split(":", 1)[1].strip()
                if value and value not in ids:
                    ids.append(value)
    return ids


def gold_sample_id(gold: str) -> str:
    """The sample id a query has to declare to be scored against the gold standard *gold*.

    AMBER pairs a query with the truth by sample id and refuses anything it cannot pair, so this is
    not cosmetic: naming a query after its bin directory makes every query unpairable
    (``Sample ID concoct_bins ... not found in the genome binning gold standard``). The gold standard
    names the sample; the queries have to agree with it, not the other way round.
    """
    ids = bioboxes_sample_ids(gold)
    if not ids:
        raise AssessError(
            "%s declares no @SampleID, so AMBER has no sample to pair a query with. Add a line "
            "'@SampleID:<name>' above the @@SEQUENCEID header." % gold
        )
    if len(ids) > 1:
        raise AssessError(
            "%s declares %d samples (%s). A variant's bins are one co-assembly and can only be "
            "scored against one sample's truth, so scoring against this file would silently count "
            "the other samples' contigs as misassigned. Pool the gold standard into a single "
            "@SampleID section first." % (gold, len(ids), ", ".join(ids))
        )
    return ids[0]


def run_amber(
    bin_dirs_: Sequence[str],
    gold: str,
    out_dir: str,
    *,
    labels: Optional[Sequence[str]] = None,
    threads: int = 1,
    env: Optional[str] = _UNSET,
) -> str:
    """Run AMBER on *bin_dirs_* against the gold standard *gold*. Returns *out_dir*.

    All of a variant's bin sets go into one AMBER invocation rather than one each, because AMBER's
    own run is the expensive part and it is happy to take several queries at once - and because its
    combined output puts the sets next to each other, which is the comparison being made. The sets
    are told apart by ``-l``, which is what AMBER labels its rows with; the sample id is the gold
    standard's and is the same in every query, because that is how AMBER finds the truth at all.
    """
    env = amber_env() if env is _UNSET else env
    os.makedirs(out_dir, exist_ok=True)
    sample_id = gold_sample_id(gold)

    names = list(labels or [os.path.basename(p) for p in bin_dirs_])
    if labels is not None and len(names) != len(bin_dirs_):
        raise AssessError(
            "run_amber got %d labels for %d bin sets; AMBER matches them by position, so a "
            "mismatch silently relabels the results" % (len(names), len(bin_dirs_))
        )
    # The converted query used to be named for the bin directory, and bin_refinement names its final
    # directory for the thresholds it was cut at - so scoring two arms whose sets were both cut at
    # -c 70 -x 10 wrote metawrap_70_10_bins.binning twice, the second over the first, and handed AMBER
    # the same file under both labels. It reported two identical rows to four decimal places for bin
    # sets of 74 and 60 bins, and nothing in its output said the queries were the same file. Named for
    # the label instead, which is unique by construction because it is what tells the rows apart.
    if len(set(names)) != len(names):
        duplicates = sorted({n for n in names if names.count(n) > 1})
        raise AssessError(
            "run_amber was given the same label more than once (%s); the labels name the output "
            "rows, so duplicates cannot be told apart" % ", ".join(duplicates)
        )

    queries = []
    for path, name in zip(bin_dirs_, names):
        query = os.path.join(out_dir, name.replace(os.sep, "_") + ".binning")
        bins_to_bioboxes(path, query, sample_id=sample_id)
        queries.append(query)
    if not queries:
        raise AssessError("no bins to score")

    cmd = [
        "amber.py",
        "-g",
        gold,
        "-l",
        ",".join(names),
        "-o",
        out_dir,
        # AMBER's HTML report wants a browser to read it and takes minutes to render; the TSV is
        # what gets parsed, and --stdout keeps the numbers in the log too.
        "--silent",
        # Don't score the gold standard against itself. It is a row of perfect scores labelled
        # "Gold standard", and the scorer reads every labelled row as a bin set - so leaving it in
        # files a 1.0 F1 under a name no directory on disk has, which would win any ranking it
        # entered.
        "--skip_gs",
    ] + queries
    _run(
        cmd,
        env=env,
        tool="amber.py",
        hint="AMBER is a benchmark-only dependency. Install it with "
        "'metawrap2 install-env %s' (conda package cami-amber)." % AMBER_ENV_YAML,
    )
    return out_dir


def run_bin_quality(
    bin_dirs_: Sequence[str],
    *,
    threads: int,
    memory_gb: int,
    checkm2: bool = False,
    env: Optional[str] = None,
    checkm2_env: Optional[str] = None,
) -> List[str]:
    """Run CheckM over each of *bin_dirs_*, returning the stats files written.

    One run per bin set, because completeness is a property of a bin set and two sets sharing a
    directory would be scored as one. The per-set tables are what the scorer reads; it is also the
    only form in which a three-binner variant's numbers can be read without silently counting the
    same genome three times.

    This is the call site where :func:`metawrap2.checkm.run_checkm_cached`'s pre-filter earns anything,
    and it took a measurement to find that out. bin_refinement scores bins that ``copy_bins`` has
    already held to 50 kb - 20 Mb, so inside refinement the filter is a no-op: of 10,168 scored bins
    there, none is over the ceiling and 56 (0.02% of the sequence) fall under the floor. Here the bin
    sets are whatever the binner emitted, and metaBAT2 emits its unbinnable remainder as a bin: of 9,019
    raw bins, 64 are over 20 Mb and they hold **11.8% of all the sequence**, which is what CheckM's cost
    follows. The cache helps wherever a bin recurs, which across a variant's own sets it does not, so one
    is created per call rather than shared.
    """
    written = []
    for path in bin_dirs_:
        written.append(
            _checkm.run_checkm_cached(
                path,
                cache=_checkm.BinScoreCache(path + ".checkm_scores.json"),
                threads=threads,
                mem_gb=memory_gb,
                env=env,
                checkm2=checkm2,
                checkm2_env=checkm2_env,
                binner=os.path.basename(path)[: -len("_bins")] or None,
                log_path=path + ".checkm.log",
            )
        )
    return written


def gunc_settings() -> Tuple[Optional[str], str]:
    """The env to run GUNC in and the reference database to run it against.

    Both come from the config for the same reason :func:`checkm_envs` does: scoring runs the module's
    own tool, and a benchmark holding its own copy of GUNC and its own copy of a multi-gigabyte
    reference database would be two more things to keep in step with what the pipeline actually uses.
    """
    settings = _config.load_settings()
    # GUNC lives in bin_refinement's own env. It used to have one of its own, and this line looked it
    # up by that key - so consolidating the envs turned scoring into a KeyError that only appeared on
    # a variant whose bins had not already been assessed. Resolved through MODULE_ENVS so there is one
    # place the answer comes from.
    env = _config.MODULE_ENVS[CHECKM_MODULE] if settings.use_conda_envs else None
    return env, settings.db("GUNC_DB") or ""


def run_chimerism(
    bin_dirs_: Sequence[str],
    *,
    threads: int,
    db: str,
    env: Optional[str],
    reassess: bool = False,
    say=print,
) -> List[str]:
    """Run GUNC over each of *bin_dirs_*, returning the reports written or already present.

    One run per bin set, like CheckM, and for the same reason: a chimera rate is a property of a set,
    and the report goes to ``<bin set>.gunc`` so the scorer can tell which set each one describes.

    A set whose report is already there is left alone unless *reassess*. That skip matters more here
    than anywhere else in this module - GUNC's diamond stage measured 3.2 to 5.7 CPU-hours per bin
    set on this project's benchmark arms, so re-reading a table to add one metric must not cost an
    afternoon of somebody's machine.
    """
    written = []
    for path in bin_dirs_:
        out_dir = path + ".gunc"
        existing = _gunc.find_report(out_dir) if os.path.isdir(out_dir) else None
        if existing and not reassess:
            say("    %s already has a GUNC report, skipping" % os.path.basename(path))
            written.append(existing)
            continue
        report = _gunc.run_gunc(
            path,
            threads=threads,
            db=db,
            env=env,
            out_dir=out_dir,
            log_path=path + ".gunc.log",
        )
        if report is not None:
            written.append(report)
    return written


def assess_bins(
    root: str,
    *,
    metrics: Sequence[str],
    gold_binning: Optional[str],
    assembly: Optional[str] = None,
    threads: int,
    memory_gb: int,
    checkm2: bool = False,
    checkm_env: Optional[str] = _UNSET,
    checkm2_env: Optional[str] = _UNSET,
    reassess: bool = False,
    say=print,
) -> Tuple[List[str], Dict[str, str]]:
    """Assess every bin set under *root* with whichever assessors *metrics* calls for.

    Returns the bin directories found and a map of assessor to what it produced or why it did not
    run. The reasons are returned rather than raised: a comparison is five variants deep by the time
    anyone scores it, and losing four variants' numbers because the first one has no bins would be
    a worse outcome than saying so per variant.

    *assembly* is the contig file the binners were given; AMBER needs it for the sequence lengths its
    bp-weighted metrics are computed over, and only when the gold standard does not carry them
    itself.
    """
    dirs = bin_dirs(root)
    outcome: Dict[str, str] = {}
    if not dirs:
        outcome["bins"] = "no non-empty *_bins directory under %s" % root
        return dirs, outcome

    if needs_amber(metrics):
        if not gold_binning:
            outcome["amber"] = "the dataset declares no gold-standard binning to score against"
        else:
            out_dir = os.path.join(root, "amber")
            os.makedirs(out_dir, exist_ok=True)
            try:
                lengths = fasta_lengths(assembly) if assembly and os.path.isfile(assembly) else None
                gold = gold_to_bioboxes(
                    gold_binning, os.path.join(out_dir, "gold.binning"), lengths=lengths
                )
                run_amber(dirs, gold, out_dir, threads=threads)
                outcome["amber"] = out_dir
            except Exception as exc:  # noqa: BLE001 - reported per variant, not fatal
                outcome["amber"] = "failed: %s" % exc
                say("    AMBER failed on %s: %s" % (root, exc))

    if needs_checkm(metrics):
        if checkm_env is _UNSET or checkm2_env is _UNSET:
            default_env, default_checkm2_env = checkm_envs()
            checkm_env = default_env if checkm_env is _UNSET else checkm_env
            checkm2_env = default_checkm2_env if checkm2_env is _UNSET else checkm2_env
        try:
            paths = run_bin_quality(
                dirs,
                threads=threads,
                memory_gb=memory_gb,
                checkm2=checkm2,
                env=checkm_env,
                checkm2_env=checkm2_env,
            )
            outcome["checkm"] = ", ".join(paths)
        except Exception as exc:  # noqa: BLE001 - reported per variant, not fatal
            outcome["checkm"] = "failed: %s" % exc
            say("    CheckM failed on %s: %s" % (root, exc))

    if needs_gunc(metrics):
        gunc_env, gunc_db = gunc_settings()
        if not (gunc_db and os.path.exists(gunc_db)):
            # Said rather than raised, and said with the same hint the module uses: a reference database
            # that was never fetched is the likeliest reason this metric is missing, and finding that
            # out from the scorer is better than from a diamond run that dies an hour in.
            outcome["gunc"] = "no GUNC database at %r. %s" % (gunc_db or "<not set>", _gunc.DB_HINT)
            say("    GUNC skipped on %s: %s" % (root, outcome["gunc"]))
        else:
            try:
                reports = run_chimerism(
                    dirs,
                    threads=threads,
                    db=gunc_db,
                    env=gunc_env,
                    reassess=reassess,
                    say=say,
                )
                outcome["gunc"] = ", ".join(reports)
            except Exception as exc:  # noqa: BLE001 - reported per variant, not fatal
                outcome["gunc"] = "failed: %s" % exc
                say("    GUNC failed on %s: %s" % (root, exc))

    return dirs, outcome
