"""Assessing a variant's bins at score time, after and outside the measured run.

The measured cost of a variant is the cost of binning. ``metawrap2 binning`` does not run CheckM
unless asked, and CheckM's cost scales with the number of bins emitted, so including it in the
measured command would add a cost the pipeline does not have at that step and would penalise
variants that emit more bins.

Assessors:

- **AMBER** (Meyer et al., *Gigascience* 2018) compares bins with the dataset's contig-to-genome
  truth. It does not use marker genes, which several binners use internally.
- **CheckM** gives the completeness/contamination tiers used in papers (NC, HQ, MQ). It does not
  use the truth and cannot detect a bin made from parts of two genomes.
- **GUNC** checks whether a bin's genes agree on one lineage (chimerism).

Each runs only when the comparison's declared metrics need it.
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

#: Key of the AMBER env in :data:`metawrap2.config.OPTIONAL_ENVS`. AMBER is benchmark-only, so it
#: has its own env instead of sharing a module's.
AMBER_ENV_KEY = "amber"

#: Argument for ``metawrap2 install-env`` (the yaml stem). It differs from the key above, so error
#: messages must print this one.
AMBER_ENV_YAML = "benchmark-amber"

#: AMBER's summary file, read by :func:`~metawrap2.commands.benchmark._collect_scores`. Written
#: under a per-variant subdirectory so bin sets do not overwrite each other.
AMBER_RESULTS = "results.tsv"

#: Sentinel for an env argument that was not passed; an explicit None means "use PATH".
_UNSET = "\0unset"

#: Bioboxes binning-format header. AMBER reads only this format, so the gold standard and every
#: query set are converted to it.
_BIOBOXES_HEADER = "@Version:0.9.1\n@SampleID:%s\n\n@@SEQUENCEID\tBINID\n"

#: Gold-standard header with a length column. AMBER's bp-weighted metrics need a length per
#: sequence and read it from the gold standard only. Without the column AMBER fails with
#: ``KeyError: "['LENGTH'] not in index"``.
_BIOBOXES_GOLD_HEADER = "@Version:0.9.1\n@SampleID:%s\n\n@@SEQUENCEID\tBINID\t_LENGTH\n"

#: Metrics provided by each assessor. An assessor runs only if the comparison declares one of its
#: metrics.
_AMBER_METRICS = ("amber_f1_bp", "amber_ari_bp", "amber_completeness_bp", "amber_purity_bp")
_CHECKM_METRICS = ("bins_near_complete", "bins_high_quality", "bins_medium_quality", "bins_total")
_GUNC_METRICS = ("gunc_pass_fraction",)


def amber_env() -> Optional[str]:
    """The conda env to run AMBER in, or None when the config says not to use conda envs.

    Resolved at call time so that edits to ``use_conda_envs`` in the config take effect.
    """
    settings = _config.load_settings()
    if not settings.use_conda_envs:
        return None
    return _config.OPTIONAL_ENVS[AMBER_ENV_KEY]


#: Per-bin-set score cache, beside the set, named after the scorer (``checkm1`` or ``checkm2``).
#: CheckM1 here uses the full genome tree, as runs before 2.1.4 did.
SCORE_CACHE = ".checkm_scores.%s.json"

#: The module whose env provides CheckM1. Scoring reuses it instead of keeping a second CheckM.
CHECKM_MODULE = "bin_refinement"


def checkm_envs() -> Tuple[Optional[str], Optional[str]]:
    """The envs to run CheckM1 and CheckM2 in, or ``(None, None)`` when conda envs are off.

    Without an env, ``checkm`` is looked up on PATH, where a standard install does not put it
    (``No such file or directory: 'checkm'``).
    """
    settings = _config.load_settings()
    if not settings.use_conda_envs:
        return None, None
    return _config.MODULE_ENVS.get(CHECKM_MODULE), _config.OPTIONAL_ENVS.get("checkm2")


class AssessError(RuntimeError):
    """Raised when a variant's bins cannot be assessed, with the reason and what to do about it."""


def needs_amber(metrics: Iterable[str]) -> bool:
    """True if any of *metrics* comes from AMBER."""
    return any(m in _AMBER_METRICS for m in metrics)


def needs_checkm(metrics: Iterable[str]) -> bool:
    """True if any of *metrics* comes from CheckM."""
    return any(m in _CHECKM_METRICS for m in metrics)


def needs_gunc(metrics: Iterable[str]) -> bool:
    """True if any of *metrics* comes from GUNC.

    GUNC detects chimeric bins that CheckM and AMBER score as nearly clean: a bin that is 95% one
    genome and 5% another has few duplicated markers and little bp-weighted impurity. In this
    project's measurements, a 1.3-point AMBER purity gap between two bin sets came with a 2.8-fold
    difference in chimeric bins.
    """
    return any(m in _GUNC_METRICS for m in metrics)


def bin_dirs(root: str) -> List[str]:
    """Every bin directory under *root*, in a stable order.

    Binning writes one ``<binner>_bins`` directory per binner; refinement writes
    ``metawrap_<c>_<x>_bins``, which has the same suffix. Choosing the variant's final set is up to
    the caller (see :func:`metawrap2.constants.is_refined_bins_dirname`).

    Refined sets one level down are also found, because a binning comparison refines each variant
    into a subdirectory; without this, multi-binner variants had no consolidated set and were left
    unranked. Only refined names are taken from subdirectories: refinement restores its input sets
    there under their original names, and taking those would count the raw sets twice.
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

    Bin headers keep the assembler's description, while the gold standard has only the contig id.
    Untrimmed ids match nothing, and AMBER reports that as zero purity without an error.
    """
    ids = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split(None, 1)[0].strip())
    return ids


def fasta_lengths(path: str) -> Dict[str, int]:
    """Sequence id to sequence length for every record in the FASTA at *path*.

    Use the assembly the binners were given. Lengths from the bins would cover only binned contigs
    and make every completeness denominator equal to what was recovered.
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

    The bin id is the FASTA file name without extension (``bin_007`` for ``bin_007.fasta``).
    Unbinned contigs are left out, which is how the format represents them; an "unbinned" bin would
    be scored as a genome.
    """
    # fasta_files drops the unbinned catch-all and is shared with the modules. The dotfile check
    # skips editor files such as ``.bin_001.fasta.swp``.
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
    """True if a bioboxes header has a length column (``LENGTH`` or ``_LENGTH``)."""
    for line in head.splitlines():
        if line.startswith("@@"):
            columns = [c.strip().upper() for c in line[2:].split("\t")]
            return "LENGTH" in columns or "_LENGTH" in columns
    return False


#: Columns of a CAMI2 contig mapping, past the two every mapping has.
_CAMI_CONTIG_ID, _CAMI_START, _CAMI_END = 3, 5, 6


def _cami_source_contig(fields: List[str]) -> Optional[str]:
    """The name a CAMI2 gold-standard assembly gives the contig described by *fields*, or None.

    ``gsa.fa`` names a contig after its source region
    (``CP013019.1_from_3_to_4351861_total_4351859``); the mapping has those values in separate
    columns. The total includes both ends: ``end - start + 1``.

    Returns None for rows of another shape, so they are reported as unmatched.
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

    Used by :func:`gold_to_bioboxes` and :mod:`metawrap2.benchmark.abundance`, so both resolve
    contig names the same way.

    Only contigs in *lengths* (the assembly the binners were given) are kept. Order follows the
    file and the first occurrence of a repeated contig wins, so CAMI2's raw mapping (one row per
    source region) and its de-duplicated mapping give the same result.

    A CAMI2 mapping's first column is an anonymous id (``PC0``) that ``gsa.fa`` does not use;
    ``gsa.fa`` names contigs from the ``contig_id``, ``start_position`` and ``end_position``
    columns. On the gastrointestinal challenge none of the 210,618 anonymous ids appears in the
    210,633-contig assembly, and the reconstructed names match all 210,618. Each row tries the first
    column, then the reconstructed name.
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
                # Not an assembly contig name; try the name built from the other columns.
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

    Accepts CAMI2's per-body-site mapping (a ``#``-commented table whose first two columns are the
    anonymous contig id and its genome, followed by taxonomy and coordinates) or a bioboxes file. A
    bioboxes file that already has a length column is returned unchanged.

    *lengths* maps contig id to length, from the assembly the binners were given. AMBER's
    bp-weighted metrics read lengths from the gold standard only, so a bioboxes file without them is
    rewritten with them; otherwise AMBER fails with an uninformative pandas ``KeyError``. Contigs
    absent from the assembly are dropped so they do not inflate the recall denominators.

    Parsing is done by :func:`contig_to_genome`, shared with the abundance scoring.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        head = fh.read(4096)
    already_bioboxes = "@@SEQUENCEID" in head
    if already_bioboxes and _declares_lengths(head):
        return path
    if lengths is None:
        raise AssessError(
            "%s has no sequence lengths and none were supplied; AMBER needs them for the "
            "bp-weighted metrics. Pass the assembly the binners were given "
            "with --assembly." % path
        )
    if already_bioboxes:
        # Keep the dataset's own sample id; the dataset's other files and AMBER's report use it.
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
            "%s has no contig-to-genome pairs for contigs in the assembly "
            "(%d contig(s) in the gold standard are not in the "
            "assembly). Expected a bioboxes binning file or a tab-separated table with "
            "the contig id and its genome in the first two columns, over the contigs of the "
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

    AMBER pairs queries with the truth by sample id and rejects queries it cannot pair
    (``Sample ID concoct_bins ... not found in the genome binning gold standard``), so every query
    uses the gold standard's sample id.
    """
    ids = bioboxes_sample_ids(gold)
    if not ids:
        raise AssessError(
            "%s declares no @SampleID, so AMBER cannot pair a query with it. Add a line "
            "'@SampleID:<name>' above the @@SEQUENCEID header." % gold
        )
    if len(ids) > 1:
        raise AssessError(
            "%s declares %d samples (%s). A variant's bins are one co-assembly and are scored "
            "against one sample's truth; with this file the other samples' contigs would count "
            "as misassigned. Pool the gold standard into a single "
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

    All of a variant's bin sets go into one AMBER run, which is cheaper and puts the sets side by
    side in one table. Sets are distinguished by the ``-l`` labels; every query uses the gold
    standard's sample id.
    """
    env = amber_env() if env is _UNSET else env
    os.makedirs(out_dir, exist_ok=True)
    sample_id = gold_sample_id(gold)

    names = list(labels or [os.path.basename(p) for p in bin_dirs_])
    if labels is not None and len(names) != len(bin_dirs_):
        raise AssessError(
            "run_amber got %d labels for %d bin sets; AMBER matches them by position, so the "
            "two counts must be equal" % (len(names), len(bin_dirs_))
        )
    # Query files are named by label, which must be unique. Naming them by bin directory let two
    # arms refined at the same thresholds (both metawrap_70_10_bins) overwrite each other's query,
    # and AMBER reported identical rows for different bin sets.
    if len(set(names)) != len(names):
        duplicates = sorted({n for n in names if names.count(n) > 1})
        raise AssessError(
            "run_amber was given the same label more than once (%s); labels name the output "
            "rows and must be unique" % ", ".join(duplicates)
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
        # The HTML report takes minutes to render and is not used; the scorer parses the TSV.
        "--silent",
        # Skip scoring the gold standard against itself. Its row of perfect scores would otherwise
        # be read as a bin set.
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
    checkm2: bool = True,
    env: Optional[str] = None,
    checkm2_env: Optional[str] = None,
) -> List[str]:
    """Run CheckM over each of *bin_dirs_*, returning the stats files written.

    One run per bin set, so a multi-binner variant does not count the same genome once per binner.

    The size pre-filter in :func:`metawrap2.checkm.run_checkm_cached` matters here. In
    bin_refinement, ``copy_bins`` has already limited bins to 50 kb - 20 Mb (of 10,168 bins there,
    none was over the ceiling and 56, 0.02% of the sequence, under the floor). Raw binner output is
    not filtered, and metaBAT2 emits its unbinnable remainder as a bin: of 9,019 raw bins, 64 were
    over 20 Mb and held 11.8% of the sequence, which drives CheckM's cost. Bins do not recur across
    a variant's sets, so each call gets its own cache. The cache is named after the scorer
    (:data:`SCORE_CACHE`), so CheckM1 scores are never read back as CheckM2 scores.
    """
    scorer = "checkm2" if checkm2 else "checkm1"
    written = []
    for path in bin_dirs_:
        written.append(
            _checkm.run_checkm_cached(
                path,
                cache=_checkm.BinScoreCache(path + SCORE_CACHE % scorer),
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

    Both come from the config, as in :func:`checkm_envs`, so scoring uses the pipeline's GUNC and
    database.
    """
    settings = _config.load_settings()
    # GUNC is in bin_refinement's env. Looked up through MODULE_ENVS so there is one source; a
    # separate GUNC env key caused a KeyError after the envs were merged.
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

    One run per bin set, with the report in ``<bin set>.gunc``.

    A set that already has a report is skipped unless *reassess*. GUNC's diamond stage took 3.2 to
    5.7 CPU-hours per bin set on this project's benchmark arms.
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
    checkm2: bool = True,
    checkm_env: Optional[str] = _UNSET,
    checkm2_env: Optional[str] = _UNSET,
    reassess: bool = False,
    say=print,
) -> Tuple[List[str], Dict[str, str]]:
    """Assess every bin set under *root* with whichever assessors *metrics* calls for.

    Returns the bin directories found and a map of assessor to its output or the reason it did not
    run. Failures are returned instead of raised so one variant does not stop the others.

    *assembly* is the contig file the binners were given. AMBER uses it for sequence lengths when
    the gold standard does not include them.
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
            # Report the missing database with the module's hint before starting a long diamond run.
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
