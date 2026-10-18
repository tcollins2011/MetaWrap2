"""`metawrap2 install-db [name ...] [--all]` - download and configure MetaWrap2's databases.

Database setup used to be the roughest part of getting MetaWrap2 working: five separate
manual downloads copy-pasted out of the docs, each with its own unpacking and indexing step,
followed by hand-editing ``metawrap2.toml`` with the right key names, plus a
``checkm data setRoot`` that is easy to forget and silently breaks bin_refinement later.

This command does all of that: it fetches each database, unpacks it, runs whatever indexing
step it needs, writes the matching key into the config file, and verifies the result. Every
step is skipped if it already looks complete, so re-running it is cheap and safe.

    metawrap2 install-db --list            # what is available, how big, who needs it
    metawrap2 install-db checkm taxdump    # just these two
    metawrap2 install-db --all             # everything (large: see --list for sizes)

Sizes are the real thing, so ``--all`` is a many-hour, several-hundred-GB commitment. The
per-database ``--small`` variants download a capped-but-real stand-in instead, which is what
you want for a smoke test of the pipeline rather than a production analysis.

Because that run is long, it is built to be interrupted: each database's config key is written as
that database finishes, every installer short-circuits on one that is already complete, and
``--list`` gives both numbers you need to plan it - what crosses the network and what is left on
disk afterwards.
"""

from __future__ import annotations

import collections
import concurrent.futures
import contextlib
import glob
import gzip
import os
import shutil
import signal
import subprocess
import sys
import tarfile
import threading
import time
from typing import Callable, Deque, Dict, Iterator, List, NamedTuple, Optional, Set, Tuple

from .. import interrupt as _interrupt
from ..cliformat import make_parser
from ..config import MODULE_ENVS, OPTIONAL_ENVS, env_prefix, find_config

# Databases are independent of each other and mostly network-bound, so fetching several at
# once is close to a linear speedup. Modest default so a slow link is not swamped.
DEFAULT_PARALLEL_DBS = 3
# Serialises the progress prints from concurrent installers so lines do not interleave.
_PRINT_LOCK = threading.Lock()

# How much of a failed installer's own output to show under the FAILED line, and how much of it to
# keep in order to be able to show that much.
#
# A failed database used to report nothing but its name. The reporter's dbcan failure was three
# lines - the header, the command, "FAILED: dbcan" - because the subprocess wrote its error to the
# shared stdout where it was interleaved with two other databases' progress meters under `-t 3`.
# Re-running the identical command by hand worked, so the failure was transient and the one piece
# of information needed to know that was the part that got lost. Each installer's child output is
# captured per database instead, and the tail of it printed where the reason belongs.
#
# 15 lines is enough for a traceback or a curl message without burying the summary. The ring buffer
# is what makes capturing safe for a 300 GB download: curl's progress meter to a pipe writes
# indefinitely, and holding all of it would be a leak proportional to download time.
FAILURE_TAIL_LINES = 15
_CAPTURE_LINES = 400
#: The ring buffer the installer running on *this* thread writes its child output into, or None
#: outside an installer. Thread-local rather than one shared buffer because up to `-t N` installers
#: run at once and the whole point is to keep their output apart.
_CAPTURED = threading.local()


def _say(*lines: str) -> None:
    with _PRINT_LOCK:
        for line in lines:
            print(line)
        sys.stdout.flush()


class Database(NamedTuple):
    key: Optional[str]  # metawrap2.toml [databases] key, None if the tool stores its own path
    used_by: str
    #: What crosses the network, and what is left on disk afterwards - two numbers, because for
    #: half of these they are not the same number and the difference decides whether the install
    #: fits. One "size" column used to carry whichever of the two was known, so a row said
    #: "6.7 GB download" (gunc, ~13 GB once unpacked) beside one that said "1.4 GB" (checkm, a
    #: bare size) and the column could not be added up. ``on_disk`` is marked ``~`` when it is an
    #: estimate rather than a measurement; see the note under ``--list``.
    download: str
    on_disk: str
    small_size: Optional[str]
    describe: str
    #: Written into the config instead of the install directory, when the tool wants a file rather
    #: than a directory - GUNC takes ``--db_file``, so GUNC_DB has to name the .dmnd itself. Also
    #: used for a subdirectory the downloader creates under the install directory: geNomad writes
    #: ``genomad_db/`` into the destination it is given and wants that inner path back.
    key_file: Optional[str] = None
    #: ``key_file`` for a ``--small`` install, when the capped variant lands somewhere else:
    #: ``bakta_db download --type light`` writes ``db-light/`` where ``--type full`` writes
    #: ``db/``, so the configured path depends on which one was asked for.
    small_key_file: Optional[str] = None
    #: Like ``key_file``, for a path whose name is not known until it is downloaded. SingleM's
    #: metapackage directory carries its own version (``S6.5.0.GTDB_r232.metapackage_...``), and
    #: that version is part of the answer - a per-cell fraction is only as current as the reference
    #: behind it - so the versioned directory is what goes into the config rather than a `latest`
    #: symlink that would hide it.
    key_glob: Optional[str] = None

    @property
    def size(self) -> str:
        """The download size, under the name ``quickstart`` reads it by.

        ``quickstart._print_database_guidance`` prints its own short table from this registry and
        reads ``db.size``; splitting the field in two would have made ``metawrap2 quickstart``
        raise AttributeError on a machine where nothing else had changed. One spelling, one
        meaning: the size column anywhere outside ``--list`` is the download.
        """
        return self.download


# Sizes below are measurements where one exists and the download size carried across where one
# does not - never a number invented to fill the column. The ones worth knowing the provenance of:
# gunc's ~13 GB is a reported on-disk measurement against a 6.7 GB download, checkm2's ~3 GB and
# pfam's 1.6 GB come from their own release notes, singlem's 9.7 GB is the unpacked metapackage,
# and the five tool-managed databases (bakta, gtdbtk, amrfinder, genomad, checkv) carry the figures
# the README's table has always quoted. A row whose on-disk figure is an estimate says so with `~`.
DATABASES: Dict[str, Database] = {
    "checkm": Database(
        key="CHECKM_DB",
        used_by="bin_refinement, reassemble_bins, binning --run-checkm",
        download="1.4 GB",
        on_disk="~1.4 GB",
        small_size=None,
        describe="CheckM1 marker-gene reference data",
    ),
    "taxdump": Database(
        key="TAXDUMP",
        used_by="blobology, classify_bins",
        download="550 MB",
        on_disk="~550 MB",
        small_size=None,
        describe="NCBI taxonomy dump (names.dmp / nodes.dmp)",
    ),
    "kraken2": Database(
        key="KRAKEN2_DB",
        used_by="kraken2",
        download="~90 GB",
        on_disk="~90 GB",
        small_size="8 GB",
        describe="Kraken2 standard database (--small: the capped 8 GB standard build)",
    ),
    "blast": Database(
        key="BLASTDB",
        used_by="blobology, classify_bins",
        download="~300 GB",
        on_disk="~300 GB",
        small_size="~2 GB",
        describe="NCBI nt (--small: complete representative prokaryote genomes)",
    ),
    "checkm2": Database(
        key="CHECKM2DB",
        used_by="bin_refinement --checkm2, binning --run-checkm --checkm2",
        download="1.7 GB",
        on_disk="~3 GB",
        small_size=None,
        describe="CheckM2's DIAMOND reference (uniref100.KO; not CheckM1's data - see checkm)",
        key_file=os.path.join("CheckM2_database", "uniref100.KO.1.dmnd"),
    ),
    "gunc": Database(
        key="GUNC_DB",
        used_by="bin_refinement (unless --skip-gunc)",
        download="6.7 GB",
        on_disk="~13 GB",
        small_size=None,
        describe="GUNC proGenomes 2.1 reference (diamond .dmnd, gunzipped to about twice the "
        "download; this one runs by default, so budget the 13 GB)",
        key_file="gunc_db_progenomes2.1.dmnd",
    ),
    "gtdbtk": Database(
        key="GTDBTK_DATA_PATH",
        used_by="classify_bins --gtdbtk",
        download="~110 GB",
        on_disk="~110 GB",
        small_size=None,
        describe="GTDB-Tk reference package, by its own download-db.sh (the tarball is deleted "
        "only after it unpacks, so the install needs about twice this free)",
    ),
    "dbcan": Database(
        key="DBCAN_DB",
        used_by="annotate --cazy",
        download="7.4 GB",
        on_disk="~7.4 GB",
        small_size=None,
        describe="dbCAN CAZyme HMMs, CAZy DIAMOND database and substrate mapping",
    ),
    "pfam": Database(
        key="PFAM_DB",
        used_by="annotate --pfam",
        download="399 MB",
        on_disk="1.6 GB",
        small_size=None,
        describe="Pfam-A HMMs with their per-family gathering thresholds, plus the release notes",
        key_file="Pfam-A.hmm",
    ),
    "eggnog": Database(
        key="EGGNOG_DB",
        used_by="annotate --eggnog",
        download="12 GB",
        on_disk="48 GB",
        small_size=None,
        # The release is not named here: EGGNOG_DB_VERSION is defined with the other URLs, below
        # this table, and a second copy of "5.0.2" would be correct until one of them was edited.
        describe="eggNOG annotation SQLite (39 GB unpacked) and the DIAMOND database of eggNOG "
        "proteins (8.7 GB)",
    ),
    "bakta": Database(
        key="BAKTA_DB",
        used_by="annotate --bakta, annotate_bins --bakta",
        download="30 GB",
        on_disk="~30 GB",
        small_size="1.3 GB",
        describe="Bakta's reference database, by bakta_db (--small: the `light` type, which drops "
        "the UPS/IPS protein sets and annotates from the HMMs alone)",
        key_file="db",
        small_key_file="db-light",
    ),
    "card": Database(
        key="CARD_DB",
        used_by="amr --rgi",
        download="4 MB",
        on_disk="~4 MB",
        small_size=None,
        describe="CARD's resistance ontology and reference sequences (card.json plus its FASTAs)",
        key_file="card.json",
    ),
    "amrfinder": Database(
        key="AMRFINDER_DB",
        used_by="amr (the default backend)",
        download="200 MB",
        on_disk="~200 MB",
        small_size=None,
        describe="AMRFinderPlus reference database, by amrfinder_update, which also builds its "
        "BLAST and HMMER indexes; the dated version directory is what gets configured",
        key_glob="20[0-9][0-9]-[0-9][0-9]-[0-9][0-9]*",
    ),
    "genomad": Database(
        key="GENOMAD_DB",
        used_by="viruses, amr --mobile-elements",
        download="1.5 GB",
        on_disk="~1.5 GB",
        small_size=None,
        describe="geNomad marker database, by `genomad download-database`; one download serves "
        "both the viral calls and the plasmid calls amr's mobility work reads",
        key_file="genomad_db",
    ),
    "checkv": Database(
        key="CHECKV_DB",
        used_by="viruses (unless --skip-checkv)",
        download="6.4 GB",
        on_disk="~6.4 GB",
        small_size=None,
        describe="CheckV reference database, by `checkv download_database`; the versioned "
        "checkv-db-vX directory it creates is what gets configured",
        key_glob="checkv-db-*",
    ),
    "singlem": Database(
        key="SINGLEM_METAPACKAGE",
        used_by="appraise",
        download="2.2 GB",
        on_disk="9.7 GB",
        small_size=None,
        describe="SingleM GTDB marker metapackage (single-copy marker genes + their taxonomy). "
        "Zenodo serves this at a sustained 150-300 KB/s, so the 2.2 GB takes 2-4 hours where "
        "every other database here arrives at 0.5-20 MB/s: start it early, it is not stuck",
        key_glob="*.smpkg.zb",
    ),
    "kofam": Database(
        key="KOFAM_DB",
        used_by="annotate --kofam",
        download="1.5 GB",
        on_disk="~1.5 GB",
        small_size=None,
        describe="KOfam HMM profiles and their per-profile thresholds (ko_list)",
    ),
    "kegg": Database(
        key="KEGG_BRITE",
        used_by="pathways (functional hierarchy and Krona plots)",
        download="4 MB",
        on_disk="~4 MB",
        small_size=None,
        describe="KEGG BRITE ko00001: KO -> pathway -> category, projected to one table",
    ),
    "host": Database(
        key="BMTAGGER_DB",
        used_by="read_qc host removal",
        download="~20 GB",
        on_disk="~20 GB",
        small_size="~9 GB",
        describe="bmtagger index of the host genome (--small: human chr21 only). The bitmask and "
        "srprism index are built locally, which is the slow step rather than the download",
    ),
}

#: Size strings in this registry are human-readable on purpose - they are read far more often than
#: they are summed - so the totals row parses them back. GB throughout; MB and TB convert.
_SIZE_UNITS: Dict[str, float] = {"MB": 1.0 / 1024, "GB": 1.0, "TB": 1024.0}


def _gigabytes(size: str) -> float:
    """Parse one of this registry's size strings ("399 MB", "~110 GB") into GB.

    The `~` is dropped rather than rejected: an estimate still has to be counted, or the totals
    row would omit the largest databases in the table and read as smaller than the smallest
    useful install.
    """
    number, _, unit = size.replace("~", "").strip().partition(" ")
    return float(number) * _SIZE_UNITS[unit.strip()]


def _total(sizes: List[str]) -> str:
    """Sum size strings for the totals row, always marked `~` because the inputs are mixed."""
    return "~%d GB" % round(sum(_gigabytes(size) for size in sizes))


def _print_list() -> None:
    """The authoritative database table: two numbers per row and a column that can be added up.

    What this replaced put four kinds of number in one "full size" column - a bare size, a
    download size, a download-and-unpacked pair, and "~20 GB + indexing" - so only 4 of 13 rows
    said how much disk the install needed, the two-number rows overflowed into the `--small`
    column and broke the alignment, and planning a disk budget from it was not possible. GUNC was
    the costly case: it read "6.7 GB download" and occupies about 13 GB, and it runs by default in
    bin_refinement, so the under-statement landed on everybody.
    """
    print("\nMetaWrap2 databases\n")
    print("  %-10s %-11s %-11s %-11s %s" % ("name", "download", "on disk", "--small", "used by"))
    print("  " + "-" * 84)
    for name, db in DATABASES.items():
        print(
            "  %-10s %-11s %-11s %-11s %s"
            % (name, db.download, db.on_disk, db.small_size or db.download, db.used_by)
        )
    print("  " + "-" * 84)
    print(
        (
            "  %-10s %-11s %-11s %-11s"
            % (
                "total",
                _total([db.download for db in DATABASES.values()]),
                _total([db.on_disk for db in DATABASES.values()]),
                _total([db.small_size or db.download for db in DATABASES.values()]),
            )
        ).rstrip()
    )
    print()
    for name, db in DATABASES.items():
        print("  %-10s %s" % (name, db.describe))
    print(
        "\n  Both numbers matter and they are not the same: `download` crosses the network, "
        "`on disk`\n  is what is left after unpacking and indexing. An on-disk figure marked `~` "
        "is an estimate;\n  where no measurement exists the download size is carried across "
        "rather than guessed upwards,\n  so read the on-disk total as a floor. `--small` is what "
        "that flag downloads - for a database\n  with no capped variant it repeats the full "
        "download, because that is what --small installs."
    )
    print("\n  Install into a directory of your choice with --dir (default ~/metawrap2_dbs).")
    print("  Each database's path is written to your config file as it finishes.\n")


# --- shared helpers ---------------------------------------------------------------------


#: The process group of every download running right now, one per child started by :func:`_run`.
#: Kept so that an interrupt can kill the *group* rather than only the process we started: curl is
#: the one writing the file, and a downloader from a conda env is several processes deep.
_DOWNLOAD_GROUPS: Set[int] = set()
_GROUPS_LOCK = threading.Lock()


@contextlib.contextmanager
def _tracked(child: subprocess.Popen) -> Iterator[None]:
    """Remember *child*'s process group for as long as it runs.

    Started with ``start_new_session``, so the child leads its own group and ``child.pid`` is that
    group's id. Both halves matter: a new session is what stops a terminal's Ctrl-C from reaching
    curl behind this command's back, and tracking the group is what lets
    :func:`_terminate_downloads` then kill it on purpose. Without the second half the first one
    turns a Ctrl-C into the very orphan it was meant to prevent.
    """
    with _GROUPS_LOCK:
        _DOWNLOAD_GROUPS.add(child.pid)
    try:
        yield
    finally:
        with _GROUPS_LOCK:
            _DOWNLOAD_GROUPS.discard(child.pid)


def _terminate_downloads() -> None:
    """SIGTERM every download process group still running, and forget it.

    This is the half of the orphan fix that acts. A killed `install-db` used to leave its curl
    running: the reporter's `ps` showed a 44-minute curl from a stopped run and a 16-minute one
    from the next, both writing the same `.part` at their own offsets, and the garbage that makes
    surfaces later as a `tar` or `gzip` failure that reads like a bad download. Killing the group
    rather than the child covers the downloaders that are a shell script around a tool around curl.

    A group that has already exited is not an error: the race being closed is a writer that
    outlives this process, not one that finished before it.
    """
    with _GROUPS_LOCK:
        groups = sorted(_DOWNLOAD_GROUPS)
        _DOWNLOAD_GROUPS.clear()
    for pgid in groups:
        try:
            os.killpg(pgid, signal.SIGTERM)
        except OSError:
            continue


def _run(cmd: List[str], **kw) -> int:
    """Run *cmd*, capturing its output into this thread's buffer when there is one.

    Outside an installer there is no buffer and the child writes straight to this process's
    stdout, which is what a direct call (a test, another command) wants. Inside one, stderr is
    folded into stdout and both are kept so that a failure can quote them: the alternative is the
    status quo, where three databases' output arrives on one terminal interleaved and the one that
    failed is unreadable.

    Every child goes into its own session - see :func:`_tracked` - and is launched through Popen
    even when its output is not captured, because a pid is what the cleanup needs and
    ``subprocess.call`` does not hand one back.
    """
    print("  $ %s" % " ".join(cmd))
    sink: Optional[Deque[str]] = getattr(_CAPTURED, "lines", None)
    if sink is None:
        child = subprocess.Popen(cmd, start_new_session=True, **kw)
        with _tracked(child):
            return child.wait()
    sink.append("$ %s" % " ".join(cmd))
    child = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
        start_new_session=True,
        **kw,
    )
    with _tracked(child):
        if child.stdout is not None:
            for line in child.stdout:
                sink.append(line.rstrip("\n"))
        return child.wait()


def _install_capturing(name: str, target: str, small: bool, jobs: int) -> Tuple[bool, List[str]]:
    """Run one installer with its child output captured, and hand back whatever it said.

    The exception text is appended to the captured lines as well as printed: an installer that
    raises is the case where the reason is most likely to be the last thing printed rather than
    the subprocess's own output.

    Only `Exception` is caught, deliberately: a Ctrl-C has to leave this function and stop the run
    rather than be recorded as "this database failed" while the next one starts downloading.
    """
    sink: Deque[str] = collections.deque(maxlen=_CAPTURE_LINES)
    _CAPTURED.lines = sink
    try:
        ok = INSTALLERS[name](target, small, jobs)
    except Exception as exc:  # noqa: BLE001 - one bad database must not stop the others
        sink.append("%s: %s" % (type(exc).__name__, exc))
        _say("  ERROR installing %s: %s" % (name, exc))
        ok = False
    finally:
        _CAPTURED.lines = None
    return ok, list(sink)


def _failure_tail(output: List[str]) -> List[str]:
    """The last few captured lines, indented under the FAILED line they explain."""
    tail = output[-FAILURE_TAIL_LINES:]
    if not tail:
        return ["    (the installer produced no output)"]
    return ["    %s" % line for line in tail]


#: How long a `.part` has to have been untouched before this command will resume it.
#:
#: Resuming is the useful default - these downloads are hours long and a stale `.part` is most of
#: the work already done - but `curl -C -` appends from wherever the file currently ends, so two
#: curls resuming the same `.part` write at their own offsets and interleave into a file that is
#: not the download. That is what the reporter hit, and it surfaces as a `tar`/`gzip` failure
#: rather than as a collision.
#:
#: 120 seconds because the quiet gaps a *live* curl can have are all much shorter: `--retry 3` with
#: no delay backs off 1s, 2s, 4s, and a stalled TCP connection goes quiet for a connection timeout,
#: not for minutes. Erring long is deliberate - resuming over a live writer corrupts silently,
#: while refusing a `.part` that is actually dead costs one message that names the file to delete.
LIVE_PART_SECONDS = 120


def _another_writer(part: str) -> bool:
    """Does *part* look like a file another `install-db` is still writing? Prints the reason.

    Mtime rather than a lock file: the writer is curl, started by whatever process, possibly
    orphaned from a run that no longer exists, so there is nobody to hold a lock and nothing to
    ask. What can be observed is that the bytes are still arriving.
    """
    try:
        idle = time.time() - os.path.getmtime(part)
    except OSError:
        return False
    if idle >= LIVE_PART_SECONDS:
        return False
    print(
        "  another install-db appears to be downloading this database "
        "(.part modified %.0fs ago); wait for it or delete %s" % (max(idle, 0.0), part)
    )
    return True


def _fetch(url: str, dest: str) -> bool:
    """Download *url* to *dest* with curl, resuming a partial file that nobody else is writing."""
    if os.path.isfile(dest) and os.path.getsize(dest):
        print("  already downloaded: %s" % os.path.basename(dest))
        return True
    tmp = dest + ".part"
    if os.path.isfile(tmp) and _another_writer(tmp):
        return False
    rc = _run(["curl", "-fL", "--retry", "3", "-C", "-", "-o", tmp, url])
    if rc != 0:
        print("  download failed: %s" % url)
        return False
    os.replace(tmp, dest)
    return True


def _extract(archive: str, into: str) -> bool:
    """Unpack *archive* into *into*, refusing members that would land outside it.

    These archives come off the public internet, and a bare `extractall` writes wherever a member's
    name says to - including `/etc/...` and `../../..`. Python 3.12's `data` filter rejects exactly
    that (and strips the ownership and mode bits a database tarball has no business setting); the
    explicit check below is the same guard for older versions, so the refusal does not depend on
    which interpreter happens to be running.

    This is the only place any archive here is unpacked. It used to be one of two: the CARD
    installer hand-rolled the same check against realpath() prefixes and got it wrong, which made
    that database impossible to install at all. One guard, used by everything, is the point -
    compression is detected from the file, so bz2 and gz both arrive here.
    """
    print("  unpacking %s" % os.path.basename(archive))
    try:
        with tarfile.open(archive) as tf:
            safe = []
            for member in tf.getmembers():
                if _escapes(member):
                    print("  skipping %s: it would be written outside %s" % (member.name, into))
                    continue
                if not (member.isfile() or member.isdir() or member.issym() or member.islnk()):
                    # A device node or fifo in a database tarball is not something to write out.
                    print("  skipping %s: not a file, directory or link" % member.name)
                    continue
                safe.append(member)
            if hasattr(tarfile, "data_filter"):
                tf.extractall(into, members=safe, filter="data")
            else:
                tf.extractall(into, members=safe)
    except (OSError, tarfile.TarError) as problem:
        print("  could not unpack %s: %s" % (os.path.basename(archive), problem))
        return False
    return True


def _untar(archive: str, into: str) -> None:
    """`_extract`, and remove the archive once it has been unpacked.

    A failed unpack keeps the archive: it is usually a truncated download, and `_fetch` resumes
    one with `curl -C -` on the next run rather than starting the download again.
    """
    if _extract(archive, into):
        os.remove(archive)


def _gunzip(archive: str, dest: str) -> bool:
    """Decompress *archive* to *dest*, removing the archive. Streamed, not read into memory."""
    print("  unpacking %s" % os.path.basename(archive))
    tmp = dest + ".part"
    with gzip.open(archive, "rb") as src, open(tmp, "wb") as out:
        shutil.copyfileobj(src, out, length=8 << 20)
    os.replace(tmp, dest)
    os.remove(archive)
    return os.path.isfile(dest) and os.path.getsize(dest) > 0


def _escapes(member: tarfile.TarInfo) -> bool:
    """Would extracting *member* write outside the destination directory?

    Checked on the link target as well as the name: a symlink to `/etc` is how the next member in
    the archive gets written there.
    """
    for name in (member.name, member.linkname or ""):
        if name.startswith(("/", "\\")):
            return True
        if ".." in name.replace("\\", "/").split("/"):
            return True
    return False


def _mamba_run(env: str, cmd: List[str]) -> int:
    # Flags come from _run_flags: --no-capture-output is conda-only and mamba 2.x chokes on it.
    from ..command import _run_flags

    return _run(["mamba", "run"] + _run_flags("mamba") + ["-n", env] + cmd)


def dry_run() -> bool:
    """True if this is a ``--dry-run``: print the plan, touch nothing.

    Read from the one place that knows, ``command.runner``, which is what the global ``--dry-run``
    flag sets in ``cli`` and what every module asks through ``modules._common.dry_run``. This
    command used to have no notion of it at all - ``_run`` went straight to subprocess, so
    ``metawrap2 --dry-run install-db gtdbtk`` downloaded 110 GB from a flag whose entire purpose
    is to not do that. A second switch of our own would have been a second answer to a question
    that already has one.
    """
    from ..command import runner

    return bool(runner.dry_run)


# --- what each installer fetches or runs -------------------------------------------------
#
# URLs and downloader argv live here, not inline in the installers, because `--dry-run` has to
# print what the real run would fetch. One definition read by both: a plan composed from a second
# copy of these strings would be correct exactly until one of them was edited.

CHECKM_URL = "https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
KRAKEN2_URL = "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz"
#: A real standard build, capped to 8 GB - right for testing, not for production.
KRAKEN2_SMALL_URL = "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz"
GUNC_URL = "https://swifter.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.gz"
PFAM_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
PFAM_FILES = ("Pfam-A.hmm.gz", "Pfam-A.hmm.dat.gz")
CARD_URL = "https://card.mcmaster.ca/latest/data"

#: eggNOG's emapper database files. Fetched here rather than through eggNOG-mapper's own
#: `download_eggnog_data.py`, which hardcodes `eggnogdb.embl.de` - a host that no longer answers.
#: Run against it the downloader printed "Finished" over three zero-byte files, because it shells
#: out to wget per file and checks neither the exit status nor the result. The live mirror is
#: eggnog5.embl.de, and fetching the three files directly also gets this project's own resume,
#: live-`.part` refusal and progress reporting.
EGGNOG_DB_VERSION = "5.0.2"
EGGNOG_BASE_URL = "http://eggnog5.embl.de/download/emapperdb-%s" % EGGNOG_DB_VERSION
#: What `-m diamond` needs: the annotation SQLite, the taxonomy, and the DIAMOND protein database.
#: Each is gzipped on the server except the taxonomy, which is a tarball.
EGGNOG_FILES: Tuple[str, ...] = (
    "eggnog.db.gz",
    "eggnog.taxa.tar.gz",
    "eggnog_proteins.dmnd.gz",
)
HOST_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
#: Real human sequence, one chromosome - enough to exercise host removal.
HOST_SMALL_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz"
#: What the --small BLAST build is assembled from; the full one is fetched by update_blastdb.pl.
BLAST_SMALL_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"

#: The databases a tool downloads for itself: the env that owns the tool, and the argv, as a
#: function of ``--small``. Both the installer and the dry run read this, so "the command it would
#: run" and "the command it runs" cannot drift apart. The env comes from config.py for the reason
#: the CheckM installer records below - a hardcoded env name survived a rename and failed quietly.
DOWNLOADERS: Dict[str, Tuple[str, Callable[[str, bool], List[str]]]] = {
    "checkm2": (
        OPTIONAL_ENVS["checkm2"],
        lambda target, small: ["checkm2", "database", "--download", "--path", target],
    ),
    "dbcan": (
        MODULE_ENVS["annotate"],
        lambda target, small: ["run_dbcan", "database", "--db_dir", target, "--aws_s3"],
    ),
    "singlem": (
        MODULE_ENVS["appraise"],
        lambda target, small: ["singlem", "data", "--output-directory", target],
    ),
    "amrfinder": (
        MODULE_ENVS["amr"],
        lambda target, small: ["amrfinder_update", "--database", target],
    ),
    "genomad": (
        MODULE_ENVS["viruses"],
        lambda target, small: ["genomad", "download-database", target],
    ),
    "checkv": (
        MODULE_ENVS["viruses"],
        lambda target, small: ["checkv", "download_database", target],
    ),
    "bakta": (
        MODULE_ENVS["annotate"],
        lambda target, small: [
            "bakta_db",
            "download",
            "--output",
            target,
            "--type",
            "light" if small else "full",
        ],
    ),
    "gtdbtk": (
        OPTIONAL_ENVS["gtdbtk"],
        lambda target, small: ["download-db.sh", target],
    ),
}


def _run_downloader(name: str, target: str, small: bool) -> int:
    """Run *name*'s own downloader in the env that owns it.

    The env is checked first because five of these databases are fetched by their own tool, and
    `quickstart` builds only the default envs - so `install-db checkm2` on a fresh machine failed
    with mamba's own words, `critical libmamba The given prefix does not exist:
    ".../envs/metawrap2-checkm2"`, followed by `FAILED: checkm2`. That names the prefix rather
    than the thing to do about it. One readable line instead, naming the command that fixes it.
    """
    env, argv = DOWNLOADERS[name]
    if not dry_run() and env_prefix(env) is None:
        print(
            "  cannot install %s: its downloader lives in the %s environment, which is not "
            "built yet.\n  Run `metawrap2 install-env %s` first, then this command again."
            % (name, env, env.replace("metawrap2-", ""))
        )
        return 1
    return _mamba_run(env, argv(target, small))


# --- is it already installed? ------------------------------------------------------------
#
# One answer per database, because two callers need it: the installer, to short-circuit (which is
# what makes re-running this command cheap), and the dry run, which has to say "already present,
# would skip" rather than announce it would download 90 GB that is already on the disk.
#
# The general rule is that a database is installed when the path that would go into the config is
# there and not empty - which is the right question by construction, since that is the path every
# module will open. Only the databases whose key names the install directory itself need a marker
# file naming here.
PRESENT_WHEN: Dict[str, Tuple[str, ...]] = {
    "checkm": ("taxon_marker_sets.tsv",),
    "taxdump": ("names.dmp",),
    "kraken2": ("hash.k2d",),
    "dbcan": ("dbCAN.hmm",),
    "kofam": ("ko_list.gz", "profiles"),
    #: GTDB-Tk's package is directories all the way down and its own download-db.sh decides it is
    #: installed the same way - by whether the destination holds a subdirectory.
    "gtdbtk": ("*/",),
    #: What read_qc actually opens. hg38.fa is an input to the indexing, not a product of it.
    "host": ("hg38.bitmask", "hg38.srprism.idx"),
    #: Both files emapper needs, so a 27 GB download interrupted between them is not "installed".
    #: EGGNOG_DB itself stays the *directory*, because emapper takes --data_dir - naming the file
    #: in the config is the shape that lost Pfam's release version in 2.1.2.
    "eggnog": ("eggnog.db", "eggnog_proteins.dmnd"),
}


def _already_installed(name: str, target: str, small: bool) -> bool:
    """Is *name* complete in *target*? Never writes, so it is safe under a dry run."""
    if name == "blast":
        # BLAST's database is a dozen files under a name the config can change (BLASTDB_NAME), so
        # a filename pattern here would be a second and worse answer to a question blastdb.find
        # already answers.
        from .. import blastdb
        from ..config import Settings

        present, _ = blastdb.find(Settings(databases={"BLASTDB": target}))
        return bool(present)
    if name == "kegg":
        from ..hierarchy import HIERARCHY_FILE

        patterns: Tuple[str, ...] = (HIERARCHY_FILE,)
    elif name in PRESENT_WHEN:
        patterns = PRESENT_WHEN[name]
    else:
        db = DATABASES[name]
        key_file = db.small_key_file if (small and db.small_key_file) else db.key_file
        patterns = (db.key_glob or key_file or "",)
    return all(_exists_nonempty(os.path.join(target, pattern)) for pattern in patterns)


def _exists_nonempty(pattern: str) -> bool:
    """Does *pattern* match a file with bytes in it, or a directory with something in it?

    Globbed because three of these databases name themselves after the release they are
    (`checkv-db-v1.5`, `*.smpkg.zb`, AMRFinderPlus's dated directories), and emptiness is checked
    because an interrupted download leaves the directory it was writing into behind.
    """
    for match in glob.glob(pattern):
        if os.path.isdir(match) and os.listdir(match):
            return True
        if os.path.isfile(match) and os.path.getsize(match):
            return True
    return False


# --- per-database installers -------------------------------------------------------------


def _install_checkm(target: str, small: bool, jobs: int = 8) -> bool:
    if _already_installed("checkm", target, small):
        print("  CheckM data already present.")
    else:
        os.makedirs(target, exist_ok=True)
        archive = os.path.join(target, "checkm_data.tar.gz")
        if not _fetch(CHECKM_URL, archive):
            return False
        _untar(archive, target)

    # CheckM stores this path inside its own env, so this is belt and braces rather than the
    # mechanism: every CheckM1 call also gets CHECKM_DATA_PATH from the configured CHECKM_DB (see
    # checkm.checkm1_data_env), which is what makes the setting survive the env being rebuilt. Doing it
    # anyway means `checkm` run by hand from the env also works.
    #
    # The env name comes from MODULE_ENVS rather than being written out: it used to be spelled
    # "metawrap2-bin_refinement" here, and consolidating the envs left this line pointing at one that
    # no longer exists - which fails quietly, since the message below is only a NOTE.
    checkm_env = MODULE_ENVS["bin_refinement"]
    print("  pointing CheckM at the data (checkm data setRoot)")
    if _mamba_run(checkm_env, ["checkm", "data", "setRoot", target]) != 0:
        print(
            "  NOTE: could not run `checkm data setRoot` automatically. CHECKM_DATA_PATH is passed\n"
            "        per call so bin_refinement will work regardless, but to use `checkm` directly:\n"
            "        mamba run -n %s checkm data setRoot %s" % (checkm_env, target)
        )
    return os.path.isfile(os.path.join(target, "taxon_marker_sets.tsv"))


def _install_taxdump(target: str, small: bool, jobs: int = 8) -> bool:
    if _already_installed("taxdump", target, small):
        print("  taxdump already present.")
        return True
    os.makedirs(target, exist_ok=True)
    archive = os.path.join(target, "taxdump.tar.gz")
    if not _fetch(TAXDUMP_URL, archive):
        return False
    _untar(archive, target)
    return _already_installed("taxdump", target, small)


def _install_kraken2(target: str, small: bool, jobs: int = 8) -> bool:
    if _already_installed("kraken2", target, small):
        print("  Kraken2 database already present.")
        return True
    os.makedirs(target, exist_ok=True)
    archive = os.path.join(target, "k2_standard_08gb.tar.gz" if small else "k2_standard.tar.gz")
    if not _fetch(KRAKEN2_SMALL_URL if small else KRAKEN2_URL, archive):
        return False
    _untar(archive, target)
    return _already_installed("kraken2", target, small)


def _install_blast(target: str, small: bool, jobs: int = 8) -> bool:
    # Every other installer short-circuits when its database is already there; without the
    # same check here, re-running install-db re-ran makeblastdb over the whole database.
    if _already_installed("blast", target, small):
        print("  BLAST database already present.")
        return True
    os.makedirs(target, exist_ok=True)
    if small:
        print("  building a representative-genome BLAST database (stand-in for nt)")
        script = os.path.join(target, "build_small_nt.sh")
        with open(script, "w") as fh:
            fh.write(_SMALL_BLAST_SCRIPT)
        os.chmod(script, 0o755)
        # Thousands of small genome downloads: run them concurrently or this dominates.
        if _run(["bash", script, str(jobs)], cwd=target) != 0:
            return False
    else:
        print("  downloading NCBI nt with update_blastdb.pl (this is very large)")
        if (
            _mamba_run(MODULE_ENVS["classify_bins"], ["update_blastdb.pl", "--decompress", "nt"])
            != 0
        ):
            return False
    return True


_SMALL_BLAST_SCRIPT = r"""#!/bin/bash
# Build a small but real BLAST nucleotide database from complete RefSeq representative
# prokaryote genomes, named "nt" so it matches MetaWrap2's default BLASTDB_NAME.
set -euo pipefail
JOBS="${1:-8}"
curl -fsSL -o assembly_summary.txt \
  https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
awk -F'\t' '$12=="Complete Genome" && $11=="latest" && ($5=="reference genome" || $5=="representative genome")' \
  assembly_summary.txt | awk -F'\t' '{print $1"\t"$6"\t"$20}' > download_list.tsv
mkdir -p genomes

# Fetch the genomes JOBS at a time - serially this is the slowest part of the whole build.
fetch_one() {
  url="$1"; base=$(basename "$url"); f="genomes/${base}_genomic.fna.gz"
  [ -s "$f" ] && return 0
  curl -fsSL --retry 2 -o "$f.part" "${url}/${base}_genomic.fna.gz" && mv "$f.part" "$f" \
    || { rm -f "$f.part"; return 0; }
}
export -f fetch_one
cut -f3 download_list.tsv | xargs -P "$JOBS" -I{} bash -c 'fetch_one "$@"' _ {}

# Build the seqid -> taxid map after the downloads, so a retry never appends twice.
: > taxid_map.txt
while IFS=$'\t' read -r acc taxid url; do
  base=$(basename "$url"); f="genomes/${base}_genomic.fna.gz"
  [ -s "$f" ] || continue
  zcat "$f" | grep '^>' | sed 's/^>//' | awk -v t="$taxid" '{print $1"\t"t}' >> taxid_map.txt
done < download_list.tsv
zcat genomes/*.fna.gz > all_genomes.fna
makeblastdb -in all_genomes.fna -dbtype nucl -out nt -parse_seqids -taxid_map taxid_map.txt \
  -title "RefSeq representative prokaryote genomes"
rm -f all_genomes.fna
"""


def _install_checkm2(target: str, small: bool, jobs: int = 8) -> bool:
    """CheckM2's DIAMOND reference, fetched by CheckM2 itself.

    ``checkm2 database --download`` is the only supported way to get this: the file lives on Zenodo
    under a versioned record that CheckM2 pins, and a hand-picked URL goes stale the next time
    CheckM2 is upgraded. So this installer runs the tool rather than fetching a URL of its own.

    It is a different database from ``checkm``, which is the one thing about CheckM2 that reliably
    catches people out - both are "the CheckM database", they share no files, and having one
    installed does nothing for the other. ``--small`` is ignored: completeness estimates come from a
    model trained against this exact reference, so a trimmed copy would not be a faster CheckM2 but
    a differently-calibrated one.
    """
    if _already_installed("checkm2", target, small):
        print("  CheckM2 database already present.")
        return True
    os.makedirs(target, exist_ok=True)
    if _run_downloader("checkm2", target, small) != 0:
        return False
    return _already_installed("checkm2", target, small)


def _install_gunc(target: str, small: bool, jobs: int = 8) -> bool:
    """The GUNC proGenomes 2.1 diamond database.

    One file, so there is no indexing step and nothing to configure inside a conda env - which makes
    this the simplest installer here. ``--small`` is ignored because GUNC's verdict is only as good
    as the breadth of what it compares against: a capped reference would report bins as
    under-represented rather than as clean, which is a worse answer than a slower one.
    """
    if _already_installed("gunc", target, small):
        print("  GUNC database already present.")
        return True
    os.makedirs(target, exist_ok=True)
    dmnd = os.path.join(target, DATABASES["gunc"].key_file or "")
    archive = dmnd + ".gz"
    if not _fetch(GUNC_URL, archive):
        return False
    return _gunzip(archive, dmnd)


def _install_dbcan(target: str, small: bool, jobs: int = 8) -> bool:
    """dbCAN's databases, fetched by dbCAN itself.

    ``--aws_s3`` rather than the default HTTP path, and that is the whole decision here: dbCAN's
    own help calls the HTTP route "db_current (moving snapshot)" while the S3 route is a pinned
    release. A CAZyme family call is only reproducible against the release that produced it, so a
    database that changes under a stable name is worse than a large one.

    ``--small`` is ignored: the HMM sets are the method, and the three-method agreement this module
    reports needs all three present.
    """
    if _already_installed("dbcan", target, small):
        print("  dbCAN database already present.")
        return True
    os.makedirs(target, exist_ok=True)
    if _run_downloader("dbcan", target, small) != 0:
        return False
    return _already_installed("dbcan", target, small)


def _install_eggnog(target: str, small: bool, jobs: int = 8) -> bool:
    """eggNOG's annotation database, fetched file by file from the live mirror.

    Not through eggNOG-mapper's own `download_eggnog_data.py`, for the reason recorded above
    EGGNOG_FILES: it points at a host that no longer answers and reports success regardless. The
    three files are gunzipped into the flat layout emapper's `--data_dir` expects. 12 GB on the
    wire and 48 GB unpacked, which is why the two figures are reported separately in `--list`.

    ``--small`` is ignored because eggNOG has no capped variant: its optional flags fetch *more*
    data (Pfam HMMs, MMseqs indexes), not less. Saying so beats offering a flag that does nothing.
    """
    if _already_installed("eggnog", target, small):
        print("  eggNOG database already present.")
        return True
    os.makedirs(target, exist_ok=True)
    for name in EGGNOG_FILES:
        archive = os.path.join(target, name)
        if not _fetch("%s/%s" % (EGGNOG_BASE_URL, name), archive):
            return False
        if name.endswith(".tar.gz"):
            if not _extract(archive, target):
                return False
            # The tarball is the only one of the three that is not a single gzipped file, and
            # emapper looks for `eggnog.taxa.db` beside the others rather than in a subdirectory.
            _flatten_into(target, os.path.splitext(os.path.splitext(name)[0])[0])
            continue
        if not _gunzip(archive, archive[: -len(".gz")]):
            return False
    print("  eggNOG database %s installed." % EGGNOG_DB_VERSION)
    return _already_installed("eggnog", target, small)


def _flatten_into(target: str, subdir: str) -> None:
    """Move anything the taxonomy tarball unpacked into *subdir* up into *target*, then remove it.

    emapper resolves `eggnog.taxa.db` beside `eggnog.db`, so a tarball that unpacks into its own
    directory leaves the tool looking in the wrong place - and the error it gives names a path the
    user never typed.
    """
    inner = os.path.join(target, subdir)
    if not os.path.isdir(inner):
        return
    for name in sorted(os.listdir(inner)):
        source = os.path.join(inner, name)
        destination = os.path.join(target, name)
        if not os.path.exists(destination):
            os.replace(source, destination)
    try:
        os.rmdir(inner)
    except OSError:
        # Something unexpected is still in there; leaving it is better than deleting blind.
        pass


def _install_pfam(target: str, small: bool, jobs: int = 8) -> bool:
    """Pfam-A, from EBI's current release.

    The HMMs carry their own gathering thresholds, which is why no separate cutoff table is needed -
    and why this is one file rather than the half-dozen InterProScan's data would be.

    `relnotes.txt` is fetched with them and is not decoration: Pfam renumbers and merges families
    between releases, so a domain call means "PF00001 as of release 37.0" and the release is what
    makes it reproducible. `--small` is ignored; the profile set is the method.
    """
    if _already_installed("pfam", target, small):
        print("  Pfam already present (release %s)." % _pfam_release(target))
        return True
    os.makedirs(target, exist_ok=True)
    for name in PFAM_FILES:
        if not _fetch("%s/%s" % (PFAM_BASE_URL, name), os.path.join(target, name)):
            return False
    _fetch("%s/relnotes.txt" % PFAM_BASE_URL, os.path.join(target, "relnotes.txt"))
    for name in PFAM_FILES:
        if _run(["gunzip", "-f", os.path.join(target, name)]) != 0:
            return False
    if not _already_installed("pfam", target, small):
        return False
    print("  Pfam release %s" % _pfam_release(target))
    return True


def _pfam_release(target: str) -> str:
    from ..pfam import release_from_notes

    return release_from_notes(target)


def _install_card(target: str, small: bool, jobs: int = 8) -> bool:
    """CARD's data, from CARD's own `latest/data` endpoint.

    Fetched by URL rather than by a tool, because RGI has no downloader: `rgi load` takes a
    card.json that is already on disk. The tarball is 4 MB, which is why this one is listed with no
    `--small` variant - there is nothing to trim.

    The version inside card.json is what dates every call made against it, so it is printed here and
    recorded per run in the AMR table. CARD releases every few months and the URL always points at
    the newest one, so two installs months apart are different databases under one path: rename the
    directory after the version if that matters to you.
    """
    marker = os.path.join(target, DATABASES["card"].key_file or "")
    if _already_installed("card", target, small):
        print("  CARD database already present (version %s)." % _card_version(marker))
        return True
    os.makedirs(target, exist_ok=True)
    archive = os.path.join(target, "card_data.tar.bz2")
    if not _fetch(CARD_URL, archive):
        return False
    # Unpacked by `_extract` like every other archive here. This installer used to hand-roll the
    # same escape check against realpath() prefixes, and that copy made CARD impossible to
    # install: CARD's tarball has 19 members and the first is `.`, the archive root, for which
    # realpath("<target>/.") normalises to "<target>" - not a match for "<target>" plus a
    # separator, so the one member that could not possibly escape was the one rejected, and the
    # installer aborted the whole run. All 18 real members passed. `_escapes` asks the question
    # that decides it - a leading `/`, or a `..` component - under which `.` is simply fine.
    if not _extract(archive, target):
        return False
    os.remove(archive)
    if not (os.path.isfile(marker) and os.path.getsize(marker)):
        return False
    print("  CARD version %s" % _card_version(marker))
    return True


def _card_version(card_json: str) -> str:
    from ..amr import card_version

    return card_version(card_json)


# --- the five that each tool downloads for itself ----------------------------------------
#
# amrfinder, genomad, checkv, bakta and gtdbtk were advertised in the README's table and in
# `usage.md` (`--genomad-db $GENOMAD_DB`, as though it already existed) but were not in this
# registry, so `install-db gtdbtk` answered "Unknown database(s): gtdbtk" and four modules -
# `amr`, `viruses`, `annotate --bakta`, `classify_bins --gtdbtk` - had no documented way to get
# their data at all.
#
# All five shell out to the tool's own downloader, the pattern `checkm2` already uses here:
# the tool knows which database version it can read, the URL moves with the release, and several
# of these build indexes after downloading that a plain fetch would leave undone. What this
# command adds is the part the tools do not do - running inside the right env, recording the
# resolved path under the key the modules read, and being idempotent on a re-run.
#
# They are checked as "the directory the downloader was told to fill has something in it" rather
# than against a filename inside somebody else's data package. Those layouts change between
# releases, and a successful 30 GB bakta download reported as FAILED because a marker file was
# renamed upstream is a worse failure than a half-finished one reported as present - which the
# tool itself catches the first time a module runs against it.


def _has_content(path: str) -> bool:
    return os.path.isdir(path) and bool(os.listdir(path))


def _subdirectories(path: str) -> List[str]:
    return sorted(p for p in glob.glob(os.path.join(path, "*")) if os.path.isdir(p))


def _install_amrfinder(target: str, small: bool, jobs: int = 8) -> bool:
    """The AMRFinderPlus reference database, fetched by amrfinder_update.

    `amrfinder_update` rather than a URL because it does two things a fetch would not: it picks the
    database version the installed AMRFinderPlus can actually read, and it runs makeblastdb and
    hmmpress over what it downloaded, without which `amr` fails at its first call.

    It writes a dated directory per version (`2024-07-22.1`) plus a `latest` symlink, and the dated
    one is what goes into the config - see AMRFINDER_DB in config.py. An ARG call is only as
    datable as its database, and a rerun pointed at a symlink that has since moved is a different
    experiment wearing the same name. `--small` is ignored: NCBI curates one per-gene threshold set
    and a subset of it would not be a faster `amr` but an incomparable one.
    """
    versions = _amrfinder_versions(target)
    if versions:
        print("  AMRFinderPlus database already present: %s" % os.path.basename(versions[-1]))
        return True
    os.makedirs(target, exist_ok=True)
    if _run_downloader("amrfinder", target, small) != 0:
        return False
    versions = _amrfinder_versions(target)
    if versions:
        print("  AMRFinderPlus database version %s" % os.path.basename(versions[-1]))
    return bool(versions)


def _amrfinder_versions(target: str) -> List[str]:
    """The dated version directories under *target*, oldest first. Excludes `latest`."""
    glob_pattern = DATABASES["amrfinder"].key_glob or ""
    return sorted(p for p in glob.glob(os.path.join(target, glob_pattern)) if os.path.isdir(p))


def _install_genomad(target: str, small: bool, jobs: int = 8) -> bool:
    """geNomad's marker database, by `genomad download-database`.

    One download, two answers: geNomad classifies a contig as viral, plasmid or chromosomal in one
    pass against this database, so `viruses` and the mobile-element work in `amr` share it rather
    than needing a second plasmid reference. That is also why it is worth installing even for a
    study that never asks about phage.

    The command creates `genomad_db/` *inside* the destination it is given and `genomad end-to-end`
    wants that inner directory, so GENOMAD_DB names it rather than the install directory.
    `--small` is ignored; the markers are the method.
    """
    if _already_installed("genomad", target, small):
        print("  geNomad database already present.")
        return True
    os.makedirs(target, exist_ok=True)
    if _run_downloader("genomad", target, small) != 0:
        return False
    return _already_installed("genomad", target, small)


def _install_checkv(target: str, small: bool, jobs: int = 8) -> bool:
    """CheckV's reference database, by `checkv download_database`.

    The command names the directory after the database release (`checkv-db-v1.5`), and that
    versioned directory is what `checkv end_to_end -d` takes and what CHECKV_DB records: a
    completeness figure is a comparison against this reference, so which release it was is part of
    the number. `--small` is ignored for the same reason it is for GUNC - a trimmed reference makes
    complete genomes look fragmentary, which is a wrong answer rather than a cheap one.
    """
    existing = sorted(glob.glob(os.path.join(target, DATABASES["checkv"].key_glob or "")))
    if existing:
        print("  CheckV database already present: %s" % os.path.basename(existing[-1]))
        return True
    os.makedirs(target, exist_ok=True)
    if _run_downloader("checkv", target, small) != 0:
        return False
    found = sorted(glob.glob(os.path.join(target, DATABASES["checkv"].key_glob or "")))
    if found:
        print("  CheckV database %s" % os.path.basename(found[-1]))
    return bool(found)


def _install_bakta(target: str, small: bool, jobs: int = 8) -> bool:
    """Bakta's database, by `bakta_db download`.

    The one database here where `--small` is a real choice rather than a stand-in: `--type light`
    is 1.3 GB against 30 GB and is a supported Bakta configuration, not a capped sample of one. It
    drops the UPS/IPS protein sets, so genes are annotated from the HMMs and the per-gene identity
    evidence is weaker - fine for a run where Bakta is one input among several, not for publishing
    Bakta's own calls.

    The two land in different directories (`db-light/` and `db/`) under the output directory, which
    is why BAKTA_DB's configured value depends on which was installed, and why installing one does
    not hide the other.
    """
    if _already_installed("bakta", target, small):
        print("  Bakta %s database already present." % ("light" if small else "full"))
        return True
    os.makedirs(target, exist_ok=True)
    if _run_downloader("bakta", target, small) != 0:
        return False
    return _already_installed("bakta", target, small)


def _install_gtdbtk(target: str, small: bool, jobs: int = 8) -> bool:
    """The GTDB-Tk reference package, by the `download-db.sh` that ships in its env.

    GTDB-Tk's own script rather than a URL of ours: the release it needs is pinned by the tool's
    version (see envs/gtdbtk.yaml, held at 2.4.x), the script carries the matching release URL, and
    a hand-picked tarball that disagrees with the installed GTDB-Tk fails after the 110 GB has
    already been downloaded. It takes the destination as its one positional argument.

    Two things about it are worth knowing before starting, because neither is recoverable halfway:
    it refuses a destination that already contains a subdirectory, so a half-finished attempt has
    to be emptied by hand rather than resumed; and it deletes the tarball only after unpacking it,
    so the install needs room for both at once - about twice what the table's figure says.

    Its last step tries `conda env config vars set` to record the path in the env, which fails
    harmlessly outside an activated conda shell: `classify_bins` reads GTDBTK_DATA_PATH out of the
    config, which is what this command writes, so the run does not depend on it having worked.
    """
    if _already_installed("gtdbtk", target, small):
        print("  GTDB-Tk reference data already present.")
        return True
    os.makedirs(target, exist_ok=True)
    if _run_downloader("gtdbtk", target, small) != 0:
        return False
    return _already_installed("gtdbtk", target, small)


def _install_singlem(target: str, small: bool, jobs: int = 8) -> bool:
    """The SingleM metapackage, fetched by SingleM itself.

    ``singlem data`` is the supported route for the same reason ``checkm2 database --download`` is:
    the package lives in a Zenodo record whose DOI SingleM pins per release, it is distributed as a
    ZenodoBackpack with checksums SingleM verifies after unpacking, and a hand-picked URL goes stale
    at the next upgrade.

    ``--small`` is ignored. The metapackage is the marker set and its taxonomy; a trimmed copy would
    not be a faster appraisal but a different reference, and the whole figure is relative to that
    reference.
    """
    existing = glob.glob(os.path.join(target, DATABASES["singlem"].key_glob or ""))
    if existing:
        print("  SingleM metapackage already present: %s" % os.path.basename(existing[0]))
        return True
    os.makedirs(target, exist_ok=True)
    if _run_downloader("singlem", target, small) != 0:
        return False
    return _already_installed("singlem", target, small)


def _install_host(target: str, small: bool, jobs: int = 8) -> bool:
    """The bmtagger index of the host genome: download the FASTA, then build two indexes from it.

    The early return matters more here than anywhere else. read_qc opens the bitmask and the
    srprism index, not hg38.fa, so the FASTA is an input to this installer rather than a product
    of it - and without this check, a user who had reclaimed the 3 GB FASTA after indexing
    re-downloaded 20 GB on every subsequent `install-db host`, to rebuild two indexes that were
    already there.
    """
    if _already_installed("host", target, small):
        print("  bmtagger index already present.")
        return True
    fasta = os.path.join(target, "hg38.fa")
    os.makedirs(target, exist_ok=True)
    if not os.path.isfile(fasta):
        gz = os.path.join(target, "chr21.fa.gz" if small else "hg38.fa.gz")
        if not _fetch(HOST_SMALL_URL if small else HOST_URL, gz):
            return False
        print("  decompressing to hg38.fa")
        with gzip.open(gz, "rb") as fin, open(fasta, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=8 << 20)
        os.remove(gz)

    bitmask = os.path.join(target, "hg38.bitmask")
    if not os.path.isfile(bitmask):
        print("  building the bmtagger bitmask (bmtool)")
        if (
            _mamba_run(
                MODULE_ENVS["read_qc"],
                ["bmtool", "-d", fasta, "-o", bitmask, "-A", "0", "-w", "18"],
            )
            != 0
        ):
            return False
    srprism = os.path.join(target, "hg38.srprism")
    if not os.path.isfile(srprism + ".idx"):
        print("  building the srprism index (this is the slow step)")
        if (
            _mamba_run(
                MODULE_ENVS["read_qc"],
                ["srprism", "mkindex", "-i", fasta, "-o", srprism, "-M", "60000"],
            )
            != 0
        ):
            return False
    return _already_installed("host", target, small)


def _install_kofam(target: str, small: bool, jobs: int = 8) -> bool:
    """KOfam: the HMM profiles and ko_list, their per-profile thresholds, from KEGG's FTP.

    The two belong together: a profile's threshold is what makes a hit a KO call, so they are
    downloaded into one directory and `annotate` reads the thresholds from beside the profiles.
    `--small` is ignored; the profile set is the method.
    """
    from ..kofam import KO_LIST_MD5, KOFAM_URL, ko_list_md5

    ko_list = os.path.join(target, "ko_list.gz")
    if _already_installed("kofam", target, small):
        print("  KOfam already present.")
        return True
    os.makedirs(target, exist_ok=True)
    if not _fetch(KOFAM_URL + "ko_list.gz", ko_list):
        return False
    if ko_list_md5(ko_list) != KO_LIST_MD5:
        print("  note: this ko_list differs from the one this release was validated on.")
    archive = os.path.join(target, "profiles.tar.gz")
    if not _fetch(KOFAM_URL + "profiles.tar.gz", archive):
        return False
    _untar(archive, target)
    return os.path.isdir(os.path.join(target, "profiles"))


def _install_kegg(target: str, small: bool, jobs: int = 8) -> bool:
    """KEGG BRITE ko00001, projected to the per-(KO, pathway) table `pathways` reads.

    The download is recorded in SNAPSHOT.json (date and md5), and every hierarchy row names it:
    BRITE is revised, so a percentage is only comparable with one over the same hierarchy.
    """
    from ..hierarchy import BRITE_URL, HIERARCHY_FILE, project_htext, write_snapshot

    if _already_installed("kegg", target, small):
        print("  KEGG hierarchy already present.")
        return True
    os.makedirs(target, exist_ok=True)
    htext = os.path.join(target, "ko00001.keg")
    if not _fetch(BRITE_URL, htext):
        return False
    rows = project_htext(htext, os.path.join(target, HIERARCHY_FILE))
    if rows == 0:
        print("  the download did not look like BRITE ko00001: no KO rows found")
        return False
    write_snapshot(target, htext, BRITE_URL, rows)
    os.remove(htext)
    print("  %d KO-pathway rows" % rows)
    return True


INSTALLERS: Dict[str, Callable[..., bool]] = {
    "checkm": _install_checkm,
    "taxdump": _install_taxdump,
    "kraken2": _install_kraken2,
    "blast": _install_blast,
    "checkm2": _install_checkm2,
    "gunc": _install_gunc,
    "gtdbtk": _install_gtdbtk,
    "card": _install_card,
    "amrfinder": _install_amrfinder,
    "genomad": _install_genomad,
    "checkv": _install_checkv,
    "bakta": _install_bakta,
    "pfam": _install_pfam,
    "eggnog": _install_eggnog,
    "dbcan": _install_dbcan,
    "singlem": _install_singlem,
    "kofam": _install_kofam,
    "kegg": _install_kegg,
    "host": _install_host,
}

#: Databases whose installer runs a tool out of a conda env instead of fetching a URL itself, so
#: without mamba there is nothing it can even start. Checked up front: finding this out after the
#: other databases have downloaded wastes the hours they took.
TOOL_INSTALLED_DBS: Tuple[str, ...] = (
    "checkm",
    "checkm2",
    "blast",
    "host",
    "dbcan",
    "singlem",
    "gtdbtk",
    "amrfinder",
    "genomad",
    "checkv",
    "bakta",
)


class Outcome(NamedTuple):
    """What one installer produced: enough to configure it, or to say why it did not."""

    name: str
    key: Optional[str]
    value: str
    ok: bool
    #: The installer's captured child output, for quoting under a FAILED line.
    output: List[str]


def _size_label(db: Database, small: bool) -> str:
    """How a database's size is announced, for the plan and the real run alike."""
    if small and db.small_size:
        return "%s, --small" % db.small_size
    return "%s download, %s on disk" % (db.download, db.on_disk)


def _planned_text(db: Database, small: bool, on_disk: bool) -> str:
    """The size string that applies to this install: the capped one under ``--small``.

    A `--small` variant is quoted as a download size and there is no separate measurement of what
    it unpacks to, so it stands for both columns rather than being inflated into a second figure
    nobody measured.
    """
    if small and db.small_size:
        return db.small_size
    return db.on_disk if on_disk else db.download


def _planned_size(db: Database, small: bool, on_disk: bool) -> float:
    """``_planned_text`` in GB, for totalling."""
    return _gigabytes(_planned_text(db, small, on_disk))


def _planned_value(db: Database, target: str, small: bool) -> str:
    """The path that would go into the config for *db* installed in *target*.

    Shared with the dry run, which is why an unresolved ``key_glob`` comes back as the pattern
    rather than as the install directory: a plan that printed the parent would be describing
    something the real run does not do. In a real run the glob always resolves, because an
    installer that leaves no match reports failure and nothing is written.
    """
    if db.key_glob:
        matches = sorted(glob.glob(os.path.join(target, db.key_glob)))
        return matches[-1] if matches else os.path.join(target, db.key_glob)
    key_file = db.small_key_file if (small and db.small_key_file) else db.key_file
    return os.path.join(target, key_file) if key_file else target


def _human_gb(size_gb: float) -> str:
    """GB as bytes-with-a-unit, so a 4 MB plan does not total to "~0 GB"."""
    from .. import limits

    return limits.human(size_gb * 1024**3)


def _disk_line(path: str, need_gb: float) -> str:
    """ "N GB needed, <filesystem> is X% full (ceiling Y%, Z of headroom)".

    One sentence, printed by both the plan and the real run, and built from ``limits.disk_state``
    like ``benchmark fetch``'s disk verdict: a download total means very little without the
    filesystem it is going onto beside it, and this command's totals are large enough that the
    answer is sometimes no.
    """
    from .. import limits

    need = int(need_gb * 1024**3)
    state = limits.disk_state(path, limits.resolve())
    if state is None:
        return "%s needed; could not measure the filesystem holding %s" % (
            limits.human(need),
            path,
        )
    return "%s needed, %s" % (limits.human(need), state.describe())


# --- the dry run -------------------------------------------------------------------------
#
# `metawrap2 --dry-run install-db <name>` prints what it would fetch, where it would put it, and
# which config key it would set - and touches nothing. It used to ignore the flag completely,
# because `_run` goes straight to subprocess and nothing here consulted `command.runner`: the
# first user to type it with `gtdbtk` got 110 GB they were explicitly trying not to commit to.
#
# The plan is composed from the same URL constants and DOWNLOADERS entries the installers use, and
# the "already present" test is the same `_already_installed` they short-circuit on - so this is a
# description of what would happen rather than a second implementation of it.


def _fetch_step(url: str) -> str:
    return "download %s" % url


def _env_step(env: str, argv: List[str]) -> str:
    return "run `%s` in conda env %s" % (" ".join(argv), env)


def _plan_steps(name: str, target: str, small: bool, present: bool) -> List[str]:
    """The steps installing *name* would run, in order."""
    if present:
        # checkm is the only one with anything left to do once its data is there: it re-points
        # CheckM at the directory on every install, because that setting lives inside the conda
        # env and does not survive the env being rebuilt.
        if name == "checkm":
            return [_env_step(MODULE_ENVS["bin_refinement"], ["checkm", "data", "setRoot", target])]
        return []
    if name in DOWNLOADERS:
        env, argv = DOWNLOADERS[name]
        return [_env_step(env, argv(target, small))]
    if name == "checkm":
        return [
            _fetch_step(CHECKM_URL),
            "unpack it into the install directory",
            _env_step(MODULE_ENVS["bin_refinement"], ["checkm", "data", "setRoot", target]),
        ]
    if name == "taxdump":
        return [_fetch_step(TAXDUMP_URL), "unpack it into the install directory"]
    if name == "kraken2":
        return [
            _fetch_step(KRAKEN2_SMALL_URL if small else KRAKEN2_URL),
            "unpack it into the install directory",
        ]
    if name == "blast":
        if small:
            return [
                _fetch_step(BLAST_SMALL_SOURCE),
                "download the complete representative prokaryote genomes it lists, %d at a time"
                % DEFAULT_PARALLEL_DBS,
                "build a BLAST database named nt from them with makeblastdb",
            ]
        return [
            _env_step(MODULE_ENVS["classify_bins"], ["update_blastdb.pl", "--decompress", "nt"])
        ]
    if name == "gunc":
        return [
            _fetch_step(GUNC_URL),
            "gunzip it to %s" % (DATABASES["gunc"].key_file or ""),
        ]
    if name == "pfam":
        return [_fetch_step("%s/%s" % (PFAM_BASE_URL, what)) for what in PFAM_FILES] + [
            _fetch_step("%s/relnotes.txt" % PFAM_BASE_URL),
            "gunzip both HMM files",
        ]
    if name == "card":
        return [_fetch_step(CARD_URL), "unpack it into the install directory"]
    if name == "eggnog":
        return [_fetch_step("%s/%s" % (EGGNOG_BASE_URL, what)) for what in EGGNOG_FILES] + [
            "gunzip the SQLite and the DIAMOND database, unpack the taxonomy beside them",
        ]
    if name == "kofam":
        from ..kofam import KOFAM_URL

        return [
            _fetch_step(KOFAM_URL + "ko_list.gz"),
            _fetch_step(KOFAM_URL + "profiles.tar.gz"),
            "unpack the profiles into the install directory",
        ]
    if name == "kegg":
        from ..hierarchy import BRITE_URL, HIERARCHY_FILE

        return [_fetch_step(BRITE_URL), "project it to %s" % HIERARCHY_FILE]
    if name == "host":
        return [
            _fetch_step(HOST_SMALL_URL if small else HOST_URL),
            "decompress it to hg38.fa",
            _env_step(MODULE_ENVS["read_qc"], ["bmtool", "-d", "hg38.fa", "-o", "hg38.bitmask"]),
            _env_step(MODULE_ENVS["read_qc"], ["srprism", "mkindex", "-i", "hg38.fa"]),
        ]
    # Unreachable while every registered database is planned above, which a test asserts.
    return ["install %s" % name]


def _print_plan(
    names: List[str], directory: str, small: bool, config_path: str, update_config: bool
) -> int:
    """Print the plan for *names* and return the exit code. Writes nothing, runs nothing."""
    print("\n[dry-run] this is the plan. Nothing is downloaded, unpacked, indexed or written.\n")
    download: List[str] = []
    on_disk: List[str] = []
    globbed = False
    for name in names:
        db = DATABASES[name]
        target = os.path.join(directory, name.upper())
        present = _already_installed(name, target, small)
        print("=== %s (%s) -> %s" % (name, _size_label(db, small), target))
        if present:
            print("  already present, would skip")
        else:
            download.append(_planned_text(db, small, on_disk=False))
            on_disk.append(_planned_text(db, small, on_disk=True))
        for step in _plan_steps(name, target, small, present):
            print("  would %s" % step)
        if db.key and update_config:
            value = _planned_value(db, target, small)
            globbed = globbed or "*" in value
            print("  would set %s = %s in %s" % (db.key, value, config_path))
        elif db.key:
            print("  would leave %s unset (--no-config-update)" % db.key)
        else:
            print("  nothing to configure: the tool stores this path itself")

    print()
    if download:
        needed = sum(_gigabytes(size) for size in on_disk)
        print(
            "  %d of %d to install: would download %s, leaving %s on disk."
            % (
                len(download),
                len(names),
                _human_gb(sum(_gigabytes(size) for size in download)),
                _human_gb(needed),
            )
        )
        print("  %s" % _disk_line(directory, needed))
    else:
        print("  Nothing to do: every database asked for is already installed.")
    if globbed:
        print("  A path with a `*` in it is named by the download; the real run resolves it.")
    print("\n  Drop --dry-run to do any of this.\n")
    return 0


# --- writing the config ------------------------------------------------------------------


def _config_path(explicit: Optional[str]) -> str:
    return explicit or find_config() or os.path.expanduser("~/.metawrap2/config.toml")


def _set_config_keys(path: str, values: Dict[str, str]) -> None:
    """Set ``key = "value"`` under [databases] in *path*, creating the file/section if needed.

    Deliberately a line edit rather than a parse-and-rewrite: TOML comments in a user's
    config are worth preserving, and the standard library has no TOML writer.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    lines = []
    if os.path.isfile(path):
        with open(path) as fh:
            lines = fh.read().splitlines()

    if not any(ln.strip() == "[databases]" for ln in lines):
        if lines and lines[-1].strip():
            lines.append("")
        lines.append("[databases]")

    start = next(i for i, ln in enumerate(lines) if ln.strip() == "[databases]")
    end = len(lines)
    for i in range(start + 1, len(lines)):
        if lines[i].strip().startswith("["):
            end = i
            break

    for key, value in values.items():
        entry = '%s = "%s"' % (key, value)
        for i in range(start + 1, end):
            stripped = lines[i].lstrip()
            if stripped.split("=")[0].strip() == key and not stripped.startswith("#"):
                lines[i] = entry
                break
        else:
            lines.insert(end, entry)
            end += 1

    with open(path, "w") as fh:
        fh.write("\n".join(lines).rstrip("\n") + "\n")


def main(argv: List[str]) -> int:
    ap = make_parser(
        prog="metawrap2 install-db",
        description="Download, index, and configure the databases MetaWrap2 modules need.",
    )
    ap.add_argument("names", nargs="*", help="databases to install (see --list)")
    ap.add_argument("--all", action="store_true", help="install every database")
    ap.add_argument("--list", action="store_true", help="show available databases and sizes")
    ap.add_argument(
        "--dir",
        default=os.path.expanduser("~/metawrap2_dbs"),
        help="parent directory to install into (default ~/metawrap2_dbs)",
    )
    ap.add_argument(
        "--small",
        action="store_true",
        help="where a capped-but-real variant exists, install that instead "
        "(much faster; good for testing, not for production analysis)",
    )
    ap.add_argument(
        "-t",
        "--threads",
        type=int,
        default=DEFAULT_PARALLEL_DBS,
        help="how many databases to fetch concurrently (default %d). Also sets "
        "the parallelism of the many small downloads inside a single "
        "database build." % DEFAULT_PARALLEL_DBS,
    )
    ap.add_argument("--config", help="config file to update (default: the one in effect)")
    ap.add_argument(
        "--no-config-update",
        action="store_true",
        help="don't write the resulting paths into the config file",
    )
    args = ap.parse_args(argv)

    if args.list or (not args.names and not args.all):
        _print_list()
        return 0 if args.list else 1

    names = list(DATABASES) if args.all else args.names
    unknown = [n for n in names if n not in DATABASES]
    if unknown:
        print("Unknown database(s): %s\nKnown: %s" % (", ".join(unknown), ", ".join(DATABASES)))
        return 1

    config_path = _config_path(args.config)
    needs_mamba = sorted(n for n in names if n in TOOL_INSTALLED_DBS)
    if not shutil.which("mamba") and needs_mamba:
        if not dry_run():
            print(
                "mamba is not on PATH, and these are downloaded by a tool from a module env: %s\n"
                "Install the envs first: metawrap2 install-env --all" % ", ".join(needs_mamba)
            )
            return 1
        # Not fatal while planning: deciding whether to commit to a 110 GB download is exactly
        # what somebody does before building the envs, and a plan is still the answer they asked
        # for. The note says what would have to happen first.
        print(
            "\n[dry-run] note: mamba is not on PATH, so these could not run yet: %s\n"
            "          `metawrap2 install-env --all` builds the envs they come from."
            % ", ".join(needs_mamba)
        )

    if dry_run():
        return _print_plan(names, args.dir, args.small, config_path, not args.no_config_update)

    installed: Dict[str, str] = {}
    #: What has already been written to the config, so the summary pass does not rewrite it.
    written: Dict[str, str] = {}
    failed: List[str] = []

    def install(name: str) -> Outcome:
        db = DATABASES[name]
        target = os.path.join(args.dir, name.upper())
        note = (
            "  (no smaller variant for this one; installing the full database)"
            if args.small and not db.small_size
            else ""
        )
        _say(
            "=== %s (%s) -> %s" % (name, _size_label(db, args.small), target),
            *([note] if note else []),
        )
        ok, output = _install_capturing(name, target, args.small, args.threads)
        value = _planned_value(db, target, args.small)
        if db.key_glob:
            # The path is not known until the download names it; take the newest match so an
            # upgraded metapackage beside an old one points at the new one.
            matches = sorted(glob.glob(os.path.join(target, db.key_glob)))
            if len(matches) > 1:
                _say(
                    "  %s holds %d versions of this database; configuring %s"
                    % (target, len(matches), os.path.basename(value))
                )
        return Outcome(name, db.key, value, ok, output)

    workers = max(1, min(args.threads, len(names)))
    # Only what is actually missing counts towards the headroom line: on a re-run after a failure
    # most of these are already on the disk being measured, and totalling them again would warn
    # about space that is in use by the very databases it is about to skip.
    missing = [
        n for n in names if not _already_installed(n, os.path.join(args.dir, n.upper()), args.small)
    ]
    needed = sum(_planned_size(DATABASES[n], args.small, on_disk=True) for n in missing)
    print(
        "\nInstalling %d database(s) into %s with %d worker(s).\n%s\n"
        % (len(names), args.dir, workers, _disk_line(args.dir, needed))
    )
    # A scheduler's SIGTERM, or a dropped SSH session's SIGHUP, now unwinds like Ctrl-C so the
    # cleanup below runs. Without it those signals killed this process outright and every download
    # it had started carried on in its own session - the orphan, by another route.
    _interrupt.install_handlers()
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        try:
            for done in pool.map(install, names):
                if done.ok:
                    _say("  OK: %s installed" % done.name)
                    if done.key:
                        installed[done.key] = done.value
                        if not args.no_config_update:
                            # Written here, as this database lands, not once after the whole pool
                            # has finished. `--all` is hundreds of gigabytes over many hours, so
                            # being interrupted is the ordinary case rather than the exception: a
                            # user who stopped a 9-database run after 4 had completed - 17 GB on
                            # disk - got a config that named none of them and a TAXDUMP still
                            # holding the `config init` placeholder. Every key that got as far as
                            # being installed is now on disk before the next download starts.
                            #
                            # Safe to do per key despite -t N: pool.map hands results back to this
                            # thread, so the read-modify-write of the config file stays serialised.
                            _set_config_keys(config_path, {done.key: done.value})
                            written[done.key] = done.value
                else:
                    _say("  FAILED: %s" % done.name, *_failure_tail(done.output))
                    failed.append(done.name)
        finally:
            # Before the pool is joined, not after: the worker threads are blocked in
            # ``child.wait()`` on children that no longer share this process's session, so a
            # Ctrl-C does not reach them and the executor's own shutdown would wait for the
            # downloads to finish - turning an interrupt into "it ignored me". Killing the groups
            # first is what lets those waits return. A normal finish leaves nothing to kill.
            _terminate_downloads()

    if installed and not args.no_config_update:
        # Only what the incremental writes above did not already cover - normally nothing, but a
        # key whose value changed is rewritten rather than left stale.
        pending = {k: v for k, v in installed.items() if written.get(k) != v}
        if pending:
            _set_config_keys(config_path, pending)
        print("\nUpdated %s:" % config_path)
        for key, value in installed.items():
            print("  %s = %s" % (key, value))

    print()
    if failed:
        # The retry advice is here because the idempotency it relies on is invisible: every
        # installer short-circuits on an already-complete database, so re-running the same command
        # retries only what is missing. Without saying so, a transient failure in one database of
        # nine reads like "start the whole thing again".
        print(
            "%d database(s) failed: %s. Re-run the same command to retry; completed databases "
            "are detected and skipped." % (len(failed), ", ".join(failed))
        )
        return 1
    print("All requested databases are installed. Verify with: metawrap2 check\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
