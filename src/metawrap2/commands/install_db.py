"""`metawrap2 install-db [name ...] [--all]`: download and configure MetaWrap2's databases.

For each database this downloads, unpacks and indexes it, writes its key into ``metawrap2.toml``,
runs any tool-side setup (such as ``checkm data setRoot``) and verifies the result. Steps that
already look complete are skipped, so re-running is safe.

    metawrap2 install-db --list            # what is available, how big, who needs it
    metawrap2 install-db checkm taxdump    # just these two
    metawrap2 install-db --all             # everything (large: see --list for sizes)

``--all`` takes many hours and several hundred GB. ``--small`` installs a capped variant where one
exists, for smoke tests.

An interrupted run can be restarted: each config key is written as its database finishes, and
installers skip databases that are already complete. ``--list`` shows the download size and the
on-disk size of each database.
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
from ..benchmark import fetch as _fetcher
from ..cliformat import make_parser
from ..config import MODULE_ENVS, OPTIONAL_ENVS, env_prefix, find_config

# Databases are independent and mostly network-bound, so several are fetched at once. The default
# is modest so a slow link is not saturated.
DEFAULT_PARALLEL_DBS = 3
# Serialises the progress prints from concurrent installers so lines do not interleave.
_PRINT_LOCK = threading.Lock()

# Lines of a failed installer's output shown under its FAILED line, and lines kept to show them.
# Child output is captured per database so concurrent installers (`-t N`) do not interleave and a
# failure can be explained. A bounded buffer is used because curl's progress meter writes for as
# long as the download runs.
FAILURE_TAIL_LINES = 15
_CAPTURE_LINES = 400
#: Ring buffer for the child output of the installer on this thread, or None outside an
#: installer. Thread-local so concurrent installers keep their output separate.
_CAPTURED = threading.local()


def _say(*lines: str) -> None:
    with _PRINT_LOCK:
        for line in lines:
            print(line)
        sys.stdout.flush()


class Database(NamedTuple):
    key: Optional[str]  # metawrap2.toml [databases] key, None if the tool stores its own path
    used_by: str
    #: Download size, and size on disk after unpacking and indexing; these often differ (gunc
    #: downloads 6.7 GB and uses about 13 GB). ``on_disk`` is marked ``~`` when it is an estimate.
    download: str
    on_disk: str
    small_size: Optional[str]
    describe: str
    #: Path under the install directory to write into the config instead of the directory itself:
    #: a file the tool takes (GUNC's ``--db_file`` .dmnd) or a subdirectory the downloader creates
    #: (geNomad's ``genomad_db/``).
    key_file: Optional[str] = None
    #: ``key_file`` for a ``--small`` install, when it differs: ``bakta_db download --type light``
    #: writes ``db-light/`` and ``--type full`` writes ``db/``.
    small_key_file: Optional[str] = None
    #: Like ``key_file``, for a path whose name is only known after download. SingleM's
    #: metapackage directory includes its version (``S6.5.0.GTDB_r232.metapackage_...``); the
    #: versioned path is configured so the reference version stays visible.
    key_glob: Optional[str] = None

    @property
    def size(self) -> str:
        """The download size. Kept for ``quickstart._print_database_guidance``, which reads it."""
        return self.download


# Sizes are measurements where available; otherwise the download size is reused for on-disk.
# Sources: gunc's ~13 GB is a reported on-disk measurement, checkm2's ~3 GB and pfam's 1.6 GB come
# from their release notes, singlem's 9.7 GB is the unpacked metapackage, and the tool-managed
# databases (bakta, gtdbtk, amrfinder, genomad, checkv) use the figures from the README table.
DATABASES: Dict[str, Database] = {
    "checkm": Database(
        key="CHECKM_DB",
        used_by="--checkm1 in bin_refinement and binning --run-checkm",
        download="1.4 GB",
        on_disk="~1.4 GB",
        small_size=None,
        describe="CheckM1 marker-gene reference data (only for --checkm1; CheckM2 is the default)",
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
        used_by="bin_refinement, reassemble_bins, binning --run-checkm (the default scorer)",
        download="1.7 GB",
        on_disk="~3 GB",
        small_size=None,
        describe="CheckM2's DIAMOND reference (uniref100.KO; separate from CheckM1's `checkm`)",
        key_file=os.path.join("CheckM2_database", "uniref100.KO.1.dmnd"),
    ),
    "gunc": Database(
        key="GUNC_DB",
        used_by="bin_refinement (unless --skip-gunc)",
        download="6.7 GB",
        on_disk="~13 GB",
        small_size=None,
        describe="GUNC proGenomes 2.1 reference (diamond .dmnd, about twice the download size "
        "once gunzipped; GUNC runs by default)",
        key_file="gunc_db_progenomes2.1.dmnd",
    ),
    "gtdbtk": Database(
        key="GTDBTK_DATA_PATH",
        used_by="classify_bins --gtdbtk",
        download="~110 GB",
        on_disk="~110 GB",
        small_size=None,
        describe="GTDB-Tk reference package, installed by its download-db.sh (the tarball is "
        "deleted after unpacking, so the install needs about twice this space free)",
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
        # The release is not named here; it is defined once, in EGGNOG_DB_VERSION below.
        describe="eggNOG annotation SQLite (39 GB unpacked) and the DIAMOND database of eggNOG "
        "proteins (8.7 GB)",
    ),
    "bakta": Database(
        key="BAKTA_DB",
        used_by="annotate --bakta, annotate_bins --bakta",
        download="30 GB",
        on_disk="~30 GB",
        small_size="1.3 GB",
        describe="Bakta's reference database, installed by bakta_db (--small: the `light` type, "
        "without the UPS/IPS protein sets; annotates from the HMMs only)",
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
        describe="AMRFinderPlus reference database, installed by amrfinder_update, which also "
        "builds its BLAST and HMMER indexes; the dated version directory is configured",
        key_glob="20[0-9][0-9]-[0-9][0-9]-[0-9][0-9]*",
    ),
    "genomad": Database(
        key="GENOMAD_DB",
        used_by="viruses, amr --mobile-elements",
        download="1.5 GB",
        on_disk="~1.5 GB",
        small_size=None,
        describe="geNomad marker database, installed by `genomad download-database`; used for "
        "the viral calls and for the plasmid calls that amr reads",
        key_file="genomad_db",
    ),
    "checkv": Database(
        key="CHECKV_DB",
        used_by="viruses (unless --skip-checkv)",
        download="6.4 GB",
        on_disk="~6.4 GB",
        small_size=None,
        describe="CheckV reference database, installed by `checkv download_database`; the "
        "versioned checkv-db-vX directory it creates is configured",
        key_glob="checkv-db-*",
    ),
    "singlem": Database(
        key="SINGLEM_METAPACKAGE",
        used_by="appraise",
        download="2.2 GB",
        on_disk="9.7 GB",
        small_size=None,
        describe="SingleM GTDB marker metapackage (single-copy marker genes + their taxonomy). "
        "Downloads at 150-300 KB/s, so the 2.2 GB takes 2-4 hours (the other databases here "
        "download at 0.5-20 MB/s); start it early",
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
        "srprism index are built locally; building takes longer than the download",
    ),
}

#: Unit factors to GB, for parsing the registry's human-readable size strings into totals.
_SIZE_UNITS: Dict[str, float] = {"MB": 1.0 / 1024, "GB": 1.0, "TB": 1024.0}


def _gigabytes(size: str) -> float:
    """Parse one of this registry's size strings ("399 MB", "~110 GB") into GB.

    A leading `~` is ignored, so estimates are included in totals.
    """
    number, _, unit = size.replace("~", "").strip().partition(" ")
    return float(number) * _SIZE_UNITS[unit.strip()]


def _total(sizes: List[str]) -> str:
    """Sum size strings for the totals row, always marked `~` because the inputs are mixed."""
    return "~%d GB" % round(sum(_gigabytes(size) for size in sizes))


def _print_list() -> None:
    """Print the database table with download, on-disk and ``--small`` sizes, and totals."""
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
        "\n  `download` is the size fetched over the network; `on disk` is the size after "
        "unpacking\n  and indexing. An on-disk figure marked `~` is an estimate; where none was "
        "measured, the\n  download size is used, so the on-disk total is a lower bound. "
        "`--small` is the download\n  size with that flag; for a database with no smaller "
        "variant it is the full download."
    )
    print("\n  Install into a directory of your choice with --dir (default ~/metawrap2_dbs).")
    print("  Each database's path is written to your config file as it finishes.\n")


# --- shared helpers ---------------------------------------------------------------------


#: Process group of each child started by :func:`_run` that is still running. An interrupt kills
#: the whole group, since the curl writing the file may be several processes below the child.
_DOWNLOAD_GROUPS: Set[int] = set()
_GROUPS_LOCK = threading.Lock()


@contextlib.contextmanager
def _tracked(child: subprocess.Popen) -> Iterator[None]:
    """Remember *child*'s process group for as long as it runs.

    The child is started with ``start_new_session``, so ``child.pid`` is its group id. The new
    session keeps a terminal Ctrl-C from reaching curl directly; tracking the group lets
    :func:`_terminate_downloads` kill it instead, so no download is orphaned.
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

    Prevents a curl from outliving an interrupted `install-db`. Two curls writing the same `.part`
    corrupt it, which later shows up as a `tar` or `gzip` failure. Groups that already exited are
    ignored.
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

    Outside an installer the child writes to this process's stdout. Inside one, stderr is merged
    into stdout and captured so a failure can be quoted.

    Every child runs in its own session (see :func:`_tracked`) and is started with Popen so its pid
    is available for cleanup.
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

    An exception's text is appended to the captured lines as well as printed. Only `Exception` is
    caught, so Ctrl-C stops the whole run.
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


#: Seconds a `.part` must be unmodified before it is resumed.
#:
#: `curl -C -` appends from the current end of the file, so two curls resuming the same `.part`
#: interleave and corrupt it, which later shows up as a `tar`/`gzip` failure. A live curl's quiet
#: gaps are much shorter than this (`--retry 3` backs off 1s, 2s, 4s). Refusing a dead `.part`
#: only costs a message naming the file to delete.
LIVE_PART_SECONDS = 120


def _another_writer(part: str) -> bool:
    """Does *part* look like a file another `install-db` is still writing? Prints the reason.

    Uses the file's mtime instead of a lock file, because the writer may be an orphaned curl
    that holds no lock.
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
    """Download *url* to *dest*, verified, resuming a partial file that no other process is writing.

    The transfer is ``benchmark.fetch.download``: several ranged connections with aria2c when it is
    installed, curl otherwise, with the same checks either way. The databases fetched by URL here
    are ``.gz``, ``.tar.gz`` or ``.tar.bz2`` files without published checksums, so each one is
    test-decompressed before it is accepted, and a corrupt one fails here instead of at unpacking.

    ``_run`` is passed in so the child still leads its own process group and writes into this
    thread's capture buffer. The live ``.part`` check below stays here: it guards against two
    copies of this command downloading the same file.
    """
    if os.path.isfile(dest) and os.path.getsize(dest):
        print("  already downloaded: %s" % os.path.basename(dest))
        return True
    tmp = _fetcher.partial_path(dest)
    if os.path.isfile(tmp) and _another_writer(tmp):
        return False
    problem = _fetcher.download(url, dest, run=_run)
    if problem:
        print("  download failed: %s\n  %s" % (url, problem))
        return False
    return True


def _extract(archive: str, into: str) -> bool:
    """Unpack *archive* into *into*, refusing members that would land outside it.

    A bare `extractall` writes wherever a member's name points, including absolute and `..`
    paths. Python 3.12's `data` filter rejects those and strips ownership and mode bits; the
    explicit check below does the same on older versions.

    All archives in this module are unpacked here. Compression (gz, bz2) is detected from the
    file. A separate check in the CARD installer once rejected valid archives, so there is only
    this one.
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
                    # Skip device nodes and fifos.
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

    A failed unpack keeps the archive, which is usually a truncated download that `_fetch` can
    resume on the next run.
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

    Checks the link target as well as the name, since a symlink out of the directory lets later
    members be written through it.
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

    Read from ``command.runner``, which the global ``--dry-run`` flag sets and every module reads
    through ``modules._common.dry_run``. Before this check, ``metawrap2 --dry-run install-db
    gtdbtk`` downloaded 110 GB.
    """
    from ..command import runner

    return bool(runner.dry_run)


# --- what each installer fetches or runs -------------------------------------------------
#
# URLs and downloader argv are defined here so the installers and `--dry-run` read the same values.

CHECKM_URL = "https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
KRAKEN2_URL = "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz"
#: The standard build capped to 8 GB, for testing.
KRAKEN2_SMALL_URL = "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz"
GUNC_URL = "https://swifter.embl.de/~fullam/gunc/gunc_db_progenomes2.1.dmnd.gz"
PFAM_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
PFAM_FILES = ("Pfam-A.hmm.gz", "Pfam-A.hmm.dat.gz")
CARD_URL = "https://card.mcmaster.ca/latest/data"

#: eggNOG's emapper database files, fetched from eggnog5.embl.de. eggNOG-mapper's own
#: `download_eggnog_data.py` uses `eggnogdb.embl.de`, which no longer responds, and reports
#: "Finished" over zero-byte files because it does not check wget's exit status. Fetching directly
#: also gives this module's resume and `.part` handling.
EGGNOG_DB_VERSION = "5.0.2"
EGGNOG_BASE_URL = "http://eggnog5.embl.de/download/emapperdb-%s" % EGGNOG_DB_VERSION
#: Files `-m diamond` needs: annotation SQLite, taxonomy and DIAMOND protein database. All are
#: gzipped on the server except the taxonomy, which is a tarball.
EGGNOG_FILES: Tuple[str, ...] = (
    "eggnog.db.gz",
    "eggnog.taxa.tar.gz",
    "eggnog_proteins.dmnd.gz",
)
HOST_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
#: One human chromosome, enough to test host removal.
HOST_SMALL_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz"
#: Source list for the --small BLAST build; the full build is fetched by update_blastdb.pl.
BLAST_SMALL_SOURCE = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"

#: Databases fetched by their own tool: the env that has the tool, and the argv as a function of
#: ``--small``. Read by both the installer and the dry run. Env names come from config.py so they
#: follow env renames.
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

    Checks first that the env exists, since `quickstart` builds only the default envs, and prints
    the `install-env` command to run if it does not. Otherwise mamba reports only a missing prefix.
    """
    env, argv = DOWNLOADERS[name]
    if not dry_run() and env_prefix(env) is None:
        print(
            "  cannot install %s: its downloader is in the %s environment, which is not "
            "built yet.\n  Run `metawrap2 install-env %s`, then run this command again."
            % (name, env, env.replace("metawrap2-", ""))
        )
        return 1
    return _mamba_run(env, argv(target, small))


# --- is it already installed? ------------------------------------------------------------
#
# Used by the installers, to skip complete databases, and by the dry run, to report them.
#
# A database is installed when the path that would go into the config exists and is not empty.
# Databases whose key is the install directory itself list marker files here.
PRESENT_WHEN: Dict[str, Tuple[str, ...]] = {
    "checkm": ("taxon_marker_sets.tsv",),
    "taxdump": ("names.dmp",),
    "kraken2": ("hash.k2d",),
    "dbcan": ("dbCAN.hmm",),
    "kofam": ("ko_list.gz", "profiles"),
    #: Installed when the destination holds a subdirectory, as GTDB-Tk's download-db.sh checks.
    "gtdbtk": ("*/",),
    #: The index files read_qc opens. hg38.fa is only an input to indexing.
    "host": ("hg38.bitmask", "hg38.srprism.idx"),
    #: Both files emapper needs, so a download interrupted between them is not complete.
    #: EGGNOG_DB is the directory, because emapper takes --data_dir.
    "eggnog": ("eggnog.db", "eggnog_proteins.dmnd"),
}


def _already_installed(name: str, target: str, small: bool) -> bool:
    """Is *name* complete in *target*? Never writes, so it is safe under a dry run."""
    if name == "blast":
        # BLAST's files are named by the configurable BLASTDB_NAME; blastdb.find handles that.
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

    Globbed because some databases are named by release (`checkv-db-v1.5`, `*.smpkg.zb`,
    AMRFinderPlus's dated directories). Empty matches are rejected because an interrupted download
    leaves its directory behind.
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

    # Pipeline calls pass CHECKM_DATA_PATH from CHECKM_DB (see checkm.checkm1_data_env), which
    # survives env rebuilds. setRoot is also run so `checkm` works when called by hand from the env.
    # The env name comes from MODULE_ENVS so it follows env renames; a failure here is only a NOTE.
    checkm_env = MODULE_ENVS["bin_refinement"]
    print("  pointing CheckM at the data (checkm data setRoot)")
    if _mamba_run(checkm_env, ["checkm", "data", "setRoot", target]) != 0:
        print(
            "  NOTE: could not run `checkm data setRoot`. bin_refinement still works because\n"
            "        CHECKM_DATA_PATH is passed on each call. To use `checkm` directly, run:\n"
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
    # Skip if present, so a re-run does not rebuild the database with makeblastdb.
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
        # Thousands of small genome downloads, run concurrently.
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

    Uses ``checkm2 database --download``, because CheckM2 pins a versioned Zenodo record and a
    fixed URL would go stale on upgrade.

    This is separate from the ``checkm`` (CheckM1) database; they share no files. ``--small`` is
    ignored, because CheckM2's model is trained against this exact reference.
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

    A single file with no indexing step. ``--small`` is ignored, because a reduced reference would
    leave more bins unassessable.
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

    Uses dbCAN's pinned-release download source (the flag set in DOWNLOADERS) instead of the default
    HTTP source, which dbCAN's help describes as a moving snapshot. CAZyme calls are only
    reproducible against a fixed release.

    ``--small`` is ignored: the three-method agreement reported by annotate needs all databases.
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

    See EGGNOG_FILES for why eggNOG-mapper's downloader is not used. The files are unpacked into
    the flat layout emapper's `--data_dir` expects (12 GB download, 48 GB unpacked).

    ``--small`` is ignored: eggNOG has no reduced variant.
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
            # emapper expects `eggnog.taxa.db` beside the other files, not in a subdirectory.
            _flatten_into(target, os.path.splitext(os.path.splitext(name)[0])[0])
            continue
        if not _gunzip(archive, archive[: -len(".gz")]):
            return False
    print("  eggNOG database %s installed." % EGGNOG_DB_VERSION)
    return _already_installed("eggnog", target, small)


def _flatten_into(target: str, subdir: str) -> None:
    """Move anything the taxonomy tarball unpacked into *subdir* up into *target*, then remove it.

    emapper looks for `eggnog.taxa.db` beside `eggnog.db`.
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
        # Not empty; leave the remaining files in place.
        pass


def _install_pfam(target: str, small: bool, jobs: int = 8) -> bool:
    """Pfam-A, from EBI's current release.

    The HMMs include their gathering thresholds, so no separate cutoff table is needed.

    `relnotes.txt` is also fetched to record the release, since Pfam renumbers and merges families
    between releases. `--small` is ignored.
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

    Fetched by URL because RGI has no downloader (`rgi load` takes a card.json already on disk).
    The tarball is 4 MB, so there is no `--small` variant.

    The version in card.json is printed here and recorded per run in the AMR table. The URL always
    serves the newest release, so installs made at different times can differ under the same path;
    rename the directory after the version if needed.
    """
    marker = os.path.join(target, DATABASES["card"].key_file or "")
    if _already_installed("card", target, small):
        print("  CARD database already present (version %s)." % _card_version(marker))
        return True
    os.makedirs(target, exist_ok=True)
    archive = os.path.join(target, "card_data.tar.bz2")
    if not _fetch(CARD_URL, archive):
        return False
    # Uses the shared `_extract`. CARD's tarball starts with a `.` member, which an earlier
    # realpath-prefix check in this installer rejected, so CARD could not be installed.
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


# --- databases fetched by their own tool's downloader ------------------------------------
#
# amrfinder, genomad, checkv, bakta and gtdbtk use the tool's own downloader, as checkm2 does: the
# tool knows which database version it reads, and several build indexes after downloading. This
# command runs the downloader in the right env, records the resolved path in the config, and skips
# databases already installed.
#
# Completeness is checked as "the target directory is not empty", not by a marker filename, because
# upstream layouts change between releases. The tool reports an incomplete database when a module
# first uses it.


def _has_content(path: str) -> bool:
    return os.path.isdir(path) and bool(os.listdir(path))


def _subdirectories(path: str) -> List[str]:
    return sorted(p for p in glob.glob(os.path.join(path, "*")) if os.path.isdir(p))


def _install_amrfinder(target: str, small: bool, jobs: int = 8) -> bool:
    """The AMRFinderPlus reference database, fetched by amrfinder_update.

    `amrfinder_update` picks the database version the installed AMRFinderPlus can read and runs
    makeblastdb and hmmpress, which `amr` needs.

    It writes a dated directory per version (`2024-07-22.1`) plus a `latest` symlink. The dated
    directory is configured (see AMRFINDER_DB in config.py) so results can be tied to a version.
    `--small` is ignored.
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

    geNomad classifies contigs as viral, plasmid or chromosomal in one pass, so `viruses` and
    `amr --mobile-elements` share this database.

    The command creates `genomad_db/` inside the destination, and GENOMAD_DB names that inner
    directory, which `genomad end-to-end` takes. `--small` is ignored.
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

    The command names the directory after the release (`checkv-db-v1.5`). CHECKV_DB records that
    versioned directory, which `checkv end_to_end -d` takes. `--small` is ignored, because a reduced
    reference would make complete genomes look fragmentary.
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

    `--small` installs `--type light` (1.3 GB instead of 30 GB), a supported Bakta configuration.
    It drops the UPS/IPS protein sets, so genes are annotated from the HMMs with weaker per-gene
    evidence. Suitable when Bakta is one input among several, not for publishing Bakta's calls.

    The two types install into `db-light/` and `db/`, so BAKTA_DB depends on which was installed
    and both can coexist.
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

    The script carries the release URL matching the installed GTDB-Tk (pinned to 2.4.x in
    envs/gtdbtk.yaml); a mismatched tarball would fail after the 110 GB download. It takes the
    destination as its only argument.

    It refuses a destination that already contains a subdirectory, so a partial attempt must be
    emptied by hand. It deletes the tarball only after unpacking, so the install needs about twice
    the listed size.

    Its last step, `conda env config vars set`, fails harmlessly outside an activated conda shell;
    `classify_bins` reads GTDBTK_DATA_PATH from the config this command writes.
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

    ``singlem data`` fetches the Zenodo record SingleM pins per release and verifies its checksums.
    A fixed URL would go stale on upgrade.

    ``--small`` is ignored: appraise figures are relative to the full reference.
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

    Completeness is judged by the two indexes read_qc opens, not by hg38.fa, so deleting the FASTA
    after indexing does not trigger a 20 GB re-download.
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

    Both are downloaded into one directory; `annotate` reads the thresholds from beside the
    profiles. `--small` is ignored.
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

    The download date and md5 are recorded in SNAPSHOT.json and referenced by every hierarchy row,
    because BRITE is revised over time.
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

#: Databases whose installer runs a tool from a conda env and so needs mamba. Checked before any
#: download starts.
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
    """Result of one installer: the config key and value, success, and captured output."""

    name: str
    key: Optional[str]
    value: str
    ok: bool
    #: The installer's captured child output, for quoting under a FAILED line.
    output: List[str]


def _size_label(db: Database, small: bool) -> str:
    """Size label printed for a database, in both the plan and the real run."""
    if small and db.small_size:
        return "%s, --small" % db.small_size
    return "%s download, %s on disk" % (db.download, db.on_disk)


def _planned_text(db: Database, small: bool, on_disk: bool) -> str:
    """The size string that applies to this install: the capped one under ``--small``.

    A `--small` size is only known as a download size, so it is used for both columns.
    """
    if small and db.small_size:
        return db.small_size
    return db.on_disk if on_disk else db.download


def _planned_size(db: Database, small: bool, on_disk: bool) -> float:
    """``_planned_text`` in GB, for totalling."""
    return _gigabytes(_planned_text(db, small, on_disk))


def _planned_value(db: Database, target: str, small: bool) -> str:
    """The path that would go into the config for *db* installed in *target*.

    Also used by the dry run, so an unresolved ``key_glob`` is returned as the pattern. In a real
    run the glob resolves, because an installer that leaves no match reports failure.
    """
    if db.key_glob:
        matches = sorted(glob.glob(os.path.join(target, db.key_glob)))
        return matches[-1] if matches else os.path.join(target, db.key_glob)
    key_file = db.small_key_file if (small and db.small_key_file) else db.key_file
    return os.path.join(target, key_file) if key_file else target


def _human_gb(size_gb: float) -> str:
    """Format GB with a suitable unit, so small totals do not print as "~0 GB"."""
    from .. import limits

    return limits.human(size_gb * 1024**3)


def _disk_line(path: str, need_gb: float) -> str:
    """ "N GB needed, <filesystem> is X% full (ceiling Y%, Z of headroom)".

    Printed by both the plan and the real run, built from ``limits.disk_state`` as in
    ``benchmark fetch``.
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
# `metawrap2 --dry-run install-db <name>` prints what would be fetched, where it would go and
# which config key would be set, without changing anything. The plan uses the same URL constants,
# DOWNLOADERS entries and `_already_installed` check as the installers.


def _fetch_step(url: str) -> str:
    return "download %s" % url


def _env_step(env: str, argv: List[str]) -> str:
    return "run `%s` in conda env %s" % (" ".join(argv), env)


def _plan_steps(name: str, target: str, small: bool, present: bool) -> List[str]:
    """The steps installing *name* would run, in order."""
    if present:
        # checkm reruns setRoot on every install, because the setting lives in the conda env and
        # is lost when the env is rebuilt.
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
    # Unreachable while every registered database is handled above (checked by a test).
    return ["install %s" % name]


def _print_plan(
    names: List[str], directory: str, small: bool, config_path: str, update_config: bool
) -> int:
    """Print the plan for *names* and return the exit code. Writes nothing, runs nothing."""
    print("\n[dry-run] Plan only: nothing is downloaded, unpacked, indexed or written.\n")
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
        print("  A `*` in a path stands for a name set by the download; the real run fills it in.")
    print("\n  Run again without --dry-run to install.\n")
    return 0


# --- writing the config ------------------------------------------------------------------


def _config_path(explicit: Optional[str]) -> str:
    return explicit or find_config() or os.path.expanduser("~/.metawrap2/config.toml")


def _set_config_keys(path: str, values: Dict[str, str]) -> None:
    """Set ``key = "value"`` under [databases] in *path*, creating the file/section if needed.

    Edits lines instead of parsing and rewriting, to keep the user's TOML comments; the standard
    library also has no TOML writer.
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
        help="install the smaller variant where one exists "
        "(much faster; for testing, not for production analysis)",
    )
    ap.add_argument(
        "-t",
        "--threads",
        type=int,
        default=DEFAULT_PARALLEL_DBS,
        help="number of databases to fetch at once (default %d). Also sets "
        "the number of parallel downloads within a single "
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
                "mamba is not on PATH; a tool in a module env downloads these databases: %s\n"
                "Install the envs first: metawrap2 install-env --all" % ", ".join(needs_mamba)
            )
            return 1
        # Not fatal for a dry run: print the plan with a note that the envs are needed first.
        print(
            "\n[dry-run] note: mamba is not on PATH, so these could not run yet: %s\n"
            "          `metawrap2 install-env --all` builds the envs they come from."
            % ", ".join(needs_mamba)
        )

    if dry_run():
        return _print_plan(names, args.dir, args.small, config_path, not args.no_config_update)

    installed: Dict[str, str] = {}
    #: Keys already written to the config, so the summary pass does not rewrite them.
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
            # The path is named by the download; use the newest match when several versions exist.
            matches = sorted(glob.glob(os.path.join(target, db.key_glob)))
            if len(matches) > 1:
                _say(
                    "  %s holds %d versions of this database; configuring %s"
                    % (target, len(matches), os.path.basename(value))
                )
        return Outcome(name, db.key, value, ok, output)

    workers = max(1, min(args.threads, len(names)))
    # Only missing databases count towards the space needed; installed ones are already on disk.
    missing = [
        n for n in names if not _already_installed(n, os.path.join(args.dir, n.upper()), args.small)
    ]
    needed = sum(_planned_size(DATABASES[n], args.small, on_disk=True) for n in missing)
    print(
        "\nInstalling %d database(s) into %s with %d worker(s).\n%s\n"
        % (len(names), args.dir, workers, _disk_line(args.dir, needed))
    )
    # Make SIGTERM and SIGHUP unwind like Ctrl-C, so the cleanup below kills running downloads.
    _interrupt.install_handlers()
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        try:
            for done in pool.map(install, names):
                if done.ok:
                    _say("  OK: %s installed" % done.name)
                    if done.key:
                        installed[done.key] = done.value
                        if not args.no_config_update:
                            # Written as each database finishes, so an interrupted run keeps
                            # the keys of completed databases. pool.map returns results on this
                            # thread, so config writes stay serialised under -t N.
                            _set_config_keys(config_path, {done.key: done.value})
                            written[done.key] = done.value
                else:
                    _say("  FAILED: %s" % done.name, *_failure_tail(done.output))
                    failed.append(done.name)
        finally:
            # Runs before the pool is joined: workers wait on children in their own sessions,
            # which Ctrl-C does not reach, so the pool shutdown would otherwise wait for the
            # downloads to finish.
            _terminate_downloads()

    if installed and not args.no_config_update:
        # Write any keys not already written above with the same value.
        pending = {k: v for k, v in installed.items() if written.get(k) != v}
        if pending:
            _set_config_keys(config_path, pending)
        print("\nUpdated %s:" % config_path)
        for key, value in installed.items():
            print("  %s = %s" % (key, value))

    print()
    if failed:
        # Tell the user a re-run retries only the failed databases.
        print(
            "%d database(s) failed: %s. Re-run the same command to retry; completed databases "
            "are detected and skipped." % (len(failed), ", ".join(failed))
        )
        return 1
    print("All requested databases are installed. Verify with: metawrap2 check\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
