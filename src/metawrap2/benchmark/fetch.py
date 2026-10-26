"""Download benchmark data: planned, checksum-verified and within the disk budget.

* :func:`plan` lists what is missing, its size, and whether the ``[limits]`` disk ceiling has
  room, without network access. A CAMI2 body site is about 95 GB.
* Files are verified by size first (catches truncated downloads quickly), then by md5, or by
  test-decompressing when there is no md5. Size alone is not enough: a resumed curl download once
  had exactly the expected size and a corrupt gzip stream. A file that fails is deleted.
* :func:`download` uses several ranged connections with aria2c when it is installed, and curl
  otherwise. ENA limits bandwidth per connection, so this is several times faster on large read
  sets. Both paths run the same checks.
* Extraction is selective. About 95% of each CAMI2 archive is per-OTU BAM alignments that no
  binning benchmark reads. Skipping them does not reduce the download (a ``.tar.gz`` is one gzip
  stream) but limits what is written to the reads and contigs. Members with absolute or escaping
  paths are rejected, since ``extractall`` writes wherever a member's path points.

Nothing downloads on import, ``list`` or ``plan``. A fetch starts only when a user names a dataset
and confirms its size.
"""

from __future__ import annotations

import bz2
import fnmatch
import gzip
import hashlib
import os
import shutil
import subprocess
import tarfile
import time
import zlib
from dataclasses import dataclass, field
from typing import Any, Callable, List, NamedTuple, Optional, Sequence, Tuple

from .. import limits as _limits
from .registry import Dataset, Remote

#: Read size for hashing: few syscalls on large files, small memory use.
_CHUNK = 8 * 1024 * 1024

#: Returned by :func:`verify` when a file is intact.
OK = ""

#: Where bytes land before they have been checked. Nothing is renamed to the real name until a
#: check has passed, so a half file can never be read by something else as a whole one.
PARTIAL_SUFFIX = ".part"

#: Suffix of the control file aria2c writes beside a download to record which byte ranges have
#: arrived. It is removed with the partial file, since it describes those bytes.
CONTROL_SUFFIX = ".aria2"

#: Connections and pieces for a ranged download. ENA limits bandwidth per connection, so one
#: connection gets 1-3 MB/s and several get 6-10x that. aria2c allows up to 16, but publishers
#: also cap connections per client, and exceeding that cap fails the transfer.
CONNECTIONS = 8

#: Attempts per file, and the first wait between them, which doubles after each failure. curl's
#: own ``--retry`` covers a connection dropped during an attempt; this loop retries a download
#: that completed but failed a check, and deletes the partial file first.
ATTEMPTS = 3
RETRY_DELAY = 5.0


class Problem(NamedTuple):
    """Which check rejected a download, and what that check saw.

    The check is one of transport, size, checksum or decompression. Each needs a different fix,
    so the final message names it.
    """

    check: str
    detail: str

    def __str__(self) -> str:
        return "%s: %s" % (self.check, self.detail)


def md5_of(path: str, progress: Optional[Callable[[int], None]] = None) -> str:
    """Streaming md5 of *path*. Never reads the whole file into memory."""
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(_CHUNK)
            if not chunk:
                break
            digest.update(chunk)
            if progress is not None:
                progress(len(chunk))
    return digest.hexdigest()


def verify(path: str, remote: Remote, *, hash_it: bool = True) -> str:
    """Why *path* is not a good copy of *remote*, or :data:`OK` if it is.

    Size is checked first because it is cheap and catches interrupted downloads, but a file of the
    right size can still be corrupt. Pass ``hash_it=False`` to skip the md5 when the hash has
    already been checked; :func:`download` never accepts a file on size alone.
    """
    if not os.path.isfile(path):
        return "missing"
    actual = os.path.getsize(path)
    if remote.size is not None and actual != remote.size:
        return "wrong size (%d bytes on disk, %d expected)" % (actual, remote.size)
    if actual == 0:
        return "empty"
    if hash_it and remote.md5:
        got = md5_of(path)
        if got != remote.md5:
            return "md5 mismatch (%s on disk, %s expected)" % (got, remote.md5)
    return OK


@dataclass
class FetchPlan:
    """What a fetch would do, before it does any of it."""

    dataset: Dataset
    dest: str
    present: List[Remote] = field(default_factory=list)
    #: Files still needed, each with the reason: "missing", or why an existing copy was rejected,
    #: so the report separates "never downloaded" from "downloaded and corrupt".
    needed: List[Tuple[Remote, str]] = field(default_factory=list)

    @property
    def bytes_needed(self) -> int:
        return sum(r.size or 0 for r, _ in self.needed)

    @property
    def complete(self) -> bool:
        return not self.needed

    def describe(self) -> str:
        lines = ["%s -> %s" % (self.dataset.name, self.dest)]
        if self.present:
            lines.append("  %d file(s) already present and verified" % len(self.present))
        if self.needed:
            lines.append(
                "  %d file(s) to fetch, %s" % (len(self.needed), _limits.human(self.bytes_needed))
            )
            for remote, reason in self.needed:
                lines.append(
                    "    %-46s %-10s %s" % (remote.path, _limits.human(remote.size), reason)
                )
        else:
            lines.append("  nothing to do")
        return "\n".join(lines)


def _landing_path(dataset: Dataset, remote: Remote, dest: str) -> str:
    """Where *remote* lands under *dest*.

    Flattened to a basename. The published layout (reads under ``gastrooral/``, truth under
    ``per_bodysite/gastro_bam/``) carries nothing the harness needs.
    """
    return os.path.join(dest, "archives" if remote.is_archive else "gold", remote.name)


def plan(dataset: Dataset, dest: str, *, hash_it: bool = False) -> FetchPlan:
    """What is missing from *dest*. Never touches the network.

    *hash_it* is off by default to avoid hashing 95 GB. The size check catches truncation, and
    :func:`fetch` hashes everything it downloads. Pass it to re-verify an existing copy in full.

    A verified archive whose contents are not extracted still counts as needed. Without this, a
    fetch interrupted between download and extraction (three verified archives, one extracted
    sample) was reported complete, and a benchmark ran on one sample while reporting three.
    """
    result = FetchPlan(dataset=dataset, dest=dest)
    for remote in list(dataset.files) + list(dataset.gold_standard):
        if dataset.is_local:
            path = dataset.local_path_for(remote)
            reason = "missing" if not os.path.isfile(path) else OK
        else:
            path = _landing_path(dataset, remote, dest)
            reason = verify(path, remote, hash_it=hash_it)
            if not reason and remote.is_archive and not _is_extracted(dataset, remote, dest):
                reason = "downloaded but not extracted"
        if reason:
            result.needed.append((remote, reason))
        else:
            result.present.append(remote)
    return result


def _is_extracted(dataset: Dataset, remote: object, dest: str) -> bool:
    """Has *remote*'s sample been unpacked under ``dest/extracted``?

    Matched on the sample number, because the archive's directory is timestamp-prefixed
    (``2017.12.04_18.56.22_sample_17``) and the keep list extracts only a subset of members. The
    test is a directory ending in ``_sample_<n>`` that holds at least one file, which needs no
    access to the archive.
    """
    sample = getattr(remote, "sample", None)
    if sample is None:
        return True
    extracted = os.path.join(dest, "extracted")
    if not os.path.isdir(extracted):
        return False
    suffix = "_sample_%s" % sample
    for name in os.listdir(extracted):
        if not name.endswith(suffix):
            continue
        for _parent, _dirs, files in os.walk(os.path.join(extracted, name)):
            if files:
                return True
    return False


def disk_verdict(plan_: FetchPlan, budget: Optional[_limits.Budget] = None) -> Tuple[bool, str]:
    """Is there room under the ``[limits]`` disk ceiling for what *plan_* still needs?

    Compares the download total, which is an upper bound: each archive is deleted after
    extraction and the kept members are compressed FASTQ, so the extracted size cannot exceed it.
    """
    budget = budget or _limits.resolve()
    state = _limits.disk_state(plan_.dest, budget)
    if state is None:
        return True, "could not measure the filesystem holding %s; proceeding" % plan_.dest
    if state.over:
        return False, "already over the ceiling: %s" % state.describe()
    need = plan_.bytes_needed
    if need > state.headroom_bytes:
        return False, "%s needed, %s" % (_limits.human(need), state.describe())
    return True, "%s needed, %s" % (_limits.human(need), state.describe())


def partial_path(path: str) -> str:
    """Where bytes for *path* land while they are still unverified."""
    return path + PARTIAL_SUFFIX


def aria2c() -> Optional[str]:
    """Path to ``aria2c``, or None when no ranged downloader is available.

    aria2c is optional. It is looked up on PATH, then in the annotate env, where it is a
    transitive dependency of bakta (``envs/locks/annotate.linux-64.lock``) and may disappear with
    a future bakta pin. Callers must work when this returns None; :func:`download` then uses curl
    with the same checks.

    ``METAWRAP2_ARIA2C`` names the aria2c binary to use; set it empty to force curl, for a server
    that throttles or blocks clients opening several connections.
    """
    override = os.environ.get("METAWRAP2_ARIA2C")
    if override is not None:
        return override or None
    found = shutil.which("aria2c")
    if found:
        return found
    from ..config import MODULE_ENVS, env_prefix

    prefix = env_prefix(MODULE_ENVS["annotate"])
    if not prefix:
        return None
    candidate = os.path.join(prefix, "bin", "aria2c")
    return candidate if os.access(candidate, os.X_OK) else None


def _aria2c_argv(binary: str, url: str, partial: str) -> List[str]:
    """A multi-connection ranged download of *url* into *partial*.

    aria2c takes the destination as ``--dir`` and ``--out``. ``--auto-file-renaming=false`` stops
    it from writing ``<name>.1`` when the destination exists.

    ``--checksum`` and ``--check-integrity`` are not used, even when an md5 is known. All checks
    run in :func:`_verified`, so the curl path gets the same ones.

    Progress comes from aria2c's summary line. Do not compute it from the file size on disk:
    aria2c allocates the whole file before the first byte arrives, so the size reads as complete
    at once.
    """
    return [
        binary,
        "--max-connection-per-server=%d" % CONNECTIONS,
        "--split=%d" % CONNECTIONS,
        # Below this, splitting a file costs more in requests than it gains in parallelism.
        "--min-split-size=20M",
        "--max-tries=3",
        "--retry-wait=5",
        # Safe with aria2c: the control file records which ranges arrived.
        "--continue=true",
        "--auto-file-renaming=false",
        "--allow-overwrite=true",
        "--summary-interval=30",
        "--console-log-level=warn",
        "--dir",
        os.path.dirname(os.path.abspath(partial)),
        "--out",
        os.path.basename(partial),
        url,
    ]


def _curl_argv(url: str, partial: str) -> List[str]:
    """The fallback: one connection, and ``-C -`` to resume a partial file.

    Resuming saves hours on large files. A resume at the wrong offset can produce a file of the
    right size with a corrupt gzip stream, so :func:`_clear` deletes the partial whenever a check
    fails and the next attempt starts from zero.
    """
    return ["curl", "-fL", "--retry", "3", "--retry-delay", "5", "-C", "-", "-o", partial, url]


def _clear(path: str) -> None:
    """Remove the partial file and any downloader control file for *path*.

    Called before every retry, because resuming onto bytes that failed a check fails again. With
    nothing left behind, the next run and :func:`plan` treat the file as missing.
    """
    partial = partial_path(path)
    for victim in (partial, partial + CONTROL_SUFFIX, path + CONTROL_SUFFIX):
        try:
            os.remove(victim)
        except OSError:
            continue


#: Decompressor for each compression suffix, used to test that a download decompresses. For a
#: ``.tar.gz`` this tests the gzip stream only; :func:`extract` reads the members later.
_OPENERS: Tuple[Tuple[str, Callable[..., Any]], ...] = (
    (".gz", gzip.open),
    (".tgz", gzip.open),
    (".bz2", bz2.open),
    (".tbz2", bz2.open),
)


def _shown(path: str) -> str:
    """The name of *path* for messages, without the ``.part`` suffix.

    A failed attempt deletes the partial file, so messages name the requested file.
    """
    return os.path.basename(path).removesuffix(PARTIAL_SUFFIX)


def _opener(path: str) -> Optional[Callable[..., Any]]:
    """How to decompress *path*, or None for a format with no integrity check of its own."""
    # Verification runs before the rename, so strip `.part` to find the compression suffix.
    lower = path.lower().removesuffix(PARTIAL_SUFFIX)
    for suffix, opener in _OPENERS:
        if lower.endswith(suffix):
            return opener
    return None


def content_problem(path: str) -> str:
    """Why *path* does not decompress, or :data:`OK` (also for uncompressed formats).

    The in-process equivalent of ``gzip -t`` and ``bzip2 -t``: read the whole stream so the
    decompressor checks its trailer. It catches a corrupt file of the right size when no checksum
    is published.
    """
    opener = _opener(path)
    if opener is None:
        return OK
    try:
        with opener(path, "rb") as fh:
            while fh.read(_CHUNK):
                pass
    except (OSError, EOFError, zlib.error) as problem:
        return "%s does not decompress: %s" % (_shown(path), problem)
    return OK


def _verified(
    path: str, *, md5: Optional[str] = None, size: Optional[int] = None
) -> Tuple[Optional[Problem], str]:
    """``(what rejected it, what proved it)`` for a freshly downloaded *path*.

    Checks run cheapest first: size, then md5, then decompression. A file is accepted on size
    alone only when there is neither a checksum nor a decompressor. A matching md5 is enough on
    its own, which saves a second full read of a large file.
    """
    if not os.path.isfile(path):
        return Problem("transport", "%s was never written" % _shown(path)), ""
    actual = os.path.getsize(path)
    if size is not None and actual != size:
        return Problem("size", "%d bytes on disk, %d expected" % (actual, size)), ""
    if actual == 0:
        return Problem("size", "%s is empty" % _shown(path)), ""
    if md5:
        got = md5_of(path)
        if got != md5:
            return Problem("checksum", "md5 %s on disk, %s expected" % (got, md5)), ""
        return None, "md5 %s" % got
    problem = content_problem(path)
    if problem:
        return Problem("decompression", problem), ""
    if _opener(path) is not None:
        return None, "decompresses cleanly"
    # No checksum and no compression: only the size was checked, and the message says so.
    return None, "size only - no checksum published and nothing to decompress"


def download(
    url: str,
    path: str,
    *,
    md5: Optional[str] = None,
    size: Optional[int] = None,
    attempts: int = ATTEMPTS,
    delay: float = RETRY_DELAY,
    say: Callable[[str], None] = print,
    run: Optional[Callable[[List[str]], int]] = None,
) -> str:
    """Fetch *url* to *path*, verified. Returns :data:`OK`, or which check failed and what it saw.

    Uses aria2c with several ranged connections when available, otherwise curl. Bytes land in
    ``.part``, are verified there, and are renamed only after a check passes. On :data:`OK` the
    file matched its md5 or decompressed; on failure nothing is left on disk.

    *run* replaces ``subprocess.call``; ``install-db`` passes its own runner
    (:func:`metawrap2.commands.install_db._run`) for process-group tracking and output capture.
    It is looked up at call time so that tests which patch ``subprocess.call`` take effect.
    """
    run = run or subprocess.call
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    partial = partial_path(path)
    binary = aria2c()
    # A partial file left by aria2c can have gaps, which curl's `-C -` would not fill. Delete it
    # before using curl.
    if binary is None and os.path.exists(partial + CONTROL_SUFFIX):
        _clear(path)

    attempts = max(1, attempts)
    # Overwritten by every attempt; the initial value keeps the type a Problem.
    last = Problem("transport", "nothing was attempted")
    for attempt in range(1, attempts + 1):
        argv = _aria2c_argv(binary, url, partial) if binary else _curl_argv(url, partial)
        try:
            rc = run(argv)
        except FileNotFoundError:
            # Not retried. Reached when curl is missing or METAWRAP2_ARIA2C names a missing file.
            return "transport check failed: %s is not installed" % argv[0]
        if rc != 0:
            last = Problem("transport", "%s exited %d" % (os.path.basename(argv[0]), rc))
        else:
            problem, proof = _verified(partial, md5=md5, size=size)
            if problem is None:
                os.replace(partial, path)
                say("  verified (%s)" % proof)
                return OK
            last = problem
        say("  attempt %d of %d failed: %s" % (attempt, attempts, last))
        _clear(path)
        if attempt < attempts:
            time.sleep(delay * 2 ** (attempt - 1))
    return "%s check failed after %d attempt(s): %s" % (last.check, attempts, last.detail)


def _member_is_safe(name: str) -> bool:
    """Would extracting *name* stay inside the destination directory?

    Rejects absolute paths and ``..`` components. Python 3.12 added a tarfile filter for this; this
    check covers every supported Python version.
    """
    if name.startswith(("/", "\\")):
        return False
    parts = name.replace("\\", "/").split("/")
    return ".." not in parts


def _wanted(name: str, keep: Sequence[str]) -> bool:
    """Does *name* match any keep-glob? An empty *keep* means keep everything."""
    if not keep:
        return True
    return any(fnmatch.fnmatch(name, pattern) for pattern in keep)


def extract(archive: str, into: str, keep: Sequence[str] = ()) -> Tuple[int, int]:
    """Extract the members of *archive* matching *keep* into *into*.

    Returns ``(members, bytes)`` written. Symlinks and device nodes are skipped: a reads archive
    does not need them, and they can also point outside the destination.
    """
    os.makedirs(into, exist_ok=True)
    members_written = 0
    bytes_written = 0
    # Use the Python 3.12+ 'data' filter where available, in addition to the explicit checks
    # (which keep behaviour the same on older versions). It also strips mode bits and ownership.
    filtered = hasattr(tarfile, "data_filter")
    with tarfile.open(archive) as tf:
        for member in tf:
            if not (member.isfile() or member.isdir()):
                continue
            if not _member_is_safe(member.name):
                continue
            if member.isfile() and not _wanted(member.name, keep):
                continue
            if filtered:
                tf.extract(member, into, filter="data")
            else:
                tf.extract(member, into)
            if member.isfile():
                members_written += 1
                bytes_written += member.size
    return members_written, bytes_written


def fetch(
    dataset: Dataset,
    dest: str,
    *,
    budget: Optional[_limits.Budget] = None,
    keep_archives: bool = False,
    dry_run: bool = False,
    say: Callable[[str], None] = print,
) -> int:
    """Fetch everything *dataset* needs into *dest*. Returns an exit code.

    Each archive is downloaded, checksum-verified, extracted, and deleted unless *keep_archives*.
    Peak disk use is then about one archive plus the extracted reads.
    """
    if dataset.is_local:
        current = plan(dataset, dest)
        if current.complete:
            say("%s is a committed fixture; nothing to fetch." % dataset.name)
            return 0
        say(
            "%s is a committed fixture but %d of its files are missing from %s:\n%s"
            % (
                dataset.name,
                len(current.needed),
                dataset.local_root,
                "\n".join("  %s" % r.path for r, _ in current.needed),
            )
        )
        return 1

    current = plan(dataset, dest)
    say(current.describe())
    if current.complete:
        return 0

    ok, message = disk_verdict(current, budget)
    say("\nDisk: %s" % message)
    if not ok:
        say(
            "\nNot starting: this would exceed the [limits] disk ceiling.\n"
            "Fetch fewer samples, choose a --dest on another filesystem, or raise\n"
            "max_disk_fraction in [limits]."
        )
        return 1
    if dry_run:
        say("\n--dry-run: stopping here without downloading anything.")
        return 0

    extracted = os.path.join(dest, "extracted")
    failures: List[str] = []
    for remote, _reason in current.needed:
        path = _landing_path(dataset, remote, dest)
        say("\n%s (%s)" % (remote.path, _limits.human(remote.size)))
        if _reason == "downloaded but not extracted":
            # Downloaded and verified; extract without downloading again.
            say("  already downloaded and verified; extracting")
            members, written = extract(path, extracted, dataset.keep)
            say("  extracted %d file(s), %s" % (members, _limits.human(written)))
            if not keep_archives:
                os.remove(path)
                say("  removed the archive (pass --keep-archives to keep it)")
            continue
        # `download` verifies against the registry's md5 and size before renaming, and deletes a
        # file that fails, so there is no second check here.
        problem = download(dataset.url_for(remote), path, md5=remote.md5, size=remote.size, say=say)
        if problem:
            failures.append("%s: %s" % (remote.path, problem))
            say("  FAILED: %s" % problem)
            continue
        if remote.is_archive:
            members, written = extract(path, extracted, dataset.keep)
            say(
                "  extracted %d file(s), %s%s"
                % (
                    members,
                    _limits.human(written),
                    " (kept %s of the archive)" % ", ".join(dataset.keep) if dataset.keep else "",
                )
            )
            if not keep_archives:
                os.remove(path)
                say("  removed the archive (pass --keep-archives to keep it)")

    if failures:
        say("\n%d file(s) failed:\n%s" % (len(failures), "\n".join("  " + f for f in failures)))
        return 1
    say("\n%s is complete under %s" % (dataset.name, dest))
    return 0
