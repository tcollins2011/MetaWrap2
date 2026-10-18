"""Getting benchmark data onto this machine: planned, checksum-verified, and inside the budget.

Three things make this more than a `curl` loop, and each of them is here because the alternative
fails in a way that costs hours:

* **It plans before it fetches.** A CAMI2 body site is ~95 GB. :func:`plan` says exactly what is
  missing, what that weighs, and whether the ``[limits]`` disk ceiling has room, and it does that
  without opening a socket - so "will this fit" is answerable before anything is committed to.
* **It verifies by size first, then by checksum.** A truncated 10 GB download looks exactly like a
  good one until something tries to read the end of it. Comparing the length catches the common
  case in milliseconds; the md5 catches the rest.
* **It extracts selectively and safely.** ~95% of each CAMI2 archive is per-OTU BAM alignments no
  binning benchmark reads. Skipping them cannot shrink the download - a ``.tar.gz`` is one gzip
  stream, so there is nothing to seek past - but it keeps what lands to the reads and contigs.
  Members are also checked for absolute and escaping paths, because these are archives from the
  internet and ``extractall`` will happily write wherever a member says to.

Nothing here downloads on import, on ``list``, or on ``plan``. A fetch happens when a user names a
dataset and confirms its size, and not otherwise.
"""

from __future__ import annotations

import fnmatch
import hashlib
import os
import subprocess
import tarfile
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Sequence, Tuple

from .. import limits as _limits
from .registry import Dataset, Remote

#: Read size for hashing. Large enough that a 10 GB file is not 10 million syscalls, small enough
#: not to matter to memory.
_CHUNK = 8 * 1024 * 1024

#: Reported by :func:`verify` when a file is intact. Kept as a constant so callers compare against
#: a name rather than an empty string.
OK = ""


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

    Size is checked first because it is free and catches the interrupted download, which is the
    failure that actually happens. *hash_it* exists for the caller that has already paid for the
    hash once and only wants the cheap check.
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
    #: Files that are absent, paired with why - "missing", or the reason an existing copy was
    #: rejected. Keeping the reason means the report can distinguish "never downloaded" from
    #: "downloaded and corrupt", which call for different responses.
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

    Flattened to a basename on purpose: the published layout puts reads under ``gastrooral/`` and
    truth under ``per_bodysite/gastro_bam/``, which is the publisher's filing system, not
    information the harness needs to carry around.
    """
    return os.path.join(dest, "archives" if remote.is_archive else "gold", remote.name)


def plan(dataset: Dataset, dest: str, *, hash_it: bool = False) -> FetchPlan:
    """What is missing from *dest*. Never touches the network.

    *hash_it* is off by default because hashing 95 GB to answer "what do I still need" is a poor
    trade; the size check catches truncation, and :func:`fetch` hashes everything it downloads.
    Pass it to re-verify an existing copy in full.

    An archive that is present *and* verified still counts as needed when its contents are not on
    disk. That is not hypothetical: a fetch interrupted after downloading and before extracting
    left three verified archives and one extracted sample, and the next run said the dataset was
    "complete". A benchmark then runs on one sample while reporting three - which is the kind of
    wrong number this project exists to avoid, and it has no error anywhere in it.
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

    Matched on the sample number rather than on a member list, because the archive's own directory
    is timestamp-prefixed (``2017.12.04_18.56.22_sample_17``) and the keep list means the extracted
    tree is deliberately a subset of what the archive holds. A directory whose name ends in
    ``_sample_<n>`` and holds at least one file is the honest test available without opening the
    archive, which is the whole point of not opening it.
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

    The figure compared is the download total, which is a genuine upper bound rather than an
    estimate: each archive is deleted once extracted, and what is kept out of it is
    already-compressed FASTQ, so the extracted footprint cannot exceed the bytes that arrived.
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


def _download(url: str, path: str) -> str:
    """Fetch *url* to *path*, resuming a partial file. Returns "" on success, else the reason.

    Downloads into ``.part`` and renames, so an interrupted fetch can never leave a short file
    that passes for complete - the size check would catch it anyway, but a half file under the
    real name invites something else to read it first.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    partial = path + ".part"
    argv = ["curl", "-fL", "--retry", "3", "--retry-delay", "5", "-C", "-", "-o", partial, url]
    try:
        rc = subprocess.call(argv)
    except FileNotFoundError:
        return "curl is not installed"
    if rc != 0:
        return "curl exited %d" % rc
    os.replace(partial, path)
    return OK


def _member_is_safe(name: str) -> bool:
    """Would extracting *name* stay inside the destination directory?

    ``tarfile.extractall`` writes wherever a member says to, including ``/etc/...`` and
    ``../../..``. Python 3.12 gained a filter for this; this check works on every version the
    project supports and is explicit about what it rejects.
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

    Returns ``(members, bytes)`` actually written. Symlinks and device nodes are skipped outright:
    nothing in a reads archive needs them, and they are the other half of how an archive escapes
    its destination.
    """
    os.makedirs(into, exist_ok=True)
    members_written = 0
    bytes_written = 0
    # Python 3.12 added a member filter for exactly this class of problem. Used where available on
    # top of the checks above rather than instead of them: the explicit checks are what make the
    # behaviour the same on every version the project supports, and 'data' also strips the mode
    # bits and ownership that a tarball has no business setting.
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

    Each archive is downloaded, checksum-verified, extracted, and then deleted unless
    *keep_archives* - deleting as it goes is what keeps peak disk at roughly one archive plus the
    extracted reads instead of the whole 95 GB twice.
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
            "\nRefusing to start: this would cross the [limits] disk ceiling.\n"
            "Fetch fewer samples, choose a --dest on another filesystem, or raise\n"
            "max_disk_fraction if you know the rest of the machine can spare it."
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
            # Already here and already verified; the only thing missing is the unpacking. Saying so
            # matters because the alternative - re-downloading 10 GB that is on disk - is what the
            # old "present and verified" shortcut was trying to avoid in the first place.
            say("  already downloaded and verified; extracting")
            members, written = extract(path, extracted, dataset.keep)
            say("  extracted %d file(s), %s" % (members, _limits.human(written)))
            if not keep_archives:
                os.remove(path)
                say("  removed the archive (pass --keep-archives to keep it)")
            continue
        problem = _download(dataset.url_for(remote), path)
        if problem:
            failures.append("%s: %s" % (remote.path, problem))
            say("  FAILED: %s" % problem)
            continue
        problem = verify(path, remote)
        if problem:
            # Leave the bad file in place rather than deleting it: a checksum failure on a file
            # this size is worth looking at, and a silent delete-and-retry loop can spin forever
            # on a host serving corrupt bytes.
            failures.append("%s: %s" % (remote.path, problem))
            say("  FAILED verification: %s" % problem)
            continue
        say("  verified (md5 %s)" % remote.md5)
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
