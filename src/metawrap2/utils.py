"""Small utility functions shared across MetaWrap2 modules.

Has no heavy dependencies, so any module or helper script can import it.
"""

from __future__ import annotations

import hashlib
import os
import tempfile

from .constants import CONTAMINATION_WEIGHT, DEREPLICATION_TIEBREAK


def quality_score(completeness: float, contamination: float, tiebreak: float = 0.0) -> float:
    """Bin quality score: completeness - 5*contamination (+ tiny tiebreaker).

    The tiebreaker (typically genome size) only separates otherwise-equal bins and does
    not meaningfully shift ranking.
    """
    return completeness - CONTAMINATION_WEIGHT * contamination + DEREPLICATION_TIEBREAK * tiebreak


def bin_name_from_filename(filename: str) -> str:
    """Return a bin's name by stripping the directory and final extension.

    e.g. ``/path/bin.3.fa`` -> ``bin.3`` and ``some/dir/maxbin.007.fasta`` -> ``maxbin.007``.
    """
    base = os.path.basename(filename)
    stem, _ext = os.path.splitext(base)
    return stem


def strip_trailing_slash(path: str) -> str:
    """Remove a single trailing slash, matching the legacy scripts' path handling."""
    return path.removesuffix("/")


#: Linux gives a unix-domain socket address 108 bytes of ``sun_path`` including the terminator, so
#: 107 usable characters. Several tools MetaWrap2 runs set their temp directory process-wide
#: (CheckM1 assigns ``tempfile.tempdir``; SemiBin2 exports ``TMPDIR``), and Python's
#: ``multiprocessing`` binds its listener socket inside it. With a long temp directory ``bind()``
#: fails in the child and the parent reports an unexplained ``EOFError`` from ``recv()``. A
#: 90-character tmpdir under a benchmark work directory failed this way; the same bins passed under
#: /tmp.
SOCKET_PATH_LIMIT = 107
#: What ``multiprocessing`` adds below the temp directory before binding: ``get_temp_dir()``
#: creates a subdirectory ``pymp-`` plus 8 random characters, then ``arbitrary_address()`` adds
#: ``listener-`` plus 8 more inside it.
#:
#: This gives 75 usable characters, checked by binding real sockets under directories of 60 to 90
#: characters: 75 succeeds, 76 raises ``AF_UNIX path too long`` in the child (a bare ``EOFError``
#: in the parent).
_MP_TEMP_DIR = len("/pymp-") + 8
_LISTENER_NAME = len("/listener-") + 8
MAX_TMPDIR_LEN = SOCKET_PATH_LIMIT - _MP_TEMP_DIR - _LISTENER_NAME


def short_tmpdir(real: str) -> str:
    """Create *real* and return a path to it short enough to hold a unix socket.

    The directory stays where the caller asked (these temp directories hold GBs and belong on the
    run's filesystem, not /tmp). When the path is too long (see :data:`MAX_TMPDIR_LEN`), a short
    symlink to it under the system temp directory is returned instead, so socket names stay short
    while data is still written to *real*. The symlink name is a hash of the real path, so
    ``--resume`` reuses it.
    """
    os.makedirs(real, exist_ok=True)  # harmless under --dry-run; the tools need it to exist
    if len(real) <= MAX_TMPDIR_LEN:
        return real
    digest = hashlib.sha1(os.path.abspath(real).encode()).hexdigest()[:8]
    link = os.path.join(tempfile.gettempdir(), "mw2-tmp-" + digest)
    if os.path.islink(link) and os.path.realpath(link) != os.path.realpath(real):
        os.remove(link)  # a stale link from a deleted run that hashed the same
    if not os.path.exists(link):
        os.symlink(real, link)
    return link
