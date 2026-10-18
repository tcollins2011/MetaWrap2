"""Small utility functions shared across MetaWrap2 modules.

Kept deliberately dependency-light so any module or helper script can import from here
without pulling in heavy libraries.
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


#: Linux gives a unix-domain socket address 108 bytes of ``sun_path`` including its terminator, so
#: 107 usable characters. That number decides where a tool is allowed to put its temp files, for a
#: reason that is not obvious from either end: several of the tools MetaWrap2 drives take a temp
#: directory and then set it process-wide (CheckM1 assigns ``tempfile.tempdir``; SemiBin2 exports
#: ``TMPDIR``), and Python's ``multiprocessing`` binds its listener socket inside whatever that
#: resolves to. The name is ``/listener-`` plus eight random characters, so a long enough temp
#: directory leaves no room, ``bind()`` fails in the child, and the parent reports an ``EOFError``
#: from a ``recv()`` that says nothing about paths.
#:
#: Measured rather than reasoned about: a 90-character tmpdir under a benchmark work directory
#: failed exactly this way on every variant of a comparison - taking its completeness metrics - and
#: the same bins passed under /tmp. 90 is not a deep path by any standard, which is why this is
#: handled rather than documented as a limitation.
SOCKET_PATH_LIMIT = 107
#: What ``multiprocessing`` adds below the temp directory before it binds anything. Two components,
#: not one, and missing the first is what made the first version of this fix too generous by 14
#: characters: ``get_temp_dir()`` first creates a *subdirectory* named ``pymp-`` plus 8 random
#: characters, and only then does ``arbitrary_address()`` put ``listener-`` plus 8 more inside it.
#:
#: The boundary this gives - 75 usable characters - was then checked by binding real sockets under
#: directories of length from 60 to 90: 75 succeeds, 76 raises ``AF_UNIX path too long`` in the
#: child and surfaces in the parent as a bare ``EOFError``. A 77-character tmpdir is what exposed
#: the gap, having been judged short enough.
_MP_TEMP_DIR = len("/pymp-") + 8
_LISTENER_NAME = len("/listener-") + 8
MAX_TMPDIR_LEN = SOCKET_PATH_LIMIT - _MP_TEMP_DIR - _LISTENER_NAME


def short_tmpdir(real: str) -> str:
    """Create *real* and return a path to it short enough to hold a unix socket.

    The directory itself stays where the caller asked for it, which is the point of passing one at
    all: these temp directories hold GBs, and they belong on the filesystem chosen for the run's
    output rather than on whatever /tmp happens to be. When that path is too long (see
    :data:`MAX_TMPDIR_LEN`) the returned path is instead a short symlink to it under the system temp
    directory - ``bind()`` resolves the symlink, so the socket name stays short while every byte
    still lands where it was going to. The symlink is named from a hash of the real path, so a
    ``--resume`` finds the same one rather than leaking a fresh link per attempt.
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
