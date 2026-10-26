"""Tests for the maximum temp-directory length that still fits a unix socket path.

CheckM1 sets ``tempfile.tempdir`` and SemiBin2 exports ``TMPDIR``, and ``multiprocessing`` binds its
listener socket under that directory. If the path is too long, ``bind()`` fails in the child and the
parent reports a bare ``EOFError``.

The socket goes in a ``pymp-XXXXXXXX`` subdirectory as ``listener-XXXXXXXX``; an earlier limit
ignored the subdirectory and was 14 characters too high. The boundary is tested against the real
socket API.
"""

from __future__ import annotations

import multiprocessing as mp
import os
import sys
import tempfile

import pytest

from metawrap2.utils import MAX_TMPDIR_LEN, SOCKET_PATH_LIMIT, short_tmpdir


def _bind_under(directory: str) -> bool:
    """Whether a multiprocessing Manager can bind a socket with *directory* as the temp dir."""
    previous = tempfile.tempdir
    tempfile.tempdir = directory
    try:
        manager = mp.Manager()
        manager.dict()
        manager.shutdown()
        return True
    except (OSError, EOFError):
        return False
    finally:
        tempfile.tempdir = previous


def _dir_of_length(root, length: int) -> str:
    """A directory whose absolute path is exactly *length* characters."""
    base = str(root)
    pad = length - len(base) - 1
    if pad < 1:
        pytest.skip("tmp_path %r is already longer than %d characters" % (base, length))
    path = os.path.join(base, "x" * pad)
    os.makedirs(path, exist_ok=True)
    assert len(path) == length
    return path


def test_one_character_over_the_limit_cannot(tmp_path):
    """A temp dir of MAX_TMPDIR_LEN binds and one character longer does not."""
    assert MAX_TMPDIR_LEN == SOCKET_PATH_LIMIT - (len("/pymp-") + 8) - (len("/listener-") + 8)
    assert _bind_under(_dir_of_length(tmp_path, MAX_TMPDIR_LEN)) is True
    if sys.version_info < (3, 14):
        # Python 3.14's multiprocessing falls back to a short temp directory itself. The limit
        # still applies to CheckM1 and SemiBin2, which run on older Pythons in their own envs.
        assert _bind_under(_dir_of_length(tmp_path, MAX_TMPDIR_LEN + 1)) is False

    # A path over the limit is replaced by a short symlink, which bind() resolves.
    real = _dir_of_length(tmp_path, MAX_TMPDIR_LEN + 20)
    link = short_tmpdir(real)
    assert len(link) <= MAX_TMPDIR_LEN
    assert os.path.realpath(link) == os.path.realpath(real)
    assert _bind_under(link) is True
