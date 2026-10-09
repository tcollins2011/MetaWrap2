"""The temp-directory length a unix socket still fits under, and the arithmetic behind it.

Several tools MetaWrap2 drives take a temp directory and then set it process-wide - CheckM1 assigns
``tempfile.tempdir``, SemiBin2 exports ``TMPDIR`` - and Python's ``multiprocessing`` binds its listener
socket inside whatever that resolves to. Too long, and ``bind()`` fails in the child while the parent
reports a bare ``EOFError`` that mentions nothing about paths.

This was fixed once and the fix was too generous by 14 characters, because ``multiprocessing`` creates
a ``pymp-XXXXXXXX`` *subdirectory* before putting ``listener-XXXXXXXX`` inside it. A 77-character temp
directory was therefore judged short enough and still failed, taking a whole comparison's metrics with
it. The boundary is pinned here against the real socket API rather than against the arithmetic.
"""

from __future__ import annotations

import multiprocessing as mp
import os
import tempfile

import pytest

from metawrap2.utils import MAX_TMPDIR_LEN, SOCKET_PATH_LIMIT, short_tmpdir


def test_the_budget_accounts_for_both_components_multiprocessing_adds():
    """107 usable bytes, minus /pymp- + 8, minus /listener- + 8."""
    assert SOCKET_PATH_LIMIT == 107
    assert MAX_TMPDIR_LEN == 107 - (len("/pymp-") + 8) - (len("/listener-") + 8)
    assert MAX_TMPDIR_LEN == 75


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


def test_a_directory_at_the_limit_can_still_hold_a_socket(tmp_path):
    assert _bind_under(_dir_of_length(tmp_path, MAX_TMPDIR_LEN)) is True


def test_one_character_over_the_limit_cannot(tmp_path):
    """The assertion that would have caught the off-by-14: 76 has to fail, not merely 90."""
    assert _bind_under(_dir_of_length(tmp_path, MAX_TMPDIR_LEN + 1)) is False


def test_short_tmpdir_leaves_a_path_within_the_budget_alone(tmp_path):
    inside = _dir_of_length(tmp_path, MAX_TMPDIR_LEN - 4)
    assert short_tmpdir(inside) == inside


def test_short_tmpdir_symlinks_a_path_over_the_budget_and_the_link_works(tmp_path):
    real = _dir_of_length(tmp_path, MAX_TMPDIR_LEN + 20)
    link = short_tmpdir(real)
    assert link != real
    assert len(link) <= MAX_TMPDIR_LEN
    assert os.path.realpath(link) == os.path.realpath(real)
    # The point of the symlink is that bind() resolves it, so a socket under the link must work.
    assert _bind_under(link) is True


def test_the_same_real_path_reuses_one_link(tmp_path):
    """A --resume must find the same link rather than leaking a fresh one per attempt."""
    real = _dir_of_length(tmp_path, MAX_TMPDIR_LEN + 20)
    assert short_tmpdir(real) == short_tmpdir(real)
