"""Tests for ``benchmark.fetch.download``: the checks a download must pass, and cleanup on failure.

A resumed ``curl -C -`` download can have exactly the expected size and a corrupt gzip stream.
Several tests build such a file. aria2c is used for ranged downloads when available and curl
otherwise; both paths run the same checks.

No test opens a socket or starts a process. The downloader is a function that writes the bytes a
test chooses to wherever the argv points.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import tarfile
from typing import Callable, List, Optional, Tuple

import pytest

from metawrap2.benchmark import fetch
from metawrap2.benchmark.registry import Dataset, Remote

#: Enough bytes that the gzip stream has a body rather than only a header and a trailer.
PAYLOAD = b"ACGT" * 4096


@pytest.fixture(autouse=True)
def waits(monkeypatch) -> List[float]:
    """Record the retry waits instead of sleeping, and return them.

    ``fetch.fetch`` has no way to pass a shorter delay, so the sleep is patched.
    """
    recorded: List[float] = []
    monkeypatch.setattr(fetch.time, "sleep", recorded.append)
    return recorded


def _corrupt(data: bytes) -> bytes:
    """*data* with its last byte flipped: same length, and the gzip trailer check fails."""
    return data[:-1] + bytes([data[-1] ^ 0xFF])


def _destination(argv: List[str]) -> str:
    """The output path in *argv*: ``-o`` for curl, ``--dir`` and ``--out`` for aria2c."""
    if "-o" in argv and "--out" not in argv:
        return argv[argv.index("-o") + 1]
    return os.path.join(argv[argv.index("--dir") + 1], argv[argv.index("--out") + 1])


def _downloader(
    *arrivals: Optional[bytes], rc: int = 0, control: bool = False
) -> Tuple[List[List[str]], List[List[str]], Callable[[List[str]], int]]:
    """A stand-in downloader that writes *arrivals* in order, one per attempt.

    Returns the recorded argvs, the directory listing each attempt found, and the callable.
    The last arrival repeats for later attempts. ``None`` writes nothing. *control* also writes an
    ``.aria2`` control file, as aria2c does.
    """
    calls: List[List[str]] = []
    found: List[List[str]] = []
    queue = list(arrivals)

    def run(argv: List[str]) -> int:
        calls.append(list(argv))
        target = _destination(argv)
        directory = os.path.dirname(target)
        found.append(sorted(os.listdir(directory)) if os.path.isdir(directory) else [])
        payload = queue.pop(0) if len(queue) > 1 else (queue[0] if queue else None)
        if payload is not None:
            with open(target, "wb") as fh:
                fh.write(payload)
            if control:
                with open(target + fetch.CONTROL_SUFFIX, "wb") as fh:
                    fh.write(b"control")
        return rc

    return calls, found, run


def _curl_only(monkeypatch) -> None:
    """Make aria2c unavailable, so curl is used."""
    monkeypatch.setattr(fetch, "aria2c", lambda: None)


def _with_aria2c(monkeypatch, path: str = "/opt/conda/bin/aria2c") -> None:
    monkeypatch.setattr(fetch, "aria2c", lambda: path)


def test_a_size_correct_but_corrupt_gzip_is_rejected(tmp_path, monkeypatch):
    """A gzip file of the right size with a corrupt stream is rejected and not left in place."""
    _curl_only(monkeypatch)
    good = gzip.compress(PAYLOAD)
    _calls, _found, run = _downloader(_corrupt(good))
    dest = tmp_path / "reads.fastq.gz"

    problem = fetch.download(
        "https://example.invalid/reads.fastq.gz",
        str(dest),
        size=len(good),
        attempts=1,
        run=run,
        say=lambda line: None,
    )

    assert problem.startswith("decompression check failed"), problem
    assert "does not decompress" in problem
    assert not dest.exists(), "a file that failed a check must not be left under the real name"


def test_a_preallocated_full_length_file_is_not_a_finished_download(tmp_path, monkeypatch):
    """A full-length file of zeros, as aria2c preallocates, fails the decompression check."""
    _with_aria2c(monkeypatch)
    good = gzip.compress(PAYLOAD)
    _calls, _found, run = _downloader(b"\x00" * len(good))
    dest = tmp_path / "reads.fastq.gz"

    problem = fetch.download(
        "https://example.invalid/reads.fastq.gz",
        str(dest),
        size=len(good),
        attempts=1,
        run=run,
        say=lambda line: None,
    )

    assert problem.startswith("decompression check failed"), problem
    assert not dest.exists()


def test_a_checksum_mismatch_is_rejected_and_names_both_digests(tmp_path, monkeypatch):
    """An md5 mismatch fails the checksum check, and the message shows both digests."""
    _curl_only(monkeypatch)
    arrived = gzip.compress(PAYLOAD)
    _calls, _found, run = _downloader(arrived)
    wrong = hashlib.md5(b"something else").hexdigest()
    dest = tmp_path / "reads.fastq.gz"

    problem = fetch.download(
        "https://example.invalid/reads.fastq.gz",
        str(dest),
        md5=wrong,
        attempts=2,
        run=run,
        say=lambda line: None,
    )

    assert problem == "checksum check failed after 2 attempt(s): md5 %s on disk, %s expected" % (
        hashlib.md5(arrived).hexdigest(),
        wrong,
    )
    assert not dest.exists()


def test_a_truncated_download_is_rejected_by_size(tmp_path, monkeypatch):
    """A short file fails the size check, which runs first."""
    _curl_only(monkeypatch)
    good = gzip.compress(PAYLOAD)
    _calls, _found, run = _downloader(good[:20])
    dest = tmp_path / "reads.fastq.gz"

    problem = fetch.download(
        "https://example.invalid/reads.fastq.gz",
        str(dest),
        size=len(good),
        attempts=1,
        run=run,
        say=lambda line: None,
    )

    assert problem == "size check failed after 1 attempt(s): 20 bytes on disk, %d expected" % len(
        good
    )


def test_a_transport_failure_names_the_transport_check(tmp_path, monkeypatch):
    """A non-zero downloader exit is reported as a transport check failure."""
    _curl_only(monkeypatch)
    _calls, _found, run = _downloader(None, rc=7)
    dest = tmp_path / "reads.fastq.gz"

    problem = fetch.download(
        "https://example.invalid/reads.fastq.gz",
        str(dest),
        attempts=3,
        run=run,
        say=lambda line: None,
    )

    assert problem == "transport check failed after 3 attempt(s): curl exited 7"


def test_a_missing_downloader_is_reported_once_and_not_retried(tmp_path, monkeypatch):
    """A missing downloader binary is reported once, without retries."""
    _curl_only(monkeypatch)
    attempts: List[int] = []

    def run(argv: List[str]) -> int:
        attempts.append(1)
        raise FileNotFoundError(argv[0])

    problem = fetch.download(
        "https://example.invalid/reads.fastq.gz",
        str(tmp_path / "reads.fastq.gz"),
        attempts=3,
        run=run,
        say=lambda line: None,
    )

    assert problem == "transport check failed: curl is not installed"
    assert len(attempts) == 1


def test_the_partial_and_control_file_are_gone_before_a_retry(tmp_path, monkeypatch):
    """Each retry starts with no partial file and no control file, and none is left at the end."""
    _with_aria2c(monkeypatch)
    good = gzip.compress(PAYLOAD)
    bad = _corrupt(good)
    calls, found, run = _downloader(bad, bad, bad, control=True)
    dest = tmp_path / "reads.fastq.gz"

    problem = fetch.download(
        "https://example.invalid/reads.fastq.gz",
        str(dest),
        size=len(good),
        attempts=3,
        run=run,
        say=lambda line: None,
    )

    assert problem.startswith("decompression check failed after 3 attempt(s)")
    assert len(calls) == 3
    assert found == [[], [], []], "every attempt has to start with an empty directory"
    assert (
        sorted(os.listdir(tmp_path)) == []
    ), "a failed fetch leaves no partial and no control file"


def test_a_retry_after_a_bad_attempt_can_succeed(tmp_path, monkeypatch):
    """A corrupt first attempt followed by a good one leaves only the final file."""
    _curl_only(monkeypatch)
    good = gzip.compress(PAYLOAD)
    calls, _found, run = _downloader(_corrupt(good), good)
    dest = tmp_path / "reads.fastq.gz"

    problem = fetch.download(
        "https://example.invalid/reads.fastq.gz",
        str(dest),
        md5=hashlib.md5(good).hexdigest(),
        size=len(good),
        attempts=3,
        run=run,
        say=lambda line: None,
    )

    assert problem == fetch.OK
    assert len(calls) == 2
    assert dest.read_bytes() == good
    assert sorted(os.listdir(tmp_path)) == ["reads.fastq.gz"]


def test_the_wait_between_attempts_doubles_and_none_follows_the_last_one(
    tmp_path, monkeypatch, waits
):
    """The wait doubles after each failed attempt, with no wait after the last one."""
    _curl_only(monkeypatch)
    _calls, _found, run = _downloader(None, rc=7)

    fetch.download(
        "https://example.invalid/reads.fastq.gz",
        str(tmp_path / "reads.fastq.gz"),
        attempts=3,
        run=run,
        say=lambda line: None,
    )

    assert waits == [fetch.RETRY_DELAY, fetch.RETRY_DELAY * 2]


def test_aria2c_is_used_when_present_and_curl_when_it_is_not(tmp_path, monkeypatch):
    """aria2c with ranged connections is used when found; otherwise curl with ``-C -``."""
    good = gzip.compress(PAYLOAD)
    digest = hashlib.md5(good).hexdigest()

    _with_aria2c(monkeypatch, "/opt/conda/bin/aria2c")
    ranged, _found, run = _downloader(good)
    assert (
        fetch.download(
            "https://example.invalid/reads.fastq.gz",
            str(tmp_path / "a" / "reads.fastq.gz"),
            md5=digest,
            run=run,
            say=lambda line: None,
        )
        == fetch.OK
    )
    assert ranged[0][0] == "/opt/conda/bin/aria2c"
    assert "--split=%d" % fetch.CONNECTIONS in ranged[0]
    assert "--max-connection-per-server=%d" % fetch.CONNECTIONS in ranged[0]

    _curl_only(monkeypatch)
    single, _found, run = _downloader(good)
    assert (
        fetch.download(
            "https://example.invalid/reads.fastq.gz",
            str(tmp_path / "b" / "reads.fastq.gz"),
            md5=digest,
            run=run,
            say=lambda line: None,
        )
        == fetch.OK
    )
    assert single[0][0] == "curl"
    assert "-C" in single[0], "curl resumes partial files"


def test_aria2c_is_asked_for_its_own_progress_rather_than_measured_by_file_size(monkeypatch):
    """aria2c is asked for periodic summaries and is not silenced, and gets ``--dir``/``--out``."""
    argv = fetch._aria2c_argv(
        "/opt/conda/bin/aria2c", "https://example.invalid/x.gz", "/d/x.gz.part"
    )

    assert any(a.startswith("--summary-interval=") for a in argv)
    assert "--dir" in argv and "/d" == argv[argv.index("--dir") + 1]
    assert "--out" in argv and "x.gz.part" == argv[argv.index("--out") + 1]
    assert not any(a.startswith("--quiet") for a in argv), "aria2c output must not be silenced"


def test_a_control_file_from_a_ranged_download_is_not_resumed_onto_by_curl(tmp_path, monkeypatch):
    """A partial file with an aria2c control file is deleted before curl runs.

    Such a file can have gaps, and curl's ``-C -`` would append after them.
    """
    _curl_only(monkeypatch)
    good = gzip.compress(PAYLOAD)
    dest = tmp_path / "reads.fastq.gz"
    partial = tmp_path / ("reads.fastq.gz" + fetch.PARTIAL_SUFFIX)
    partial.write_bytes(b"\x00" * len(good))
    control = tmp_path / ("reads.fastq.gz" + fetch.PARTIAL_SUFFIX + fetch.CONTROL_SUFFIX)
    control.write_bytes(b"control")

    _calls, found, run = _downloader(good)
    assert (
        fetch.download(
            "https://example.invalid/reads.fastq.gz",
            str(dest),
            md5=hashlib.md5(good).hexdigest(),
            run=run,
            say=lambda line: None,
        )
        == fetch.OK
    )
    assert found == [[]], "the ranged partial and its control file go before curl is started"


def test_a_file_with_no_checksum_and_no_compression_says_it_was_checked_by_size_only(
    tmp_path, monkeypatch
):
    """A file with no md5 and no compression is reported as checked by size only."""
    _curl_only(monkeypatch)
    said: List[str] = []
    _calls, _found, run = _downloader(b"ko00001 hierarchy\n")

    assert (
        fetch.download(
            "https://example.invalid/ko00001.keg",
            str(tmp_path / "ko00001.keg"),
            run=run,
            say=said.append,
        )
        == fetch.OK
    )
    assert said == ["  verified (size only - no checksum published and nothing to decompress)"]


def test_a_checksummed_download_is_not_also_decompressed(tmp_path, monkeypatch):
    """A matching md5 is accepted without a decompression check.

    The bytes match their md5 but are not valid gzip, so a decompression check would fail.
    """
    _curl_only(monkeypatch)
    arrived = b"\x1f\x8b" + b"\x00" * 200
    _calls, _found, run = _downloader(arrived)
    said: List[str] = []
    digest = hashlib.md5(arrived).hexdigest()

    assert (
        fetch.download(
            "https://example.invalid/reads.fastq.gz",
            str(tmp_path / "reads.fastq.gz"),
            md5=digest,
            run=run,
            say=said.append,
        )
        == fetch.OK
    )
    assert said == ["  verified (md5 %s)" % digest]


def test_the_content_check_covers_every_compressed_format_the_fetcher_lands(tmp_path):
    """``content_problem`` passes good gzip, bzip2 and plain files, and flags corrupt ones."""
    import bz2

    cases = {
        "reads.fastq.gz": gzip.compress(PAYLOAD),
        "archive.tgz": gzip.compress(PAYLOAD),
        "card.tar.bz2": bz2.compress(PAYLOAD),
        "db.dmnd.bz2": bz2.compress(PAYLOAD),
    }
    for name, data in cases.items():
        good = tmp_path / name
        good.write_bytes(data)
        assert fetch.content_problem(str(good)) == fetch.OK, name
        broken = tmp_path / ("broken_" + name)
        broken.write_bytes(_corrupt(data))
        assert fetch.content_problem(str(broken)), "%s: corruption not caught" % name
        assert broken.stat().st_size == good.stat().st_size, "corrupt file has the same size"

    plain = tmp_path / "card.json"
    plain.write_bytes(b"{}")
    assert fetch.content_problem(str(plain)) == fetch.OK


def _dataset(tmp_path, data: bytes, md5: str) -> Dataset:
    """A one-file remote dataset, as `benchmarks/datasets.toml` would declare it."""
    return Dataset(
        name="probe",
        title="a one-file stand-in for a CAMI2 body site",
        base_url="https://example.invalid/probe",
        samples=(1,),
        keep=("*/reads/*",),
        files=(Remote(path="sample_1.tar.gz", md5=md5, size=len(data), sample=1),),
    )


def _tarball(tmp_path) -> bytes:
    """A .tar.gz with one member under the keep glob, so `extract` has something to do."""
    reads = tmp_path / "staging" / "2017.12.04_18.56.22_sample_1" / "reads"
    reads.mkdir(parents=True)
    (reads / "anonymous_reads.fq.gz").write_bytes(gzip.compress(b"@r\nACGT\n+\nIIII\n"))
    archive = tmp_path / "staging.tar.gz"
    with tarfile.open(archive, "w:gz") as tf:
        tf.add(reads.parent, arcname="2017.12.04_18.56.22_sample_1")
    return archive.read_bytes()


def test_a_dataset_fetch_verifies_against_the_registry_md5_before_extracting(tmp_path, monkeypatch):
    """``fetch`` rejects an archive that fails the registry md5 and extracts one that passes.

    The size matches in both cases, so the md5 decides.
    """
    _curl_only(monkeypatch)
    data = _tarball(tmp_path)
    dataset = _dataset(tmp_path, data, hashlib.md5(data).hexdigest())
    dest = tmp_path / "data"

    _calls, _found, bad = _downloader(_corrupt(data))
    monkeypatch.setattr(fetch.subprocess, "call", bad)
    said: List[str] = []
    assert fetch.fetch(dataset, str(dest), dry_run=False, say=said.append) == 1
    report = "\n".join(said)
    assert "checksum check failed" in report, report
    assert not (dest / "archives" / "sample_1.tar.gz").exists()

    _calls, _found, good = _downloader(data)
    monkeypatch.setattr(fetch.subprocess, "call", good)
    said = []
    assert fetch.fetch(dataset, str(dest), dry_run=False, say=said.append) == 0
    report = "\n".join(said)
    assert "verified (md5 %s)" % hashlib.md5(data).hexdigest() in report, report
    assert (dest / "extracted" / "2017.12.04_18.56.22_sample_1").is_dir()
    assert not (dest / "archives" / "sample_1.tar.gz").exists(), "the archive is deleted once read"


def test_aria2c_is_found_without_being_a_dependency(monkeypatch, tmp_path):
    """aria2c is found on PATH or in the annotate env, else None; ``METAWRAP2_ARIA2C`` overrides."""
    monkeypatch.delenv("METAWRAP2_ARIA2C", raising=False)
    monkeypatch.setattr(fetch.shutil, "which", lambda name: "/usr/bin/aria2c")
    assert fetch.aria2c() == "/usr/bin/aria2c"

    monkeypatch.setattr(fetch.shutil, "which", lambda name: None)
    monkeypatch.setattr("metawrap2.config.env_prefix", lambda env: None)
    assert fetch.aria2c() is None

    binary = tmp_path / "envs" / "metawrap2-annotate" / "bin" / "aria2c"
    binary.parent.mkdir(parents=True)
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)
    monkeypatch.setattr("metawrap2.config.env_prefix", lambda env: str(binary.parent.parent))
    assert fetch.aria2c() == str(binary)

    monkeypatch.setenv("METAWRAP2_ARIA2C", "")
    assert fetch.aria2c() is None, "an empty override forces the curl path"
    monkeypatch.setenv("METAWRAP2_ARIA2C", "/elsewhere/aria2c")
    assert fetch.aria2c() == "/elsewhere/aria2c"
