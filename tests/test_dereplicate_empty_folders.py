"""dereplicate skips empty bin folders with a warning; only an all-empty input stops it."""

import pytest

from metawrap2.modules import dereplicate


def _folder(path, n):
    path.mkdir()
    for i in range(n):
        (path / ("bin.%d.fa" % i)).write_text(">c%d\nACGT\n" % i)
    return str(path)


class _Stop(Exception):
    pass


def test_an_empty_folder_is_dropped(tmp_path, monkeypatch, capsys):
    full = _folder(tmp_path / "s1", 2)
    empty = _folder(tmp_path / "s2", 0)
    seen = {}

    def stop(_name, checks):
        seen["paths"] = [path for _check, path, _what in checks]
        raise _Stop

    monkeypatch.setattr(dereplicate, "validate_inputs", stop)
    with pytest.raises(_Stop):
        dereplicate.main(["-o", str(tmp_path / "out"), "-b", full, "-b", empty])
    assert [p.rstrip("/").split("/")[-1] for p in seen["paths"]] == ["s1"]
    out = capsys.readouterr()
    assert "no bins" in out.out + out.err


def test_every_folder_empty_is_an_error(tmp_path):
    a = _folder(tmp_path / "s1", 0)
    b = _folder(tmp_path / "s2", 0)
    with pytest.raises(SystemExit):
        dereplicate.main(["-o", str(tmp_path / "out"), "-b", a, "-b", b])
