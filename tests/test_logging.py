from metawrap2 import logging as mwlog


def test_banner_borders_and_width():
    out = mwlog.format_comment("hello world", "-")
    lines = out.split("\n")
    # leading blank line + top border + content + bottom border + trailing blank
    assert lines[0] == ""
    assert lines[1] == "-" * 120
    assert lines[-2] == "-" * 120
    assert lines[-1] == ""
    # every content line is exactly 120 chars wide and framed by 5 delimiters
    for line in lines[2:-2]:
        assert len(line) == 120
        assert line.startswith("-----") and line.endswith("-----")
        assert "hello world" in line


def test_delimiters_differ_by_level():
    assert mwlog.format_comment("x", "#").split("\n")[1] == "#" * 120
    assert mwlog.format_comment("x", "*").split("\n")[1] == "*" * 120


def test_long_message_wraps_to_multiple_lines():
    msg = "word " * 40  # far exceeds the 90-char wrap width
    content = [ln for ln in mwlog.format_comment(msg.strip(), "-").split("\n") if set(ln) != {"-"} and ln]
    assert len(content) >= 2
    for line in content:
        assert len(line) == 120


def test_error_exits():
    import pytest

    with pytest.raises(SystemExit) as exc:
        mwlog.error("bad")
    assert exc.value.code == 1
