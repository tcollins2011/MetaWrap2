"""Tests for the two scripts rewritten out of their Python-2-era form.

Both rewrites were validated against the originals on real data before these tests were
written: ``choose_best_bin`` produced the same winners, and
``filter_reads_for_bin_reassembly`` produced byte-for-byte identical FASTQ output for all 12
files of a real bwa-mem stream. These tests pin the behaviour those checks confirmed.
"""

import os

from metawrap2.scripts import choose_best_bin as cbb
from metawrap2.scripts import filter_reads_for_bin_reassembly as frb

# --- choose_best_bin --------------------------------------------------------------------

STATS = (
    "bin\tcompleteness\tcontamination\tGC\tlineage\tN50\tsize\n"
    "bin.1.orig\t80.0\t2.0\t0.51\tBacteria\t20000\t3000000\n"
    "bin.1.strict\t88.0\t1.5\t0.51\tBacteria\t35000\t3200000\n"
    "bin.1.permissive\t90.0\t3.0\t0.51\tBacteria\t40000\t3300000\n"
    "bin.2.orig\t75.0\t4.0\t0.40\tBacteroides\t15000\t2500000\n"
    "bin.2.strict\t79.0\t3.0\t0.40\tBacteroides\t21000\t2600000\n"
    "bin.2.permissive\t77.0\t5.0\t0.40\tBacteroides\t19000\t2600000\n"
)


def _write(tmp_path, text=STATS):
    path = tmp_path / "reassembled_bins.stats"
    path.write_text(text)
    return str(path)


def test_contamination_is_weighted_five_times_completeness():
    # 88.0 + 5*(100-1.5) = 580.5  beats  90.0 + 5*(100-3.0) = 575.0
    strict = cbb.BinVersion("bin.1", "strict", 88.0, 1.5, 35000)
    permissive = cbb.BinVersion("bin.1", "permissive", 90.0, 3.0, 40000)
    assert strict.score > permissive.score
    assert strict.beats(permissive)


def test_ties_are_broken_by_n50():
    a = cbb.BinVersion("bin.1", "strict", 80.0, 2.0, 30000)
    b = cbb.BinVersion("bin.1", "permissive", 80.0, 2.0, 10000)
    assert a.score == b.score
    assert a.beats(b) and not b.beats(a)


def test_parse_stats_splits_bin_and_style():
    versions = cbb.parse_stats(STATS.splitlines())
    assert len(versions) == 6
    assert {v.bin_name for v in versions} == {"bin.1", "bin.2"}
    assert {v.style for v in versions} == {"orig", "strict", "permissive"}
    assert versions[0].full_name == "bin.1.orig"


def test_chooses_the_expected_winners(tmp_path):
    # matches what the original script produced for this input
    assert cbb.best_bin_names(_write(tmp_path), 5, 30) == ["bin.1.strict", "bin.2.strict"]


def test_versions_failing_thresholds_are_not_candidates(tmp_path):
    # only bin.1.permissive clears 89% completeness
    assert cbb.best_bin_names(_write(tmp_path), 89, 30) == ["bin.1.permissive"]
    # nothing clears 99%
    assert cbb.best_bin_names(_write(tmp_path), 99, 30) == []
    # a strict contamination ceiling excludes the otherwise-best versions
    assert cbb.best_bin_names(_write(tmp_path), 5, 1.6) == ["bin.1.strict"]


def test_summarize_counts_winning_styles(tmp_path):
    assert cbb.summarize(_write(tmp_path), 5, 30) == (0, 2, 0)  # (orig, strict, permissive)


def test_unreadable_rows_are_skipped(tmp_path):
    bad = STATS + "bin.3.orig\tnot-a-number\t1.0\t0.4\tBacteria\t100\t200\n" + "short\trow\n"
    assert cbb.best_bin_names(_write(tmp_path, bad), 5, 30) == ["bin.1.strict", "bin.2.strict"]


def test_choose_best_bin_cli_usage_error(capsys):
    assert cbb.main([]) == 2
    assert "usage:" in capsys.readouterr().err


# --- filter_reads_for_bin_reassembly ------------------------------------------------------


def _sam(name, flag, ref, seq="ACGT", qual="IIII", nm=0):
    tags = ["NM:i:%d" % nm] if nm is not None else []
    return "\t".join([name, str(flag), ref, "1", "60", "4M", "=", "1", "0", seq, qual] + tags)


def _bins(tmp_path):
    folder = tmp_path / "bins"
    folder.mkdir()
    (folder / "bin.1.fa").write_text(">c1 total_depth=41.98 sample_depths=1,2\nACGT\n>c2\nTTTT\n")
    (folder / "bin.2.fa").write_text(">c3\nGGGG\n")
    (folder / "bins.stats").write_text("not a genome\n")  # must be ignored
    return str(folder)


def test_contig_bins_keyed_on_contig_id_not_full_header(tmp_path):
    """Bins carrying metaBAT2's header annotations must still match their own alignments."""
    mapping = frb.load_contig_bins(_bins(tmp_path))
    assert mapping == {"c1": "bin.1", "c2": "bin.1", "c3": "bin.2"}


def test_reverse_complement():
    assert frb.reverse_complement("ACGTN") == "NACGT"
    assert frb.reverse_complement("acgt") == "acgt"


def test_iter_pairs_uses_bitwise_flags_not_string_indexing():
    """bin(flag) indexing raised IndexError for small flags; masks do not."""
    stream = iter(
        [
            "@HD\tVN:1.6",  # header, skipped
            _sam("r1", 0x40, "c1"),  # first mate
            _sam("r1", 0x80, "c1"),  # second mate
            _sam("r2", 0, "c1"),  # flag 0: used to crash, now simply not a pair
            "truncated\tline",  # too few fields, skipped
            _sam("r3", 0x40, "c2"),
            _sam("r3", 0x80, "c2"),
        ]
    )
    pairs = list(frb.iter_pairs(stream))
    assert [p[0][0] for p in pairs] == ["r1", "r3"]


def test_mismatches_sums_nm_tags():
    assert frb.mismatches(_sam("r", 0x40, "c1", nm=3).split("\t")) == 3
    assert frb.mismatches(_sam("r", 0x40, "c1", nm=None).split("\t")) == 0
    assert frb.has_nm_tag(_sam("r", 0x40, "c1", nm=0).split("\t")) is True
    assert frb.has_nm_tag(_sam("r", 0x40, "c1", nm=None).split("\t")) is False


def test_pair_must_agree_on_one_bin(tmp_path):
    mapping = frb.load_contig_bins(_bins(tmp_path))
    same = (_sam("r", 0x40, "c1").split("\t"), _sam("r", 0x80, "c2").split("\t"))
    assert frb.bin_for_pair(*same, contig_bins=mapping) == "bin.1"  # c1+c2 are both bin.1
    across = (_sam("r", 0x40, "c1").split("\t"), _sam("r", 0x80, "c3").split("\t"))
    assert frb.bin_for_pair(*across, contig_bins=mapping) is None  # different bins
    unmapped = (_sam("r", 0x40, "*").split("\t"), _sam("r", 0x80, "*").split("\t"))
    assert frb.bin_for_pair(*unmapped, contig_bins=mapping) is None
    unbinned = (_sam("r", 0x40, "cX").split("\t"), _sam("r", 0x80, "cX").split("\t"))
    assert frb.bin_for_pair(*unbinned, contig_bins=mapping) is None


def test_reverse_aligned_reads_are_restored_to_original_orientation():
    forward = _sam("r", 0x40 | frb.FLAG_REVERSE, "c1", seq="ACGT", qual="ABCD").split("\t")
    record = frb.as_fastq(forward, 1)
    assert record == "@r/1\nACGT\n+\nDCBA\n".replace("ACGT", "ACGT")
    # sequence reverse-complemented, quality reversed
    assert record.splitlines()[1] == frb.reverse_complement("ACGT")
    assert record.splitlines()[3] == "ABCD"[::-1]


def test_strict_and_permissive_cutoffs(tmp_path):
    bins = _bins(tmp_path)
    out = tmp_path / "out"
    stream = iter(
        [
            # 0 mismatches -> both sets
            _sam("clean", 0x40, "c1", nm=0),
            _sam("clean", 0x80, "c1", nm=0),
            # 3 mismatches total -> permissive only (strict cutoff is 2)
            _sam("mid", 0x40, "c1", nm=1),
            _sam("mid", 0x80, "c1", nm=2),
            # 9 mismatches -> neither
            _sam("noisy", 0x40, "c1", nm=4),
            _sam("noisy", 0x80, "c1", nm=5),
        ]
    )
    counts = frb.filter_reads(stream, bins, str(out), strict=2, permissive=5)
    assert counts == {"pairs": 3, "recruited": 3, "strict": 1, "permissive": 2}

    strict_1 = (out / "bin.1.strict_1.fastq").read_text()
    permissive_1 = (out / "bin.1.permissive_1.fastq").read_text()
    assert "@clean/1" in strict_1 and "@mid/1" not in strict_1
    assert "@clean/1" in permissive_1 and "@mid/1" in permissive_1
    assert "@noisy/1" not in permissive_1
    # mates land in the _2 file with a /2 suffix
    assert "@clean/2" in (out / "bin.1.strict_2.fastq").read_text()


def test_only_bins_that_recruited_reads_get_files(tmp_path):
    bins = _bins(tmp_path)
    out = tmp_path / "out"
    stream = iter([_sam("r", 0x40, "c1", nm=0), _sam("r", 0x80, "c1", nm=0)])
    frb.filter_reads(stream, bins, str(out), strict=2, permissive=5)
    produced = sorted(os.listdir(out))
    assert produced == [
        "bin.1.permissive_1.fastq",
        "bin.1.permissive_2.fastq",
        "bin.1.strict_1.fastq",
        "bin.1.strict_2.fastq",
    ]


def test_pair_with_no_nm_tag_on_either_mate_is_not_recruited(tmp_path):
    bins = _bins(tmp_path)
    out = tmp_path / "out"
    stream = iter([_sam("r", 0x40, "c1", nm=None), _sam("r", 0x80, "c1", nm=None)])
    counts = frb.filter_reads(stream, bins, str(out), strict=2, permissive=5)
    assert counts["recruited"] == 0


def test_filter_reads_cli_usage_error(capsys):
    assert frb.main([]) == 2
    assert "usage:" in capsys.readouterr().err
