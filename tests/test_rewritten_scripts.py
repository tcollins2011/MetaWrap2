"""Tests for the two scripts rewritten out of their Python-2-era form.

Both rewrites were validated against the originals on real data before these tests were
written: ``choose_best_bin`` produced the same winners, and
``filter_reads_for_bin_reassembly`` produced byte-for-byte identical FASTQ output for all 12
files of a real bwa-mem stream. This file pins choose_best_bin; the read filter is pinned in
test_reassembly_read_filters.py.
"""

from metawrap2.scripts import choose_best_bin as cbb

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


def test_chooses_the_expected_winners(tmp_path):
    # end-to-end: parse, score, threshold - matches what the original script produced
    assert cbb.best_bin_names(_write(tmp_path), 5, 30) == ["bin.1.strict", "bin.2.strict"]

    # Contamination is weighted five times completeness:
    # 88.0 + 5*(100-1.5) = 580.5  beats  90.0 + 5*(100-3.0) = 575.0
    strict = cbb.BinVersion("bin.1", "strict", 88.0, 1.5, 35000)
    permissive = cbb.BinVersion("bin.1", "permissive", 90.0, 3.0, 40000)
    assert strict.score > permissive.score
    assert strict.beats(permissive)
