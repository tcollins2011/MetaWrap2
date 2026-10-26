"""GUNC: report parsing, the stats annotation, and the wiring into bin_refinement."""

import pytest

from metawrap2 import gunc

REPORT_HEADER = (
    "genome\tn_genes_called\tn_genes_mapped\ttaxonomic_level\t"
    "clade_separation_score\tpass.GUNC\n"
)


def _report(tmp_path, rows, name="GUNC.progenomes_2.1.maxCSS_level.tsv"):
    out = tmp_path / "gunc"
    out.mkdir(exist_ok=True)
    (out / name).write_text(REPORT_HEADER + "".join(rows))
    return str(out)


def test_parse_report_prefers_guncs_own_verdict_over_the_threshold(tmp_path):
    # Row two is below the CSS threshold but GUNC failed it because too few of its genes are in
    # the reference. Recomputing from CSS alone would pass it.
    out = _report(
        tmp_path,
        [
            "bin_001\t900\t850\tkingdom\t0.9100\tFalse\n",
            "bin_002\t900\t20\tphylum\t0.0100\tFalse\n",
            "bin_003\t900\t880\tphylum\t0.0200\tTrue\n",
        ],
    )
    verdicts = gunc.parse_report(gunc.find_report(out))
    assert [v.passed for v in (verdicts["bin_001"], verdicts["bin_002"], verdicts["bin_003"])] == [
        False,
        False,
        True,
    ]
    assert verdicts["bin_001"].css == pytest.approx(0.91)
    assert verdicts["bin_001"].level == "kingdom"

    # GUNC derives genome names with basename.split(file_suffix)[0], so "fasta" without the dot
    # would name bins "bin_001." and match nothing in the stats table.
    cmd = gunc.GUNC_RUN.format(bins="b", ext=".fasta", db="db.dmnd", out="o", threads=4)
    assert "--file_suffix .fasta" in cmd
