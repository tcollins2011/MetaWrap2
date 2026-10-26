"""Tests for the appraise module's sample naming, bin extension and read-length handling.

These cases produce wrong results without an error, so they are tested apart from the parsing in
test_recovery.py:

* SingleM labels rows with the read file's basename; its ``ERR011347_1`` must map to the
  pipeline's ``ERR011347`` or the appraisal joins to nothing.
* SingleM's default genome extension is ``fna`` and MetaWrap2 writes ``.fasta``; a wrong guess
  appraises zero genomes.
* A read shorter than SingleM's 60 bp window carries no marker, so trimmed reads appraise as an
  empty community.
"""

from __future__ import annotations

from metawrap2.modules import appraise


def test_singlem_label_and_pipeline_sample_name_are_both_derived(tmp_path):
    for path, label, sample in (
        ("/x/ERR011347_1.fastq.gz", "ERR011347_1", "ERR011347"),
        ("/x/gut_R1.fq", "gut_R1", "gut"),
        ("/x/sample.fastq", "sample", "sample"),
    ):
        assert appraise._singlem_label(path) == label
        assert appraise._sample_name(path) == sample
    # The extension is inferred from the bin files, since SingleM defaults to `fna`.
    for name in ("bin.1.fasta", "bin.2.fasta", "notes.txt"):
        (tmp_path / name).write_text(">c\nACGT\n")
    assert appraise._infer_extension(str(tmp_path), None) == "fasta"
