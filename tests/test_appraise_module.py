"""The appraise module's own logic: sample naming, bin extensions, and the read-length guard.

Each of these is a way the module can produce a confident wrong answer rather than an error, which
is why they are tested separately from the parsing in test_recovery.py:

* SingleM labels rows with the read file's basename, so its ``ERR011347_1`` has to become the
  pipeline's ``ERR011347`` or the appraisal joins to nothing;
* SingleM's default genome extension is ``fna`` and MetaWrap2 writes ``.fasta``, so a wrong guess
  appraises zero genomes and reports that none of the community is represented;
* a read shorter than SingleM's 60 bp window carries no marker, so trimmed reads appraise as an
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
    # SingleM's default extension is `fna`; a wrong guess appraises zero genomes.
    for name in ("bin.1.fasta", "bin.2.fasta", "notes.txt"):
        (tmp_path / name).write_text(">c\nACGT\n")
    assert appraise._infer_extension(str(tmp_path), None) == "fasta"
