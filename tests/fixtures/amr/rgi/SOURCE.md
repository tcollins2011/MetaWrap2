# RGI output, as it really looks

`bin_007.txt` is a real RGI 6.0.8 run against CARD 4.0.2 over one CAMI2 bin, with three columns
removed (`Predicted_DNA`, `Predicted_Protein`, `CARD_Protein_Sequence`) because they are whole
sequences and nothing here parses them. Every other column and value is unedited.

It is the fixture for what matters about RGI: `Cut_Off` tiers (one `Perfect`, four `Strict`), an
`ORF_ID` that is a Prodigal header rather than a gene id, `Best_Identities` as low as 27.84% on a
Strict hit, and bare ARO numbers where every other tool writes `ARO:`-prefixed accessions.
