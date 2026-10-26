# RGI output

`bin_007.txt` is RGI 6.0.8 output against CARD 4.0.2 for one CAMI2 bin. Three columns are removed
(`Predicted_DNA`, `Predicted_Protein`, `CARD_Protein_Sequence`) because they hold full sequences and
are not parsed. All other columns and values are unedited.

The file covers:

- `Cut_Off` tiers: one `Perfect` hit and four `Strict` hits;
- `ORF_ID` values that are Prodigal headers, not gene ids;
- `Best_Identities` down to 27.84% on a Strict hit;
- bare ARO numbers in place of the `ARO:`-prefixed accessions other tools write.
