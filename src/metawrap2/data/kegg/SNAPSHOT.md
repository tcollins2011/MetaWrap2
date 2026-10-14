# KEGG reference data

| File | Contents | Source |
|---|---|---|
| `modules_table.tsv` | 573 KEGG module definitions (`module`, `definition`, `name`, `class`), taken 2026-10-12 | [kegg-pathways-completeness](https://github.com/EBI-Metagenomics/kegg-pathways-completeness-tool) 1.4.4 |

`metawrap2.kegg` reads it and writes its date into every row. A module's completeness changes when
its definition changes, so results are only comparable over the same definitions.

Two more KEGG files are downloaded at install time rather than shipped:

- `metawrap2 install-db kofam` fetches the KOfam profiles and `ko_list.gz`, their per-profile
  thresholds. This release was validated on the `ko_list.gz` with md5
  `e037ff1135d820c775da78e6d5f256c7` (KOfam, 2026-09-24). `annotate` warns when it finds a different
  one, and every KO row records which list it used.
- `metawrap2 install-db kegg` fetches BRITE `ko00001` and projects it to `ko_hierarchy.tsv.gz`, one
  row per KO and pathway. `SNAPSHOT.json` beside it records the download date and md5, and every
  hierarchy row names it.

## Updating the module definitions

Do it deliberately, never mid-study:

1. Replace `modules_table.tsv` and update the version and date above and `kegg.DEFINITION_SNAPSHOT`.
2. Run `tests/test_kegg.py`. It checks properties of the definitions that the output flags depend
   on, such as 47 single-step modules and 14 with a reaction that has no KO.
3. Note in the changelog that numbers from before and after are not comparable.
