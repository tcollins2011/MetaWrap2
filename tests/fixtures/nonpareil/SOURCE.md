# Nonpareil example curves

`HumanGut.npo`, `IowaSoil.npo` and `LakeLanier.npo` are copied unmodified from Nonpareil 3.5.5,
`utils/Nonpareil/inst/extdata/`, where they ship as the R package's example data (Artistic-2.0).

`metawrap2.nonpareil` is a Python port of the R package's model fit. The expected values in
`tests/test_nonpareil.py` come from the R package, run on these three files in a separate conda env
built from bioconda's `nonpareil` 3.5.4 package:

    np <- Nonpareil.curve(file, plot = FALSE)
    c(np$C, np$kappa, np$diversity, np$LRstar, np$modelR, coef(np$model))

The Python port matches all six quantities on all three datasets to six significant figures.
