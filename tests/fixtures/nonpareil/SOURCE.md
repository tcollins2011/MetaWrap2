# Nonpareil example curves

`HumanGut.npo`, `IowaSoil.npo` and `LakeLanier.npo` are copied unmodified from Nonpareil 3.5.5,
`utils/Nonpareil/inst/extdata/`, where they ship as the R package's example data (Artistic-2.0).

They are here because `metawrap2.nonpareil` ports the R package's model fit into Python, and a port
of numerical code is only worth having if it agrees with the original. The expected values in
`tests/test_nonpareil.py` are not this implementation's own output: they were produced by running
the R package itself on these three files, in a throwaway conda env built from bioconda's
`nonpareil` 3.5.4 package, with

    np <- Nonpareil.curve(file, plot = FALSE)
    c(np$C, np$kappa, np$diversity, np$LRstar, np$modelR, coef(np$model))

The Python port reproduces all six quantities on all three datasets to six significant figures.
