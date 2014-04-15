deseq_varstab = function(physeq, sampleConditions = rep("A", nsamples(physeq)), 
    ...) {
    require("DESeq")
    # Enforce orientation.
    if (!taxa_are_rows(physeq)) {
        physeq <- t(physeq)
    }
    x = as(otu_table(physeq), "matrix")
    # The same tweak as for edgeR to avoid NaN problems that cause the workflow
    # to stall/crash.
    x = x + 1
    # Create annotated data.frame with the taxonomy table, in case it is useful
    # later
    taxADF = as(data.frame(as(tax_table(physeq), "matrix"), stringsAsFactors = FALSE), 
        "AnnotatedDataFrame")
    cds = newCountDataSet(x, sampleConditions, featureData = taxADF)
    # First estimate library size factors
    cds = estimateSizeFactors(cds)
    # Variance estimation, passing along additional options
    cds = estimateDispersions(cds, ...)
    # Determine which column(s) have the dispersion estimates
    dispcol = grep("disp\\_", colnames(fData(cds)))
    # Enforce that there are no infinite values in the dispersion estimates
    if (any(!is.finite(fData(cds)[, dispcol]))) {
        fData(cds)[which(!is.finite(fData(cds)[, dispcol])), dispcol] <- 0
    }
    vsmat = exprs(varianceStabilizingTransformation(cds))
    otu_table(physeq) <- otu_table(vsmat, taxa_are_rows = TRUE)
    return(physeq)
}