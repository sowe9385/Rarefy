CSS_varstab = function(x) {
    require("metagenomeSeq")
    obj = biom2MRexperiment(x)
    p = cumNormStatFast(obj)
	data = cumNorm(obj, p = p)
	mat = MRcounts(data, norm = TRUE, log = TRUE)
	exportStats(obj, file = file.path("~/Desktop", "CSS_stats.tsv"))
	head(read.csv(file = file.path("~/Desktop", "CSS_stats.tsv"), sep = "\t"))
    return(mat)
}