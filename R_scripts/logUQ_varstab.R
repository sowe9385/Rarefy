logUQ_varstab = function(x) {
    require("metagenomeSeq")
    obj = biom2MRexperiment(x)
	counts = MRcounts(obj)
	xx=counts
	xx[xx==0]=NA
	s75 = colQuantiles(xx,p=.75, na.rm=T)
	normalizedCounts = log2( sweep(counts, 2, s75, "/") + 1)
    return(normalizedCounts)
}