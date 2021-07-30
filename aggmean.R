library(ggplot2)
dat = subset(read.table("features.txt"), V6 == "test")
# Aggregate the data.

# Mean-aggregate scores per barcode, also return the 
# length, which will be useful for computing the
# standard deviation.
aggmean = aggregate(dat$V2, FUN=mean, by=list(brcd=dat$V1,
           chrom=dat$V10, strand=dat$V11, pos=dat$V12,
           MM=dat$V3, GC1=dat$V8, GC2=dat$V9))

write.table(aggmean, file="aggfeat.txt",
   sep="\t", quote=FALSE, row.names=FALSE)
