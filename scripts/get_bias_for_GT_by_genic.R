dat = subset(read.table("barcode_view_GT_with_genic.txt"), V6 == "test")

aggmean = aggregate(dat$V2, FUN=mean, by=list(MM=dat$V3,
   time=dat$V4, method=dat$V5, genic=dat$V14, expt=dat$V13))

write.table(aggmean, file="bias_for_GT_by_genic.txt",
   sep="\t", quote=FALSE, row.names=FALSE)
