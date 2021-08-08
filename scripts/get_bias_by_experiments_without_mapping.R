dat = subset(read.table("barcode_view_all_events_without_mapping.txt"), V6 == "test")

aggmean = aggregate(dat$V2, FUN=mean, by=list(MM=dat$V3,
   time=dat$V4, method=dat$V5, expt=dat$V7))

write.table(aggmean, file="bias_by_experiments_without_mapping.txt",
   sep="\t", quote=FALSE, row.names=FALSE)
