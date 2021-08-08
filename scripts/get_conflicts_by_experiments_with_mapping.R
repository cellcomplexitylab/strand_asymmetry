dat = subset(read.table("barcode_view_all_events_with_mapping.txt"), V6 == "test")

# 0 or 1 means no conflict.
conflict = (dat$V2 != dat$V2^2)
dat$V2 = conflict

aggmean = aggregate(dat$V2, FUN=mean, by=list(MM=dat$V3,
   time=dat$V4, method=dat$V5, expt=dat$V13))

write.table(aggmean, file="conflicts_by_experiments_with_mapping.txt",
   sep="\t", quote=FALSE, row.names=FALSE)
