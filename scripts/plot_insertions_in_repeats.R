library(GenomicRanges)
repeats = read.delim("misc/repeats-mm9.txt.bz2", comm="#")

INDEX = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
   "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
   "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
   "chr19", "chrX", "chrY")

AG1 = subset(read.table("mapping/AG1.ins.gz", as.is=TRUE), V2 %in% INDEX)
AG2 = subset(read.table("mapping/AG2.ins.gz", as.is=TRUE), V2 %in% INDEX)
TG1 = subset(read.table("mapping/TG1.ins.gz", as.is=TRUE), V2 %in% INDEX)
TG2 = subset(read.table("mapping/TG2.ins.gz", as.is=TRUE), V2 %in% INDEX)
AC1 = subset(read.table("mapping/AC1.ins.gz", as.is=TRUE), V2 %in% INDEX)
AC2 = subset(read.table("mapping/AC2.ins.gz", as.is=TRUE), V2 %in% INDEX)
TC1 = subset(read.table("mapping/TC1.ins.gz", as.is=TRUE), V2 %in% INDEX)
TC2 = subset(read.table("mapping/TC2.ins.gz", as.is=TRUE), V2 %in% INDEX)
GT1 = subset(read.table("mapping/GT1.ins.gz", as.is=TRUE), V2 %in% INDEX)
GT2 = subset(read.table("mapping/GT2.ins.gz", as.is=TRUE), V2 %in% INDEX)

AG = rbind(AG1, AG2)
TG = rbind(TG1, TG2)
AC = rbind(AC1, AC2)
TC = rbind(TC1, TC2)
GT = rbind(GT1, GT2)

gAG = GRanges(Rle(AG$V2), IRanges(start=AG$V4, width=1))
gTG = GRanges(Rle(TG$V2), IRanges(start=TG$V4, width=1))
gAC = GRanges(Rle(AC$V2), IRanges(start=AC$V4, width=1))
gTC = GRanges(Rle(TC$V2), IRanges(start=TC$V4, width=1))
gGT = GRanges(Rle(GT$V2), IRanges(start=GT$V4, width=1))

grepeats = makeGRangesFromDataFrame(repeats,
   seqnames.field="genoName", start.field="genoStart",
   end.field="genoEnd")

total_AG = sum(countOverlaps(gAG, grepeats))
total_TG = sum(countOverlaps(gTG, grepeats))
total_AC = sum(countOverlaps(gAC, grepeats))
total_TC = sum(countOverlaps(gTC, grepeats))
total_GT = sum(countOverlaps(gGT, grepeats))

ins = data.frame(
   total    = c(total_AG, total_TG, total_AC, total_TC, total_GT),
   fraction = c(total_AG, total_TG, total_AC, total_TC, total_GT) /
              c(nrow(AG), nrow(TG), nrow(AC), nrow(TC), nrow(GT)))
print(ins)

COL = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B")

pdf("figures/insertions_in_repeats.pdf", width=4, useDingbats=FALSE)
par(mar=c(1,3.5,.5,.5))
barplot(ins$fraction, col=COL, ylim=c(0,.5),
   bty="n", xaxt="n", yaxt="n", col.lab="grey25",
   line=2.2, ylab="Insertions in repeats")
abline(h=(1:4)/10, lwd=2, col="white")
axis(side=2, col="grey50", cex.axis=.8, col.axis="grey25",
   at=(0:5)/10)
dev.off()
