library(GenomicRanges)

GT = subset(read.table("barcode_view_all_events_with_mapping.txt"), V3 == "GT")
gGTins = GRanges(Rle(GT$V10), IRanges(start=GT$V12, width=1))

chrom = read.delim("misc/CrPs_and_hM_cM_peaks_in_0-95_probability_segments.txt.gz")

gH3K9me3 = GRanges(Rle(chrom$chr[chrom$H3K9me3 == 1]),
      IRanges(start=chrom$start_1bases[chrom$H3K9me3 == 1], width=200))
gH3K36me3 = GRanges(Rle(chrom$chr[chrom$H3K36me3 == 1]),
      IRanges(start=chrom$start_1bases[chrom$H3K36me3 == 1], width=200))

idxH3K9me3 = findOverlaps(gGTins, gH3K9me3, select="arbitrary")
idxH3K36me3 = findOverlaps(gGTins, gH3K36me3, select="arbitrary")

GT$category = "rest"
GT$category[!is.na(idxH3K36me3)] = "H3K36me3"
GT$category[!is.na(idxH3K9me3)] = "H3K9me3"

write.table(GT, file="barcode_view_GT_with_histones.txt",
   sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
