library(GenomicRanges)

ins = subset(read.table("barcode_view_all_events_with_mapping.txt"), V6 == "test")
gins = GRanges(Rle(ins$V10), IRanges(start=ins$V12, width=1))

chrom = read.delim("misc/CrPs_and_hM_cM_peaks_in_0-95_probability_segments.txt.gz")
gchrom = GRanges(Rle(chrom$chr), IRanges(start=chrom$start_1bases, width=200))


idx = findOverlaps(gins, gchrom, select="arbitrary")

ins = ins[!is.na(idx),]
feat = chrom[idx[!is.na(idx)],9:ncol(chrom)]

# One-hot encoding of the construct
cstr = with(ins, model.matrix(~0 + V3))
colnames(cstr) = sub("V3", "", colnames(cstr))
# Time point (0 for 24 h, 1 for 48 h)
time = ifelse(ins$V4 == 48, 1, 0)
# Time point (0 for LA, 1 for 6xPCR)
method = ifelse(ins$V5 == "6xPCR", 1, 0)
# Replicate (0-based).
rep = ins$V7 - 1
# GC1 and GC2
GC1 = ins$V8
GC2 = ins$V9

feature_table = data.frame(
   y = ins$V2,                    # output variable
   cstr[,1:4], time, method, rep, # 7 features
   GC1, GC2, feat,                # 79 features
   chrom=ins$V10, strand=ins$V11, pos=ins$V12, expt=ins$V13,
   bcd = ins$V1)                  # barcode

# Remove a few NAs for the G+C content.
feature_table = subset(feature_table, complete.cases(feature_table))

write.table(feature_table, file="barcode_view_all_features_with_mapping.txt",
   sep="\t", quote=FALSE, row.names=FALSE)
