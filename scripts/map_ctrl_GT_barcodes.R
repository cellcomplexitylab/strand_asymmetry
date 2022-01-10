library(GenomicRanges)
barcodes_with_errror = read.table("ctrl_GT_barcodes_with_errors.txt")
barcodes_with_meCpG_like_error = subset(barcodes_with_errror, V5 == 1)
all_barcodes = scan("all_ctrl_GT_barcodes.txt", what="")


GT1 = read.table("mapping/GT1.ins.gz")
GT2 = read.table("mapping/GT2.ins.gz")
ins = subset(rbind(GT1, GT2), V1 %in% all_barcodes)

# Get expression data in mouse ES cells.
exprs = read.delim("misc/GSE93238_gene.fpkm.txt.gz")
head(exprs)
# Remove version number from expression data set (ENSMUSG00000000125.5)
exprs$gene = sub("\\..*", "", exprs$gene)
score = exprs$ES1 + exprs$ES2
cutoff = 0
mean(score > cutoff)
higenes = exprs$gene[score > cutoff]
logenes = exprs$gene[score <= cutoff]
head(higenes)
head(logenes)

# Get mm9 genes and their coordinates.
GTF = read.table("misc/Mus_musculus.NCBIM37.67.gtf.gz", se="\t", as.is=TRUE)
GTF = subset(GTF, !grepl("NT_", V1))
# Parse GTF to get coordinates.
chrom = paste("chr", GTF$V1, sep="")
start = GTF$V4
end = GTF$V5
strand_ = GTF$V7

gene = gsub("^ gene_id (ENSMUSG[0-9]{11});.*", "\\1", GTF$V9)
gene_chrom = tapply(X=chrom, INDEX=gene, unique)
gene_start = tapply(X=start, INDEX=gene, min)
gene_end = tapply(X=end, INDEX=gene, max)
gene_strand = tapply(X=strand_, INDEX=gene, unique)

ENSG = data.frame(gene=names(gene_chrom), chrom=gene_chrom,
      start=gene_start, end=gene_end, strand=gene_strand, row.names=NULL)

# Subset to genes expressed in ES cells.
ENSG = subset(ENSG, gene %in% exprs$gene)
head(ENSG)

gENSG = GRanges(seqnames=Rle(ENSG$chrom),
   ranges=IRanges(start=ENSG$start, end=ENSG$end),
   strand=Rle(ENSG$strand))

ENSG_hi = subset(ENSG, gene %in% higenes)
ENSG_lo = subset(ENSG, gene %in% logenes)
gENSG_hi = GRanges(seqnames=Rle(ENSG_hi$chrom),
   ranges=IRanges(start=ENSG_hi$start, end=ENSG_hi$end),
   strand=Rle(ENSG_hi$strand))
gENSG_lo = GRanges(seqnames=Rle(ENSG_lo$chrom),
   ranges=IRanges(start=ENSG_lo$start, end=ENSG_lo$end),
   strand=Rle(ENSG_lo$strand))

gprom = promoters(gENSG_hi, upstream=1500, downstream=500)

gins = GRanges(Rle(ins$V2), IRanges(start=ins$V4, width=1))

idxlo = findOverlaps(gins, gENSG_lo, select="arbitrary")
idxhi = findOverlaps(gins, gENSG_hi, select="arbitrary")
idxprom = findOverlaps(gins, gprom, select="arbitrary")

ins$gtype = "intergenic"
ins$gtype[!is.na(idxlo)] = "gene/low"
ins$gtype[!is.na(idxhi)] = "gene/high"

ins$is_active_prom = FALSE
ins$is_active_prom[!is.na(idxprom)] = TRUE

chrom = read.delim("misc/CrPs_and_hM_cM_peaks_in_0-95_probability_segments.txt.gz")
gchrom = GRanges(Rle(chrom$chr), IRanges(start=chrom$start_1bases, width=200))

idx = findOverlaps(gins, gchrom, select="arbitrary")

ins = ins[!is.na(idx),]
feat = chrom[idx[!is.na(idx)],9:ncol(chrom)]

feature_table = data.frame(
   bcd=ins$V1, GC1=ins$V6, GC2=ins$V7,
   feat,  # 79 features
   chrom=ins$V2, strand=ins$V3, pos=ins$V4,
   gtype=ins$gtype, is_active_prom=ins$is_active_prom
)

feature_table$meCpG = feature_table$bcd %in% barcodes_with_meCpG_like_error$V1

write.table(feature_table, file="features_ctrl_GT_barcodes.txt",
   sep="\t", quote=FALSE, row.names=FALSE)
