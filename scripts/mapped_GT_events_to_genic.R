library(GenomicRanges)

GT = subset(read.table("barcode_view_all_events_with_mapping.txt"), V3 == "GT")

gGTins = GRanges(Rle(GT$V10), IRanges(start=GT$V12, width=1))

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
GTF = read.table("misc/Mus_musculus.NCBIM37.67.gtf.gz", se="\t")
GTF = subset(GTF, !grepl("NT_", V1))
# Parse GTF to get coordinates.
chrom = paste("chr", GTF$V1, sep="")
start = GTF$V4
end = GTF$V5
gene = gsub("^ gene_id (ENSMUSG[0-9]{11});.*", "\\1", GTF$V9)
gene_chrom = tapply(X=chrom, INDEX=gene, unique)
gene_start = tapply(X=start, INDEX=gene, min)
gene_end = tapply(X=end, INDEX=gene, max)
ENSG = data.frame(gene=names(gene_chrom), chrom=gene_chrom,
      start=gene_start, end=gene_end, row.names=NULL)
# Subset to genes expressed in ES cells.
ENSG = subset(ENSG, gene %in% exprs$gene)
head(ENSG)
gENSG = GRanges(Rle(ENSG$chrom),
      IRanges(start=ENSG$start, end=ENSG$end))

ENSG_hi = subset(ENSG, gene %in% higenes)
ENSG_lo = subset(ENSG, gene %in% logenes)
gENSG_hi = GRanges(Rle(ENSG_hi$chrom),
      IRanges(start=ENSG_hi$start, end=ENSG_hi$end))
gENSG_lo = GRanges(Rle(ENSG_lo$chrom),
      IRanges(start=ENSG_lo$start, end=ENSG_lo$end))

idxlo = findOverlaps(gGTins, gENSG_lo, select="arbitrary")
idxhi = findOverlaps(gGTins, gENSG_hi, select="arbitrary")

GT$category = "intergenic"
GT$category[!is.na(idxlo)] = "gene/low"
GT$category[!is.na(idxhi)] = "gene/high"

write.table(GT, file="barcode_view_GT_with_genic.txt",
   sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
