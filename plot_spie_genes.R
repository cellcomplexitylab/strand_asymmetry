library(GenomicRanges)

source("scripts/spie.R")


AG1 = subset(read.table("mapping/AG1.ins"), V2 != "pT2")
AG2 = subset(read.table("mapping/AG2.ins"), V2 != "pT2")
TG1 = subset(read.table("mapping/TG1.ins"), V2 != "pT2")
TG2 = subset(read.table("mapping/TG2.ins"), V2 != "pT2")
AC1 = subset(read.table("mapping/AC1.ins"), V2 != "pT2")
AC2 = subset(read.table("mapping/AC2.ins"), V2 != "pT2")
TC1 = subset(read.table("mapping/TC1.ins"), V2 != "pT2")
TC2 = subset(read.table("mapping/TC2.ins"), V2 != "pT2")
CA1 = subset(read.table("mapping/CA1.ins"), V2 != "pT2")
CA2 = subset(read.table("mapping/CA2.ins"), V2 != "pT2")

allins = rbind(AG1, AG2, TG1, TG2, AC1, AC2, TC1, TC2, CA1, CA2)
gallins = GRanges(Rle(allins$V2), IRanges(start=allins$V4, width=1))

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

Olo = sum(countOverlaps(gallins, gENSG_lo))
Ohi = sum(countOverlaps(gallins, gENSG_hi))
Orest = nrow(allins) - Olo - Ohi

# Total size of mm9: 2780285931.
Elo = sum(sapply(coverage(gENSG_lo), sum))
Ehi = sum(sapply(coverage(gENSG_hi), sum))
Erest = 2780285931. - Elo - Ehi

obsv = c(Ohi, Olo, Orest)
expt = c(Ehi, Elo, Erest)

print(obsv)
print(expt / sum(expt) * sum(obsv))

COL=c("#0B2639", "#598078", "#E1DBC0")

pdf("figures/spie_genes.pdf", useDingbats=FALSE)
par(mar=c(0,0,0,0))
spiechart(expected=expt, observed=obsv, col=COL)
dev.off()
