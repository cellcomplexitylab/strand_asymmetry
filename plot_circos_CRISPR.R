library(RCircos)
library(GenomicRanges)

# Overwrite RCircos default fuctions.
source("scripts/improvedRCircos.R")

# Keep only the first 4 columns.
CA1 = subset(read.table("mapping/CA1.ins.gz"), V2 != "pT2")[,1:4]
CA2 = subset(read.table("mapping/CA2.ins.gz"), V2 != "pT2")[,1:4]

gRNA_CA1 = subset(read.table("misc/gRNA_counts_CA1.txt"), V1 %in% CA1$V1)
gRNA_CA2 = subset(read.table("misc/gRNA_counts_CA2.txt"), V1 %in% CA2$V1)

colnames(CA1) = c("barcode", "chrom", "strand", "pos")
colnames(CA2) = c("barcode", "chrom", "strand", "pos")
colnames(gRNA_CA1) = c("barcode", "g1", "g2")
colnames(gRNA_CA2) = c("barcode", "g1", "g2")

CA1 = merge(CA1, gRNA_CA1)
CA2 = merge(CA2, gRNA_CA2)

head(CA1)
head(CA2)

gCA1 = GRanges(Rle(CA1$chrom), IRanges(start=CA1$pos, width=1))
gCA2 = GRanges(Rle(CA2$chrom), IRanges(start=CA2$pos, width=1))


data(UCSC.Mouse.GRCm38.CytoBandIdeogram)
chr.exclude = NULL
cyto.info = UCSC.Mouse.GRCm38.CytoBandIdeogram
tracks.inside = 5
tracks.outside = 0
RCircos.Set.Core.Components(cyto.info, chr.exclude,
      tracks.inside, tracks.outside)

# Keep only insertions mapping to annotations.
gGRCm38 = makeGRangesFromDataFrame(UCSC.Mouse.GRCm38.CytoBandIdeogram,
   seqnames.field="Chromosome", start.field="ChromStart",
   end.field="ChromEnd")

CA1 = subset(CA1, countOverlaps(gCA1, gGRCm38) > 0)
CA2 = subset(CA2, countOverlaps(gCA2, gGRCm38) > 0)

print(nrow(CA1))
print(nrow(CA2))


trackCA1 = CA1[,c("chrom", "pos", "pos")]
trackCA2 = CA2[,c("chrom", "pos", "pos")]


rcircos.params = RCircos.Get.Plot.Parameters()
rcircos.params$track.height = 0.08
rcircos.params$track.in.start = 1.2
rcircos.params$chrom.paddings = 500

RCircos.Reset.Plot.Parameters(rcircos.params)

colmap = function(x, y) {
   # Construct 2D color map.
   NCOL = 16
   left = colorRampPalette(c("#00047d", "#fa017d"))(NCOL)
   right = colorRampPalette(c("#00fb89", "#f6fb89"))(NCOL)
   map = matrix(NA, nrow=NCOL, ncol=NCOL)
   for (i in 1:NCOL) {
      map[i,] = colorRampPalette(c(left[i], right[i]))(NCOL)
   }
   # Return colors as indexed.
   idx = matrix(c(ceiling(.5 + (NCOL-1) * x),
         ceiling(.5 + (NCOL-1) * y)), ncol=2)
   return(map[idx])
}

colCA1 = colmap(CA1$g1, CA1$g2)
colCA2 = colmap(CA2$g1, CA2$g2)

pdf("figures/CRISPR_circos.pdf", useDingbats=FALSE)
par(mar=c(0,0,0,0))
RCircos.Set.Plot.Area()
RCircos.Draw.Chromosome.Ideogram()
RCircos.Label.Chromosome.Names()
Tile.Plot(trackCA1, 1, "in", col=colCA1, lwd=.5)
Tile.Plot(trackCA2, 2, "in", col=colCA2, lwd=.5)
dev.off()

# Plot color map.
pdf("figures/CRISPR_circos_color_map.pdf", useDingbats=FALSE)
par(mar=c(0,0,0,0))
plot(c(0,16), c(0,16), type="n", bty="n", xaxt="n", yaxt="n",
     xlab="", ylab="")
all.x = rep(1:16, each=16)
all.y = rep(1:16, times=16)
for (i in 1:256) {
   x = all.x[i]
   y = all.y[i]
   col = colmap((x-.5)/15, (y-.5)/15)
   rect(xleft=x-1, xright=x, ybottom=y-1, ytop=y, col=col, border=col)
}
dev.off()

#pdf("figures/CA_CRISPR.pdf", useDingbats=FALSE, height=2, width=6)
#par(mar=c(0,0,0,0))
#plot(c(0, max(CA1$pos, CA2$pos)), c(0,3), type="n")
#points(CA1$pos, rep(1, nrow(CA1)), pch=19, col=colmap(CA1$g1, CA1$g2))
#points(CA2$pos, rep(2, nrow(CA2)), pch=19, col=colmap(CA2$g1, CA2$g2))
##plot(c(0, max(CA1$pos, CA2$pos)), c(0,3), type="n")
##points(CA1$pos, CA1$g1, pch=19, cex=.8, col="#fa017d80")
##points(CA1$pos, CA1$g2, pch=19, cex=.8, col="#00047d80")
##points(CA2$pos, CA2$g1 + 2, pch=19, cex=.8, col="#fa017d80")
##points(CA2$pos, CA2$g2 + 2, pch=19, cex=.8, col="#00047d80")
##abline(h=c(.5, 2.5), col="gray80")
#dev.off()
