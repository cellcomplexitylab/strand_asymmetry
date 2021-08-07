library(RCircos)
library(GenomicRanges)

# Overwrite RCircos default fuctions.
source("scripts/improvedRCircos.R")

# Keep only the first 4 columns.
GT1 = subset(read.table("mapping/GT1.ins.gz"), V2 != "pT2")[,1:4]
GT2 = subset(read.table("mapping/GT2.ins.gz"), V2 != "pT2")[,1:4]

gRNA_GT1 = subset(read.table("misc/gRNA_counts_GT1.txt"), V1 %in% GT1$V1)
gRNA_GT2 = subset(read.table("misc/gRNA_counts_GT2.txt"), V1 %in% GT2$V1)

colnames(GT1) = c("barcode", "chrom", "strand", "pos")
colnames(GT2) = c("barcode", "chrom", "strand", "pos")
colnames(gRNA_GT1) = c("barcode", "g1", "g2")
colnames(gRNA_GT2) = c("barcode", "g1", "g2")

GT1 = merge(GT1, gRNA_GT1)
GT2 = merge(GT2, gRNA_GT2)

head(GT1)
head(GT2)

gGT1 = GRanges(Rle(GT1$chrom), IRanges(start=GT1$pos, width=1))
gGT2 = GRanges(Rle(GT2$chrom), IRanges(start=GT2$pos, width=1))


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

GT1 = subset(GT1, countOverlaps(gGT1, gGRCm38) > 0)
GT2 = subset(GT2, countOverlaps(gGT2, gGRCm38) > 0)

print(nrow(GT1))
print(nrow(GT2))


trackGT1 = GT1[,c("chrom", "pos", "pos")]
trackGT2 = GT2[,c("chrom", "pos", "pos")]


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

colGT1 = colmap(GT1$g1, GT1$g2)
colGT2 = colmap(GT2$g1, GT2$g2)

pdf("figures/CRISPR_circos.pdf", useDingbats=FALSE)
par(mar=c(0,0,0,0))
RCircos.Set.Plot.Area()
RCircos.Draw.Chromosome.Ideogram()
RCircos.Label.Chromosome.Names()
Tile.Plot(trackGT1, 1, "in", col=colGT1, lwd=.5)
Tile.Plot(trackGT2, 2, "in", col=colGT2, lwd=.5)
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
