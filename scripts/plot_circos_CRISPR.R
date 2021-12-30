library(RCircos)
library(GenomicRanges)

# Overwrite RCircos default fuctions.
source("scripts/improvedRCircos.R")

# Keep only the first 4 columns.
GT1 = subset(read.table("mapping/GT1.ins.gz"), V2 != "pT2")[,1:4]
GT2 = subset(read.table("mapping/GT2.ins.gz"), V2 != "pT2")[,1:4]
GT = rbind(GT1, GT2)

gRNA_GT1 = subset(read.table("misc/gRNA_counts_GT1.txt"), V1 %in% GT1$V1)
gRNA_GT2 = subset(read.table("misc/gRNA_counts_GT2.txt"), V1 %in% GT2$V1)
gRNA_GT = rbind(gRNA_GT1, gRNA_GT2)

colnames(GT) = c("barcode", "chrom", "strand", "pos")
colnames(gRNA_GT) = c("barcode", "g1", "g2")

GT = merge(GT, gRNA_GT)

head(GT)

gGT = GRanges(Rle(GT$chrom), IRanges(start=GT$pos, width=1))


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

GT = subset(GT, countOverlaps(gGT, gGRCm38) > 0)

print(nrow(GT))
set.seed(123)
GT = GT[sample(nrow(GT), 1000),]

trackGT = GT[,c("chrom", "pos", "pos")]

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

colGT = colmap(GT$g1, GT$g2)

pdf("figures/CRISPR_circos.pdf", useDingbats=FALSE)
par(mar=c(0,0,0,0))
RCircos.Set.Plot.Area()
RCircos.Draw.Chromosome.Ideogram()
RCircos.Label.Chromosome.Names()
Tile.Plot(trackGT, 1, "in", col=colGT, lwd=.5)
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
