library(GenomicRanges)
library(RCircos)

# Overwrite RCircos default fuctions.
source("scripts/improvedRCircos.R")

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

AG = rbind(AG1, AG2)
TG = rbind(TG1, TG2)
AC = rbind(AC1, AC2)
TC = rbind(TC1, TC2)
CA = rbind(CA1, CA2)

gAG = GRanges(Rle(AG$V2), IRanges(start=AG$V4, width=1))
gTG = GRanges(Rle(TG$V2), IRanges(start=TG$V4, width=1))
gAC = GRanges(Rle(AC$V2), IRanges(start=AC$V4, width=1))
gTC = GRanges(Rle(TC$V2), IRanges(start=TC$V4, width=1))
gCA = GRanges(Rle(CA$V2), IRanges(start=CA$V4, width=1))

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

AG = subset(AG, countOverlaps(gAG, gGRCm38) > 0)
TG = subset(TG, countOverlaps(gTG, gGRCm38) > 0)
AC = subset(AC, countOverlaps(gAC, gGRCm38) > 0)
TC = subset(TC, countOverlaps(gTC, gGRCm38) > 0)
CA = subset(CA, countOverlaps(gCA, gGRCm38) > 0)

print(nrow(AG))
print(nrow(TG))
print(nrow(AC))
print(nrow(TC))
print(nrow(CA))


# Take 1000 at random.
set.seed(123)
trackAG = AG[sample(nrow(AG),1000),c(2,4,4)]
trackTG = TG[sample(nrow(TG),1000),c(2,4,4)]
trackAC = AC[sample(nrow(AC),1000),c(2,4,4)]
trackTC = TC[sample(nrow(TC),1000),c(2,4,4)]
trackCA = CA[sample(nrow(CA),1000),c(2,4,4)]


rcircos.params = RCircos.Get.Plot.Parameters()
rcircos.params$track.height = 0.08
rcircos.params$track.in.start = 1.2
rcircos.params$chrom.paddings = 500

RCircos.Reset.Plot.Parameters(rcircos.params)

COL = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B")

pdf("figures/integ_circos.pdf", useDingbats=FALSE)
par(mar=c(0,0,0,0))
RCircos.Set.Plot.Area()
RCircos.Draw.Chromosome.Ideogram()
RCircos.Label.Chromosome.Names()
Tile.Plot(trackAG, 1, "in", col=COL[1])
Tile.Plot(trackTG, 2, "in", col=COL[2])
Tile.Plot(trackAC, 3, "in", col=COL[3])
Tile.Plot(trackTC, 4, "in", col=COL[4])
Tile.Plot(trackCA, 5, "in", col=COL[5])
dev.off()
