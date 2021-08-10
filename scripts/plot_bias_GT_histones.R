biasdat = read.delim("bias_for_GT_by_histones.txt")

POS1 = c(1, 3.5, 6)
names(POS1) = c("H3K9me3", "H3K36me3", "rest")
POS2 = c(0,1)
names(POS2) = c("24", "48")
COL = c("#30110D80", "#72262080", "#F2BC9480")
names(COL) = c("H3K9me3", "H3K36me3", "rest")
SYM = c(1,2)
names(SYM) = c("LA", "6xPCR")

pos = POS1[as.character(biasdat$histone)] + POS2[as.character(biasdat$time)]
col = COL[as.character(biasdat$histone)]
pch1 = SYM[as.character(biasdat$method)]
pch2 = pch1 + 15

set.seed(123)
jitter = rnorm(n=length(pos), sd=.1)
pdf("figures/bias_GT_histones.pdf", useDingbats=FALSE, height=5, width=2.8)
par(mar=c(0,3.5,0,1))
plot(pos + jitter, biasdat$x,
   ylim=c(0,1), xlim=c(.5,7.5),
   pch=pch1, col="gray50", panel.first=grid(nx=NA, ny=NULL),
   bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
points(pos + jitter, biasdat$x, col=col, pch=pch2)
legend(x="topright", pch=c(1,2), col="grey50", text.col="grey40",
   legend=c("Lin. Amp.", "6x PCR"), bg="white", box.col="white",
   cex=.85, pt.cex=1.1, inset=.04)
axis(side=2, col.axis="grey50", col="grey50")
dev.off()
