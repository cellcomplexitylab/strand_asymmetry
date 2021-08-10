biasdat = read.delim("bias_for_GT_by_genic.txt")

POS1 = c(1, 3.5, 6)
names(POS1) = c("gene/high", "gene/low", "intergenic")
POS2 = c(0,1)
names(POS2) = c("24", "48")
COL = c("#0B263980", "#59807880", "#E1DBC080")
names(COL) = c("gene/high", "gene/low", "intergenic")
SYM = c(1,2)
names(SYM) = c("LA", "6xPCR")

pos = POS1[as.character(biasdat$genic)] + POS2[as.character(biasdat$time)]
col = COL[as.character(biasdat$genic)]
pch1 = SYM[as.character(biasdat$method)]
pch2 = pch1 + 15

set.seed(123)
jitter = rnorm(n=length(pos), sd=.1)
pdf("figures/bias_GT_genic.pdf", useDingbats=FALSE, height=5, width=2.8)
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
