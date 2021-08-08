biasdat = read.delim("bias_by_experiments_without_mapping.txt")

POS1 = c(1, 3.5, 6, 8.5, 11)
names(POS1) = c("AG", "TG", "AC", "TC", "GT")
POS2 = c(0,1)
names(POS2) = c("24", "48")
COL = c("#44228880", "#6CA2EA80", "#B5D33D80", "#FED23F80", "#EB7D5B80")
names(COL) = c("AG", "TG", "AC", "TC", "GT")
SYM = c(1,2)
names(SYM) = c("LA", "6xPCR")

pos = POS1[as.character(biasdat$MM)] + POS2[as.character(biasdat$time)]
col = COL[as.character(biasdat$MM)]
pch1 = SYM[as.character(biasdat$method)]
pch2 = pch1 + 15

set.seed(123)
jitter = rnorm(n=length(pos), sd=.1)
pdf("figures/bias.pdf", useDingbats=FALSE, height=5, width=4)
par(mar=c(0,3.5,0,1))
plot(pos + jitter, biasdat$x,
   ylim=c(0,1), xlim=c(.5,12.5),
   pch=pch1, col="gray50", panel.first=grid(nx=NA, ny=NULL),
   bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
points(pos + jitter, biasdat$x, col=col, pch=pch2)
legend(x=9.5, y=.98, pch=c(1,2), col="grey50", text.col="grey40",
   legend=c("Lin. Amp.", "6x PCR"), bg="white", box.col="white",
   cex=.85, pt.cex=1.1, inset=.04)
axis(side=2, col.axis="grey50", col="grey50")
dev.off()
