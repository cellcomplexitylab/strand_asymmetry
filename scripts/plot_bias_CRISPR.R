biasdat = read.delim("bias_by_CRISPR_without_mapping.txt")

POS1 = c(1, 3.5)
names(POS1) = c("gRNA1", "gRNA2")
POS2 = c(0,1)
names(POS2) = c("24", "48")

pos = POS1[as.character(biasdat$gRNA)] + POS2[as.character(biasdat$time)]

set.seed(123)
jitter = rnorm(n=length(pos), sd=.1)
pdf("figures/bias_CRISPR.pdf", useDingbats=FALSE, height=5, width=2)
par(mar=c(0,3.5,0,1))
plot(pos + jitter, biasdat$x,
   ylim=c(0,1), xlim=c(.5,5),
   pch=2, col="gray50", panel.first=grid(nx=NA, ny=NULL),
   bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
points(pos + jitter, biasdat$x, col="#EB7D5B80", pch=17)
axis(side=2, col.axis="grey50", col="grey50")
dev.off()
