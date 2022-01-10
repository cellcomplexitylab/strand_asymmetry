mut = read.table("mutations.txt")

# 52 C's, of which 19 CG's.
ratio = (mut$V2/19) / (mut$V3/33)

set.seed(123)
jitter = rnorm(n=length(ratio), sd=.1)
pdf("figures/C_to_T_transitions.pdf", useDingbats=FALSE, height=5, width=2)
par(mar=c(1.5,3.5,0,1))
plot(jitter, ratio, ylim=c(0,12), xlim=c(-1,1),
   pch=2, col="gray50", panel.first=grid(nx=NA, ny=NULL),
   bty="n", xaxt="n", yaxt="n", xlab="",
   line=2.2, col.lab="gray25", ylab="CG>TG / C>T ratio")
points(jitter, ratio, col="#EB7D5B80", pch=17)
axis(side=2, col.axis="grey50", col="grey50")
title(xlab="Replicates", line=0)
dev.off()
