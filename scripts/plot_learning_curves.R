# The performance on the test set is in column V3.
a = read.table("learning_curves_full.txt")
b = read.table("learning_curves_null.txt")

am = tapply(X=a$V3, INDEX=a$V1, median)
aup = tapply(X=a$V3, INDEX=a$V1, quantile, .985)
adn = tapply(X=a$V3, INDEX=a$V1, quantile, .015)

bm = tapply(X=b$V3, INDEX=b$V1, median)[1:100]
bup = tapply(X=b$V3, INDEX=b$V1, quantile, .985)
bdn = tapply(X=b$V3, INDEX=b$V1, quantile, .015)

#ymax = max(aup, bup)
#ymin = min(adn, bdn)

tCOL=c("#9CD1C7B0", "#CAAE53B0")
COL=c("#4CA2A1", "#807C6F")

pdf("figures/learning_curves.pdf", useDingbats=FALSE, height=3, width=5)
par(mar=c(3,3,0.5,1))
plot(1:100, am, lwd=3, type="n", ylim=c(-3.3,-2.5),
   panel.first=grid(nx=NA, ny=NULL),
   bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
polygon(c(1:100,100:1), c(bup, rev(bdn)), lty=0, col=tCOL[2])
polygon(c(1:100,100:1), c(aup, rev(adn)), lty=0, col=tCOL[1])
lines(1:100, am, lwd=2, col=COL[1])
lines(1:100, bm, lwd=2, col=COL[2])
legend(x="topright", col=COL, text.col="grey40", lwd=2,
   legend=c("w/ chromatin", "w/o chromatin"), bg="white",
   box.col="white", cex=.85, pt.cex=1.1, inset=.04)
axis(side=2, col.axis="grey50", col="grey50", cex.axis=.8, padj=.8)
axis(side=1, col.axis="grey50", col="grey50", cex.axis=.8, padj=-.8)
title(xlab="Epoch", col="grey50", line=1.8, col.lab="grey30")
title(ylab="Test loss", col="grey50", line=1.8, col.lab="grey30")
dev.off()
