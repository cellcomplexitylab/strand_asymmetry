biasdat = read.table("stats.txt")

idx = biasdat$V7

num1 = tapply(INDEX=idx, X=biasdat$V3 * biasdat$V6, FUN=sum)
num2 = tapply(INDEX=idx, X=biasdat$V3^2 * biasdat$V6, FUN=sum)
den1 = tapply(INDEX=idx, X=biasdat$V6, FUN=sum)

COL = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B")

means = num1/den1
sds = sqrt(num2/den1 - means^2)
print(sds)


pdf("bias.pdf", useDingbats=FALSE)
par(mar=c(0,3.5,0,1))
plot(c(0,20), c(-.4,.4), type="n",
   bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
abline(h=seq(from=-.3, to=.3, by=.1), lwd=1, col="grey90")
# Error bars.
segments(x0=(1:20)-.416, y0=means-.5, x1=(1:20)-.416,
   y1=means-.5+sign(means-.5)*sds, lwd=2, col="grey40")
barplot(means -.5, col=rep(COL, each=4), add=TRUE,
   width=1/1.2,
   xaxt="n", yaxt="n")
axis(side=2, col.axis="grey25", col="grey50",
   at=c(-.4,-.2,0,.2,.4), labels=c(.1,.3,.5,.7,.9))
dev.off()


pos = c(rep(1,8),rep(2,6), rep(3.5,8),rep(4.5,8),
   rep(6,8),rep(7,8), rep(8.5,8),rep(9.5,8), rep(11,8),rep(12,8))
col = c(rep(1,14), rep(2,16), rep(3,16), rep(4,16), rep(5,16))
pch = c(1,2)[1 + (idx %% 2 == 0)]
tpch = c(16,17)[1 + (idx %% 2 == 0)]
TCOL = c("#44228880", "#6CA2EA80", "#B5D33D80", "#FED23F80", "#EB7D5B80")

set.seed(123)
jitter = rnorm(n=length(pos), sd=.1)
pdf("bias2.pdf", useDingbats=FALSE, height=5, width=4)
par(mar=c(0,3.5,0,1))
plot(pos + jitter, biasdat$V3,
   ylim=c(0,1), xlim=c(.5,12.5),
   pch=pch, col=COL[col], panel.first=grid(nx=NA, ny=NULL),
   bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
points(pos + jitter, biasdat$V3, col=TCOL[col], pch=tpch)
legend(x=9.5, y=.98, pch=c(1,2), col="grey50", text.col="grey40",
   legend=c("Lin. Amp.", "6x PCR"), bg="white", box.col="white",
   cex=.85, pt.cex=1.1, inset=.04)
axis(side=2, col.axis="grey50", col="grey50")
dev.off()

pdf("zo.pdf", useDingbats=FALSE, height=5, width=4)
par(mar=c(0,3.5,0,1))
plot(pos + jitter, 1-biasdat$V2,
   ylim=c(0,1), xlim=c(.5,12.5),
   pch=pch, col=COL[col], panel.first=grid(nx=NA, ny=NULL),
   bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
points(pos + jitter, 1-biasdat$V2, col=TCOL[col], pch=tpch)
legend(x=9.5, y=.98, pch=c(1,2), col="grey50", text.col="grey40",
   legend=c("Lin. Amp.", "6x PCR"), bg="white", box.col="white",
   cex=.85, pt.cex=1.1, inset=.04)
axis(side=2, col.axis="grey50", col="grey50")
dev.off()

zo = read.table("zo_48.txt")
x = zo[,2] / rowSums(zo)
jitter = rnorm(length(x), sd=.1) * rnorm(length(x), sd=.1) * 10

pdf("ex.pdf", useDingbats=FALSE, height=5, width=2)
par(mar=c(0,3.5,0,1))
TCOL = c("#44228820", "#6CA2EA20", "#B5D33D20", "#FED23F20", "#EB7D5B20")
plot(jitter, x, cex=.7,
   xlim=c(-.7,.7), 
   pch=1, col=TCOL[2], panel.first=grid(nx=NA, ny=NULL),
   bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
TCOL = c("#44228820", "#6CA2EA20", "#B5D33D20", "#FED23F20", "#EB7D5B20")
points(jitter, x, col=TCOL[2], pch=16, cex=.7)
axis(side=2, col.axis="grey50", col="grey50")
dev.off()
