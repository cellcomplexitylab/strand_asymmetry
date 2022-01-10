f = read.delim("features_ctrl_GT_barcodes.txt")

means_w__feature = rep(NA, ncol(f))
means_wo_feature = rep(NA, ncol(f))
c_low = rep(NA, ncol(f))
c_hi = rep(NA, ncol(f))
n_w_feature = rep(NA, ncol(f))
pvals = rep(NA, ncol(f))
for (i in 68:79) {
   means_w__feature[i] = mean(f$meCpG[f[,i] == 1])
   means_wo_feature[i] = mean(f$meCpG[f[,i] == 0])
   n_w_feature[i] = sum(f[,i] == 1)
   test = t.test(f$meCpG[f[,i] == 0], f$meCpG[f[,i] == 1], conf.level=.99)
   pvals[i] = test$p.value
   c_low[i] = test$conf.int[1]
   c_hi[i] = test$conf.int[2]
}

features = .2021 - (means_wo_feature[!is.na(means_wo_feature)] - means_w__feature[!is.na(means_w__feature)])
n_w_feature = n_w_feature[!is.na(n_w_feature)]
mid = 100 * (c_hi[!is.na(c_hi)] - c_low[!is.na(c_low)]) / 2
sem = sqrt(features * (1-features) / n_w_feature)

N = length(features)
print(pvals[!is.na(pvals)])

COL = c(
   "#491212", "#491212",
   "#831a19", "#831a19",
   "#c24229", "#c24229",
   "#dd8437", "#dd8437",
   "#eeb440", "#eeb440",
   "#d8cea7", "#d8cea7")


pdf("figures/meCpG_in_chromatin.pdf", useDingbats=FALSE, height=5, width=4)
par(mar=c(5,3.5,.5,.5))
plot(100 * features, ylim=c(0,25), pch=19, col=COL,
   panel.first=grid(), bty="n", xaxt="n", yaxt="n", col.lab="grey25",
   line=2.2, ylab="Reporters with CG>TG errors (%)", xlab="")
abline(h=20.21, col="gray50")
segments(x0=1:N, y0=100*features - mid, y1=100*features + mid,
   col=COL, lwd=2)
#abline(h=.05 * 1:5, col="gray50", lty=2)
axis(side=2, col="grey50", cex.axis=.8, col.axis="grey25")
axis(side=1, col="grey50", cex.axis=.8, col.axis="grey25", at=1:12,
   labels=colnames(f)[68:79], las=3)
dev.off()
