MIdat = read.table("mutual_information.txt")

COL = rep(c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"), each=3)

pdf("figures/mutual_information.pdf", width=4, useDingbats=FALSE)
par(mar=c(1,3.5,.5,.5))
barplot(MIdat$V7, col=COL, ylim=c(0,1),
   bty="n", xaxt="n", yaxt="n", col.lab="grey25",
   line=2.2, ylab="Information between replicates (bits)")
axis(side=2, col="grey50", cex.axis=.8, col.axis="grey25",
   at=(0:5)/5)
dev.off()
