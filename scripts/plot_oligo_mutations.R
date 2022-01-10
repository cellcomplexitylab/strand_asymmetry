freq = matrix(byrow=TRUE, ncol=2, c(
    0 / 52309 / 2,    0 / 44735 / 2, # A > C
   21 / 52309 / 2,   19 / 44735 / 2, # A > G
    3 / 52309 / 2,    2 / 44735 / 2, # A > T
               NA,               NA,
   70 / 52309 / 2,   47 / 44735 / 2, # G > A
    0 / 52309 / 2,    0 / 44735 / 2, # G > C
   85 / 52309 / 2,   73 / 44735 / 2, # G > T
               NA,               NA,
    2 / 52309 / 2,    1 / 44735 / 2, # T > A
   23 / 52309 / 2,   28 / 44735 / 2, # T > C
    2 / 52309 / 2,    4 / 44735 / 2, # T > G
               NA,               NA,
    2 / 52309 / 2,    1 / 44735 / 2, # C > A
    1 / 52309 / 2,    1 / 44735 / 2, # C > G
   13 / 52309    ,   13 / 44735    , # C > T
               NA,               NA,
    2 / 52309 / 2,    1 / 44735 / 2, # CG > AG
    1 / 52309 / 2,    1 / 44735 / 2, # CG > GG
   30 / 52309    ,   46 / 44735    , # CG > TG
               NA,               NA,
    4 / 72718    ,   10 / 71488    , # mCG > AG
    0 / 72718    ,    0 / 71488    , # mCG > GG
  169 / 72718    ,  183 / 71488      # mCG > TG
))

# BrBG palette.
COL=c("#8c510aff", "#8c510ae8", "#8c510adf", "#8c510ac8", "#8c510abf", "#8c510aa8", NA, NA,
      "#bf812dff", "#bf812de8", "#bf812ddf", "#bf812dc8", "#bf812dbf", "#bf812da8", NA, NA,
      "#dfc27dff", "#dfc27de8", "#dfc27ddf", "#dfc27dc8", "#dfc27dbf", "#dfc27da8", NA, NA,
      "#f6e8c3ff", "#f6e8c3e8", "#f6e8c3df", "#f6e8c3c8", "#f6e8c3bf", "#f6e8c3a8", NA, NA,
      "#c7eae5ff", "#c7eae5e8", "#c7eae5df", "#c7eae5c8", "#c7eae5bf", "#c7eae5a8", NA, NA,
      "#01665eff", "#01665ee8", "#01665edf", "#01665ec8", "#01665ebf", "#01665ea8")

av = mean(freq[1:12,], na.rm=TRUE)

pdf("figures/oligo_mutations.pdf", height=5, width=6, useDingbats=FALSE)
par(mar=c(4.4,3.5,.5,.5))
barplot(100 * t(freq), beside=TRUE, col=COL, border="gray50", space=c(0,0.3),
   bty="n", xaxt="n", yaxt="n", col.lab="grey25",
   line=2.2, ylab="Mutation frequency (%)")
abline(h=.05 * 1:5, col="gray50", lty=2)
axis(side=2, col="grey50", cex.axis=.8, col.axis="grey25")
axis(side=1, col="grey50", cex.axis=.8, col.axis="grey25",
     at=c(1.3,3.6,5.9) + rep(9.2*(0:5), each=3), las=3,
     labels=c("A>C", "A>G", "A>T",  "G>A", "G>C", "G>T", "T>A", "T>C", "T>G",
              "C>A", "C>G", "C>T", "CG>AG", "CG>GG", "CG>TG", "mCG>AG", "mCG>GG", "mCG>TG"))
   
dev.off()
