#!/usr/bin/env Rscript

png("recreated_boxplot.png",height=1000,width=800)
combined <- read.csv("combined.txt")
boxplot(combined, log = 'y', las = 2, par(mar = c(7, 4.25, 0.20, 0.20)+ 0.1), col = c("#f2b57c", "#ef7d7c", "#b67eec", "#7f7cf1", "#7cb7f1", "#7af1ef"))
mtext("Number of predicted non-reference insertions per sample", 
      side = 2, line = 3.5, cex = .9)
dev.off()
