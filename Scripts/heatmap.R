#!/usr/bin/env Rscript

library(gplots)

args = commandArgs(trailingOnly = TRUE)

heatmap = function(df, filename) {

  # merge CpG coordinates
  df[, "CpG_sites"] = paste(df$chr, ":", df$start, "-", df$end, sep = "")

  # change DataFrame names to CpG coordinates
  row.names(df) = df[, "CpG_sites"]
  df = df[, 4:(ncol(df)-1)]
  
  # plot histogram
  pdf(filename, width = 13, height = 12)
  heatmap.2(as.matrix(df), srtCol = 360, adjCol = c(0.5, 0), cexCol = 0.86, cexRow = 0.7, margins = c(3, 12.5), hclustfun = hclust, trace = "none", col = colorpanel(100, "red3", "gray", "blue4"))
  dev.off()

}

df_methylation = read.csv(args[1], sep = "\t")

heatmap(df_methylation, "Methylation_heatmap.pdf")
