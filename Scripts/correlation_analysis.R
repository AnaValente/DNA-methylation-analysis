#!/usr/bin/env Rscript

library("corrgram")
library("ggplot2")
library("RColorBrewer")
library("factoextra")
library("ggrepel")

args = commandArgs(trailingOnly = TRUE)

correlation = function(df) {

  df = read.csv(df, sep = "\t")
  df = df[, 4:ncol(df)] # remove string collumns
  png("correlation_analysis.png", width = 1500, height = 1400, res = 180)
  corrgram(df, upper.panel = panel.cor)
  dev.off()

  # principal component analysis
  pca = prcomp(df, scale = TRUE)

  summ = summary(pca)

  scores = as.data.frame(pca$rotation)
  scores = scores[1:3] # choose only the first principal components

  colours = c("#f8766D", "#f8766D", "#7cae00", "#7cae00", "#00bfc4", "#00bfc4", "#C77CFF", "#C77CFF", "#00B0F6", "#00B0F6", "#D89000", "#D89000")

  ggplot(scores, aes(PC1, PC2)) + geom_point(colour = head(colours, nrow(scores)), size = 2) +
    geom_text_repel(aes(label = rownames(scores)), size = 3) + 
    labs(x = sprintf("PC1 (%s%% of the variance explained)", round(summ$importance[2, 1], digits = 2) * 100), y = sprintf("PC2 (%s%% of the variance explained)", round(summ$importance[2, 2], digits = 2) * 100)) + 
    theme_linedraw() + theme(axis.line = element_line(color = "black"), plot.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) # plot of the PCA scores of PC1 and PC2
  ggsave("PCA_analysis_PC1_PC2.png", width = 5, height = 4, dpi = 400)

  ggplot(scores, aes(PC2, PC3)) + geom_point(colour = head(colours, nrow(scores)), size = 2) +
    geom_text_repel(aes(label = rownames(scores)), size = 3) + 
    labs(x = sprintf("PC2 (%s%% of the variance explained)", round(summ$importance[2, 3], digits = 2) * 100), y=sprintf("PC3 (%s%% of the variance explained)", round(summ$importance[2, 4], digits = 2) * 100)) + 
    theme_linedraw() + theme(axis.line = element_line(color = "black"), plot.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank())# plot of the PCA scores of PC2 and PC3
  ggsave("PCA_analysis_PC2_PC3.png", width = 5, height = 4, dpi = 400)

  fviz_eig(pca, addlabels = TRUE, ylim = c(0, 100)) + theme_bw() + theme(plot.background = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + scale_y_continuous(expand = c(0, 0)) # scree plot
  ggsave("PCA_scree_plot.png", width = 6, height = 5, dpi = 400)

}

correlation(args[1])
