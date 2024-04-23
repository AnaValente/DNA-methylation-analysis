#!/usr/bin/env Rscript

library("corrgram")
library("ggplot2")
library("RColorBrewer")
library("factoextra")
library("ggrepel")
library("scales")

args = commandArgs(trailingOnly = TRUE)

correlation = function(df) {

  df = df[, 4:ncol(df)] # Remove string columns
  png("correlation_analysis.png", width = 1500, height = 1400, res = 180)
  corrgram(df, upper.panel = panel.cor) # Correlation matrix
  dev.off()

}

PCA = function(df, n_samples, n_replicates, control) {
  # Principal component analysis
  df = df[, 4:ncol(df)]
  pca = prcomp(df, scale = TRUE)

  summ = summary(pca)

  scores = as.data.frame(pca$rotation)
  
  if (n_replicates <= 1 & n_samples <= 2) {
    scores = scores[1:2] # Choose only the first two principal components
  }else{
    scores = scores[1:3] # Choose only the first three principal components
  }

  colors = c()
  if (n_replicates > 1){
    for (i in hue_pal()(as.integer(n_samples)+1)) {
      colors =  c(colors,rep(i,n_replicates))
    }
  }else{
    for (i in hue_pal()(as.integer(n_samples)+1)) {
      colors =  c(colors,i)
  }}
  
  if (!(paste(control,'_rep2',sep = '') %in% colnames(df)) & any(grepl("rep2", colnames(df)))) {
    colors = colors[-1]
  }else{
    colors = colors
  }

  # Plot PCA
  ggplot(scores, aes(PC1, PC2)) + geom_point(colour = colors, size = 2) +
    geom_text_repel(aes(label = rownames(scores)), size = 3) + 
    labs(x = sprintf("PC1 (%s%% of the variance explained)", round(summ$importance[2, 1], digits = 2) * 100),
    y = sprintf("PC2 (%s%% of the variance explained)", round(summ$importance[2, 2], digits = 2) * 100)) +
    theme_linedraw() + theme(axis.line = element_line(color = "black"), plot.background = element_blank(),
    panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    xlim(min(scores$PC1)-0.01, max(scores$PC1)+0.01) # Plot of the PCA scores of PC1 and PC2
  ggsave("PCA_analysis_PC1_PC2.png", width = 5, height = 4, dpi = 400)

  if (n_replicates > 1 & n_samples > 2){
    ggplot(scores, aes(PC2, PC3)) + geom_point(colour = colors, size = 2) +
      geom_text_repel(aes(label = rownames(scores)), size = 3) + 
      labs(x = sprintf("PC2 (%s%% of the variance explained)", round(summ$importance[2, 2], digits = 2) * 100),
      y=sprintf("PC3 (%s%% of the variance explained)", round(summ$importance[2, 3], digits = 2) * 100)) +
      theme_linedraw() + theme(axis.line = element_line(color = "black"), plot.background = element_blank(),
      panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
      xlim(min(scores$PC2)-0.01, max(scores$PC3)+0.01) # Plot of the PCA scores of PC2 and PC3
    ggsave("PCA_analysis_PC2_PC3.png", width = 5, height = 4, dpi = 400)
    }

  fviz_eig(pca, addlabels = TRUE, ylim = c(0, 100)) + theme_bw() + theme(plot.background = element_blank(), 
  panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + scale_y_continuous(expand = c(0, 0)) # Scree plot
  ggsave("PCA_scree_plot.png", width = 6, height = 5, dpi = 400)

}

df = read.csv(args[1], sep = "\t")
correlation(df)
PCA(df,args[2],args[3],args[4])
