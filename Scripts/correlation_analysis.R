#!/usr/bin/env Rscript

library("corrgram")
library("ggplot2")
library("RColorBrewer")
library("factoextra")
library("ggrepel")
library("scales")
library("gtools")
library("data.table")
library("scales")

args <- commandArgs(trailingOnly = TRUE)

table_trimming <- function(df, n_samples, n_replicates) {

  var_df <- data.frame(random = runif(nrow(df)))

  # Get sequence with column numbers
  vect <- seq(4, ((as.numeric(n_samples) + 1) * as.numeric(n_replicates) + 3))
  fst <- 1
  lst <- as.numeric(n_replicates)

  if (fst == lst) {
    lst <- lst + 1
  }

  if ((ncol(df[, -c(1,2,3)]) %% as.numeric(n_replicates)) != 0) {
    fst <- 2
    lst <- as.numeric(n_replicates) + 1
    var_df[, paste0("sample", 1)] <- df[, 4]

    for (sample in 2:(as.numeric(n_samples) + 1)) {
      mean_rep <- rowMeans(df[, vect[fst:lst]]) # Row mean of column replicates
      var_df[, paste0("sample", sample)] <- mean_rep # Add mean of replicate to the variance dataframe
      fst <- fst + as.numeric(n_replicates)
      lst <- lst + as.numeric(n_replicates)
    }
  } else {

    if (as.numeric(n_replicates) == 1) {
      mean_rep <- rowMeans(df[, vect[fst:lst]])
      var_df[, paste0("sample")] <- mean_rep

    } else {

      for (sample in 1:(as.numeric(n_samples) + 1)) {
        mean_rep <- rowMeans(df[, vect[fst:lst]])
        var_df[, paste0("sample", sample)] <- mean_rep
        fst <- fst + as.numeric(n_replicates)
        lst <- lst + as.numeric(n_replicates)
      }
    }
  }

  var_df <- var_df[-1]

  data_matx <- data.matrix(var_df, rownames.force = NA)

  variance <- apply(data_matx, 1, var) # Calculate variance

  var_df[, "var"] <- variance

  var_df <- var_df[order(var_df$var, decreasing = TRUE), ]

  # Row names of the 10% most variable CpGs
  var_names <- rownames(var_df[0:(nrow(var_df) * 0.1), ])

  df <- df[as.numeric(var_names), ] # Select the 10% most variable CpGs

  return(df)
}

correlation <- function(df) {

  df <- df[, 4:(ncol(df))] # Remove string columns

  png("correlation_analysis.png", width = 1500, height = 1400, res = 180)
  corrgram(df, upper.panel = panel.cor) # Correlation matrix
  dev.off()

}

pca <- function(df, n_samples, n_replicates, control) {
  # Principal component analysis
  df <- df[, 4:ncol(df)]
  pca <- prcomp(df, scale = TRUE, center = TRUE)

  summ <- summary(pca)

  scores <- as.data.frame(pca$rotation)

  if (n_replicates <= 1 && n_samples <= 2) {
    scores <- scores[1:2] # Choose only the first two principal components
  }else {
    scores <- scores[1:3] # Choose only the first three principal components
  }

  colors <- c()
  if (n_replicates > 1) {
    for (i in hue_pal()(as.integer(n_samples) + 1)) {
      colors <-  c(colors, rep(i, n_replicates))
    }
  }else {
    for (i in hue_pal()(as.integer(n_samples) + 1)) {
      colors <-  c(colors, i)
    }
  }

  if (!(paste(control,"_rep2",sep = "") %in% colnames(df)) && any(grepl("rep2", colnames(df)))) {
    colors <- colors[-1]
  }else {
    colors <- colors
  }

  # Plot PCA
  ggplot(scores, aes(PC1, PC2)) + geom_point(colour = colors, size = 2) +
    geom_text_repel(aes(label = rownames(scores)), size = 3, max.overlaps=Inf) + 
    labs(x = sprintf("PC1 (%s%% of the variance explained)", round(summ$importance[2, 1], digits = 2) * 100),
    y = sprintf("PC2 (%s%% of the variance explained)", round(summ$importance[2, 2], digits = 2) * 100)) +
    theme_linedraw() + theme(axis.line = element_line(color = "black"), plot.background = element_blank(),
    panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    xlim(min(scores$PC1)-0.01, max(scores$PC1)+0.01) # Plot of the PCA scores of PC1 and PC2
  ggsave("PCA_analysis_PC1_PC2.png", width = 5, height = 4, dpi = 400)

  if (n_replicates > 1 && n_samples > 1) {
    ggplot(scores, aes(PC2, PC3)) + geom_point(colour = colors, size = 2) +
      geom_text_repel(aes(label = rownames(scores)), size = 3, max.overlaps=Inf) + 
      labs(x = sprintf("PC2 (%s%% of the variance explained)", round(summ$importance[2, 2], digits = 2) * 100),
      y = sprintf("PC3 (%s%% of the variance explained)", round(summ$importance[2, 3], digits = 2) * 100)) +
      theme_linedraw() + theme(axis.line = element_line(color = "black"), plot.background = element_blank(),
      panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
      xlim(min(scores$PC2)-0.01, max(scores$PC2)+0.01) # Plot of the PCA scores of PC2 and PC3
    ggsave("PCA_analysis_PC2_PC3.png", width = 5, height = 4, dpi = 400)
  }

  fviz_eig(pca, addlabels = TRUE, ylim = c(0, 100)) + theme_bw() + theme(plot.background = element_blank(), 
  panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + scale_y_continuous(expand = c(0, 0)) # scree plot
  ggsave("PCA_scree_plot.png", width = 6, height = 5, dpi = 400)

}

df <- read.csv(args[1], sep = "\t")
correlation(df)
final_df <- table_trimming(df, args[2], args[3])
pca(final_df, args[2], args[3], args[4])