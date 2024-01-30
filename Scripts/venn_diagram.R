#!/usr/bin/env Rscript

# libraries
library("ggVennDiagram")
library(RColorBrewer)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)

diagram = function() {

  # list all files with *.txt prefix
  gene_files = list.files(pattern = "*.txt", full.names = TRUE)

  # read list files
  genes_matrix = lapply(gene_files, function(i) {as.character(as.matrix(read.csv(i, header = FALSE)))})

  # trim sample names
  sample_names = lapply(gene_files, function(i) {sub("_.*", "", as.character(i))})

  sample_names = lapply(sample_names, function(i) {sub("./", "", as.character(i))})

  # plot veen diagram
  diagram = ggVennDiagram(genes_matrix, label = "count", label_alpha = 0, label_size = 6, set_size = 5, category.names = sample_names, edge_size = 0) + scale_fill_gradient(low = "#F4FAFE", high = "#3f648e")

  ggsave("venn_diagram.png", diagram, scale = 2, width = 15, height = 10, dpi = 200, units = c("cm"), bg = "White")
}


diagram()