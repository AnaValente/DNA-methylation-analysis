#!/usr/bin/env Rscript

library(ggplot2)
library("annotatr")
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(forcats)
library(reshape2)
library(dplyr)
library(tidyr)

annotations <- function() {

  annots <- c("hg38_cpgs", "hg38_genes_1to5kb", "hg38_genes_promoters",
              "hg38_genes_cds", "hg38_genes_5UTRs", "hg38_genes_exons",
              "hg38_genes_firstexons", "hg38_genes_introns",
              "hg38_genes_intronexonboundaries",
              "hg38_genes_exonintronboundaries",
              "hg38_genes_3UTRs", "hg38_genes_intergenic", "hg38_basicgenes")

  # Build annotations for hg38 reference genome
  annotations <- build_annotations(genome = "hg38", annotations = annots)

  return(annotations)
}

change_annotation <- function(df) {
  # Replace annotations names for new annotation names
  annots_types_old <- c("hg38_cpg_islands", "hg38_genes_1to5kb",
                        "hg38_genes_3UTRs", "hg38_genes_5UTRs",
                        "hg38_genes_cds", "hg38_genes_exonintronboundaries",
                        "hg38_genes_exons", "hg38_genes_firstexons",
                        "hg38_genes_intergenic",
                        "hg38_genes_intronexonboundaries",
                        "hg38_genes_introns", "hg38_genes_promoters")

  annots_types_new <- c("CpG islands", "1to5kb", "3UTRs",
                        "5UTRs", "cds", "exon/intron boundaries",
                        "exons", "first exons", "intergenic",
                        "intron/exon boundaries", "introns", "promoters")

  for (n in seq(1, nrow(df))) {
    for (j in 1:12) {
      if (df[n, "Regions"] == annots_types_old[j]) {
        df[n, "Regions"] <- annots_types_new[j]
      }
    }
  }

  return(df)
}

genomic_regions_meth <- function(annotations) {

  # List all files with *_methylation.csv prefix
  annotatr_files <- list.files(pattern = "*__methylation.csv", full.names = TRUE)

  # Read files regions
  sample_regions <- lapply(annotatr_files, function(i) {read_regions(i, genome = "hg38", format = "bedgraph")})

  # Get samples annotations
  sample_annotation <- lapply(sample_regions, function(i) {annotate_regions(regions = i, annotations = annotations, ignore.strand = TRUE, quiet = FALSE)})

  # Summarize samples annotations
  annots_sum <- lapply(sample_annotation, function(i) {summarize_annotations(annotated_regions = i, quiet = TRUE)})

  for (i in seq(1, length(annots_sum))) {
    # Remove unwanted characters from sample name
    name <- sub("__.*", "", as.character(annotatr_files[[i]]))
    name  <- sub("./", "", name)

    colnames(annots_sum[[i]]) <- c("Regions", name)
  }

  # Merge list annotations into one DataFrame
  df_list <- Reduce(function(x, y, all = TRUE) merge(x, y, all = TRUE), annots_sum, accumulate = FALSE)

  # Replace annotations names for new annotation names
  df_list <- change_annotation(df_list)

  df_list <- as.data.frame(df_list)
  df_list <- na.omit(df_list)
  row.names(df_list) <- df_list$Regions
  df_list <- df_list[, 2:ncol(df_list), drop = FALSE] # Remove annotation names

  # Melt DataFrame columns into rows
  matrix <- melt(as.matrix(df_list), value.name = "Count", varnames = c("Annotations", "Nanomaterial"))

  # Plot matrix
  ggplot(matrix, aes(fill = Nanomaterial, y = Count, x = Annotations)) +
    geom_bar(stat = "identity", position = position_dodge(0.92), width = 0.9) +
    theme_light() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Nº of CpG sites") + scale_x_discrete() +
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  ggsave("genomic_regions_CpG.png", width = 20, height = 15, units = "cm")

  return(matrix)
}

genomic_regions_dmr <- function(annotations) {

  ids <- c("hg38_cpg_islands", "hg38_genes_1to5kb",
           "hg38_genes_3UTRs", "hg38_genes_5UTRs",
           "hg38_genes_cds", "hg38_genes_exonintronboundaries",
           "hg38_genes_exons", "hg38_genes_firstexons",
           "hg38_genes_intergenic",
           "hg38_genes_intronexonboundaries",
           "hg38_genes_introns", "hg38_genes_promoters")

  # List all files with *0.05.bedgraph prefix
  dmr_files <- list.files(pattern = "*0.05.bedgraph", full.names = TRUE)

  # Read all CSVs in list
  read_dmr <- lapply(dmr_files, function(i) {read.csv(i, sep = "\t", header = FALSE)})

  # Select the hypomethylated DMRs and hypermethylated DMRs
  hypo_dmr <- lapply(read_dmr, function(i) {i[i[, "V4"] <= 0, ]})
  hyper_dmr <- lapply(read_dmr, function(i) {i[i[, "V4"] >= 0, ]})

  # Genomic ranges of the hypomethylated and hypermethylated DMRs
  regions_hypo <- lapply(hypo_dmr, function(i) {GRanges(i[, 1], IRanges(i[, 2], i[, 3]))})
  regions_hyper <- lapply(hyper_dmr, function(i) {GRanges(i[, 1], IRanges(i[, 2], i[, 3]))})

  # Hypomethlylated and hypermethylated DMR annotations
  hypo_annotation <- lapply(regions_hypo, function(i) {annotate_regions(regions = i, annotations = annotations, ignore.strand = TRUE, quiet = FALSE)})
  hyper_annotation <- lapply(regions_hyper, function(i) {annotate_regions(regions = i, annotations = annotations, ignore.strand = TRUE, quiet = FALSE)})

  # Summarize hypomethylated annotations
  hypo_annsum <- lapply(hypo_annotation, function(i) {summarize_annotations(annotated_regions = i, quiet = TRUE)})
  hypo_col <- lapply(hypo_annsum, setNames, c("Regions", "Nº of DMR"))
  hypo_col <- lapply(hypo_col, function(i) {i %>% complete(Regions = ids, fill = list(`Nº of DMR` = 0))})
  hypo_col <- lapply(hypo_col, change_annotation)

  # Summarize hypermethylated annotations
  hyper_annsum <- lapply(hyper_annotation, function(i) {summarize_annotations(annotated_regions = i, quiet = TRUE)})
  hyper_col <- lapply(hyper_annsum, setNames, c("Regions", "Nº of DMR"))
  hyper_col <- lapply(hyper_col, function(i) {i %>% complete(Regions = ids, fill = list(`Nº of DMR` = 0))})
  hyper_col <- lapply(hyper_col, change_annotation)
  

  for (i in seq(1, length(dmr_files))) {
    name <- sub("./metilene_", "", as.character(dmr_files[[i]]))
    name  <- sub("__.*", "", name)
    # Add hypomethylation collumn
    hypo_col[[i]][, paste(name, "DMR", sep = " ")] <- rep("Hypomethylation", each = nrow(hypo_col[[i]]))
    # Add hypermethylation collumn
    hyper_col[[i]][, paste(name, "DMR", sep = " ")] <- rep("Hypermethylation", each = nrow(hyper_col[[i]]))
    sample_list <- list(hypo_col[[i]], hyper_col[[i]])
    # Bind the rows of the hypermethylated DMR dataframe and the hypomethylated DMR dataframe
    sample_list <- do.call("rbind", sample_list)

    fill_name <- colnames(sample_list)[3]

    ggplot(sample_list, aes(fill = .data[[fill_name]], y = `Nº of DMR`, x = Regions)) +
      geom_bar(stat = "identity", position = position_dodge(0.92),
               width = 0.9) +
      theme_light() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), legend.position = "top",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12)) +
      scale_y_continuous(expand = expansion(mult = c(0, .1)))

    ggsave(paste(name, "_genomic_regions_DMR.png", sep = ""), width = 20, height = 15, units = "cm")
  }
}

annotations <- annotations()
genomic_regions_meth(annotations)
genomic_regions_dmr(annotations)
