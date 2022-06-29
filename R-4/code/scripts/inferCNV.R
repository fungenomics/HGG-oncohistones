# Author: Selin Jessa <selin.jessa@mail.mcgill.ca>
# Date: 2020-10-20
#
# Purpose:
# A script to run inferCNV on a sample, given the Seurat object, normal reference,
# and gene annotation. This script produces a custom heatmap where the cells are
# labeled by the transcriptional cluster, # counts, and expressioon of some
# genes which mark common normal cell types.
#
# Usage:
# This should be run in the scRNAseq/data/pipeline/10X/sample/inferCNV folders,
# and expects the following config file to be present:
#
# - info.experiment.tsv: TSV indicating parameters used in the analysis, with two columns, "param" and "value"
#   NOTE: currently, missing fields from the below are not handled, so a complete info.experiment.tsv
#   must be supplied, with the defaults below if desired.
#
#       cluster_by_group      Logical, whether to cluster the observations (tumor sample) within the groups
#                             defined in the annotations. NOTE: This is different than the clustering/normalization
#                             of different reference group cell types, which is handled automatically
#                             by the script. Setting to TRUE will basically allow visualization of CNV profiles within
#                             transcriptional clusters, but will not allow recovery of new subclusters based on CNVs.
#                             Default: FALSE
#
#       window_length         Numeric, must be an odd number, length of genomic windows
#                             used for CNV calling, see "window_length" argument in infercnv::run for reference.
# 				                  	Default: 101
#
#       expression_threshold  Numeric, use 0.1 for 10X and 1 for Smartseq2, see
#                             "cutoff" argument in infercnv::run for reference. Default: 0.1
#
#       analysis_mode         String, allowed values: "samples", "subclusters", "cells",
#                             specifies the level in inferCNV,
#                             see "analysis_mode" argument in infercnv::run for reference.
#                             Using "subclusters" increases running time. Default: samples
#
#       gene_annotation       String, path to annotation file giving the positions of
#                             each gene in each chromosome. InferCNV will restrict
#                             the analysis to genes present in this file.
#
#       reference             String, path to a Seurat object containing the normal
#                             reference dataset.
#
#       reference_key         String, a short codeword to summarize the reference used,
#                             will be printed in the name of the output directory to
#                             allow for runs with different references.
#      
#       HMM			            	String, allowed values: "none", "i3", "i6". Use "none' to
#                             disable running of an HMM to predict gross CNVs. To run the HMM,
#                             use "i3" or "i6" to specify which type to use; see "HMM_type" argument
#                             in infercnv::run for reference, and inferCNV wiki for more details.
#                             Note that this substantially increases running time. Default: none.
#
#

# Config -----------------------------------------------------------------------

# load libraries
library(here)
library(readr)
library(tibble)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(glue)
library(ggplot2)
library(Seurat)
library(infercnv)

set.seed(100)

here::here()

message("@ config...")

# convert to a named list where names are parameters and elements are values
info_experiment_tsv <- read_tsv("info.experiment.tsv")
info_experiment <- as.list(tibble::deframe(info_experiment_tsv[, c(1, 2)]))
info_experiment$cluster_by_groups    <- as.logical(info_experiment$cluster_by_groups)
info_experiment$intermediates        <- as.logical(info_experiment$intermediates)
info_experiment$window_length        <- as.numeric(info_experiment$window_length)
info_experiment$expression_threshold <- as.numeric(info_experiment$expression_threshold)

(out <- glue("window{info_experiment$window_length}_exp{info_experiment$expression_threshold}_ref{info_experiment$reference_key}_HMM{info_experiment$HMM}"))

dir.create(out)

# copy the config for reproducibility, and save the last commit in the repository for version control
info_experiment_tsv %>% add_case(param = "git_last_commit",
					   value = git2r::last_commit()$sha %>% stringr::str_sub(1, 7)) %>%
	write_tsv(glue("{out}/info.experiment.copy.tsv"))

# Prep data --------------------------------------------------------------------
message("@ loading data...")

# load the normal reference dataset
reference         <- readRDS(info_experiment$reference)
Idents(reference) <- reference@meta.data$orig.ident

# load the sample dataset
load("../seurat.Rda")

# merge into one object
all_data <- merge(reference, seurat)

# write the cell annotations to file
annot <- Idents(all_data)
if (!file.exists("annotation.tsv")) write.table(annot,
            file      = "annotation.tsv",
            sep       = '\t',
            row.names = TRUE,
            col.names = FALSE,
            quote     = FALSE)

# write the expression matrix to file
expression <- as.data.frame(GetAssayData(all_data, slot = "counts"))
if (!file.exists("expression_matrix.tsv.gz")) write.table(expression,
            file          = gzfile("expression_matrix.tsv.gz"),
            sep           = '\t',
            row.names     = TRUE,
            col.names     = TRUE,
            quote         = FALSE)

# create the inferCNV object
infercnv_obj = infercnv::CreateInfercnvObject(
  raw_counts_matrix = "expression_matrix.tsv.gz",
  annotations_file  = "annotation.tsv",
  delim             = "\t",
  gene_order_file   = info_experiment$gene_annotation,
  ref_group_names   = levels(Idents(reference)))



# Run inferCNV -----------------------------------------------------------------
message("@ running inferCNV...")
infercnv_obj_default <- infercnv::run(
  infercnv_obj,
  cutoff                = info_experiment$expression_threshold,
  window_length         = info_experiment$window_length,
  analysis_mode         = info_experiment$analysis_mode,
  out_dir               = out,
  HMM                   = ifelse(info_experiment$HMM == "none", FALSE, TRUE),
  HMM_type              = ifelse(info_experiment$HMM == "none", "i3", info_experiment$HMM), # if HMM == "none", this won't have an effect
  cluster_by_groups     = info_experiment$cluster_by_groups,
  tumor_subcluster_pval = 0.1, # default
  sd_amplifier          = TRUE,
  num_threads           = 4)



# Update Seurat object ---------------------------------------------------------
message("@ post-processing...")

# add CNV results to Seurat object
# all_data <- add_to_seurat(all_data, infercnv_output_path = out)

# Generate custom plots --------------------------------------------------------

gene_annotation <- read_tsv(info_experiment$gene_annotation, col_names = c("gene", "contig", "start", "end"))

#' Plot CNV profile based on inferCNV
#'
#' This was adapted/refactored from code from Samantha Worme's function plotSampleCNVs()
#'
#' @param out Character, name of the output directory
plot_profile <- function(out,
                         cluster_rows     = TRUE,
                         return_hm        = FALSE) {

  # prep and loading -----------------------------------------------------------
  message("@ custom heatmap: preparing data...")

  # load gene x cell matrix
  obs <- read.table(glue("{out}/infercnv.observations.txt"))

  # set NA values to 1 (i.e. no CNV)
  obs[is.na(obs)] <- 1

  # transpose to get cell x gene matrix
  obs_t <- t(obs)

  # order genes based on gene annotation
  contigs <- gene_annotation %>% filter(gene %in% colnames(obs_t))
  obs_t   <- obs_t[, order(match(colnames(obs_t), contigs$gene))]

  message("@ custom heatmap: producing annotations...")
  # produce the column (contig) annotation -------------------------------------
  # 1. for each contig (chromosome), get the index of the last gene, which will be
  # used to put gaps at those points in the heatmap
  indices <- contigs %>%
    group_by(contig) %>%
    mutate(last = which(gene %in% last(gene))) %>%
    filter(gene %in% last(gene)) %>%
    pull(last)

  # 2. now that we have the index of the last gene in each contig, we need to sum
  # these up to get the index of the last gene in each contig, relative
  # to the whole gene annotation
  last_gene_idx <- cumsum(indices)

  # select only the contig column (= chromosome number) for annotating the heatmap
  contigs <- contigs %>%
    select(contig, gene) %>%
    tibble::column_to_rownames(var = "gene")

  # produce row annotations ----------------------------------------------------
  # we'll colour by transcriptional cluster and # counts (i.e. # UMIs for 10X data),
  # as well as genes marking a few common normal populations:
  # DOCK8 (microglia/macrophages), VWF (endothelial), KCNJ8 (pericytes), MBP (myelinating oligo)
  row_annotations <- data.frame("cell"       = colnames(seurat),
                                "nCount"     = as.numeric(seurat$nCount_RNA),
                                "mito"       = as.numeric(seurat$percent.mito),
                                "cluster"    = as.character(Idents(seurat)),
                                stringsAsFactors = FALSE)
  
  # TODO: Make this an argument, where the info_experiment specifies genes to plot
  # and the counts are added here in a loop
  if ("DOCK8" %in% rownames(seurat)) row_annotations$DOCK8_expr <- as.numeric(seurat[["RNA"]]["DOCK8", ]) else warning("DOCK8 not detected")
  if ("VWF"   %in% rownames(seurat)) row_annotations$VWF_expr   <- as.numeric(seurat[["RNA"]]["VWF",   ]) else warning("VWF not detected")
  if ("KCNJ8" %in% rownames(seurat)) row_annotations$KCNJ8_expr <- as.numeric(seurat[["RNA"]]["KCNJ8", ]) else warning("KCNJ8 not detected")
  if ("MBP"   %in% rownames(seurat)) row_annotations$MBP_expr   <- as.numeric(seurat[["RNA"]]["MBP",   ]) else warning("MBP not detected")
  
  row_annotations <- row_annotations %>%
    # infercnv puts a "." in place of a "-" in the cell names
    mutate(cell = gsub("-", ".", cell, fixed = TRUE)) %>% 
    column_to_rownames(var = "cell") %>% 
    .[rownames(obs_t), ]

  # create palettes ------------------------------------------------------------

  # for contigs, we'll use a repeating gray palette:
  palette_contig <- rep(c("gray90", "gray70"), length.out = length(unique(contigs$contig)))
  names(palette_contig) <- unique(contigs$contig)

  # for clusters, we'll take the palette from the Seurat object
  palette_clusters <- seurat@misc$colours
  
  # for expression of specific genes, we'll use a red-gray palette
  expr_palette <- grDevices::colorRampPalette(c("gray90", "#E09797", "red"))(n = 200)

  anno_palettes = list(contig     = palette_contig,
                       cluster    = palette_clusters,
                       nCount     = viridisLite::magma(100),
                       mito       = viridisLite::magma(100),
                       DOCK8_expr = expr_palette,
                       VWF_expr   = expr_palette,
                       KCNJ8_expr = expr_palette,
                       MBP_expr   = expr_palette)

  # for the heatmap itself, generate a red-blue heatmap (this returns a function)
  hm_palette_generator <- infercnv:::color.palette(c("darkblue", "white", "darkred"), c(2, 2))

  message("@ custom heatmap: generating heatmap...")

  # label chromosomes by repurposing the column labels: generate a vector which will 
  # mostly be empty, but for one index in each contig, we'll put the chromosome labels
  # for the first few chromosomes, shift the label slightly to the left so it's more centered...
  labels_col <- rep("", ncol(obs_t))
  labels_col[c(last_gene_idx[1:7] - 200,
               last_gene_idx[8:length(last_gene_idx)])] <- unique(contigs$contig)

  # generate the heatmap -------------------------------------------------------
  hm <- pheatmap(obs_t,
                 scale             = "none",
                 border_color      = NA,
                 cluster_rows      = cluster_rows,
                 cluster_cols      = FALSE,
                 clustering_method = "ward.D",
                 labels_col        = labels_col,
                 show_colnames     = TRUE,
                 show_rownames     = FALSE,
                 color             = hm_palette_generator(100),
                 annotation_colors = anno_palettes,
                 annotation_col    = contigs,
                 annotation_row    = row_annotations,
                 gaps_col          = last_gene_idx,
                 treeheight_row    = 70, # make the dendrogram a little bit taller
                 filename          = glue("{out}/infercnv_custom.png"),
                 main              = glue("inferCNV - ref:{info_experiment$reference_key}"),
                 annotation_legend = FALSE,
                 fontsize_col      = 7,
                 width             = 11,
                 height            = 8)

  if (return_hm) return(hm)

}

plot_profile(out)

message("@ done.")

sessionInfo()
