# Author: Selin Jessa <selin.jessa@mail.mcgill.ca>
# Date: 2020-11-05
#
# Purpose:
# This script runs correlation-based cell type prediction on a sample of interest, given
# the mean gene expression per cluster for a reference. This code is adapted
# from  Marie Coutelier, specifically the function label_correlation(). The results
# (best cell type and the correlation value) will be saved in the Seurat object,
# and some diagnostic plots will be generated in the output directory.
#
# Usage:
#   
#   Rscript predict_celltype_cor.R
#
# This script expects the following config files to be present:
# 
# - info.experiment.tsv
#
#	    reference_mat	  String, path to reference dataset expression matrix, an 
#                       .Rds file of mean gene experssion per cell type (gene x cell type)
#                       The rownames will be coerced to uppercase, no need to do so
#                       in advance.
#
#     	annotation		  String, path to cell type annotation file for the reference.
#                       Expects at least two two columns, "Cell_type" and "Color"
#

# Set up -----------------------------------------------------------------------

# load libraries
library(here)
library(readr)
library(glue)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(Seurat)

here::here()

source(here("include/style.R"))

# load config
message("@ config...")
info_experiment_tsv <- read_tsv("info.experiment.tsv")
info_experiment     <- as.list(tibble::deframe(info_experiment_tsv[, c(1, 2)])) # make it a list
info_experiment$threshold_common_genes <- as.numeric(info_experiment$threshold_common_genes)

if (length(info_experiment$threshold_common_genes) == 0) info_experiment$threshold_common_genes <- 0.5

(out <- glue("ref.{info_experiment$reference_key}"))
dir.create(out, showWarnings = FALSE)

file.copy("info.experiment.tsv", file.path(out, "info.experiment.tsv"))

# load data
load("../seurat.Rda")

# Correlations -----------------------------------------------------------------

# load reference
ref_mean_expression <- readRDS(info_experiment$reference_mat)

#' adapted from Marie Coutelier
#' 
#' @param test_expr_mat Matrix, gene x cell
#' @param ref_expr_mat Matrix, gene x cell
#' @param threshold_common_genes numeric, percentage of test dataset genes required
#' to be in the reference dataset in order to proceed
label_correlation <- function(test_expr_mat,
                              ref_expr_mat,
                              threshold_common_genes = 0.5) {
  
  rownames(test_expr_mat) <- toupper(rownames(test_expr_mat))
  rownames(ref_expr_mat) <- toupper(rownames(ref_expr_mat))
  
  # Testing how many genes are in common and stopping if not enough
  common_genes <- intersect(rownames(test_expr_mat), rownames(ref_expr_mat))
  prop_common <- length(common_genes) / length(rownames(test_expr_mat))
  message("@@ ", round(prop_common*100, digits = 2), "% of test dataset genes are in the reference dataset")
  
  if (prop_common < threshold_common_genes) stop("Proportion of common genes below threshold.")
  
  # Reducing matrices to common subset
  mat1 <- as.matrix(test_expr_mat[common_genes, ])
  mat2 <- as.matrix(ref_expr_mat[common_genes, ])
  
  # sanity check
  nrow(mat1) == nrow(mat2)
  
  # Computing correlations
  cor_matrix <- cor(mat1, mat2, method = "spearman", use = "complete.obs")
  
  # Getting the best one
  cor_label <- as.data.frame(cor_matrix) %>%
    mutate("cell" = rownames(cor_matrix)) %>%
    gather("celltype", "correlation", -cell) %>%
    group_by(cell) %>%
    top_n(1, correlation) %>%
    dplyr::select(cell, celltype, correlation) %>% 
    arrange(cell)
  
  # Returning the results
  return(cor_label)
  
}

message("@ Computing correlations...")
labels <- label_correlation(test_expr_mat = GetAssayData(seurat),
                            ref_expr_mat = ref_mean_expression,
                            threshold_common_genes = info_experiment$threshold_common_genes) 

write_tsv(labels, glue("{out}/correlation.predicted_labels.tsv"))

# Post-processing --------------------------------------------------------------

anno    <- read_tsv(info_experiment$annotation)
palette <- anno %>% select(Cell_type, Color) %>% tibble::deframe()

message("@ Generating outputs...")

# add into Seurat object
seurat <- AddMetaData(seurat, metadata = labels$celltype,    col.name = glue("COR_ref.{info_experiment$reference_key}"))
seurat <- AddMetaData(seurat, metadata = labels$correlation, col.name = glue("COR_ref.{info_experiment$reference_key}_score"))

# save 
save(seurat, file = "../seurat.Rda")

# generate tSNE/UMAP plots by projection
if ("tsne" %in% names(seurat@reductions)) {

    p1 <- TSNEPlot(seurat,
                   group.by = glue("COR_ref.{info_experiment$reference_key}"), cols = palette) +
        NoLegend() + ggtitle("Cell type")
    
    p2 <- FeaturePlot(seurat, reduction = "tsne",
                      glue("COR_ref.{info_experiment$reference_key}_score")) +
        scale_colour_gradientn(colours = brewer.pal(9, "YlOrRd")) +
        ggtitle("Correlation")
    
    p3 <- plot_grid(p1, p2, rel_widths = c(0.45, 0.55))
    ggsave(plot = p3, glue("{out}/tsne.png"), width = 11, height = 5)
        
}


p1 <- UMAPPlot(seurat,
               group.by = glue("COR_ref.{info_experiment$reference_key}"), cols = palette) +
  NoLegend() + ggtitle("Cell type")

p2 <- FeaturePlot(seurat, reduction = "umap",
                  glue("COR_ref.{info_experiment$reference_key}_score")) +
  scale_colour_gradientn(colours = brewer.pal(9, "YlOrRd")) +
  ggtitle("Correlation")

p3 <- plot_grid(p1, p2, rel_widths = c(0.45, 0.55))
ggsave(plot = p3, glue("{out}/umap.png"), width = 11, height = 5)

# get data to generate summary plots
df <- data.frame(Cluster = as.character(Idents(seurat)),
                 Cell_type = labels$celltype)

p4 <- df %>% 
  ggplot(aes(x = Cluster)) +
  geom_bar(aes(fill = Cell_type)) +
  scale_fill_manual(values = palette) +
  theme_min()

ggsave(plot = p4, glue("{out}/predictions_per_cluster.png"), width = 10, height = 5)

p5 <- df %>% 
  ggplot(aes(x = Cell_type)) +
  geom_bar(aes(fill = Cluster)) +
  scale_fill_manual(values = seurat@misc$colours) +
  theme_min() +
  rotate_x()

ggsave(plot = p5, glue("{out}/clusters_per_celltype.png"), width = 12, height = 4)

message("@ done.")
