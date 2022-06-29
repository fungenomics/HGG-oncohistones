# Selin Jessa
# 2020-11-28
#
# This script is run for an individual Seurat object, on which inferCNV
# and cell type projections with correlations have been performed. The script
# integrates these calls to generate a consensus call for each cell of whether
# it's normal or malignant. The resulting call is contained in a new column
# in the metadata, "Malignant_normal_consensus", with a string saved at
# `seurat@misc$malignant_normal_consensus` encoding the method by which the malignant/
# normal call is made.

library(tidyr)
library(purrr)
library(ggplot2)
library(Seurat)

message("@ loading data...")
load("seurat.Rda")

#' @param seurat Seurat object with columns COR_ref.joint_mouse_extended and inferCNV
#' containing 
#'
#' @return Seurat object, with new column Malignant_normal_consensus, consisting
#' of three values: Malignant, Normal, or Likely normal
call_normal_malignant_consensus1 <- function(seurat) {
  
  # Iterate over each cell, comparing CNV/celltype projections
  # We make an assumption based on the biology of pediatric brain tumors
  # about which cell types are highly likely to be normal.
  # If a cell is called normal based on inferCNV *OR* has a projection
  # to microglia, macrophages, endothelial, pericyte, vascular smooth muscle,
  # vascular leptomeningeal, meninges --> then call it Normal
  call_1 <- map2_chr(seurat@meta.data$COR_ref.joint_mouse_extended, seurat@meta.data$inferCNV,
                     function(proj, cnv) {
                       
                       if (cnv == "Normal" | grepl("MGL|MAC|ENDO|PERI|VSMC|VLM|MNG", proj)) "Normal"
                       else "Malignant"
                       
                     })
  
  # Get the proportion in each cluster assigned as normal or malignant
  prop_category <- prop.table(table(Idents(seurat), call_1), margin = 1)
  
  # If more than 50% of cells in one cluster are called Normal, call the cells
  # in the cluster which weren't as Likely normal
  # If more than 50% of cells in one cluster are called Malignant, call the cell
  # in the cluster which weren't as Likely malignant
  call_2 <- map2_chr(Idents(seurat), call_1,
                     function(cluster, call1) {
                       
                       if      (call1 == "Malignant" & prop_category[cluster, "Normal"] > 0.5) "Likely normal"
                       else if (call1 == "Normal"    & prop_category[cluster, "Malignant"] > 0.5) "Likely malignant"
                       else call1
                       
                     })
  
  # Add the call to the object metadata
  seurat <- AddMetaData(seurat, metadata = call_2, col.name = "Malignant_normal_consensus")
  
  seurat@misc$malignant_normal_consensus <- "call_normal_malignant_consensus1: COR_ref.joint_mouse_extended & inferCNV"
  
  return(seurat)
  
}

message("@ calling...")
seurat <- call_normal_malignant_consensus1(seurat)

message("@ saving...")
save(seurat, file = "seurat.Rda")

message("@ done.")

UMAPPlot(seurat, group.by = "Malignant_normal_consensus",
         cols = c("Normal"           = "gray80",
                  "Malignant"        = "red",
                  "Likely normal"    = "gray60",
                  "Likely malignant" = "darksalmon"))

ggsave(".malignant_normal_consensus.png", width = 6, height = 5)
