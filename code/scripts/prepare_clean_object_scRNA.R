library(here)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(readr)
library(Seurat)


# ------------------------------------------------------------------------------
message("@ loading data...")
load("seurat.Rda")


# ------------------------------------------------------------------------------
message("@ cleaning up...")
columns_keep <- c(
    # QC cols for RNA
    "nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo",
    # Cell cycle
    "G1.S.score_Whitfield", "G2.M.score_Whitfield", "S.Score", "G2M.Score", "Phase",
    # cNMF
    colnames(seurat@meta.data)[grepl("cNMF_program_malignant", colnames(seurat@meta.data))],
    # Cell type assignment / normal malignant calls
    "Cell_type_granular_mouse_correlations",
    "Cell_type_mouse_correlations",
    "Cell_type_consensus_Jessa2022",
    "Malignant_normal_consensus_Jessa2022")


seurat@meta.data <- seurat@meta.data[, which(colnames(seurat@meta.data) %in% columns_keep)]
seurat$ID_paper <- basename(getwd())
seurat$Technology <- "sc/snRNAseq"
seurat$Cell_barcode <- colnames(seurat)
metadata <- seurat@meta.data %>% 
    relocate(ID_paper, Cell_barcode, Technology, .before = 1)
rownames(metadata) <- NULL

seurat@assays$SCENIC <- NULL
seurat@assays$SCT    <- NULL


# ------------------------------------------------------------------------------
message("@ saving...")
save(seurat, file = paste0(basename(getwd()), "__seurat.clean.Rda"))
write_tsv(metadata, file = paste0(basename(getwd()), "__metadata.clean.tsv"))

message("@ done.")
