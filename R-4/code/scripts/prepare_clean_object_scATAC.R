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
seurat <- seurat_atac
rm(seurat_atac)

# ------------------------------------------------------------------------------
message("@ cleaning up...")
columns_keep <- c(
    # ATAC
    "nCount_peaks", "nFeature_peaks", "nCount_promoters", "nFeature_promoters",
    "nucleosome_signal", "TSS.enrichment",
    "total"                       ,     "duplicate"                 ,       "chimeric"                        ,
    "unmapped"                    ,     "lowmapq"                   ,       "mitochondrial"                   ,
    "nonprimary"                  ,     "passed_filters"            ,       "is__cell_barcode"                ,
    "excluded_reason"            ,      "TSS_fragments"            ,        "DNase_sensitive_region_fragments",
    "enhancer_region_fragments"  ,      "promoter_region_fragments",        "on_target_fragments"             ,
    "blacklist_region_fragments" ,      "peak_region_fragments"    ,        "peak_region_cutsites"       ,
    # Cell type assignment / normal malignant calls
    "Cell_type_mouse_correlations",
    "Malignant_normal_consensus_Jessa2022")


seurat@meta.data <- seurat@meta.data[, which(colnames(seurat@meta.data) %in% columns_keep)]
seurat$ID_paper <- basename(getwd())
seurat$Technology <- "scATAC"
seurat$Cell_barcode <- colnames(seurat)
metadata <- seurat@meta.data %>% 
    relocate(ID_paper, Cell_barcode, Technology, .before = 1)
rownames(metadata) <- NULL

seurat@assays$TF       <- NULL
seurat@assays$ATAC     <- NULL

# ------------------------------------------------------------------------------
message("@ saving...")
save(seurat, file = paste0(basename(getwd()), "__seurat.clean.Rda"))
write_tsv(metadata, file = paste0(basename(getwd()), "__metadata.clean.tsv"))

message("@ done.")
