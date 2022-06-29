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
    # Multiome
    "nCount_peaks", "nFeature_peaks", "nCount_promoters", "nFeature_promoters",
    "nucleosome_signal", "TSS.enrichment",
    "gex_raw_reads"                     , "gex_mapped_reads"                  ,"gex_conf_intergenic_reads"        ,
    "gex_conf_exonic_reads"             , "gex_conf_intronic_reads"           ,"gex_conf_exonic_unique_reads"     ,
    "gex_conf_exonic_antisense_reads"   , "gex_conf_exonic_dup_reads"         ,"gex_exonic_umis"                  ,
    "gex_conf_intronic_unique_reads"    , "gex_conf_intronic_antisense_reads" ,"gex_conf_intronic_dup_reads"      ,
    "gex_intronic_umis"                 , "gex_conf_txomic_unique_reads"      ,"gex_umis_count"                   ,
    "gex_genes_count"                   , "atac_raw_reads"                    ,"atac_unmapped_reads"              ,
    "atac_lowmapq"                      , "atac_dup_reads"                    ,"atac_chimeric_reads"              ,
    "atac_mitochondrial_reads"          , "atac_fragments"                    ,"atac_TSS_fragments"               ,
    "atac_peak_region_fragments"        , "atac_peak_region_cutsites",
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
seurat$Technology <- "scMultiome"
seurat$Cell_barcode <- colnames(seurat)
metadata <- seurat@meta.data %>% 
    relocate(ID_paper, Cell_barcode, Technology, .before = 1)
rownames(metadata) <- NULL

seurat@assays$chromVAR <- NULL
seurat@assays$SCENIC   <- NULL
seurat@assays$TF       <- NULL
seurat@assays$SCT      <- NULL
seurat@assays$ATAC     <- NULL

# ------------------------------------------------------------------------------
message("@ saving...")
save(seurat, file = paste0(basename(getwd()), "__seurat.clean.Rda"))
write_tsv(metadata, file = paste0(basename(getwd()), "__metadata.clean.tsv"))

message("@ done.")
