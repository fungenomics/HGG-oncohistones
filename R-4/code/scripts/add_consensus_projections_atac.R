library(here)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(glue)
library(readr)
library(Seurat)
library(cowplot)

source(here("code/functions/scRNAseq.R"))
source(here("include/style.R"))


# ------------------------------------------------------------------------------
message("@ loading data...")
load("seurat.Rda")

meta_sc <- read_tsv(here("output/00B/metadata_sc.tsv"))


# ------------------------------------------------------------------------------
message("@ making new calls...")

# clean up if this script has been run previously
if ("Cell_type_mouse_correlations" %in% colnames(seurat_atac@meta.data))  seurat_atac$Cell_type_mouse_correlations <- NULL
if ("Malignant_normal_consensus_Jessa2022" %in% colnames(seurat_atac@meta.data))  seurat_atac$Malignant_normal_consensus_Jessa2022 <- NULL

new_cols <- seurat_atac@meta.data %>%
    mutate(Cell_type_mouse_correlations = cluster_predicted.id) %>% 
    mutate(Malignant_normal_consensus_Jessa2022 = case_when(
        Cell_type_mouse_correlations %in% c("Neurons", "Microglia/macrophages", "Immune") ~ "Normal",
        TRUE ~ "Malignant"
    )) %>% 
    select(Cell_type_mouse_correlations, Malignant_normal_consensus_Jessa2022)

# sanity check
all(rownames(new_cols) == rownames(seurat_atac@meta.data))

seurat_atac@meta.data$Type <- NULL
seurat_atac <- AddMetaData(seurat_atac, new_cols)


# ------------------------------------------------------------------------------
message("@ saving...")
save(seurat_atac, file = "seurat.Rda")


# ------------------------------------------------------------------------------
p1 <- umap(seurat_atac, color_by = "Malignant_normal_consensus_Jessa2022", colors = c("Malignant" = "red", "Normal" = "gray50"),
     label = FALSE, legend = TRUE, point_size = 0.3) +
    theme(legend.position = "bottom")

p2 <- umap(seurat_atac, color_by = "Cell_type_mouse_correlations", colors = c(palette_type, "Uncertain" = "gray90"),
     label = FALSE, legend = TRUE, point_size = 0.3) +
    theme(legend.position = "bottom")

plot_grid(p1, p2, ncol = 2, align = "h", axis = "tb")

ggsave(".Jessa2022_calls.png", width = 9.5, height = 5.5)

message("@ done.")
