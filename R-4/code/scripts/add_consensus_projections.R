# Selin Jessa
# 2022-03-31
#
# This script is run for an individual Seurat object, on which inferCNV
# and cell type projections with correlations have been performed. The script
# integrates these calls along with 1) consensus projections made in the
# document R-4/code/04... and 2) removes cells that form distinct neuronal clusters
# in the integrated space, also identified in document R-4/code/04...
#
# Thus, this script will load each object, and add a new column:
# - Cell_type_consensus_Jessa2022
# - Malignant_normal_consensus_Jessa2022
#
# We will load each individual Seurat object and the consensus cell type projections
# at R-4/output/projections_agg_consensus_final.Rda
#
# We assign Malignant_normal_consensus_Jessa2022 as follows:
# - if in the neuronal clusters ---> Normal
# - if Malignant_normal_consensus is Normal or Likely normal ---> Normal
# - if Malignant_normal_consensnsus is Malignant or Likely malignant ---> Malignant
#
# We assign Cell_type_consensus_Jessa2022 as follows:
# - if in the neuronal clusters ---> Uncertain
# - otherwise, assign the consensus projection from R-4/output/projections_agg_consensus_final.Rda


# ------------------------------------------------------------------------------
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
load(here("R-4/output/04A/projections_agg_consensus_final.Rda"))

neurons_remove <- bind_rows(readRDS(here("R-4/output/04A/neurons_remove_H3.12.Rds")),
                            readRDS(here("R-4/output/04A/neurons_remove_PFA.Rds")))
meta_sc <- read_tsv(here("output/00B/metadata_sc.tsv"))


# ------------------------------------------------------------------------------
message("@ making new calls...")

# clean up if this script has been run previously
if ("Cell_type_granular_mouse_correlations" %in% colnames(seurat@meta.data)) seurat$Cell_type_granular_mouse_correlations <- NULL
if ("Cell_type_mouse_correlations" %in% colnames(seurat@meta.data)) seurat$Cell_type_mouse_correlations <- NULL
if ("Cell_type_consensus_Jessa2022" %in% colnames(seurat@meta.data))  seurat$Cell_type_consensus_Jessa2022 <- NULL
if ("Malignant_normal_consensus_Jessa2022" %in% colnames(seurat@meta.data))  seurat$Malignant_normal_consensus_Jessa2022 <- NULL

seurat <- summarize_cell_types.Seurat(seurat, "COR_ref.joint_mouse_extended")

seurat_meta <- seurat@meta.data %>%
    dplyr::select(-gex_barcode) %>%
    tibble::rownames_to_column(var = "cell.barcode") %>%
    mutate(rownames = cell.barcode) %>% 
    separate(cell.barcode, into = c("cellname", "drop"), sep = "_") %>%
    dplyr::select(rownames, cellname, sample.id = orig.ident, Malignant_normal_consensus, Cell_type_mouse_correlations = Type) %>%
    left_join(meta_sc %>% separate(scMultiome_path, sep = "/", into = c("1", "2", "3", "sample.id", "4")) %>% select(sample.id, ID_paper)) %>% 
    mutate(Data = case_when(
        grepl("Multiome", getwd()) ~ "scMultiome",
        TRUE ~ "scRNAseq"
    )) %>%
    inner_join(projections_agg_consensus_final %>% dplyr::select(ID_paper = Sample, Data, cellname, Consensus_class),
               by = c("ID_paper", "Data", "cellname"))

# implement rules described above
new_cols <- seurat_meta %>% 
    mutate(
        Cell_type_consensus_Jessa2022 = case_when(
            Consensus_class == "Normal" & cellname %in% neurons_remove$cellname ~ "Uncertain",
            TRUE ~ Consensus_class),
        Malignant_normal_consensus_Jessa2022 = case_when(
            Consensus_class == "Normal" & cellname %in% neurons_remove$cellname ~ "Normal",
            Malignant_normal_consensus %in% c("Normal", "Likely normal") ~ "Normal",
            Malignant_normal_consensus %in% c("Malignant", "Likely malignant") ~ "Malignant"
        )) %>%
    dplyr::select(rownames, Cell_type_mouse_correlations, Cell_type_consensus_Jessa2022, Malignant_normal_consensus_Jessa2022) %>% 
    tibble::column_to_rownames(var = "rownames")

# sanity check
all(rownames(new_cols) == rownames(seurat@meta.data))

seurat@meta.data$Type <- NULL
seurat@meta.data$Cell_type_granular_mouse_correlations <- seurat@meta.data$COR_ref.joint_mouse_extended
seurat <- AddMetaData(seurat, new_cols)


# ------------------------------------------------------------------------------
message("@ saving...")
save(seurat, file = "seurat.Rda")


# ------------------------------------------------------------------------------
p1 <- umap(seurat, color_by = "Malignant_normal_consensus_Jessa2022", colors = c("Malignant" = "red", "Normal" = "gray50"),
     label = FALSE, legend = TRUE, point_size = 0.3) +
    theme(legend.position = "bottom")

p2 <- umap(seurat, color_by = "Cell_type_consensus_Jessa2022", colors = c(palette_type, "Uncertain" = "gray90"),
     label = FALSE, legend = TRUE, point_size = 0.3) +
    theme(legend.position = "bottom")

plot_grid(p1, p2, ncol = 2, align = "h", axis = "tb")

ggsave(".Jessa2022_calls.png", width = 9.5, height = 5.5)

message("@ done.")
