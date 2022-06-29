# Author: Selin Jessa <selin.jessa@mail.mcgill.ca>
# Date: 2020-10-15
#
# Purpose:
# A script to merge several Seurat objects into one, produce "naively joined"
# objects, and then integrate them using the Harmony method
# 
# Usage:
# This should be run in the scRNAseq/integrations folders, and expects three config
# files to be present:
#
# - info.samples.tsv: TSV file with at least two columns, "Path", which has the path
#   to the Seurat object, and
#   "Sample" which contains the name of the associated sample. Additional columns
#   are optional, but can provide covariates, which can be integrated out, or simply
#   used to colour the low-dimensional embeddings. Must be discrete covariates.
#   The only expectation on Seurat objects is that they have counts for the assay
#   specified in info.experiment.tsv
#
# - info.groups.tsv: TSV file which specifies the groups within each covariate, and
#   the color associated with each level. Three columns: "Covariate", whose values should
#   match the optional columns in info.samples.tsv; "Levels", one row per unique value
#   of each covariate; "Color", indicating the color to use for each level of each
#   covariate, e.g.
#
#       Covariate   Level     Color
#       Location    Cortex    red
#       Location    Pons      blue
#       Mutation    H3.1K27M  green
#       Mutation    H3.3K27M  yellow
#
# - info.experiment.tsv: TSV with two columns, param, and value, indicating parameters
#   used in the analysis:
#
#   vars_integrate   String specifying which covariates to integrate out (should match
#                    columns in info.samples.tsv other than Path), use commas to
#                    separate multiple variables. All variables indicated here will be
#                    integrated out together/at once (not one at a time).
#
#   filters          String specifing expression to use to filter cells based on
#                    covariates in info.samples.tsv. Leave NA if none. NOT YET IMPLEMENTED.
#
#   verbose          Logical, whether to print lots of info, passed to functions.
#                    Default: FALSE
#
#   assay            String specifying which assay in the seurat objects to use. Default: RNA
#
#   n_pcs            Numeric, how many PCs to use for naive join, and integrated
#                    analysis.
#
#   integrate        Logical, whether to run the harmony-integration after naive
#                    join. Default: TRUE. Set to FALSE if only the naive join
#                    is needed.
#
#   malignant_only   Logical, whether to subset all seurat objects to malignant
#                    cells only for naive/integrated analysis. Default: FALSE.
#
#   vars_drop        String, for cosmetics/tidyness, a list of columns in the Seurat
#                    metadata to drop from the joint objects, if they don't make
#                    sense to use in a joint space.
#

# Config -----------------------------------------------------------------------

# load libraries
library(here)
library(readr)
library(dplyr)
library(tibble)
library(glue)
library(purrr)
library(ggplot2)
library(Seurat)
library(harmony)

set.seed(100)

here::here()

dir.create("figures")
dir.create("output")

message("@ loading data...")

# load the info files, which provide the configuration for the analysis
info_samples    <- read_tsv("info.samples.tsv")

# convert to a list of colour palettes
info_groups     <- read_tsv("info.groups.tsv") %>% split(f = .$Covariate) %>% map(~ .x %>% select(-Covariate) %>% tibble::deframe())
# make the sample palette a covariate
info_groups$Sample <- info_samples %>% select(Sample, Color) %>% tibble::deframe()

# convert to a named list where names are parameters and elements are values
info_experiment <- as.list(tibble::deframe(read_tsv("info.experiment.tsv")[, c(1, 2)]))
info_experiment$vars_integrate <- strsplit(info_experiment$vars_integrate, split = ",") %>% unlist()
info_experiment$vars_drop      <- strsplit(info_experiment$vars_drop, split = ",") %>% unlist()
info_experiment                <- lapply(info_experiment, type.convert, as.is = TRUE)


# Prepare data -----------------------------------------------------------------
# 1. load the individual objects
# 2. apply the sample-group association to the metadata
# 3. apply specified filters as they're loaded

message("@ loading Seurat objects...")
# desired behaviour: if the processed joint object already exists, do not
# regenerate it, to allow for cases where the job might fail in the second half
# of the script due to mem/time limits, and allow for restarting part way
if (file.exists("output/seurat_joint.Rda")) {
    
    message("@ found joint object at output/seurat_joint.Rda\n",
            "@ loading instead of recalculating...\n",
            "@ to recalculate, clear old results by deleting output/figures folders")
    
    load("output/seurat_joint.Rda")
    
} else {
    
    # initialize empty list
    seurat_indiv <- list()
    
    for (row in 1:nrow(info_samples)) {
        
        # since these are .Rda objects (as opposed to .Rds), we need to use this
        # paradigm of get(load(.)) in order to save the contents of the .Rda into the
        # variable called `seurat`. (Used defensively in case the object in the .Rda
        # is not named seurat)
        seurat <- get(load(info_samples[row, ]$Path))
        
        seurat@project.name <- info_samples[row, ]$Sample
        seurat$Sample <- seurat@project.name
        print(seurat@project.name)
        
        # remove certain columns from the individual seurat objects,
        # for cleanliness, when columns from the individual space would not apply
        # in the joint space
        seurat@meta.data <- seurat@meta.data[, ! colnames(seurat@meta.data) %in% info_experiment$vars_drop]
        
        # put the covariate/group info in the metadata; it will be the same for all cells
        # since this is sample-level group info
        for (covariate in names(info_groups)) seurat <- AddMetaData(seurat,
                                                                    metadata = unlist(info_samples[row, covariate]),
                                                                    col.name = covariate)
        
        # subset to malignant cells only
        if (!is.null(info_experiment$malignant_only) & info_experiment$malignant_only) {
            
            message("@ subsetting to malignant cells only...")
            
            seurat <- subset(seurat,
                             subset = Malignant_normal_consensus %in% c("Malignant", "Likely malignant"))
            
        }
        
        # populate the list
        seurat_indiv[[row]] <- seurat
        
        # clean up
        rm(seurat)
        
    }
    
    # Preprocess data ------------------------------------------------------------
    
    message("@ preprocessing data...")
    
    # merge into a single Seurat object, normalize, scale, and run PCA
    seurat_joint <- merge(x = seurat_indiv[[1]],
                          y = seurat_indiv[2:length(seurat_indiv)],
                          merge.data = FALSE)
    
    # clean up
    rm(seurat_indiv)
    
    seurat_joint <- seurat_joint %>% 
        # all samples should have been normalized the same way, but just in case, 
        # re-run it here
        Seurat::NormalizeData(verbose = info_experiment$verbose) %>% 
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
        ScaleData(verbose = FALSE) %>% 
        RunPCA(pc.genes = .@var.genes, npcs = info_experiment$n_pcs, verbose = info_experiment$verbose) %>% 
        RunTSNE(dims = 1:info_experiment$n_pcs, verbose = info_experiment$verbose, seed.use = 100, check_duplicates = FALSE) %>%
        RunUMAP(dims = 1:info_experiment$n_pcs, verbose = info_experiment$verbose, seed.use = 100)
    
    # perform clustering
    seurat_joint <- seurat_joint %>% 
        FindNeighbors(seurat,
                      reduction = "pca",
                      dims      = 1:info_experiment$n_pcs,
                      verbose   = TRUE,
                      nn.eps    = 0.5) %>% 
        FindClusters(verbose = info_experiment$verbose,
                     n.start = 10,
                     random.seed = 100, 
                     resolution = 0.5)
    
    # Save -----------------------------------------------------------------------
    
    message("@ saving joined data...")
    
    seurat_joint_dr <- list(
        "pca"  = seurat_joint@reductions$pca@cell.embeddings[, c(1, 2)],
        "tsne" = seurat_joint@reductions$tsne@cell.embeddings[, c(1, 2)],
        "umap" = seurat_joint@reductions$umap@cell.embeddings[, c(1, 2)])
    saveRDS(seurat_joint_dr, file = "output/dimred.Rds")
    
    seurat_joint_meta <- seurat_joint@meta.data
    saveRDS(seurat_joint_meta, file = "output/metadata.Rds")
    
    save(seurat_joint, file = "output/seurat_joint.Rda")
    
    # Covariates ----
    
    # colour the low-dimensional embeddings by the covariates
    iwalk(info_groups, function(palette, covariate) {
        
        plot_fun <- purrr::partial(DimPlot,
                                   object = seurat_joint,
                                   group.by = covariate,
                                   cols = palette)
        
        ggsave(plot = plot_fun(reduction = "pca"),  filename = glue("figures/PCA_{covariate}.png"),  width = 10, height = 8)
        ggsave(plot = plot_fun(reduction = "tsne"), filename = glue("figures/tSNE_{covariate}.png"), width = 10, height = 8)
        ggsave(plot = plot_fun(reduction = "umap"), filename = glue("figures/UMAP_{covariate}.png"), width = 10, height = 8)
        
    })
    
    # colour by clusters
    ggsave(plot = DimPlot(seurat_joint, reduction = "pca"),  filename = "figures/PCA_clusters.png",  width = 10, height = 8)
    ggsave(plot = DimPlot(seurat_joint, reduction = "tsne"), filename = "figures/tSNE_clusters.png", width = 10, height = 8)
    ggsave(plot = DimPlot(seurat_joint, reduction = "umap"), filename = "figures/UMAP_clusters.png", width = 10, height = 8)
    
}

# Integrate with Harmony -------------------------------------------------------

if (info_experiment$integrate) {
    
    message("@ integrating data with harmony...")
    
    # run the integration across the selected variables using Harmony; print the 
    # convergence plot for inspection
    png(filename = glue("figures/convergence.png"), width = 500, height = 400)
    seurat_joint_harmony <- seurat_joint %>% 
        RunHarmony(group.by.vars    = info_experiment$vars_integrate,
                   assay.use        = info_experiment$assay,
                   reduction        = "pca",
                   dims.use         = 1:info_experiment$n_pcs,
                   plot_convergence = TRUE,
                   verbose          = info_experiment$verbose)
    dev.off()
    
    # clean up
    rm(seurat_joint)
    
    # Downstream analysis ----------------------------------------------------------
    
    message("@ performing downstream analysis...")
    
    # run tSNE, UMAP, clustering
    seurat_joint_harmony <- seurat_joint_harmony %>% 
        RunTSNE(dims = 1:info_experiment$n_pcs, verbose = info_experiment$verbose, seed.use = 100, reduction = "harmony") %>%
        RunUMAP(dims = 1:info_experiment$n_pcs, verbose = info_experiment$verbose, seed.use = 100, reduction = "harmony")
    
    seurat_joint_harmony <- seurat_joint_harmony %>% 
        FindNeighbors(seurat,
                      reduction = "harmony",
                      dims      = 1:info_experiment$n_pcs,
                      verbose   = TRUE,
                      nn.eps    = 0.5) %>% 
        FindClusters(verbose = info_experiment$verbose,
                     n.start = 10,
                     random.seed = 100, 
                     resolution = 0.5)
    
    # Save ----
    
    message("@ saving integrated data...")
    
    seurat_joint_harmony_dr <- list(
        "pca"  = seurat_joint_harmony@reductions$pca@cell.embeddings[, c(1, 2)],
        "tsne" = seurat_joint_harmony@reductions$tsne@cell.embeddings[, c(1, 2)],
        "umap" = seurat_joint_harmony@reductions$umap@cell.embeddings[, c(1, 2)])
    saveRDS(seurat_joint_harmony_dr, file = "output/dimred.harmony.Rds")
    
    save(seurat_joint_harmony, file = "output/seurat_joint.harmony.Rda")
    
    # Covariates -------------------------------------------------------------------
    
    iwalk(info_groups, function(palette, covariate) {
        
        plot_fun <- purrr::partial(DimPlot,
                                   object = seurat_joint_harmony,
                                   group.by = covariate,
                                   cols = palette)
        
        ggsave(plot = plot_fun(reduction = "pca"),  filename = glue("figures/PCA_{covariate}.harmony.png"),  width = 10, height = 8)
        ggsave(plot = plot_fun(reduction = "tsne"), filename = glue("figures/tSNE_{covariate}.harmony.png"), width = 10, height = 8)
        ggsave(plot = plot_fun(reduction = "umap"), filename = glue("figures/UMAP_{covariate}.harmony.png"), width = 10, height = 8)
        
    })
    
    # colour by clusters
    ggsave(plot = DimPlot(seurat_joint_harmony, reduction = "pca"),  filename = "figures/PCA_clusters.harmony.png",  width = 10, height = 8)
    ggsave(plot = DimPlot(seurat_joint_harmony, reduction = "tsne"), filename = "figures/tSNE_clusters.harmony.png", width = 10, height = 8)
    ggsave(plot = DimPlot(seurat_joint_harmony, reduction = "umap"), filename = "figures/UMAP_clusters.harmony.png", width = 10, height = 8)
    
} else {
    
    message("@ Skipping harmony integration.")
    
}



# Session info -----------------------------------------------------------------
sessionInfo()
