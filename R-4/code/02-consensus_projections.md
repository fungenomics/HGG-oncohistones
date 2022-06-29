---
title: "02 - Consensus tumor cell type projections"
author: "Selin Jessa [[selin.jessa@mail.mcgill.ca](mailto:selin.jessa@mail.mcgill.ca)]"
date: "29 June, 2022"
params:
  resources: "NOT SPECIFIED"
output:
  html_document:
    keep_md: yes
    code_folding: show
    theme: flatly
    css: ../../include/style.css
    toc: yes
    toc_depth: 3
    number_sections: true
    df_print: paged
    includes:
      before_body: ../../include/header.html
      after_body:  ../../include/footer.R4.html
---

<!-- FRONT MATTER, insert configuration info -->


<!-- Load custom CSS/JS for code folding -->
<link rel="stylesheet" type="text/css" href="../../include/hideOutput.css">
<script src="../../include/hideOutput.js"></script>

***

# Configuration

Configuration of project directory & analysis outputs:

<details><summary>Show full config</summary>

```r
source(here("rr_helpers.R"))

# Set up outputs
message("Document index: ", doc_id)
```

```
## Document index: 02
```

```r
# Specify where to save outputs
out        <- here("R-4/output", doc_id); dir.create(out, recursive = TRUE)
figout     <- here("R-4/figures", doc_id); dir.create(figout, recursive = TRUE)
cache      <- paste0(readLines(here("include/project_root.txt")), "R-4.1.2/", basename(here()), "/", doc_id, "/")
```

</details>


Outputs and figures will be saved at these paths, relative to project root:

```
## public/R-4/output/02
```

```
## public/R-4/figures/02
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

Derivation of consensus cell type projections for tumors, as shown in Figure 1
and used throughout the paper.

# Libraries


```r
# Load libraries here
library(here)
library(tidyr)
library(dplyr)
library(ggrepel)
library(ggrastr)
library(data.table)
library(readr)
library(readxl)
library(glue)
library(tibble)
library(ggplot2)
library(GenomicRanges)
library(purrr)
library(feather)
library(cowplot)
library(Signac)
library(Seurat)

source(here("include/style.R"))
source(here("code/functions/scRNAseq.R"))
ggplot2::theme_set(theme_min())
```

# Load metadata


```r
meta_sc <- data.table::fread(here("data/metadata/metadata_sc.tsv"), data.table = FALSE)
```


# Post-clustering QC of integrated data {.tabset}

For the three tumor types where we performed Harmony integration on malignant cells,
we'll perform a final step of QC using the clustering in the integrated space to
identify distinct clusters of neurons, and clusters of proliferating cells.

## H3.3K27M


```r
seurat_H3.3 <- get(load(here("R-4/data/integrations/H3.3K27M_malignant/output/seurat_joint.harmony.Rda")))
seurat_H3.3 %<>% summarize_cell_types.Seurat(cluster_col = "COR_ref.joint_mouse_extended") %>%
    DietSeurat(counts = TRUE, dimreducs = "umap")
rm(seurat_joint_harmony)
```


```r
umap(seurat_H3.3, label = TRUE, legend = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.3_neurons-1.png)<!-- -->

```r
umap(seurat_H3.3, color_by = "Type", colors = palette_type, label = FALSE, legend = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.3_neurons-2.png)<!-- -->

```r
# feature(seurat_H3.3, "STMN2", point_size = 0.2)
```



```r
umap(seurat_H3.3, color_by = "S.Score", color_by_type = "continuous", colors = ylrd, label = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.3_cc-1.png)<!-- -->

```r
umap(seurat_H3.3, color_by = "G2M.Score", color_by_type = "continuous", colors = ylrd, label = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.3_cc-2.png)<!-- -->

```r
umap(seurat_H3.3, color_by = "Phase", label = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.3_cc-3.png)<!-- -->

```r
rm(seurat_H3.3)
```

No distinct cluster of neurons here to remove.

## H3.1/2K27M


```r
seurat_H3.12 <- get(load(here("R-4/data/integrations/H3.12K27M_malignant/output/seurat_joint.harmony.Rda")))
seurat_H3.12 %<>% summarize_cell_types.Seurat(cluster_col = "COR_ref.joint_mouse_extended") %>%
    DietSeurat(counts = TRUE, dimreducs = "umap")
rm(seurat_joint_harmony)
```


```r
umap(seurat_H3.12, label = TRUE, legend = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.12_neurons-1.png)<!-- -->

```r
umap(seurat_H3.12, color_by = "Type", colors = palette_type, label = FALSE, legend = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.12_neurons-2.png)<!-- -->

```r
feature(seurat_H3.12, "STMN2", point_size = 0.2)
```

```
## @ Searching assay: RNA
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.12_neurons-3.png)<!-- -->

```r
seurat_neurons_remove_H3.12 <- subset(seurat_H3.12, idents = "10")
neurons_remove_H3.12 <- seurat_neurons_remove_H3.12@meta.data %>%
    select(Sample, Technology, Label = Type) %>%
    tibble::rownames_to_column(var = "cell.barcode") %>%
    separate(cell.barcode, into = c("cellname", "drop"), sep = "_") %>%
    mutate(Data = case_when(
        grepl("Multiome", Technology) ~ "scMultiome",
        TRUE ~ "scRNAseq"
    )) %>%
    select(-drop, -Technology)

# prop of cell types in the cluster
sort(table(neurons_remove_H3.12$Label)/nrow(neurons_remove_H3.12), decreasing = TRUE)
```

```
## 
##              Neurons                  OPC     Oligodendrocytes 
##           0.67374005           0.27586207           0.01591512 
##            Ependymal     Vascular & other    Proliferating OPC 
##           0.01326260           0.00795756           0.00530504 
##           Astrocytes                  RGC    Glial progenitors 
##           0.00530504           0.00265252           0.00000000 
## Neuronal progenitors               Immune               Normal 
##           0.00000000           0.00000000           0.00000000
```

```r
saveRDS(neurons_remove_H3.12, file = glue("{out}/neurons_remove_H3.12.Rds"))
```

Here, cluster 10 expresses the neuronal marker STMN2, so we will exclude it from malignant cells.


```r
umap(seurat_H3.12, color_by = "S.Score", color_by_type = "continuous", colors = ylrd, label = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.12_cc-1.png)<!-- -->

```r
umap(seurat_H3.12, color_by = "G2M.Score", color_by_type = "continuous", colors = ylrd, label = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.12_cc-2.png)<!-- -->

```r
umap(seurat_H3.12, color_by = "Phase", label = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.12_cc-3.png)<!-- -->

```r
rm(seurat_H3.12)
```

## PFA


```r
seurat_PFA <- get(load(here("data/scRNAseq/integrations/PFA_malignant/output/seurat_joint.harmony.Rda")))
seurat_PFA %<>% summarize_cell_types.Seurat(cluster_col = "COR_ref.joint_mouse_extended") %>%
    DietSeurat(counts = TRUE, dimreducs = "umap")
rm(seurat_joint_harmony)
```


```r
umap(seurat_PFA, label = TRUE, legend = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_PFA_neurons-1.png)<!-- -->

```r
umap(seurat_PFA, color_by = "Type", colors = palette_type, label = FALSE, legend = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_PFA_neurons-2.png)<!-- -->

```r
feature(seurat_PFA, "STMN2", point_size = 0.2)
```

```
## @ Searching assay: RNA
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_PFA_neurons-3.png)<!-- -->

```r
seurat_neurons_remove_PFA <- subset(seurat_PFA, idents = "12")
neurons_remove_PFA <- seurat_neurons_remove_PFA@meta.data %>%
    select(Sample, Technology, Label = Type) %>%
    tibble::rownames_to_column(var = "cell.barcode") %>%
    separate(cell.barcode, into = c("cellname", "drop"), sep = "_") %>%
    mutate(Data = case_when(
        grepl("Multiome", Technology) ~ "scMultiome",
        TRUE ~ "scRNAseq"
    )) %>%
    select(-drop, -Technology)

# prop of cell types in the cluster
sort(table(neurons_remove_PFA$Label)/nrow(neurons_remove_PFA), decreasing = TRUE)
```

```
## 
##            Ependymal    Glial progenitors              Neurons 
##          0.699248120          0.135338346          0.110275689 
##           Astrocytes                  OPC     Vascular & other 
##          0.032581454          0.015037594          0.005012531 
##                  RGC    Proliferating OPC     Oligodendrocytes 
##          0.002506266          0.000000000          0.000000000 
## Neuronal progenitors               Immune               Normal 
##          0.000000000          0.000000000          0.000000000
```

```r
saveRDS(neurons_remove_PFA, file = glue("{out}/neurons_remove_PFA.Rds"))
```

Here, a subset of cells in cluster 12 expresses the neuronal marker STMN2, so we will exclude it from malignant cells.


```r
umap(seurat_PFA, color_by = "S.Score", color_by_type = "continuous", colors = ylrd, label = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_PFA_cc-1.png)<!-- -->

```r
umap(seurat_PFA, color_by = "G2M.Score", color_by_type = "continuous", colors = ylrd, label = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_PFA_cc-2.png)<!-- -->

```r
umap(seurat_PFA, color_by = "Phase", label = FALSE, point_size = 0.2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_PFA_cc-3.png)<!-- -->

```r
rm(seurat_PFA)
```

# Extract samples


```r
# get all single-cell transcriptome data from scRNAseq or scMultiome
rna_samples <- c(list.files(here("data/scRNAseq/pipeline_10X/"), full.names = TRUE))
rna_samples <- rna_samples[!grepl("Makefile", rna_samples)]

rna_seurat_paths <- file.path(rna_samples, "seurat.Rda")
all(file.exists(rna_seurat_paths))
```

```
## [1] TRUE
```

```r
write_lines(rna_seurat_paths, glue("{out}/rna_samples.txt"))

multi_samples <- list.files(here("R-4/data/scMultiome/pipeline_10X_Multiome/"), full.names = TRUE)
multi_samples <- multi_samples[!grepl("Makefile", multi_samples)]

multi_seurat_paths <- file.path(multi_samples, "seurat.Rda")
all(file.exists(multi_seurat_paths))
```

```
## [1] TRUE
```

```r
write_lines(multi_seurat_paths, glue("{out}/multiome_samples.txt"))

sc_samples <- c(rna_samples, multi_samples)
length(sc_samples)
```

```
## [1] 47
```

```r
datatypes_samples <- c(rep("scRNAseq", length(rna_samples)),
                       rep("scMultiome", length(multi_samples)))
```

# Projection pipeline

Each sample was projected to the mouse brain reference atlas using an in-house
Snakemake-based projection workflow which performs classification using SVM,
SciBet, and Spearman correlations. The results of that workflow are loaded here.


# Analyse projections

## Load results

First, load the cell type predictions, based on multiple methods, from the projections
pipeline:


```r
projections <- map2_dfr(sc_samples, datatypes_samples, function(i, datatype) {
    
    # message("@ ", basename(i))
    id <- basename(i)
    path_hydra <- gsub("from_narval", "from_hydra", i)
    path_hydra <- gsub("public", "stable", path_hydra)
    
    # load gene scores file from per-sample cNMF output
    proj <- data.table::fread(file.path(path_hydra, "prediction_pipeline/Prediction_Summary.tsv"),
                              data.table = FALSE, header = TRUE) %>% 
        tibble::add_column("Sample" = basename(i), .before = 1) %>% 
        tibble::add_column("Data" = datatype, .after = 1)
    
    corr <- data.table::fread(file.path(i, "correlations/ref.joint_mouse_extended/correlation.predicted_labels.tsv"),
                              data.table = FALSE, header = TRUE)
    
    # sanity check
    all(proj$cellnames == corr$cell)
    
    proj$Correlation1_Prediction <- corr$celltype
    
    return(proj)
    
})

dim(projections)
```

```
## [1] 224238     14
```

Next, we aggregate predictions from each method into the broader cell class
as used throughout the paper:


```r
projections_agg <- projections %>% 
    # Correlations performed in the first pass of cell type projection
    summarize_cell_types(cluster_col = "Correlation1_Prediction") %>% 
    dplyr::rename(Correlation1_Type = Type) %>% 
    # Correlations performed in the with the projection pipeline
    summarize_cell_types(cluster_col = "Correlation_Prediction") %>% 
    dplyr::rename(Correlation2_Type = Type) %>% 
    summarize_cell_types(cluster_col = "SciBet_Pred") %>% 
    dplyr::rename(SciBet_Type = Type) %>%
    summarize_cell_types(cluster_col = "SVM_Predictions") %>% 
    dplyr::rename(SVM_Type = Type)
```

## Projections data-wrangling

Initally, tumor cells were projected to the mouse atlas using the full
mouse atlas. However, in the prediction pipeline, the reference was downsampled
for computational tractability. Let's first confirm the correlations-based predictions
from the pipeline (using the downsampled reference) are consistent with the initial predictions using the same method.

_**NOTE:**_ if samples are profiled by two technologies, the projections from both replicates
will be combined here.


```r
p1 <- projections_agg %>% 
    select(Sample, Correlation1_Prediction, Correlation_Prediction) %>% 
    rowwise() %>% 
    mutate(Correlations_equal = ifelse(Correlation1_Prediction == Correlation_Prediction, "TRUE", "FALSE")) %>% 
    ggplot(aes(x = Sample)) +
    geom_bar(aes(fill = Correlations_equal), width = 0.5) +
    scale_fill_manual(values = c("TRUE" = "gray90", "FALSE" = "red")) +
    rotate_x() +
    ggtitle("Concordance on granular cell types") +
    ylab("# cells") +
    coord_flip() +
    theme(legend.position = "bottom")

p2 <- projections_agg %>% 
    select(Sample, Correlation1_Type, Correlation2_Type) %>% 
    rowwise() %>% 
    mutate(Correlations_equal = ifelse(Correlation1_Type == Correlation2_Type, "TRUE", "FALSE")) %>% 
    ggplot(aes(x = Sample)) +
    geom_bar(aes(fill = Correlations_equal), width = 0.5) +
    scale_fill_manual(values = c("TRUE" = "gray90", "FALSE" = "red")) +
    rotate_x() +
    ggtitle("Concordance on cell classes") +
    ylab("# cells") +
    coord_flip() +
    no_legend()

plot_grid(p1, p2, align = "h", axis = "tb")
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/correlations_concordance-1.png)<!-- -->

```r
# count how many neurons predicted by correlations are supported by both versions
projections_agg %>% 
    select(Sample, Correlation1_Type, Correlation2_Type) %>% 
    rowwise() %>% 
    mutate(Correlations_equal = ifelse(Correlation1_Type == Correlation2_Type, "TRUE", "FALSE")) %>%
    filter(Correlation1_Type == "Neurons") %>%
    pull(Correlations_equal) %>%
    table()
```

```
## .
## FALSE  TRUE 
##   749  6227
```

What we want to do is count how many methods support the predictions.


```r
projections_agg_long <- projections_agg %>%
    select(-Correlation1_Type) %>% 
    pivot_longer(matches("_Type"), names_to = "Method", values_to = "Prediction")

# for each cell, get most common prediction and the number of methods supporting it
projections_agg_counts <- projections_agg_long %>% 
    group_by(Sample, Data, cellname, Prediction) %>% 
    count() %>% 
    # peel back one layer
    ungroup(Prediction) %>% 
    # with_ties = TRUE means that for cells with two predictions equally common will be dropped,
    # both predictions will be retained
    slice_max(n, n = 1, with_ties = TRUE)

# if there are cases with multiple equally common predictions, drop them and replace by "Uncertain"2
n_preds_per_cell <- projections_agg_counts %>% 
    group_by(Sample, Data, cellname) %>%
    count()

# these cases will have 3 (a different prediction from each method)
proj_uncertain <- projections_agg_counts %>% 
    # right join will filter the cells in projections_agg_counts to those with
    # more than one equally common prediction
    right_join(n_preds_per_cell %>% filter(n == 3), by = c("Sample", "Data", "cellname")) %>% 
    select(-n.x, -n.y, -Prediction) %>% 
    mutate(Consensus_class = "Uncertain",
           N_methods = 0, # Uncertain
           Methods = "Uncertain") %>% 
    distinct(Sample, Data, cellname, Consensus_class, N_methods, Methods)

proj_certain <- projections_agg_counts %>% 
    right_join(n_preds_per_cell %>% filter(n == 1), by = c("Sample", "Data", "cellname")) %>% 
    select(-n.y) %>% 
    dplyr::rename(Consensus_class = Prediction,
                  N_methods = n.x)

# for each cell with a confident consensus prediction, get the methods that support the consensus
proj_certain2 <- projections_agg_long %>%
    left_join(proj_certain, by = c("Sample", "Data", "cellname")) %>%
    distinct() %>% 
    filter(Prediction == Consensus_class) %>%
    group_by(Sample, Data, cellname, Consensus_class, N_methods) %>%
    mutate(Method = gsub("_Type", "", Method)) %>% 
    summarize(Methods = glue_collapse(Method, sep = ","))
```

```
## `summarise()` has grouped output by 'Sample', 'Data', 'cellname', 'Consensus_class'. You can override using the `.groups` argument.
```

```r
projections_agg_consensus <- bind_rows(proj_uncertain,
                                       proj_certain2) %>% 
    arrange(Sample, Data, cellname)

# sanity check we have the same # cells/rows
nrow(projections_agg)
```

```
## [1] 224238
```

```r
nrow(projections_agg) == nrow(projections_agg_consensus)
```

```
## [1] TRUE
```

```r
save(projections_agg_long, projections_agg_consensus,
     file = glue("{out}/projections_aggregated_consensus.Rda"))

table(projections_agg_consensus$Consensus_class)
```

```
## 
##           Astrocytes            Ependymal    Glial progenitors 
##                34276                17066                 8891 
##               Immune Neuronal progenitors              Neurons 
##                24290                  585                 5834 
##     Oligodendrocytes                  OPC    Proliferating OPC 
##                14809                48123                  459 
##                  RGC            Uncertain     Vascular & other 
##                 4809                58332                 6764
```

## Evaluate predictions across methods

In this section, we'll evaluate the consensus predictions across _all_ cells (normal,
malignant).

Which combinations of methods support each predictions?

NOTE: In this figure, all cells are independent, and each column/bar describes
a specific combination of methods that support a consensus. The combinations are 
not subsets / unions of the others.


```r
# what combos?
unique(projections_agg_consensus$Methods)
```

```
## [1] "Correlation2,SciBet,SVM" "Correlation2,SciBet"    
## [3] "Uncertain"               "SciBet,SVM"             
## [5] "Correlation2,SVM"
```

```r
# get order of frequency
methods_order <- c(setdiff(
    names(sort(table(projections_agg_consensus$Methods), decreasing = TRUE)), "Uncertain"),
    "Uncertain")

palette_type2 <- c(palette_type, "Uncertain" = "red")

projections_agg_consensus %>%
    group_by(Methods) %>%
    summarize(Prop_of_all_cells = n()/nrow(projections_agg_consensus)) %>% 
    arrange(desc(Prop_of_all_cells))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Methods"],"name":[1],"type":["glue"],"align":["right"]},{"label":["Prop_of_all_cells"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"Correlation2,SciBet","2":"0.37808043"},{"1":"Uncertain","2":"0.26013432"},{"1":"Correlation2,SciBet,SVM","2":"0.23519207"},{"1":"Correlation2,SVM","2":"0.07848358"},{"1":"SciBet,SVM","2":"0.04810960"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
projections_agg_consensus %>% 
    mutate(Methods = factor(Methods, levels = methods_order)) %>% 
    ggplot(aes(x = Methods)) +
    geom_bar(aes(fill = Consensus_class), width = 0.5) +
    scale_fill_manual(values = palette_type2) +
    rotate_x() +
    ggtitle("Combinations of methods supporting consensus predictions") +
    ylab("# cells")
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/barplot_method_combos-1.png)<!-- -->


## Define consensus

Finally, we define the consensus predictions as those where the prediction from the correlations method
is supported by at least one other. Cells which belong to neuronal clusters identified
in the Harmony-integrated data are excluded.


```r
neurons_remove <- bind_rows(neurons_remove_H3.12, neurons_remove_PFA)

projections_agg_consensus_final <- projections_agg %>% 
    left_join(neurons_remove %>% mutate(Consensus_class = "Normal") %>% select(-Label), by = c("cellname", "Sample", "Data")) %>% 
    mutate_at(vars(matches("_Type")), as.character) %>% 
    mutate(
        Consensus_class = case_when(
            # Case 1: cells are in the neuronal clusters
            Consensus_class == "Normal" & cellname %in% neurons_remove$cellname ~ "Normal",
            # Case 2: prediction from correlations is supported by at least one other method
            (Correlation1_Type == SciBet_Type) | (Correlation1_Type == SVM_Type) ~ Correlation1_Type,
            TRUE ~ "Uncertain"
        ),
        N_methods = case_when(
            Consensus_class %in% c("Uncertain", "Normal") ~ 0,
            (Correlation1_Type == SciBet_Type) & (Correlation1_Type == SVM_Type) ~ 3,
            TRUE ~ 2
        )) %>% 
    select(Sample, Data, cellname, Correlation1_Type, SciBet_Type, SVM_Type, Consensus_class, N_methods)

dim(projections_agg_consensus_final)
```

```
## [1] 224238      8
```

```r
save(projections_agg_consensus_final, file = glue("{out}/projections_agg_consensus_final.Rda"))
```



# Cell-type projections {.tabset}

Next, let's quantify cell type projections in bar-plots, for malignant cells only.
Uncertain cells, lacking a consensus projection will be shown in gray, while
cells lacking a consensus and with a high G2/M phase score will be shown in a separate color.

First, we need two functions for (1) summarizing the cell-type projection
at a higher level, and (2) generating a bar plot to represent cell-type proportions
within each sample:

<details>


```r
# put uncertain cells in gray
palette_type3 <- c("High G2/M" = "#e65c12", palette_type, "Uncertain" = "gray90")

load_integration <- function(path, key, prediction_filters = quo(TRUE)) {
    
    dim  <- readRDS(here(path, glue("{key}/output/dimred.harmony.Rds"))) %>% .$umap
    
    meta1 <- readRDS(here(path, glue("{key}/output/metadata.Rds")))
    
    if (!("COR_ref.human_fetal_thalamus3" %in% colnames(meta1))) meta1$COR_ref.human_fetal_thalamus3 <- NA
    if (!("gex_barcode" %in% colnames(meta1))) meta1$gex_barcode <- NA
    
    meta <- meta1 %>%
        dplyr::rename(cellname_object = cell.barcode) %>% 
        mutate(cellname_object = case_when(
            is.na(cellname_object) & Technology == "10X Multiome" ~ paste0(orig.ident, "_", gex_barcode),
            TRUE ~ cellname_object
        )) %>% 
        tibble::rownames_to_column(var = "cell.barcode") %>%
        separate(cell.barcode, into = c("cellname", "drop"), sep = "_") %>%
        dplyr::select(cellname, cellname_object, Sample, GrowthFactorReceptor, Location, Molecular,
                      Technology, Malignant_normal_consensus, COR_ref.human_fetal_thalamus3, G2M.Score) %>%
        cbind(dim) %>%
        mutate(Data = case_when(
            grepl("Multiome", Technology) ~ "scMultiome",
            TRUE ~ "scRNAseq"
        )) %>%
        # simplify growth factor mutations
        mutate(GrowthFactorReceptor = case_when(
            grepl("ACVR1",  GrowthFactorReceptor) ~ "ACVR1",
            grepl("PDGFRA", GrowthFactorReceptor) ~ "PDGFRA",
            grepl("BRAF",   GrowthFactorReceptor) ~ "BRAF",
            TRUE ~ GrowthFactorReceptor
        )) %>%
        inner_join(projections_agg_consensus_final,
                   by = c("Sample", "Data", "cellname")) %>%
        mutate(Consensus_class = case_when(
            Consensus_class == "Uncertain" & G2M.Score > 0.5 ~ "High G2/M",
            TRUE ~ Consensus_class
        )) %>% 
        mutate(Location = ifelse(Location == "Thalamus", "Thal.", Location),
               Consensus_class = factor(Consensus_class, levels = unique(names(palette_type3)))) %>%
        filter(!!prediction_filters)
    
    return(meta)
    
}

plot_umaps <- function(df, size = 0.2) {
    
    p1 <- df %>%
        ggplot(aes(x = UMAP_1, y = UMAP_2)) +
        rasterize(geom_point(aes(colour = Consensus_class), size = size, alpha = 0.5), dpi = 600) +
        scale_colour_manual(values = palette_type3) +
        theme_min() +
        theme(legend.position = "bottom") +
        no_ticks()
    
    p1
    
}


plot_cell_type_prop <- function(df) {
    
    # calculate # of cells per sample
    sample_n <- df %>%
        select(Sample, Malignant_normal_consensus, Consensus_class, Location) %>%
        filter(Malignant_normal_consensus %in% c("Malignant", "Likely malignant")) %>%
        group_by(Sample) %>%
        summarise(n = n())
    
    # calculate the proportion of each major cell type among malignant cells of
    # each sample
    cell_type_prop <- df %>%
        select(Sample, Malignant_normal_consensus, Consensus_class, Location) %>%
        filter(Malignant_normal_consensus %in% c("Malignant", "Likely malignant")) %>%
        group_by(Sample, Location, Consensus_class) %>%
        summarise(n = n()) %>%
        mutate(freq = n / sum(n))
    
    # order by certain cell Consensus_classs
    samples_ordered_epen <- cell_type_prop %>%
        select(Sample, Consensus_class, freq, Location) %>%
        spread(Consensus_class, freq) %>%
        rowwise() %>% 
        mutate(Ependymal = ifelse(is.na(Ependymal), 0, Ependymal)) %>%
        arrange(Location, Ependymal, Astrocytes, OPC) %>%
        pull(Sample)
    
    # plot bar chart, ordered by cell Consensus_class
    p1 <- sample_n %>%
        mutate(Sample = factor(Sample, levels = samples_ordered_epen)) %>%
        ggplot(aes(x = Sample, y = n)) +
        geom_bar(fill = "gray60", stat = "identity", width = 0.8) +
        rotate_x() +
        theme(panel.grid = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank()) +
        xlab(NULL)
    
    p2 <- cell_type_prop %>%
        mutate(Sample = factor(Sample, levels = samples_ordered_epen)) %>%
        ggplot(aes(x = Sample, y = freq)) +
        geom_bar(aes(fill = Consensus_class), stat = "identity", width = 0.8) +
        scale_fill_manual(values = palette_type3) +
        rotate_x() +
        ylab("Prop. (malignant cells)") +
        theme(panel.grid = element_blank(),
              axis.ticks.x = element_blank())
    
    p3 <- df %>%
        distinct(Sample, GrowthFactorReceptor) %>%
        mutate(Sample = factor(Sample, levels = samples_ordered_epen)) %>%
        ggplot(aes(x = Sample, y = 1)) +
        geom_tile(aes(fill = GrowthFactorReceptor), colour = "white", width = 0.9, height = 0.9) +
        scale_fill_manual(values = palette_gfr, na.value = "gray90", guide = guide_legend(ncol = 3)) +
        theme_row() +
        no_legend()
    
    p4 <- df %>%
        distinct(Sample, Location) %>%
        mutate(Sample = factor(Sample, levels = samples_ordered_epen)) %>%
        ggplot(aes(x = Sample, y = 1)) +
        geom_tile(aes(fill = Location), colour = "white", width = 0.9, height = 0.9) +
        scale_fill_manual(values = palette_location, na.value = "gray90", guide = guide_legend(ncol = 3)) +
        theme_row() +
        no_legend()
    
    p_all <- plot_grid(p1, p2, p3, p4, ncol = 1, align = "v", axis = "rl",
                       rel_heights = c(0.15, 0.55, 0.15, 0.15))
    
    return(p_all)
    
}
```

</details>


### H3.3K27M (malignant cells)


```r
meta_H3.3_m <- load_integration("R-4/data/integrations", "H3.3K27M_malignant",
                                prediction_filters = quo(Consensus_class != "Normal"))

# stats
length(unique(meta_H3.3_m$Sample))
```

```
## [1] 16
```

```r
table(meta_H3.3_m$Technology)
```

```
## 
##         10X Multiome   10X Single Cell 3' 10X Single Nuclei 3' 
##                29164                19585                43658
```

```r
nrow(meta_H3.3_m)
```

```
## [1] 92407
```

```r
table(meta_H3.3_m$N_methods)
```

```
## 
##     0     2     3 
## 35172 39915 17320
```

```r
table(meta_H3.3_m$Consensus_class)
```

```
## 
##            High G2/M                  RGC    Glial progenitors 
##                 1612                  219                 1913 
##                  OPC    Proliferating OPC     Oligodendrocytes 
##                35395                  355                 8682 
##           Astrocytes            Ependymal Neuronal progenitors 
##                 8727                  313                   39 
##              Neurons               Immune     Vascular & other 
##                  774                  562                  256 
##               Normal            Uncertain 
##                    0                33560
```

```r
save(meta_H3.3_m, file = glue("{out}/metadata_joint_H3.3_m.Rda"))
```


```r
plot_umaps(meta_H3.3_m)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.3K27M_m-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/02umap_H3.3K27M_m...*]~</span>


```r
plot_cell_type_prop(meta_H3.3_m %>% filter(Consensus_class != "Uncertain"))
```

```
## `summarise()` has grouped output by 'Sample', 'Location'. You can override using the `.groups` argument.
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/bar_H3.3K27M_m-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/02bar_H3.3K27M_m...*]~</span>


### H3.3K27M thalamus (malignant cells)


```r
meta_H3.3th_m <- load_integration("R-4/data/integrations", "H3.3K27M_thalamus_malignant",
                                  prediction_filters = quo(Consensus_class != "Normal"))

# stats
length(unique(meta_H3.3th_m$Sample))
```

```
## [1] 4
```

```r
table(meta_H3.3th_m$Technology)
```

```
## 
##         10X Multiome   10X Single Cell 3' 10X Single Nuclei 3' 
##                 5603                11892                22867
```

```r
nrow(meta_H3.3th_m)
```

```
## [1] 40362
```

```r
table(meta_H3.3th_m$N_methods)
```

```
## 
##     0     2     3 
## 13888 18311  8163
```

```r
table(meta_H3.3th_m$Consensus_class)
```

```
## 
##            High G2/M                  RGC    Glial progenitors 
##                  744                  117                 1178 
##                  OPC    Proliferating OPC     Oligodendrocytes 
##                16925                  181                 4846 
##           Astrocytes            Ependymal Neuronal progenitors 
##                 2675                   47                   17 
##              Neurons               Immune     Vascular & other 
##                   50                  339                   99 
##               Normal            Uncertain 
##                    0                13144
```

```r
save(meta_H3.3th_m, file = glue("{out}/metadata_joint_H3.3th_m.Rda"))
```


```r
plot_umaps(meta_H3.3th_m)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.3thK27M_m-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/02umap_H3.3thK27M_m...*]~</span>


```r
plot_cell_type_prop(meta_H3.3th_m %>% filter(Consensus_class != "Uncertain"))
```

```
## `summarise()` has grouped output by 'Sample', 'Location'. You can override using the `.groups` argument.
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/bar_H3.3thK27M_m-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/02bar_H3.3thK27M_m...*]~</span>



### H3.1/2K27M (malignant cells)


```r
meta_H3.12_m <- load_integration("R-4/data/integrations", "H3.12K27M_malignant",
                                 prediction_filters = quo(Consensus_class != "Normal"))

# stats
length(unique(meta_H3.12_m$Sample))
```

```
## [1] 8
```

```r
table(meta_H3.12_m$Technology)
```

```
## 
##         10X Multiome 10X Single Nuclei 3' 
##                 5958                 9853
```

```r
nrow(meta_H3.12_m)
```

```
## [1] 15811
```

```r
table(meta_H3.12_m$N_methods)
```

```
## 
##    0    2    3 
## 5265 7084 3462
```

```r
table(meta_H3.12_m$Consensus_class)
```

```
## 
##            High G2/M                  RGC    Glial progenitors 
##                  519                   39                  126 
##                  OPC    Proliferating OPC     Oligodendrocytes 
##                 4722                   68                  772 
##           Astrocytes            Ependymal Neuronal progenitors 
##                 2975                 1571                   10 
##              Neurons               Immune     Vascular & other 
##                  174                   43                   46 
##               Normal            Uncertain 
##                    0                 4746
```

```r
save(meta_H3.12_m, file = glue("{out}/metadata_joint_H3.12_m.Rda"))
```


```r
plot_umaps(meta_H3.12_m)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.12K27M_m-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/02umap_H3.12K27M_m...*]~</span>


```r
plot_cell_type_prop(meta_H3.12_m %>% filter(Consensus_class != "Uncertain"))
```

```
## `summarise()` has grouped output by 'Sample', 'Location'. You can override using the `.groups` argument.
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/bar_H3.12K27M_m-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/02bar_H3.12K27M_m...*]~</span>


### PFA-EP (malignant cells)


```r
meta_pfa_m <- load_integration("data/scRNAseq/integrations", "PFA_malignant",
                               prediction_filters = quo(Consensus_class != "Normal"))

# stats
length(unique(meta_pfa_m$Sample))
```

```
## [1] 5
```

```r
table(meta_pfa_m$Technology)
```

```
## 
##   10X Single Cell 3' 10X Single Nuclei 3' 
##                14827                 8692
```

```r
nrow(meta_pfa_m)
```

```
## [1] 23519
```

```r
table(meta_pfa_m$N_methods)
```

```
## 
##     0     2     3 
##  4221 13391  5907
```

```r
table(meta_pfa_m$Consensus_class)
```

```
## 
##            High G2/M                  RGC    Glial progenitors 
##                  191                  535                 2773 
##                  OPC    Proliferating OPC     Oligodendrocytes 
##                  236                    1                   25 
##           Astrocytes            Ependymal Neuronal progenitors 
##                  451                14945                    3 
##              Neurons               Immune     Vascular & other 
##                   18                   73                  238 
##               Normal            Uncertain 
##                    0                 4030
```

```r
save(meta_pfa_m, file = glue("{out}/metadata_joint_PFA_m.Rda"))
```


```r
plot_umaps(meta_pfa_m)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_PFA_m-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/02umap_PFA_m...*]~</span>


```r
plot_cell_type_prop(meta_pfa_m %>% filter(Consensus_class != "Uncertain"))
```

```
## `summarise()` has grouped output by 'Sample', 'Location'. You can override using the `.groups` argument.
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/bar_PFA_m-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/02bar_PFA_m...*]~</span>

Clean up


```r
rm(meta_H3.3_m)
rm(meta_H3.12_m)
rm(meta_pfa_m)
```


# Gene expression levels in tumor scRNAseq

## Load data

Load joint object of malignant cells:


```r
seurat_h31 <- get(load(here("R-4/data/integrations/H3.12K27M_malignant/output/seurat_joint.Rda")))
rm(seurat_joint)

seurat_h33 <- get(load(here("R-4/data/integrations/H3.3K27M_malignant/output/seurat_joint.Rda")))
rm(seurat_joint)

seurat_h31_small <- DietSeurat(seurat_h31, data = TRUE, counts = FALSE)
seurat_h33_small <- DietSeurat(seurat_h33, data = TRUE, counts = FALSE)

rm(seurat_h31)
rm(seurat_h33)
```

## Expression

A function to plot levels of each gene per sample:


```r
expr_data <- bind_rows(FetchData(seurat_h31_small,
                                 vars = c("NKX6-1", "PAX3", "Sample",
                                          "Molecular", "Location", "GrowthFactorReceptor", "Technology")),
                       FetchData(seurat_h33_small,
                                 vars = c("NKX6-1", "PAX3", "Sample",
                                          "Molecular", "Location", "GrowthFactorReceptor", "Technology")))

# restrict cells to those with confident predictions
expr_data_anno <- expr_data %>%
    tibble::rownames_to_column(var = "cell.barcode") %>%
    separate(cell.barcode, into = c("cellname", "drop"), sep = "_") %>%
    dplyr::select(cellname, Sample, GrowthFactorReceptor, Location, Molecular,
                  Technology, `NKX6-1`, PAX3) %>%
    mutate(Data = case_when(
        grepl("Multiome", Technology) ~ "scMultiome",
        TRUE ~ "scRNAseq"
    )) %>%
    inner_join(projections_agg_consensus_final,
               by = c("Sample", "Data", "cellname")) %>%
    # filter to cells with projections supported by confident combos
    filter(!(Consensus_class %in% c("Normal", "Uncertain"))) %>%
    select(Sample = Sample, Consensus_class, everything())

dim(expr_data_anno)
```

```
## [1] 67781    14
```

```r
table(expr_data_anno$Consensus_class)
```

```
## 
##           Astrocytes            Ependymal    Glial progenitors 
##                11702                 1884                 2039 
##               Immune Neuronal progenitors              Neurons 
##                  605                   49                  948 
##     Oligodendrocytes                  OPC    Proliferating OPC 
##                 9454                40117                  423 
##                  RGC     Vascular & other 
##                  258                  302
```

```r
saveRDS(expr_data_anno, file = glue("{out}/tumor_scRNA_expr_data.Rds"))

# split by tumor type
celltype_freq <- expr_data_anno %>%
    mutate(Molecular = ifelse(Molecular == "H3.2K27M", "H3.1K27M", Molecular)) %>%
    group_by(Molecular) %>%
    mutate(N_per_molecular = n()) %>%
    group_by(Molecular, N_per_molecular, Consensus_class) %>%
    count() %>%
    mutate(Freq_per_molecular = n/N_per_molecular) %>%
    ungroup() %>%
    select(Molecular, Consensus_class, N = n, Freq = Freq_per_molecular) %>%
    pivot_wider(names_from = Molecular, values_from = c(N, Freq))

celltype_freq
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Consensus_class"],"name":[1],"type":["chr"],"align":["left"]},{"label":["N_H3.1K27M"],"name":[2],"type":["int"],"align":["right"]},{"label":["N_H3.3K27M"],"name":[3],"type":["int"],"align":["right"]},{"label":["Freq_H3.1K27M"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Freq_H3.3K27M"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"Astrocytes","2":"2975","3":"8727","4":"0.2820974777","5":"0.1524766314"},{"1":"Ependymal","2":"1571","3":"313","4":"0.1489664328","5":"0.0054686818"},{"1":"Glial progenitors","2":"126","3":"1913","4":"0.0119476579","5":"0.0334236044"},{"1":"Immune","2":"43","3":"562","4":"0.0040773753","5":"0.0098191666"},{"1":"Neuronal progenitors","2":"10","3":"39","4":"0.0009482268","5":"0.0006814012"},{"1":"Neurons","2":"174","3":"774","4":"0.0164991466","5":"0.0135231938"},{"1":"Oligodendrocytes","2":"772","3":"8682","4":"0.0732031102","5":"0.1516903992"},{"1":"OPC","2":"4722","3":"35395","4":"0.4477527024","5":"0.6184153053"},{"1":"Proliferating OPC","2":"68","3":"355","4":"0.0064479423","5":"0.0062024985"},{"1":"RGC","2":"39","3":"219","4":"0.0036980846","5":"0.0038263300"},{"1":"Vascular & other","2":"46","3":"256","4":"0.0043618434","5":"0.0044727876"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# filter out cells which represent less than 5% of the dataset
(celltype_keep_h31 <- celltype_freq %>% filter(Freq_H3.1K27M > 0.05) %>% pull(Consensus_class))
```

```
## [1] "Astrocytes"       "Ependymal"        "Oligodendrocytes" "OPC"
```

```r
(celltype_keep_h33 <- celltype_freq %>% filter(Freq_H3.3K27M > 0.05) %>% pull(Consensus_class))
```

```
## [1] "Astrocytes"       "Oligodendrocytes" "OPC"
```

```r
pct1_joint_sample <- expr_data_anno %>%
    group_by(Molecular, Sample, Consensus_class) %>%
    # calculate the number of cells of each cell type per sample
    mutate(n_total = n()) %>%
    ungroup() %>%
    gather(Gene, Expression, `NKX6-1`, PAX3) %>%
    group_by(n_total, Molecular, GrowthFactorReceptor, Location, Sample, Consensus_class, Gene) %>%
    # count the number of cells of each cell type where each gene is detected
    summarize(N_detected = sum(Expression > 0)) %>%
    # divide by the total # of cells in the cell type
    mutate(pct1 = round(N_detected / n_total, 2)) %>%
    ungroup() %>%
    select(-n_total) %>%
    mutate(Consensus_class = factor(Consensus_class, levels = rev(names(palette_type))))
```

```
## `summarise()` has grouped output by 'n_total', 'Molecular', 'GrowthFactorReceptor', 'Location', 'Sample', 'Consensus_class'. You can override using the `.groups` argument.
```

```r
plot_boxplot_per_sample <- function(df, title, color, label_ACVR1 = FALSE) {
    
    if (label_ACVR1) df <- df %>%
            mutate(ACVR1_simple = case_when(
                grepl("ACVR1-", GrowthFactorReceptor) ~ "Mutant",
                TRUE ~ "WT"
            ))
    
    p1 <- df %>%
        mutate(Gene = factor(Gene, levels = c("NKX6-1", "PAX3"))) %>%
        complete(Gene, nesting(Consensus_class), fill = list(pct1 = 0)) %>%
        ggplot(aes(x = Consensus_class, y = pct1)) +
        geom_boxplot(lwd = 0.25, fill = color, outlier.shape = NA)
    
    if (label_ACVR1) {
        
        p1 <- p1 +
            geom_point(aes(shape = ACVR1_simple), stat = "identity", alpha = 0.7) +
            scale_shape_manual(values = palette_acvr1_simple)
        
    } else p1 <- p1 + geom_point(stat = "identity", colour = "black", alpha = 0.7)
    
    p1 <- p1 +
        scale_alpha(guide = FALSE) +
        facet_wrap(~ Gene, nrow = 1) +
        coord_flip() +
        rotate_x() +
        # no_legend() +
        ylab("Proportion of cells in which the gene is detected") +
        ylim(c(0, 1)) +
        ggtitle(title) +
        theme(legend.position = "bottom")
    
    return(p1)
    
}
```


Encoding ACVR1 status:


```r
p1 <- plot_boxplot_per_sample(pct1_joint_sample %>%
                                  filter(Molecular %in% c("H3.1K27M", "H3.2K27M") &
                                             Location == "Pons" &
                                             Consensus_class %in% celltype_keep_h31),
                              title = "H3.1/2K27M Pons",
                              "orange",
                              label_ACVR1 = TRUE)

p2 <- plot_boxplot_per_sample(pct1_joint_sample %>%
                                  filter(Molecular == "H3.3K27M" &
                                             Location == "Pons" &
                                             Consensus_class %in% celltype_keep_h33),
                              title = "H3.3K27M Pons",
                              "red3",
                              label_ACVR1 = TRUE)

plot_grid(p1, p2, ncol = 1, align = "h", axis = "tb")
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/tumor_pct1_boxplots_ACVR1-1.png)<!-- -->



# Validation of tumor cell type projections across species

## Load human fetal thalamus data


```r
hf_thalamus_annotation <- data.table::fread(here("output/01A/hf_thalamus_annotation.tsv"),
                                            data.table = FALSE)

summarize_cell_types_bhaduri <- function(df, cluster_col) {
    
    cc_quo <- rlang::sym(cluster_col)
    
    df %>%
        mutate(
            # Define some broader cell type classes
            Type = case_when(
                grepl("GLIP|OPAS|Glial", !!cc_quo) ~ "Glial progenitors",
                grepl("RG|[Rr]adial|NSC|prog|Dividing", !!cc_quo) ~ "RGC",
                grepl("EXIP|INIP|NEURP|IP", !!cc_quo) & !grepl("VIP", !!cc_quo) ~ "Neuronal progenitors",
                grepl("MGE|CGE|LGE|SST|PV|INH|CEX|PEX|GABAN|NEU|SPN|NRGN|UBCN|MFN|SERN|[Nn]euron", !!cc_quo) ~ "Neurons",
                grepl("OPC-P|Proliferating OPC", !!cc_quo) ~ "Proliferating OPC",
                grepl("OPC", !!cc_quo) ~ "OPC",
                grepl("NFOL|MOL|Oligo", !!cc_quo) ~ "Oligodendrocytes",
                grepl("ASTR|Astro", !!cc_quo) ~ "Astrocytes",
                grepl("EPEN|ASEP|Epen", !!cc_quo) ~ "Ependymal",
                grepl("MGL|MAC|Micro", !!cc_quo) ~ "Immune",
                grepl("T|Fibr|Mur|Micro|Mixed|UNR|Other|ENDO|PERI|MNG|Endo", !!cc_quo) ~ "Vascular & other",
                TRUE ~ "Vascular & other")) %>%
        mutate(Type = factor(Type, levels = c("RGC",
                                              "Glial progenitors",
                                              "Proliferating OPC",
                                              "OPC",
                                              "Oligodendrocytes",
                                              "Astrocytes",
                                              "Ependymal",
                                              "Neuronal progenitors",
                                              "Neurons",
                                              "Immune",
                                              "Vascular & other",
                                              "Normal")))
    
}

hf_thalamus_info_clusters <- data.table::fread(here("output/01A/info_clusters3.tsv"),
                                               data.table = FALSE) %>%
    filter(!grepl("EXCLUDE", ID_20220321)) %>% 
    summarize_cell_types_bhaduri("ID_20220321") %>%
    mutate(Label = paste0(Sample, " #", Cluster_number))
```

## Load tumor data

We'll load the combined metadata for all H3.3K27M thalamic samples, restricted to malignant cells:


```r
meta_H3.3thal_m <- load_integration("R-4/data/integrations", "H3.3K27M_thalamus_malignant",
                                    prediction_filters = quo(!(Consensus_class %in% c("Normal", "Uncertain"))))

length(unique(meta_H3.3thal_m$Sample))
```

```
## [1] 4
```

```r
table(meta_H3.3thal_m$Technology)
```

```
## 
##         10X Multiome   10X Single Cell 3' 10X Single Nuclei 3' 
##                 4483                 7760                14975
```

```r
nrow(meta_H3.3thal_m)
```

```
## [1] 27218
```

```r
table(meta_H3.3thal_m$N_methods)
```

```
## 
##     0     2     3 
##   744 18311  8163
```

```r
table(meta_H3.3thal_m$Consensus_class)
```

```
## 
##            High G2/M                  RGC    Glial progenitors 
##                  744                  117                 1178 
##                  OPC    Proliferating OPC     Oligodendrocytes 
##                16925                  181                 4846 
##           Astrocytes            Ependymal Neuronal progenitors 
##                 2675                   47                   17 
##              Neurons               Immune     Vascular & other 
##                   50                  339                   99 
##               Normal            Uncertain 
##                    0                    0
```

```r
meta_H3.3thal_m <- meta_H3.3thal_m %>%
    left_join(hf_thalamus_info_clusters %>% select(Cluster, Type), by = c("COR_ref.human_fetal_thalamus3" = "Cluster")) %>%
    dplyr::rename(Type_human = Type)

# total agreement
meta_H3.3thal_m %>%
    mutate(Type_human = as.character(Type_human),
           Consensus_class = as.character(Consensus_class)) %>%
    filter(Type_human == Consensus_class) %>% { nrow(.)/nrow(meta_H3.3thal_m) }
```

```
## [1] 0.8056801
```


## Mouse vs. human projections

UMAP coloured by human projections:


```r
meta_H3.3thal_m %>% 
    ggplot(aes(x = UMAP_1, y = UMAP_2)) +
    rasterize(geom_point(aes(colour = Type_human), size = 0.2, alpha = 0.5), dpi = 600) +
    scale_colour_manual(values = palette_type3) +
    theme_min() +
    theme(legend.position = "bottom") +
    no_ticks()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/umap_H3.3thK27M_m_labelled_thalamus-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/02umap_H3.3thK27M_m_labelled_thalamus...*]~</span>

Plot the agreement as a confusion matrix representing the proportion of each mouse cell
types that are predicted to human cell types:


```r
# encode number of cells as well as proportion of agreement
p3 <- meta_H3.3thal_m %>%
    group_by(Consensus_class) %>%
    mutate(Total_per_class_mouse = n()) %>%
    group_by(Consensus_class, Total_per_class_mouse, Type_human) %>%
    count() %>%
    mutate(Prop = n/Total_per_class_mouse) %>%
    mutate(Type_human = factor(Type_human, levels = setdiff(names(palette_type), "Normal")),
           Consensus_class = factor(Consensus_class, levels = rev(setdiff(names(palette_type), "Normal")))) %>%
    ungroup() %>%
    complete(Consensus_class, Type_human, fill = list(Prop = 0)) %>%
    filter(Consensus_class != "Uncertain" & Type_human != "Uncertain") %>%
    ggplot(aes(x = Type_human, y = Consensus_class)) +
    geom_point(aes(colour = Prop, size = n)) +
    scale_colour_gradientn(colors = colorRampPalette(c("gray90", "red"))(n=100), limits = c(0, 1)) +
    scale_size_area(breaks = seq(2000, 10000, by = 2000), max_size = 10) +
    geom_text(aes(label = round(Prop, 2)), colour = "black", size = 3) +
    scale_x_discrete(position = "top") +
    theme_min() +
    theme(panel.border = element_blank()) +
    rotate_x() +
    theme(axis.text.x = element_text(hjust = 0)) +
    ggtitle("Proportion of cells (labelled by mouse) \nprojected to each human cell class") +
    xlab("Human") + ylab("Mouse")

p4 <- meta_H3.3thal_m %>%
    group_by(Type_human) %>%
    mutate(Total_per_class_human = n()) %>%
    group_by(Type_human, Total_per_class_human, Consensus_class) %>%
    count() %>%
    mutate(Prop = n/Total_per_class_human) %>%
    mutate(Type_human = factor(Type_human, levels = rev(setdiff(names(palette_type), "Normal"))),
           Consensus_class = factor(Consensus_class, levels = setdiff(names(palette_type), "Normal"))) %>%
    ungroup() %>%
    complete(Consensus_class, Type_human, fill = list(Prop = 0)) %>%
    filter(Consensus_class != "Uncertain" & Type_human != "Uncertain") %>%
    ggplot(aes(x = Consensus_class, y = Type_human)) +
    geom_point(aes(colour = Prop, size = n)) +
    scale_colour_gradientn(colors = colorRampPalette(c("gray90", "red"))(n=100), limits = c(0, 1)) +
    scale_size_area(breaks = seq(2000, 10000, by = 2000), max_size = 10) +
    geom_text(aes(label = round(Prop, 2)), colour = "black", size = 3) +
    scale_x_discrete(position = "top") +
    theme_min() +
    theme(panel.border = element_blank()) +
    rotate_x() +
    theme(axis.text.x = element_text(hjust = 0)) +
    ggtitle("Proportion of cells (labelled by human) \nprojected to each mouse cell class") +
    xlab("Mouse") + ylab("Human")

plot_grid(p3, p4, ncol = 2)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/02/confusion_matrix_mouse_human_thalamus-1.png)<!-- -->


<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2022-06-29 11:08:00
```

The git repository and last commit:

```
## Local:    master /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public
## Remote:   master @ origin (git@github.com:fungenomics/HGG-oncohistones.git)
## Head:     [6e4c415] 2022-06-29: Add functions for R 3.6
```

The random seed was set with `set.seed(100)`

The R session info:
<details>

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Rocky Linux 8.5 (Green Obsidian)
## 
## Matrix products: default
## BLAS/LAPACK: /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/flexiblas/3.0.4/lib64/libflexiblas.so.3.0
## 
## locale:
##  [1] LC_CTYPE=en_CA.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_CA.UTF-8        LC_COLLATE=en_CA.UTF-8    
##  [5] LC_MONETARY=en_CA.UTF-8    LC_MESSAGES=en_CA.UTF-8   
##  [7] LC_PAPER=en_CA.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices datasets  utils    
## [8] methods   base     
## 
## other attached packages:
##  [1] magrittr_2.0.1       viridis_0.5.1        viridisLite_0.3.0   
##  [4] RColorBrewer_1.1-2   SeuratObject_4.0.4   Seurat_4.0.0        
##  [7] Signac_1.3.0         cowplot_1.1.1        feather_0.3.5       
## [10] purrr_0.3.4          GenomicRanges_1.42.0 GenomeInfoDb_1.26.7 
## [13] IRanges_2.24.1       S4Vectors_0.28.1     BiocGenerics_0.36.1 
## [16] tibble_3.1.6         glue_1.6.1           readxl_1.3.1        
## [19] readr_2.1.1          data.table_1.14.2    ggrastr_0.2.3       
## [22] ggrepel_0.9.1        ggplot2_3.3.5        dplyr_1.0.7         
## [25] tidyr_1.1.4          here_1.0.1          
## 
## loaded via a namespace (and not attached):
##   [1] fastmatch_1.1-3        plyr_1.8.6             igraph_1.2.11         
##   [4] lazyeval_0.2.2         splines_4.1.2          BiocParallel_1.24.1   
##   [7] listenv_0.8.0          SnowballC_0.7.0        scattermore_0.7       
##  [10] digest_0.6.29          htmltools_0.5.2        fansi_1.0.2           
##  [13] tensor_1.5             cluster_2.1.2          ROCR_1.0-11           
##  [16] tzdb_0.2.0             globals_0.14.0         Biostrings_2.58.0     
##  [19] matrixStats_0.61.0     docopt_0.7.1           colorspace_2.0-2      
##  [22] xfun_0.29              sparsesvd_0.2          crayon_1.4.2          
##  [25] RCurl_1.98-1.5         jsonlite_1.7.3         spatstat_1.64-1       
##  [28] spatstat.data_2.1-2    survival_3.2-13        zoo_1.8-9             
##  [31] polyclip_1.10-0        gtable_0.3.0           zlibbioc_1.36.0       
##  [34] XVector_0.30.0         leiden_0.3.9           future.apply_1.8.1    
##  [37] abind_1.4-5            scales_1.1.1           DBI_1.1.2             
##  [40] miniUI_0.1.1.1         Rcpp_1.0.8             xtable_1.8-4          
##  [43] reticulate_1.23        htmlwidgets_1.5.4      httr_1.4.2            
##  [46] ellipsis_0.3.2         ica_1.0-2              farver_2.1.0          
##  [49] pkgconfig_2.0.3        ggseqlogo_0.1          sass_0.4.0            
##  [52] uwot_0.1.11            deldir_1.0-6           utf8_1.2.2            
##  [55] labeling_0.4.2         tidyselect_1.1.1       rlang_0.4.12          
##  [58] reshape2_1.4.4         later_1.3.0            munsell_0.5.0         
##  [61] cellranger_1.1.0       tools_4.1.2            generics_0.1.1        
##  [64] ggridges_0.5.3         evaluate_0.14          stringr_1.4.0         
##  [67] fastmap_1.1.0          yaml_2.2.1             goftest_1.2-3         
##  [70] knitr_1.37             fitdistrplus_1.1-6     RANN_2.6.1            
##  [73] pbapply_1.5-0          future_1.23.0          nlme_3.1-153          
##  [76] mime_0.12              slam_0.1-50            RcppRoll_0.3.0        
##  [79] compiler_4.1.2         beeswarm_0.4.0         plotly_4.10.0         
##  [82] png_0.1-7              spatstat.utils_2.3-0   tweenr_1.0.2          
##  [85] bslib_0.3.1            stringi_1.7.6          highr_0.9             
##  [88] lattice_0.20-45        Matrix_1.3-4           vctrs_0.3.8           
##  [91] pillar_1.6.4           lifecycle_1.0.1        BiocManager_1.30.15   
##  [94] lmtest_0.9-39          jquerylib_0.1.4        RcppAnnoy_0.0.19      
##  [97] bitops_1.0-7           irlba_2.3.5            httpuv_1.6.5          
## [100] patchwork_1.1.1        R6_2.5.1               promises_1.2.0.1      
## [103] renv_0.15.5            lsa_0.73.2             KernSmooth_2.23-20    
## [106] gridExtra_2.3          vipor_0.4.5            parallelly_1.30.0     
## [109] codetools_0.2-18       MASS_7.3-54            assertthat_0.2.1      
## [112] rprojroot_2.0.2        withr_2.4.3            qlcMatrix_0.9.7       
## [115] sctransform_0.3.3      Rsamtools_2.6.0        GenomeInfoDbData_1.2.4
## [118] mgcv_1.8-38            hms_1.1.1              grid_4.1.2            
## [121] rpart_4.1-15           rmarkdown_2.11         Cairo_1.5-14          
## [124] Rtsne_0.15             git2r_0.29.0           ggforce_0.3.3         
## [127] shiny_1.7.1            ggbeeswarm_0.6.0
```

</details>

The resources requested when this document was last rendered:

```
## #SBATCH --time=00:20:00
## #SBATCH --cpus-per-task=1
## #SBATCH --mem=60G
```


***

<!-- END OF END MATTER -->
