---
title: "01C - Confirmation of malignant ependymal cells"
author: "Selin Jessa [[selin.jessa@mail.mcgill.ca](mailto:selin.jessa@mail.mcgill.ca)]"
date: "08 September, 2022"
params:
  resources: "NOT SPECIFIED"
output:
  html_document:
    keep_md: yes
    code_folding: show
    theme: flatly
    css: ../include/style.css
    toc: yes
    toc_depth: 4
    number_sections: true
    df_print: paged
    includes:
      before_body: ../include/header.html
      after_body:  ../include/footer.html
---

<!-- FRONT MATTER, insert configuration info -->




<!-- Load custom CSS/JS for code folding -->
<link rel="stylesheet" type="text/css" href="../include/hideOutput.css">
<script src="../include/hideOutput.js"></script>

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
## Document index: 01C
```

```r
# Specify where to save outputs
out        <- here("output", doc_id); dir.create(out, recursive = TRUE)
figout     <- here("figures", doc_id); dir.create(figout, recursive = TRUE)
cache      <- paste0(readLines(here("include/project_root.txt")), basename(here()), "/", doc_id, "/")
```

</details>

Outputs and figures will be saved at these paths, relative to project root:


```
## public/output/01C
```

```
## public/figures/01C
```



Setting a random seed:



```r
set.seed(100)
```



***



<!-- END OF FRONT MATTER -->


# Overview

In this document, we confirm the identity of the cells projected to be
ependymal-like based on the automated cell-type projections. This is done by scoring
cells for other features (e.g. expression of the ependymal transcription FOXJ1 and its targets, re-examining the cNMF programs, and re-examining the inferCNV output).

# Libraries



```r
# Load libraries here
library(biomaRt)
library(here)
library(tidyr)
library(dplyr)
library(ggrepel)
library(readr)
library(readxl)
library(glue)
library(tibble)
library(ggplot2)
library(purrr)
library(cowplot)
library(Seurat)
library(ggrastr)
library(pheatmap)

source(here("include/style.R"))
source(here("code/functions/scRNAseq.R"))
source(here("code/functions/ssGSEA.R"))
source(here("code/functions/RNAseq.R"))
ggplot2::theme_set(theme_min())
```




# Load data

In this analysis, we specifically load samples with high ependymal content:



```r
# select samples to interrogate further
high_epen <- c("P-6292_S-8579", # PFA
               "P-6251_S-8496", # H3.1K27M
               "P-1713_S-1713", # H3.1K27M
               "P-6255_S-8500") # H3.1K27M
```



Load data and add cell-level labels for visualization:

* `Type_malig`: column labels cells as "Normal" if normal based on malignant/normal consensus calling, or their most similar cell type if Malignant
* `Type_malig_inferCNV`: column labels cells as "Normal" if normal based on inferCNV calling, or their most similar cell type if Malignant
* `Ependymal_binary` label: contains "Ependymal" or "Other" based on correlation cell type annotation




```r
seurat_list <- list()

for (i in seq_along(high_epen)) {
    
    seurat <- get(load(here(
        file.path("data/scRNAseq/pipeline_10X/", high_epen[i], "seurat.Rda")
    )))
    
    seurat$Type_malig <- seurat@meta.data %>%
        mutate(Type_malig = case_when(
            Malignant_normal_consensus_Jessa2022 == "Normal" ~ "Normal",
            TRUE ~ as.character(Cell_type_consensus_Jessa2022)
        ),
        Type_malig = factor(Type_malig, levels = c(names(palette_type), "Uncertain"))) %>% 
        pull(Type_malig)
    
    seurat$Type_malig_inferCNV <- seurat@meta.data %>%
        dplyr::select(COR_ref.joint_mouse_extended, inferCNV) %>%
        summarize_cell_types("COR_ref.joint_mouse_extended") %>% 
        mutate(Type_malig_inferCNV = case_when(
            inferCNV == "Normal" ~ "Normal",
            TRUE ~ as.character(Type)
        ),
        Type_malig_inferCNV = factor(Type_malig_inferCNV, levels = names(palette_type))) %>% 
        pull(Type_malig_inferCNV)
    
    seurat$Ependymal_binary <- seurat@meta.data %>%
        dplyr::select(COR_ref.joint_mouse_extended, inferCNV) %>%
        summarize_cell_types("COR_ref.joint_mouse_extended") %>% 
        mutate(Epen = ifelse(Type == "Ependymal", "Ependymal", "Other")) %>% 
        pull(Epen)
    
    seurat_list[[i]] <- seurat
    
}

names(seurat_list) <- high_epen
```




# Cell-type projections

Examine the cell-type projections in the space of each individual sample:

## By correlations

First, plot the cell type projections for malignant cells within each sample:



```r
imap(seurat_list, function(object, id) {
    
    p1 <- FetchData(subset(object, subset = Malignant_normal_consensus_Jessa2022 == "Malignant"),
                    c("UMAP_1", "UMAP_2", "Type_malig")) %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2)) +
        rasterize(geom_point(aes(colour = Type_malig), size = 0.5, alpha = 0.5), dpi = 500) +
        scale_colour_manual(values = c(palette_type, "Uncertain" = "gray90")) +
        theme_min() +
        no_legend() +
        no_ticks() +
        ggtitle(id)
    
    return(p1)
    
}) %>%
    {plot_grid(plotlist = ., ncol = 4, align = "h", axis = "tb")}
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/01C/projections-1.png)<!-- -->

```r
# based on the inferCNV call only
imap(seurat_list, function(object, id) {
    
    p1 <- FetchData(object, c("UMAP_1", "UMAP_2", "Type_malig_inferCNV")) %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2)) +
        rasterize(geom_point(aes(colour = Type_malig_inferCNV), size = 0.5, alpha = 0.5), dpi = 500) +
        scale_colour_manual(values = palette_type) +
        theme_min() +
        no_legend() +
        no_ticks() +
        ggtitle(id)
    
    return(p1)
    
}) %>%
    {plot_grid(plotlist = ., ncol = 4, align = "h", axis = "tb")}
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/01C/projections-2.png)<!-- -->



# Gene expression

Next, let's examine genes canonically or functionally associated with ependymal cells:



```r
feature_plots_foxj1 <- imap(seurat_list, function(object, id) {
    
    p1 <- FetchData(object, c("UMAP_1", "UMAP_2", "FOXJ1")) %>% 
        arrange(FOXJ1) %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2)) +
        rasterize(geom_point(aes(colour = FOXJ1), size = 0.5, alpha = 0.5), dpi = 500) +
        scale_colour_gradientn(colours = gryrd) +
        theme_min() +
        theme(legend.position = "bottom") +
        no_ticks() +
        ggtitle(id)
    
    return(p1)
    
})

plot_grid(plotlist = feature_plots_foxj1, ncol = 4, align = "h", axis = "tb")
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/01C/feature_plots-1.png)<!-- -->

```r
feature_plots_dnah <- imap(seurat_list, function(object, id) {
    
    p1 <- FetchData(object, c("UMAP_1", "UMAP_2", "DNAH12")) %>% 
        arrange(DNAH12) %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2)) +
        rasterize(geom_point(aes(colour = DNAH12), size = 0.5, alpha = 0.5), dpi = 500) +
        scale_colour_gradientn(colours = gryrd) +
        theme_min() +
        theme(legend.position = "bottom") +
        no_ticks() +
        ggtitle(id)
    
    return(p1)
    
})

plot_grid(plotlist = feature_plots_dnah, ncol = 4, align = "h", axis = "tb")
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/01C/feature_plots-2.png)<!-- -->


# FOXJ1 targets

We can also directly assess the expression/enrichment of FOXJ1 target genes in these tumors.

Load a list of Foxj1 targets from [Jacquet et al, 2009](https://pubmed.ncbi.nlm.nih.gov/19906869/). We convert these to human genes using [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html).



```r
foxj1_targets <- read_xlsx(here("data/misc/Jacquet2009_TableS1_Foxj1_targets.xlsx"))
head(foxj1_targets)

table(foxj1_targets$Category)

human  <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse  <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes2 <- getLDS(attributes = c("mgi_symbol"),
                 filters = "mgi_symbol",
                 values = foxj1_targets$Gene ,
                 mart = mouse,
                 attributesL = c("hgnc_symbol"),
                 martL = human,
                 uniqueRows = TRUE)

foxj1_targets_hg <- unique(genes2[, 2])
saveRDS(foxj1_targets_hg, file = glue("{out}/FOXJ1_targets_hg.Rds"))
```





```r
foxj1_targets_hg <- readRDS(glue("{out}/FOXJ1_targets_hg.Rds"))
```



## Gene set enrichment

For each sample, compute the ssGSEA enrichment of this target in the single cells:



```r
seurat_list <- map(seurat_list, function(seurat) {
    
    # run ssGSEA, then extract the vector of scores
    scores <- compute_scores_sc(seurat, genelists = list("FOXJ1_targets" = foxj1_targets_hg)) %>% 
        .[, -1] %>% t() %>% .[, 1]
    
    seurat$FOXJ1_targets <- scores
    return(seurat)
    
})
```



Plot results:



```r
feature_plots_foxj1_targets <- imap(seurat_list, function(object, id) {
    
    p1 <- FetchData(object, c("UMAP_1", "UMAP_2", "FOXJ1_targets")) %>% 
        arrange(FOXJ1_targets) %>% 
        ggplot(aes(x = UMAP_1, y = UMAP_2)) +
        rasterize(geom_point(aes(colour = FOXJ1_targets), size = 0.5, alpha = 0.5), dpi = 500) +
        scale_colour_gradientn(colours = rdbu) +
        theme_min() +
        theme(legend.position = "bottom") +
        no_ticks() +
        ggtitle(id)
    
    return(p1)
    
})

plot_grid(plotlist = feature_plots_foxj1_targets, ncol = 4, align = "h", axis = "tb")
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/01C/foxj1_ssgsea_plot-1.png)<!-- -->


# cNMF

In the two H3.1K27M samples with high ependymal content, cNMF independently
identified a gene program among malignant cells which is most highly
activated in these ependymal-like cells:



```r
p1 <- FetchData(seurat_list[[2]], c("UMAP_1", "UMAP_2", "cNMF_program_malignant_1")) %>%
    dplyr::rename(NMF_score = cNMF_program_malignant_1) %>%
    arrange(NMF_score) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2)) +
    rasterize(geom_point(aes(colour = NMF_score), size = 0.5, alpha = 0.5), dpi = 500) +
    scale_colour_gradientn(colours = rdbu, na.value = "gray90") +
    theme_min() +
    theme(legend.position = "bottom") +
    no_ticks()

p2 <- FetchData(seurat_list[[3]], c("UMAP_1", "UMAP_2", "cNMF_program_malignant_3")) %>%
    dplyr::rename(NMF_score = cNMF_program_malignant_3) %>%
    arrange(NMF_score) %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2)) +
    rasterize(geom_point(aes(colour = NMF_score), size = 0.5, alpha = 0.5), dpi = 500) +
    scale_colour_gradientn(colours = rdbu, na.value = "gray90") +
    theme_min() +
    theme(legend.position = "bottom") +
    no_ticks()

plot_grid(p1, p2, ncol = 2, align = "h", axis = "tb")
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/01C/cnmf-1.png)<!-- -->

For each of those, we can also load the top genes for that program and plot their scores:



```r
cnmf_top_genes <- readRDS(here("R-4/output/01/cNMF_top_genes.malignant.Rds"))

plot_gene_scores <- function(scores, ...) {

    scores %>%
        tibble::column_to_rownames(var = "Program") %>%
        t() %>%
        set_colnames(1:ncol(.)) %>%
        pheatmap(color = head(rdbu, 80),
                 border_color = NA,
                 scale = "row",
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 show_colnames = TRUE,
                 show_rownames = TRUE,
                 fontsize_row = 5,
                 fontsize_col = 10,
                 cellwidth = 15,
                 cellheight = 5,
                 ...)

}
```



<details>

For the first sample:



```r
gene_scores1 <- get_cnmf_gene_scores(
    path_to_sample = here(file.path("data/scRNAseq/pipeline_10X/", "P-6251_S-8496")),
    output_dir = "output_ngenes2000_niter100_malignant",
    program_string = "_cNMF_program_malignant_",
    genes_keep = cnmf_top_genes[["P-6251_S-8496_cNMF_program_malignant_1"]])

# large version with gene names; display that one
plot_gene_scores(gene_scores1, filename = glue("{figout}/cNMF_P-6251_S-8496_hm.png"))

plot_gene_scores(gene_scores1, filename = glue("{figout}/cNMF_P-6251_S-8496_hm.pdf"))

knitr::include_graphics(glue("{figout}/cNMF_P-6251_S-8496_hm.png"))
```

<img src="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/01C/cNMF_P-6251_S-8496_hm.png" width="725" />

For the second sample:



```r
gene_scores2 <- get_cnmf_gene_scores(
    path_to_sample = here(file.path("data/scRNAseq/pipeline_10X/", "P-1713_S-1713")),
    output_dir = "output_ngenes2000_niter100_malignant",
    program_string = "_cNMF_program_malignant_",
    genes_keep = cnmf_top_genes[["P-1713_S-1713_cNMF_program_malignant_3"]])

plot_gene_scores(gene_scores2, filename = glue("{figout}/cNMF_P-1713_S-1713_hm.png"))

plot_gene_scores(gene_scores2, filename = glue("{figout}/cNMF_P-1713_S-1713_hm.pdf"))

knitr::include_graphics(glue("{figout}/cNMF_P-1713_S-1713_hm.png"))
```

<img src="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/01C/cNMF_P-1713_S-1713_hm.png" width="740" />

</details>

# Copy-number profiles

Load gene annotation used for CNV calling:



```r
gene_annotation <- read_tsv(
    here("data/scRNAseq/references/inferCNV/gencode_v19_gen_pos_noHLAs_noMitoRibo.txt"),
    col_names = c("gene", "contig", "start", "end"))
```



Generate a function that can generate the inferCNV heatmap:

<details>



```r
#' Plot CNV profile based on inferCNV
#'
#' This was adapted/refactored from code from Samantha Worme's function plotSampleCNVs()
#'
#' @param out Character, name of the output directory
plot_profile <- function(seurat, infercnv_out, outfile, cells, cluster_rows = TRUE) {
    
    # prep and loading -----------------------------------------------------------
    message("@ custom heatmap: preparing data...")
    
    # load gene x cell matrix
    obs <- read.table(glue("{infercnv_out}/infercnv.observations.txt"))
    
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
    # as well as genes marking a ependymal populations:
    # FOXJ1, DNAH6
    row_annotations <- data.frame("cell"      = colnames(seurat),
                                  "celltype"  = as.character(seurat$Type_malig_inferCNV),
                                  "ependymal" = as.character(seurat$Ependymal_binary),
                                  stringsAsFactors = FALSE)
    
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
    # palette_clusters <- seurat@misc$colours
    
    # for expression of specific genes, we'll use a red-gray palette
    # expr_palette <- grDevices::colorRampPalette(c("gray90", "#E09797", "red"))(n = 200)
    
    anno_palettes = list(contig    = palette_contig,
                         celltype  = palette_type,
                         ependymal = c("Ependymal" = "red", "Other" = "gray90"))
    
    # for the heatmap itself, generate a red-blue heatmap (this line returns a function)
    hm_palette_generator <- infercnv:::color.palette(c("darkblue", "white", "darkred"), c(2, 2))
    
    message("@ custom heatmap: generating heatmap...")
    
    # label chromosomes by repurposing the column labels: generate a vector which will
    # mostly be empty, but for one index in each contig, we'll put the chromosome labels
    # for the first few chromosomes, shift the label slightly to the left so it's more centered...
    labels_col <- rep("", ncol(obs_t))
    labels_col[c(last_gene_idx[1:7] - 200,
                 last_gene_idx[8:length(last_gene_idx)])] <- unique(contigs$contig)
    
    # generate the heatmap -------------------------------------------------------
    hm_fun <- purrr::partial(pheatmap,
                             mat = obs_t[cells, ],
                             scale             = "none",
                             border_color      = NA,
                             cluster_rows      = cluster_rows,
                             cluster_cols      = FALSE,
                             clustering_method = "ward.D",
                             labels_col        = labels_col,
                             show_colnames     = TRUE,
                             show_rownames     = FALSE,
                             color             = hm_palette_generator(100),
                             breaks            = seq(0.8, 1.2, length.out = 101),
                             annotation_colors = anno_palettes,
                             annotation_col    = contigs,
                             annotation_row    = row_annotations,
                             gaps_col          = last_gene_idx,
                             treeheight_row    = 70, # make the dendrogram a little bit taller
                             main              = glue("inferCNV"),
                             annotation_legend = FALSE,
                             fontsize_col      = 7,
                             width             = 8,
                             height            = 5)
    
    hm_fun(filename = glue("{outfile}.png"))
    hm_fun(filename = glue("{outfile}.pdf"))
    
}
```



</details>

For this analysis, we focus on the two H3.1K27M samples with high proportions
of malignant cells. We'll plot separately the CNV profiles of normal cells,
malignant cells projected to the ependymal class, and malginant cells with
other identities.



```r
s1_infercnv_out <- here("data/scRNAseq/pipeline_10X/P-6251_S-8496/inferCNV/window101_exp0.1_refG34normcleaned_HMMnone/")

# correct the cell names to match the transformation by inferCNV
fix_names <- function(x) gsub("-", ".", x, fixed = TRUE)

plot_profile(seurat       = seurat_list[[2]],
             infercnv_out = s1_infercnv_out,
             cells        = WhichCells(seurat_list[[2]],
                                       expression = Cell_type_consensus_Jessa2022 != "Uncertain") %>% fix_names(),
             outfile      = glue("{figout}/P-6251_S-8496_CNV_all"))
```



```
## @ custom heatmap: preparing data...
```



```
## @ custom heatmap: producing annotations...
```



```
## @ custom heatmap: generating heatmap...
```



<img src="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/01C/P-6251_S-8496_CNV_all.png" width="2400" />

Repeat for the second sample:



```r
s2_infercnv_out <- here("data/scRNAseq/pipeline_10X/P-1713_S-1713/inferCNV/window101_exp0.1_refG34normcleaned_HMMnone/")

plot_profile(seurat       = seurat_list[[3]],
             infercnv_out = s2_infercnv_out,
             cells        = WhichCells(seurat_list[[3]],
                                       expression = Cell_type_consensus_Jessa2022 != "Uncertain") %>% fix_names(),
             outfile      = glue("{figout}/P-1713_S-1713_CNV_all"))
```



```
## @ custom heatmap: preparing data...
```



```
## @ custom heatmap: producing annotations...
```



```
## @ custom heatmap: generating heatmap...
```



<img src="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/01C/P-1713_S-1713_CNV_all.png" width="2400" />


<!-- END MATTER, insert reproducibility info -->




***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:



```
## 2022-09-08 15:01:00
```



The git repository and last commit:



```
## Local:    master /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public
## Remote:   master @ origin (git@github.com:fungenomics/HGG-oncohistones.git)
## Head:     [4101e76] 2022-09-08: Update README.md
```



The random seed was set with `set.seed(100)`

The R session info:
<details>



```
## Registered S3 method overwritten by 'cli':
##   method     from    
##   print.boxx spatstat
```



```
## Error in get(genname, envir = envir) : object 'testthat_print' not found
```



```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value                           
##  version  R version 3.6.1 (2019-07-05)    
##  os       Rocky Linux 8.6 (Green Obsidian)
##  system   x86_64, linux-gnu               
##  ui       X11                             
##  language (EN)                            
##  collate  en_CA.UTF-8                     
##  ctype    en_CA.UTF-8                     
##  tz       EST5EDT                         
##  date     2022-09-08                      
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  ! package        * version    date       lib source                           
##  P abind            1.4-5      2016-07-21 [?] CRAN (R 3.6.1)                   
##  P AnnotationDbi    1.48.0     2019-10-29 [?] Bioconductor                     
##  P askpass          1.1        2019-01-13 [?] CRAN (R 3.6.1)                   
##  P assertthat       0.2.1      2019-03-21 [?] CRAN (R 3.6.1)                   
##  P beeswarm         0.3.1      2021-03-07 [?] CRAN (R 3.6.1)                   
##  P Biobase          2.46.0     2019-10-29 [?] Bioconductor                     
##  P BiocFileCache    1.10.2     2019-11-08 [?] Bioconductor                     
##  P BiocGenerics     0.32.0     2019-10-29 [?] Bioconductor                     
##  P biomaRt        * 2.42.1     2020-03-26 [?] Bioconductor                     
##  P bit              4.0.4      2020-08-04 [?] CRAN (R 3.6.1)                   
##  P bit64            4.0.5      2020-08-30 [?] CRAN (R 3.6.1)                   
##  P blob             1.2.1      2020-01-20 [?] CRAN (R 3.6.1)                   
##  P bslib            0.2.5      2021-05-12 [?] CRAN (R 3.6.1)                   
##  P callr            3.7.0      2021-04-20 [?] CRAN (R 3.6.1)                   
##  P cellranger       1.1.0      2016-07-27 [?] CRAN (R 3.6.1)                   
##  P cli              2.5.0      2021-04-26 [?] CRAN (R 3.6.1)                   
##  P cluster          2.1.0      2019-06-19 [?] CRAN (R 3.6.1)                   
##  P codetools        0.2-16     2018-12-24 [?] CRAN (R 3.6.1)                   
##  P colorspace       2.0-1      2021-05-04 [?] CRAN (R 3.6.1)                   
##  P cowplot        * 1.1.1      2020-12-30 [?] CRAN (R 3.6.1)                   
##  P crayon           1.4.1      2021-02-08 [?] CRAN (R 3.6.1)                   
##  P curl             4.3.1      2021-04-30 [?] CRAN (R 3.6.1)                   
##  P data.table       1.14.0     2021-02-21 [?] CRAN (R 3.6.1)                   
##  P DBI              1.1.1      2021-01-15 [?] CRAN (R 3.6.1)                   
##  P dbplyr           2.1.1      2021-04-06 [?] CRAN (R 3.6.1)                   
##  P deldir           0.2-10     2021-02-16 [?] CRAN (R 3.6.1)                   
##  P desc             1.2.0      2018-05-01 [?] CRAN (R 3.6.1)                   
##  P devtools         2.3.0      2020-04-10 [?] CRAN (R 3.6.1)                   
##  P digest           0.6.27     2020-10-24 [?] CRAN (R 3.6.1)                   
##  P dplyr          * 1.0.6      2021-05-05 [?] CRAN (R 3.6.1)                   
##  P ellipsis         0.3.2      2021-04-29 [?] CRAN (R 3.6.1)                   
##  P evaluate         0.14       2019-05-28 [?] CRAN (R 3.6.1)                   
##  P fansi            0.4.2      2021-01-15 [?] CRAN (R 3.6.1)                   
##  P fastmap          1.1.0      2021-01-25 [?] CRAN (R 3.6.1)                   
##  P fitdistrplus     1.1-3      2020-12-05 [?] CRAN (R 3.6.1)                   
##  P fs               1.5.0      2020-07-31 [?] CRAN (R 3.6.1)                   
##  P future           1.21.0     2020-12-10 [?] CRAN (R 3.6.1)                   
##  P future.apply     1.7.0      2021-01-04 [?] CRAN (R 3.6.1)                   
##  P generics         0.1.0      2020-10-31 [?] CRAN (R 3.6.1)                   
##  P ggbeeswarm       0.6.0      2017-08-07 [?] CRAN (R 3.6.1)                   
##  P ggplot2        * 3.3.3      2020-12-30 [?] CRAN (R 3.6.1)                   
##  P ggrastr        * 0.2.3      2021-03-01 [?] CRAN (R 3.6.1)                   
##  P ggrepel        * 0.9.1      2021-01-15 [?] CRAN (R 3.6.1)                   
##  P ggridges         0.5.3      2021-01-08 [?] CRAN (R 3.6.1)                   
##  P git2r            0.27.1     2020-05-03 [?] CRAN (R 3.6.1)                   
##  P globals          0.14.0     2020-11-22 [?] CRAN (R 3.6.1)                   
##  P glue           * 1.4.2      2020-08-27 [?] CRAN (R 3.6.1)                   
##  P goftest          1.2-2      2019-12-02 [?] CRAN (R 3.6.1)                   
##  P gridExtra        2.3        2017-09-09 [?] CRAN (R 3.6.1)                   
##  P gtable           0.3.0      2019-03-25 [?] CRAN (R 3.6.1)                   
##  P here           * 0.1        2017-05-28 [?] CRAN (R 3.6.1)                   
##  P highr            0.9        2021-04-16 [?] CRAN (R 3.6.1)                   
##  P hms              1.0.0      2021-01-13 [?] CRAN (R 3.6.1)                   
##  P htmltools        0.5.1.1    2021-01-22 [?] CRAN (R 3.6.1)                   
##  P htmlwidgets      1.5.3      2020-12-10 [?] CRAN (R 3.6.1)                   
##  P httpuv           1.6.1      2021-05-07 [?] CRAN (R 3.6.1)                   
##  P httr             1.4.2      2020-07-20 [?] CRAN (R 3.6.1)                   
##  P ica              1.0-2      2018-05-24 [?] CRAN (R 3.6.1)                   
##  P igraph           1.2.6      2020-10-06 [?] CRAN (R 3.6.1)                   
##  P IRanges          2.20.2     2020-01-13 [?] Bioconductor                     
##  P irlba            2.3.3      2019-02-05 [?] CRAN (R 3.6.1)                   
##  P jquerylib        0.1.4      2021-04-26 [?] CRAN (R 3.6.1)                   
##  P jsonlite         1.7.2      2020-12-09 [?] CRAN (R 3.6.1)                   
##  P KernSmooth       2.23-15    2015-06-29 [?] CRAN (R 3.6.1)                   
##  P knitr            1.33       2021-04-24 [?] CRAN (R 3.6.1)                   
##  P later            1.0.0      2019-10-04 [?] CRAN (R 3.6.1)                   
##  P lattice          0.20-44    2021-05-02 [?] CRAN (R 3.6.1)                   
##  P lazyeval         0.2.2      2019-03-15 [?] CRAN (R 3.6.1)                   
##  P leiden           0.3.7      2021-01-26 [?] CRAN (R 3.6.1)                   
##  P lifecycle        1.0.0      2021-02-15 [?] CRAN (R 3.6.1)                   
##  P listenv          0.8.0      2019-12-05 [?] CRAN (R 3.6.1)                   
##  P lmtest           0.9-38     2020-09-09 [?] CRAN (R 3.6.1)                   
##  P magrittr       * 2.0.1      2020-11-17 [?] CRAN (R 3.6.1)                   
##  P MASS             7.3-54     2021-05-03 [?] CRAN (R 3.6.1)                   
##  P Matrix           1.2-18     2019-11-27 [?] CRAN (R 3.6.1)                   
##  P matrixStats      0.58.0     2021-01-29 [?] CRAN (R 3.6.1)                   
##  P memoise          1.1.0      2017-04-21 [?] CRAN (R 3.6.1)                   
##  P mgcv             1.8-35     2021-04-18 [?] CRAN (R 3.6.1)                   
##  P mime             0.10       2021-02-13 [?] CRAN (R 3.6.1)                   
##  P miniUI           0.1.1.1    2018-05-18 [?] CRAN (R 3.6.1)                   
##  P munsell          0.5.0      2018-06-12 [?] CRAN (R 3.6.1)                   
##  P nlme             3.1-152    2021-02-04 [?] CRAN (R 3.6.1)                   
##  P openssl          1.4.4      2021-04-30 [?] CRAN (R 3.6.1)                   
##  P parallelly       1.25.0     2021-04-30 [?] CRAN (R 3.6.1)                   
##  P patchwork        1.1.1      2020-12-17 [?] CRAN (R 3.6.1)                   
##  P pbapply        * 1.4-3      2020-08-18 [?] CRAN (R 3.6.1)                   
##  P pheatmap       * 1.0.12     2019-01-04 [?] CRAN (R 3.6.1)                   
##  P pillar           1.6.0      2021-04-13 [?] CRAN (R 3.6.1)                   
##  P pkgbuild         1.0.8      2020-05-07 [?] CRAN (R 3.6.1)                   
##  P pkgconfig        2.0.3      2019-09-22 [?] CRAN (R 3.6.1)                   
##  P pkgload          1.0.2      2018-10-29 [?] CRAN (R 3.6.1)                   
##  P plotly           4.9.3      2021-01-10 [?] CRAN (R 3.6.1)                   
##  P plyr             1.8.6      2020-03-03 [?] CRAN (R 3.6.1)                   
##  P png              0.1-7      2013-12-03 [?] CRAN (R 3.6.1)                   
##  P polyclip         1.10-0     2019-03-14 [?] CRAN (R 3.6.1)                   
##  P prettyunits      1.1.1      2020-01-24 [?] CRAN (R 3.6.1)                   
##  P processx         3.5.2      2021-04-30 [?] CRAN (R 3.6.1)                   
##  P progress         1.2.2      2019-05-16 [?] CRAN (R 3.6.1)                   
##  P promises         1.1.0      2019-10-04 [?] CRAN (R 3.6.1)                   
##  P ps               1.6.0      2021-02-28 [?] CRAN (R 3.6.1)                   
##  P purrr          * 0.3.4      2020-04-17 [?] CRAN (R 3.6.1)                   
##  P R6               2.5.0      2020-10-28 [?] CRAN (R 3.6.1)                   
##  P RANN             2.6.1      2019-01-08 [?] CRAN (R 3.6.1)                   
##  P rappdirs         0.3.3      2021-01-31 [?] CRAN (R 3.6.1)                   
##  P RColorBrewer   * 1.1-2      2014-12-07 [?] CRAN (R 3.6.1)                   
##  P Rcpp             1.0.6      2021-01-15 [?] CRAN (R 3.6.1)                   
##  P RcppAnnoy        0.0.18     2020-12-15 [?] CRAN (R 3.6.1)                   
##  P readr          * 1.4.0      2020-10-05 [?] CRAN (R 3.6.1)                   
##  P readxl         * 1.3.1      2019-03-13 [?] CRAN (R 3.6.1)                   
##  P remotes          2.1.1      2020-02-15 [?] CRAN (R 3.6.1)                   
##    renv             0.14.0     2021-07-21 [1] CRAN (R 3.6.1)                   
##  P reshape2         1.4.4      2020-04-09 [?] CRAN (R 3.6.1)                   
##  P reticulate       1.20       2021-05-03 [?] CRAN (R 3.6.1)                   
##  P rlang            0.4.11     2021-04-30 [?] CRAN (R 3.6.1)                   
##  P rmarkdown        2.8        2021-05-07 [?] CRAN (R 3.6.1)                   
##  P ROCR             1.0-11     2020-05-02 [?] CRAN (R 3.6.1)                   
##  P rpart            4.1-15     2019-04-12 [?] CRAN (R 3.6.1)                   
##  P rprojroot        2.0.2      2020-11-15 [?] CRAN (R 3.6.1)                   
##  P RSQLite          2.2.1      2020-09-30 [?] CRAN (R 3.6.1)                   
##  P rsvd             1.0.3      2020-02-17 [?] CRAN (R 3.6.1)                   
##  P Rtsne            0.15       2018-11-10 [?] CRAN (R 3.6.1)                   
##  P S4Vectors        0.24.4     2020-04-09 [?] Bioconductor                     
##  P sass             0.4.0      2021-05-12 [?] CRAN (R 3.6.1)                   
##  P scales           1.1.1      2020-05-11 [?] CRAN (R 3.6.1)                   
##  P sctransform      0.3.2      2020-12-16 [?] CRAN (R 3.6.1)                   
##  P sessioninfo      1.1.1      2018-11-05 [?] CRAN (R 3.6.1)                   
##  P Seurat         * 3.2.1      2020-09-07 [?] CRAN (R 3.6.1)                   
##  P shiny            1.6.0      2021-01-25 [?] CRAN (R 3.6.1)                   
##  P spatstat         1.64-1     2020-05-12 [?] CRAN (R 3.6.1)                   
##  P spatstat.data    2.1-0      2021-03-21 [?] CRAN (R 3.6.1)                   
##  P spatstat.utils   2.1-0      2021-03-15 [?] CRAN (R 3.6.1)                   
##  P stringi          1.6.1      2021-05-10 [?] CRAN (R 3.6.1)                   
##  P stringr          1.4.0      2019-02-10 [?] CRAN (R 3.6.1)                   
##  P survival         3.2-11     2021-04-26 [?] CRAN (R 3.6.1)                   
##  P tensor           1.5        2012-05-05 [?] CRAN (R 3.6.1)                   
##  P testrmd          0.0.1.9000 2021-12-06 [?] Github (rmflight/testrmd@0735c20)
##  P testthat         2.3.2      2020-03-02 [?] CRAN (R 3.6.1)                   
##  P tibble         * 3.1.1      2021-04-18 [?] CRAN (R 3.6.1)                   
##  P tidyr          * 1.1.3      2021-03-03 [?] CRAN (R 3.6.1)                   
##  P tidyselect       1.1.1      2021-04-30 [?] CRAN (R 3.6.1)                   
##  P usethis          1.6.1      2020-04-29 [?] CRAN (R 3.6.1)                   
##  P utf8             1.2.1      2021-03-12 [?] CRAN (R 3.6.1)                   
##  P uwot             0.1.10     2020-12-15 [?] CRAN (R 3.6.1)                   
##  P vctrs            0.3.8      2021-04-29 [?] CRAN (R 3.6.1)                   
##  P vipor            0.4.5      2017-03-22 [?] CRAN (R 3.6.1)                   
##  P viridis        * 0.5.1      2018-03-29 [?] CRAN (R 3.6.1)                   
##  P viridisLite    * 0.4.0      2021-04-13 [?] CRAN (R 3.6.1)                   
##  P withr            2.4.2      2021-04-18 [?] CRAN (R 3.6.1)                   
##  P xfun             0.22       2021-03-11 [?] CRAN (R 3.6.1)                   
##  P XML              3.99-0.3   2020-01-20 [?] CRAN (R 3.6.1)                   
##  P xtable           1.8-4      2019-04-21 [?] CRAN (R 3.6.1)                   
##  P yaml             2.2.1      2020-02-01 [?] CRAN (R 3.6.1)                   
##  P zoo              1.8-9      2021-03-09 [?] CRAN (R 3.6.1)                   
## 
## [1] /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/renv/library/R-3.6/x86_64-pc-linux-gnu
## [2] /tmp/Rtmp361k2q/renv-system-library
## 
##  P ── Loaded and on-disk path mismatch.
```



</details>

The resources requested when this document was last rendered:



```
## #SBATCH --time=02:00:00
## #SBATCH --cpus-per-task=1
## #SBATCH --mem=50G
```




***



<!-- END OF END MATTER -->
