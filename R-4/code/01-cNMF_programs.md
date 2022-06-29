---
title: "01 - cNMF"
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
## Document index: 01
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
## public/R-4/output/01
```

```
## public/R-4/figures/01
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

This document focuses on analysis of cell-type identity across samples,
primarily using scRNAseq data. The main analyses are consensus non-negative
matrix factorization (cNMF) to extract & annotate gene programs in an unbiased
manner, and visualization and quantification of cell-type projections, as shown
in Figure 1.

# Libraries


```r
# Load libraries here
library(biomaRt)
library(here)
library(tidyr)
library(dplyr)
library(ggrepel)
library(readr)
library(glue)
library(tibble)
library(ggplot2)
library(purrr)
library(pheatmap)
library(ape)
library(dendextend)
library(ggrastr)
library(cowplot)
library(Seurat)

source(here("include/style.R")) # contains palettes & plotting utils
source(here("code/functions/scRNAseq.R"))
ggplot2::theme_set(theme_min())
```

# Load metadata

Load the sample metadata for the project:



```r
meta      <- read_tsv(here("data/metadata/metadata_patient_samples_NGS.tsv"))
```

```
## Registered S3 method overwritten by 'cli':
##   method     from    
##   print.boxx spatstat
## Rows: 138 Columns: 54
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (48): BioID, IC_Sample, IC_Patient, ID_paper, SC_QC, Type, Source, Sex, ...
## dbl  (2): Age, RING1B
## lgl  (4): Exclude_entirely, Smartseq2, Smartseq2_path, Smartseq2_ID2
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

Single-cell metadata:


```r
meta_sc <- data.table::fread(here("data/metadata/metadata_sc.tsv"), data.table = FALSE)
```


# cNMF

To identify recurrent sources of intratumor variability in gene expression in an unsupervised manner,
we performed non-negative matrix factorization using the [consensus NMF (cNMF)](https://github.com/dylkot/cNMF)
method described in [Kotliar et al, eLife, 2019](https://elifesciences.org/articles/43803).

The main steps are:

1. Run cNMF to identify gene programs which demonstrate intra-tumor variability
within each sample
2. Extract the top genes associated with each program
3. Identify sets of similar programs across samples. These sets are referred to as "modules" in the paper, but "metaprograms" in the code.
4. Annotate programs based on overlap with other gene sets
5. Extract gene signatures for each metaprogram
 

## Running cNMF

Briefly, for each value of _k_, the number of components, this method runs 100 iterations of
NMF with different random seeds, clusters the components resulting from each replicate,
filters outlier components, and takes the median of each cluster of components as a
consensus estimate for that component. cNMF is applied to raw UMI counts for malignant 
cells and run with values of _k_ from 5-9. For each value of _k_, the Silhouette score,
measuring the stability of the components, and the Frobenius error are computed,
and the _k_ maximizing the Silhouette score and minimizing the Frobenius error was
selected for each sample. Outlier components are filtered by retaining only
components with mean distance to most similar components of 0.02 (`density_threshold = 0.02`),
resulting in a program activity matrix (the activity of each program in each cell), and
a gene scores matrix (reflecting the expected increase in transcripts per million of a
given gene for a unit increase of a given program), which is z-scored across genes.

cNMF is run for each individual sample on malignant cells only. This is performed in the
scRNA pipeline and scMultiome pipelines at `data/scRNAseq/pipeline_10X/<sample>/cNMF`
and `R-4/data/scMultiome/pipeline_10X_Multiome/<sample>/cNMF`. This does three steps:

1. Extract counts matrix, normalize and factorize the matrix, and combine
results from different iterations. Example script at `code/scripts/run_cNMF.sh`.
2. Select k (# of components) by maximizing stability and minimizing error, inspecting
the file `output_ngenes2000_niter100_malignant.k_selection.png` in the cNMF output folder.
This is done interactively (commands are described in the code at the bottom of the `code/scripts/run_cNMF.sh`
script)/
3. Add cNMF scores to Seurat object. Example script at `code/scripts/explore_cNMF.Rmd`.

Since cNMF has been run for each sample within the pipeline, this document and
the following sections mainly load the cNMF output and aggregate these to identify
and visualize recurrent sources of gene expression variation across tumors.

## Program/QC correlations

Here, we calcuate the correlation between the usage score for each program in each cell,
and QC metrics in each cell.

We loop over Seurat objects and save only the dataframe containing QC/cNMF program correlations:


```r
# paths to pipelines
rna_samples <- c(list.files(here("data/scRNAseq/pipeline_10X/"), full.names = TRUE))
rna_samples <- rna_samples[!grepl("Makefile", rna_samples)]

multi_samples <- list.files(here("R-4/data/scMultiome/pipeline_10X_Multiome/"), full.names = TRUE)
# exclude scMultiome samples that have already been profiled by scRNAseq to avoid
# breaking assumptions about independence of samples
multi_samples <- multi_samples[!grepl("Makefile|P-6253_S-8498|P-6640_S-9581|P-1764_S-1766|P-6337_S-8821|P-1709_S-1709", multi_samples)]

# combine
sc_samples <- c(rna_samples, multi_samples)

# SLOW, since it requires loading each Seurat object individually
get_qc_cnmf_correlations <- function(program_string, malignant_only) {
    
    purrr::map_dfr(sc_samples, function(i) {
        
        message("@ ", basename(i))
        
        id <- basename(i)
        load(glue("{i}/seurat.Rda"))
        
        # qc_cols will be the same across samples, cnmf_cols will vary with the selcted K
        qc_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo")
        cc_cols <- c("S.Score", "G2M.Score")
        cnmf_cols <- colnames(seurat@meta.data)[grepl(program_string, colnames(seurat@meta.data))]
        
        if (malignant_only) cells_keep <- seurat@meta.data$Malignant_normal_consensus == "Malignant"
        else cells_keep <- colnames(seurat)
        
        cor_df <- cor(seurat@meta.data[cells_keep, c(cc_cols, qc_cols)], seurat@meta.data[cells_keep, cnmf_cols]) %>%
            t() %>% 
            as.data.frame() %>% 
            tibble::rownames_to_column(var = "Program") %>% 
            tibble::add_column(.before = 1, "Sample" = id) %>% 
            mutate(Program = paste0(Sample, "_", Program))
        
        rm(seurat)
        
        return(cor_df)
        
    })
    
}

# save correlations
qc_cnmf_correlations <- get_qc_cnmf_correlations("cNMF_program_malignant", malignant_only = TRUE)
length(unique(qc_cnmf_correlations$Sample))
```

```
## [1] 42
```

```r
rr_write_tsv(qc_cnmf_correlations,
             glue("{out}/QC_cNMF_correlations.malignant.tsv"),
             "Table with correlation between cNMF programs (from malignant cells) and QC metrics across single cells, for all 10X samples (including Multiome).")
```



## Extract top genes

We next extract, for each sample, the top genes associated with each program.

The output we have from cNMF includes:

- Usage score of each program in each cell
- Contribution of each gene to each program (both TPM and the Z-scored GEPs
"which reflect how enriched a gene is in each GEP relative to all of the others")

The data that's comparable across sample in the contribution of each gene to each program.
For simplicity and to avoid loading all that data, I'll take the top 100 genes
per program and cluster programs based on overlap between those (this method follows
the NMF analysis performed in [Kinker et al,
Nature Genetics, 2020](https://www.nature.com/articles/s41588-020-00726-6)),
thanks to their code provided on [GitHub](https://github.com/gabrielakinker/CCLE_heterogeneity).



```r
get_cnmf_top_genes <- function(output_dir, program_string) {
    
    top_genes <- map(sc_samples, function(i) {
        
        message("@ ", basename(i))
        
        id <- basename(i)
        
         # load gene scores file from per-sample cNMF output
        gene_score_file <- list.files(glue("{i}/cNMF/{output_dir}/"),
                                      pattern = glob2rx(paste0(output_dir, ".gene_spectra_score.k_*.dt_0_02.txt")),
                                      full.names = TRUE)
        
        # this is a program x gene table, so we transpose
        gene_scores <- data.table::fread(gene_score_file, data.table = FALSE, sep = "\t", header = TRUE) %>%
            tibble::column_to_rownames(var = "V1") %>% t()
        
        # for each program:
        #    sort genes by scores, and take the top 100
        top_genes <- map(seq_along(colnames(gene_scores)), ~ gene_scores[, .x] %>%
                             sort(decreasing = TRUE) %>%
                             head(100) %>%
                             names())
        
        names(top_genes) <- paste0(id, program_string, 1:ncol(gene_scores))
        
        rm(gene_scores)
        
        return(top_genes)
        
    })
    
    flatten(top_genes)
    
}

cnmf_top_genes <- get_cnmf_top_genes("output_ngenes2000_niter100_malignant", "_cNMF_program_malignant_")
length(cnmf_top_genes)
```

```
## [1] 261
```

```r
saveRDS(cnmf_top_genes, file = glue("{out}/cNMF_top_genes.malignant.Rds"))
```



## Make similarity matrix

Next, we look for sets of similar programs, which we could then define as recurrent
across samples. We compute the pairwise overlap (i.e., # of genes in common)
between programs to make a similarity matrix.

This is based on [code from Kinker et al](https://github.com/gabrielakinker/CCLE_heterogeneity/blob/3c40f5fbd4b84a81094fb6caae255785506d3645/module2_rhp.R), calculating "recurrent heterogeneous
programs".



```r
# similarity matrix: pairwise overlap between programs
cnmf_intersect <- sapply(cnmf_top_genes, function(x) sapply(cnmf_top_genes, function(y) length(intersect(x, y)))) 
dim(cnmf_intersect)  
```

```
## [1] 261 261
```

```r
saveRDS(cnmf_intersect, file = glue("{out}/cNMF_intersect.malignant.Rds"))
```


## Filter rare/noisy programs

cNMF can be quite sensitive, picking up programs that are highly specific to a few
cells. To focus on programs used in the main populations of cells, we can take
advantage of the normalized program usages for each cell.

For each cell, the usages of all the programs identified by cNMF sums to 1,
so for each cell in each sample, we can calculate which program it uses most highly.
Then each program will be associated with a proportion of cells in the sample in which
that program is the most highly used.


```r
# SLOW, since it requires loading each Seurat object individually
calc_relative_usages <- function(program_string) {
    
    purrr::map_dfr(sc_samples, function(i) {
        
        message("@ ", basename(i))
        
        id <- basename(i)
        load(glue("{i}/seurat.Rda"))
        
        get_malig_cells <- function(seurat) colnames(seurat)[seurat$Malignant_normal_consensus %in%
                                                                 c("Malignant", "Likely malignant")]
        
        malig_cells <- get_malig_cells(seurat)
        cnmf_malig_cols <- colnames(seurat@meta.data)[grepl(program_string, colnames(seurat@meta.data))]
        
        # for each program, calculate the proportion of cells where that
        # program is used most highly
        props <- map_dbl(seq_along(cnmf_malig_cols),
                         ~ sum(apply(seurat@meta.data[malig_cells, cnmf_malig_cols], 1, which.max) == .x) /
                             length(malig_cells))
        
        data.frame(Sample  = id,
                   Program = cnmf_malig_cols,
                   Prop    = props) %>%
            mutate(Program = paste0(Sample, "_", Program))
        
    })
}

cnmf_program_usages <- calc_relative_usages(program_string = "cNMF_program_malignant")
saveRDS(cnmf_program_usages, file = glue("{out}/cNMF_program_usages.Rds"))
```

Plot the usages of each program:


```r
cnmf_program_usages %>%
    # mutate(Program = factor(Program, levels = hm_programs)) %>% 
    arrange(Program) %>% 
    ggplot(aes(x = Program, y = Prop)) +
    geom_bar(alpha = 0.7, stat = "identity", colour = "gray70") +
    geom_hline(yintercept = 0.05, colour = "red", alpha = 0.7) +
    scale_colour_gradientn(colours = ylrd) +
    theme_min() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    no_legend()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/01/cnmf_usages-1.png)<!-- -->


## Create heatmap

Next, we can perform hierarchical clustering over programs on the similarity matrix,
and visualize the resulting heatmap. Hierarchical clustering is performed
using `pheatmap::pheatmap()` default parameters, i.e. Euclidean distance and 
complete linkage.

Helper function:


```r
# make the column annotation for the gene programs, related to the covariates
# of each sample
program_anno <- qc_cnmf_correlations %>% 
    left_join(meta_sc, by = "Sample") %>% 
    select(Program, Sample, PRC2_group, Location) %>% 
    tibble::column_to_rownames(var = "Program")

generate_cnmf_heatmap <- purrr::partial(pheatmap,
                                        border_color = NA,
                                        show_rownames = FALSE,
                                        show_colnames = TRUE,
                                        color = custom_magma,
                                        annotation_col = program_anno %>%
                                            select(PRC2_group, Location),
                                        # annotate columns of the heatmap based on the tumor
                                        # that each program comes from
                                        annotation_colors = list(
                                            "PRC2_group" = palette_molecular,
                                            "Location"   = palette_location
                                        ),
                                        treehight_col = 0,
                                        fontsize_col = 3,
                                        cellwidth = 2,
                                        cellheight = 2)
```

Filter programs to those that are used most highly in >= 5% of cells,
and then generate the heatmap.


```r
# filter
programs_keep <- cnmf_program_usages %>% filter(Prop >= 0.05) %>% pull(Program)
cnmf_intersect_filt <- cnmf_intersect[programs_keep, programs_keep]
dim(cnmf_intersect_filt)
```

```
## [1] 146 146
```

```r
saveRDS(cnmf_intersect_filt, file = glue("{out}/cNMF_intersect.malignant_filt.Rds"))

# make heatmap
# save the heatmap output which also contains the hiearchical clustering dendrogram
cnmf_hm_filt <- generate_cnmf_heatmap(mat = cnmf_intersect_filt, silent = TRUE)
hm_programs_filt <- cnmf_hm_filt$tree_col$labels[cnmf_hm_filt$tree_col$order]
saveRDS(cnmf_hm_filt, file = glue("{out}/cNMF_heatmap.malignant_filt.Rds"))

# save both file types
generate_cnmf_heatmap(mat = cnmf_intersect_filt,
                      main = "All samples",
                      filename = glue("{figout}/cNMF_heatmap.malignant_filt.png"))

generate_cnmf_heatmap(mat = cnmf_intersect_filt,
                      main = "All samples",
                      filename = glue("{figout}/cNMF_heatmap.malignant_filt.pdf"))

knitr::include_graphics(glue("{figout}/cNMF_heatmap.malignant_filt.png"))
```

<img src="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/01/cNMF_heatmap.malignant_filt.png" width="1993" />

Therfore, after filtering, 146 programs remain.

## Define metaprograms

Gene programs identified from each tumor were then used to identify modules, i.e.
sets of programs identified recurrently across multiple samples. We designed a recursive
depth-first algorithm to traverse the hierarchical clustering dendrogram (extracted from the
heatmap above) to define discrete modules.

```
1. A set of subtrees _S_ was arbitrarily initialized by cutting the dendrogram into 5 subtrees.
2. For each subtree _t_ in _S_:
  - if there were fewer than 4 programs in _t_, it was dropped from _S_.
  - If the average inter-program similarity of the programs in _t_ was greater than 10, then t was considered a module.
  - Otherwise, _t_ was cut into 2, and each resulting subtree appended to _S_.
3. To identify the genes characterizing each module, we selected the 50 genes most frequently associated with programs belonging to the module.

```

The implementation of this algorithm is in the section in the dropdown below.

<details>


```r
#' Extract a list of k subdendrograms from a dendrogram object
#'
#' Adapted from dendextend::get_subdendrograms with bug fix
#' https://rdrr.io/cran/dendextend/src/R/get_subdendrograms.R
get_subdendrograms2 <- function(dend, k, ...) {
    clusters <- cutree(dend, k, ...)
    dend_list <- lapply(unique(clusters), function(cluster.id) {
        # bugfix: Added `names(clusters)[]` here
        find_dendrogram(dend, names(clusters)[which(clusters == cluster.id)])
    })
    class(dend_list) <- "dendlist"
    dend_list
}

#' Get the average inter-program similarity (since intra-program similarity = 100%)
avg_similarity <- function(dend) {
    
    # subset the similarity matrix to programs in the provided dendrogram
    x <- cnmf_intersect_filt[labels(dend), labels(dend)]
    
    # set the diagonal to NA to not count intra-program similarity
    diag(x) <- NA
    
    # calculate the mean similarity in the rest of the matrix
    mean(x, na.rm = TRUE)
    
}

#' Given the tree produced by pheatmap::pheatmap(), a function to extract
#' all the metaprograms from the tree
#'
#' We have arbitrarily initialized the thresholds for defining metaprograms.
#' NOTE: some programs will *not* be successfully identified within metaprograms,
#' if they don't meet the criteria. Thus, the total number of programs in the output
#' will be fewer than in number of programs in the input.
#'
#' @param tree Dendrogram
#' @param K Numeric, number of subtrees to cut \code{tree} into at the first cut
#' @param min_similarity Numeric, minimum average similarity of programs within
#' a subtree to define it as a metaprogram
#' @param min_programs Numeric, minimum number of programs within a subtree to
#' define it as a metaprogram
#'
#' @return A list, with one element per metaprogram identified. Each element is a
#' character vector containing the names of the programs in the metaprogram.
define_metaprograms <- function(tree, K = 5, min_similarity = 10, min_programs = 4) {
    
    define_metaprograms_recursive <- function(subtree,
                                              min_programs,
                                              min_similarity,
                                              debug = FALSE) {
        
        # 1. if there are fewer leaves in the tree than the minimum number of
        # programs required to define a metaprogram, drop this subtree
        if (attr(subtree, "members") < min_programs) {
            
            return(NULL)
            
            # 2. if the average similarity within this subtree is greater than the
            # minimum similarity, define this subtree as a metaprogram,
            # and add it to the list
        } else if (avg_similarity(subtree) >= min_similarity) {
            
            metaprograms[[i]] <<- labels(subtree)
            # increment the counter outside the sub-function (scoping assignment)
            i <<- i+1
            
            return(labels(subtree))
            
            # 3. if the subtree is large enough but not similar enough, cut the sutree
            # in 2, and recurse down each child/subsubtree
        } else {
            
            subsubtrees <- get_subdendrograms2(subtree, 2)
            
            lapply(subsubtrees, define_metaprograms_recursive,
                   min_programs = min_programs,
                   min_similarity = min_similarity)
            
        }
    }
    
    # initialize counter & list
    i <- 1
    metaprograms <- list()
    
    # get the first set of subtrees to initialize S, by cutting it into K subtrees
    S <- get_subdendrograms2(tree, K)
    
    # recurse!
    x <- lapply(S, define_metaprograms_recursive,
                min_programs   = min_programs,
                min_similarity = min_similarity)
    
    # enumerate metaprograms
    names(metaprograms) <- seq_along(metaprograms)
    return(metaprograms)
    
}
```

</details>

Using the helper function and the column dendrogram computed using `pheatmap::pheatmap()`,
traverse the tree and identify metaprograms:


```r
cnmf_metaprograms_filt <- define_metaprograms(as.dendrogram(cnmf_hm_filt$tree_col))
```

Do some data wrangling to get the column indices (in the heatmap) for programs
in each metaprogram, which we'll need for later visualization:


```r
# get column indices
cnmf_metaprograms_filt_idx <- map(cnmf_metaprograms_filt, ~ which(hm_programs_filt %in% .x))

# sort so they're from left to right
cnmf_metaprograms_filt_order <- names(sort(map_dbl(cnmf_metaprograms_filt_idx, 1)))

# put in the right order & rename
cnmf_metaprograms_filt_idx <- cnmf_metaprograms_filt_idx[cnmf_metaprograms_filt_order]

# rename metaprograms
names(cnmf_metaprograms_filt) <- plyr::mapvalues(names(cnmf_metaprograms_filt),
                                                 from = names(cnmf_metaprograms_filt_idx),
                                                 to = seq_along(cnmf_metaprograms_filt_idx))

# rename idx
names(cnmf_metaprograms_filt_idx) <- seq_along(cnmf_metaprograms_filt_idx)

# convert to long data frame
cnmf_metaprograms_filt_df <- imap_dfr(cnmf_metaprograms_filt,
                                      ~ data.frame(Metaprogram = as.numeric(.y), Program = .x,
                                                   stringsAsFactors = FALSE)) %>% 
    arrange(Metaprogram)

hm_metaprograms_filt <- hm_programs_filt[unname(unlist(cnmf_metaprograms_filt_idx))]
length(hm_metaprograms_filt)
```

```
## [1] 124
```

```r
save(cnmf_metaprograms_filt, cnmf_metaprograms_filt_idx, hm_metaprograms_filt,
     file = glue("{out}/cnmf_metaprograms.Rda"))

# put a gap before and after each program, and get the unique
# values, for the cases where the beginning of one program coincides
# with the end of another
# metaprograms_gaps <- unique(unlist(map(cnmf_metaprograms_filt_idx,
#                                        ~ c(.x[1] - 1, .x[length(.x)]))))
```

Re-do the heatmap for _only_ metaprogram programs:


```r
generate_cnmf_heatmap(mat = cnmf_intersect_filt[hm_metaprograms_filt, hm_metaprograms_filt],
                      # we want to put a gap in the heatmap between each metaprogram
                      # to aid visualization
                      # we can get the new gaps by taking the length of each metaprogram,
                      # and then running sum
                      gaps_row = cumsum(map(cnmf_metaprograms_filt_idx, length)),
                      gaps_col = cumsum(map(cnmf_metaprograms_filt_idx, length)),
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      filename = glue("{figout}/cNMF_heatmap_meta_only.malignant_filt_meta.png"))

generate_cnmf_heatmap(mat = cnmf_intersect_filt[hm_metaprograms_filt, hm_metaprograms_filt],
                      gaps_row = cumsum(map(cnmf_metaprograms_filt_idx, length)),
                      gaps_col = cumsum(map(cnmf_metaprograms_filt_idx, length)),
                      cluster_rows = FALSE,
                      cluster_cols = FALSE,
                      filename = glue("{figout}/cNMF_heatmap_meta_only.malignant_filt_meta.pdf"))

knitr::include_graphics(glue("{figout}/cNMF_heatmap_meta_only.malignant_filt_meta.png"))
```

<img src="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/01/cNMF_heatmap_meta_only.malignant_filt_meta.png" width="1785" />

## Program annotations {.tabset}

To annotate programs, we:

1. Calculate the correlation between each program's score across cells and QC
metrics in those cells.

2. Compute the overlap between the genes associated with
each program and reference gene signatures. Reference gene signatures were obtained
from the MSigDB collections, KEGG, PID, and Hallmark,
as well our scRNAseq mouse brain developmental dataset, restricted to non-proliferating
cell types. Since reference gene signatures differ in length, we use the
**percentage of each reference signature** overlapping program-associated genes. 

### Quality control metrics


```r
# dot plot to display correlations
qc_cnmf_correlations %>%
    filter(Program %in% hm_metaprograms_filt) %>% 
    mutate(Program = factor(Program,
                            levels = hm_metaprograms_filt)) %>%
    arrange(Program) %>%
    gather(Stat, Value, 3:ncol(.)) %>%
    rr_ggplot(aes(x = Program, y = Value), plot_num = 1) +
    geom_hline(yintercept = 0, colour = "gray80") +
    geom_point(alpha = 0.9, aes(colour = Value), size = 0.7) +
    facet_wrap(~ Stat, ncol = 1) +
    scale_colour_gradientn(colours = ylrd) +
    theme_min() +
    rotate_x() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    no_legend()
```

```
## ...writing source data of ggplot to public/R-4/figures/01/cnmf_qc_stats_filt_meta-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/01/cnmf_qc_stats_filt_meta-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/01cnmf_qc_stats_filt_meta...*]~</span>

Metaprogram 11 seems mainly technical here, so we also generate a version of
the annotations without metaprogram 11.


```r
# black and white version without M11
qc_cnmf_correlations %>%
    filter(Program %in% hm_metaprograms_filt) %>% 
    filter(!(Program %in% cnmf_metaprograms_filt$`11`)) %>% 
    mutate(Program = factor(Program,
                            levels = hm_metaprograms_filt)) %>%
    arrange(Program) %>%
    gather(Stat, Value, 3:ncol(.)) %>%
    ggplot(aes(x = Program, y = Value)) +
    geom_hline(yintercept = 0, colour = "gray80") +
    geom_point(alpha = 0.9, colour = "black", size = 0.7) +
    facet_wrap(~ Stat, ncol = 1) +
    theme_min() +
    rotate_x() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    no_legend()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/01/cnmf_qc_stats_filt_meta_no11-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/01cnmf_qc_stats_filt_meta_no11...*]~</span>

### Developmental atlas

Next, calculate the overlap between cNMF genes with the mouse atlas gene signatures:


```r
mouse_atlas_signatures <- readRDS(here("data/scRNAseq/references/mouse_atlas_extended/joint_mouse_extended.signatures_ID_20201028.Rds"))
```

Take the overlap between each program and each atlas signature:



```r
cnmf_intersect_atlas <- sapply(cnmf_top_genes,
                               function(x) sapply(mouse_atlas_signatures$hg_sym, function(y) length(intersect(x, y)) / length(y)))
dim(cnmf_intersect_atlas)
```

```
## [1] 290 261
```

```r
cnmf_intersect_atlas_long <- cnmf_intersect_atlas %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Signature") %>%
    gather(Program, Overlap, 2:ncol(.))

# filter out signatures of proliferating cells, since we do not want to confound
# cell cycle with cell type
cnmf_intersect_atlas_long_input_filt_meta <- cnmf_intersect_atlas_long %>%
    filter(!grepl("-P$|RGC|NEURP", Signature)) %>%
    filter(Program %in% hm_metaprograms_filt) %>%
    summarize_cell_types("Signature")

length(unique(cnmf_intersect_atlas_long_input_filt_meta$Signature))
```

```
## [1] 251
```

#### Empirical p-values

To assess if this overlap is significant, we can compute empirical p-values by
randomly drawing 100-gene signatures from genes detected in the mouse atlas, and
computing their overlap with the cNMF program genes.


```r
# get detected genes
atlas_path <- "/lustre03/project/6004736/sjessa/from_hydra/atlas/"
mean_expression_profile <- readRDS(file.path(atlas_path, "data/joint_mouse_extended/mean_expression_per_cluster.Rds"))
atlas_genes <- mean_expression_profile %>% select(-1, -2) %>% colnames()
length(atlas_genes)

# convert to human
ensembl <- useEnsembl(biomart = "genes")
human   <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
mouse   <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)
genes2  <- getLDS(attributes = c("mgi_symbol"),
                 filters = "mgi_symbol",
                 values = atlas_genes,
                 mart = mouse,
                 attributesL = c("hgnc_symbol"),
                 martL = human,
                 uniqueRows = TRUE)

atlas_genes_hg <- unique(genes2[, 2])
length(atlas_genes_hg)
saveRDS(atlas_genes_hg, file = glue("{out}/atlas_genes_hg.Rds"))
```

Compute null distributions: for each reference signature, generate 1,000 random signatures of 
the same length and compute overlaps with each tumour program signature. This function
then will be re-used for the MSigDB collections below.


```r
#' @param olaps_df data frame, containing at least three columns: Program, Signature, Overlap
#' @param ref_signature character, name of reference signature (should be present in olaps_df$Signature)
#' @param L numeric, length of random signatures to be sampled
#' @param N numeric, number of random iterations
#' @param universe character, vector of genes from which to sample random signatures
#'
#' @return Data frame with 3 columns: Program (from cNMF), Signature (name
#' reference signature, same as \code{ref_signature}), and p_value containing
#' the computed empirical p-value
calc_overlap_pvalues <- function(olaps_df, ref_signature, L, N, universe) {
    
    message("@ ", ref_signature)
    
    # repeat N times: sample random signature of same length as ref_signature,
    # and compute overlaps with each tumour program
    null_overlaps <- replicate(N, {
        
        # generate random signatures
        random_sig <- sample(universe, L, replace = FALSE)
        
        # compute overlap with each tumour program signature
        sapply(
            cnmf_top_genes,
            function(x) length(intersect(x, random_sig)) / length(random_sig))
        
    })
    
    # convert to tidy format
    null_overlaps_tidy <- null_overlaps %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Program") %>%
        reshape2::melt() %>% 
        suppressMessages()
    
    # compute p-values for each tumour program
    map_dfr(names(cnmf_top_genes), function(program) {
        
        # get the overlap between the signature and the tumour program
        overlap <- olaps_df %>% 
            filter(Program == program & Signature == ref_signature) %>% 
            .$Overlap
        
        # p_value := probability of observing an equal or greater overlap by chance
        # i.e. using the N random gene signatures
        data.frame(Program = program,
                   Signature = ref_signature,
                   p_value = sum(null_overlaps[names(cnmf_top_genes)[1], ] >= overlap)/N)
        
    })
    
}

# set some parameters
N <- 1000
p_value_threshold <- 1e-3
```


```r
# get genes detected in mouse atlas, computed above
atlas_genes_hg <- readRDS(glue("{out}/atlas_genes_hg.Rds"))

atlas_p_values <- imap_dfr(mouse_atlas_signatures$hg_sym,
                           ~ calc_overlap_pvalues(olaps_df = cnmf_intersect_atlas_long_input_filt_meta,
                                                  ref_signature = .y,
                                                  N = N,
                                                  L = length(.x),
                                                  universe = atlas_genes_hg))

cnmf_intersect_atlas_long_input_filt_meta <- cnmf_intersect_atlas_long_input_filt_meta %>%
    left_join(atlas_p_values, by = c("Signature", "Program")) %T>% 
    write_tsv(glue("{out}/cNMF_program_atlas_overlap_with_pvalues.tsv"))
```

Now plot only significantly overlapping signatures:


```r
# display signatures which, in at least one overlap, have:
# 1. p_value < threshold
# 2. gene overlap >= 10%
atlas_sigs_signif <- cnmf_intersect_atlas_long_input_filt_meta %>%
    mutate(Signif = ifelse(p_value < p_value_threshold & Overlap >= 0.1, TRUE, FALSE)) %>% 
    group_by(Signature) %>% 
    summarise(Signif_in_any_comparison = any(Signif)) %>% 
    filter(Signif_in_any_comparison) %>%
    pull(Signature)

cnmf_intersect_atlas_long_input_filt_meta %>%
    filter(Signature %in% atlas_sigs_signif) %>% 
    mutate(Overlap = ifelse(p_value < p_value_threshold & Overlap >= 0.1, Overlap, 0)) %>% 
    mutate(Type = factor(Type, levels = rev(names(palette_type)))) %>% 
    mutate(Program = factor(Program, levels = hm_metaprograms_filt)) %>%
    arrange(Type) %>%
    ggplot(aes(x = Program, y = Overlap, colour = Type, group = Signature)) +
    geom_line(size = 0.8, alpha = 0.5) +
    scale_colour_manual(values = palette_type) +
    theme_min() +
    guides(colour = guide_legend(ncol = 2, title = NULL)) +
    theme(legend.position = "bottom") +
    rotate_x() +
    ggtitle(paste0("N=", length(atlas_sigs_signif)))
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/01/mouse_atlas_signature_overlap_signif-1.png)<!-- -->

```r
# a version without M11 (associated with technical factors) for the main figure
cnmf_intersect_atlas_long_input_filt_meta %>%
    # show the same signatures as above
    filter(Signature %in% atlas_sigs_signif) %>% 
    mutate(Overlap = ifelse(p_value < p_value_threshold & Overlap >= 0.1, Overlap, 0)) %>% 
    # THEN, filter out M11
    filter(!(Program %in% cnmf_metaprograms_filt$`11`)) %>% 
    mutate(Type = factor(Type, levels = rev(names(palette_type)))) %>% 
    mutate(Program = factor(Program, levels = hm_metaprograms_filt)) %>%
    arrange(Type) %>% 
    ggplot(aes(x = Program, y = Overlap, colour = Type, group = Signature)) +
    geom_line(size = 0.4, alpha = 0.6) +
    scale_colour_manual(values = palette_type) +
    theme_min() +
    theme(legend.position = "bottom") +
    rotate_x()    
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/01/mouse_atlas_signature_overlap_signif-2.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/01mouse_atlas_signature_overlap_signif...*]~</span>


### Hallmark

Second, overlap cNMF genes with MSigDB [Hallmark](http://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#H) gene sets:


```r
hallmark_gmt <- readLines(here("data/misc/MSigDb_h.all.v7.4.symbols.gmt.txt")) %>%
    sapply(stringr::str_split, "\t")

# extract genes
hallmark <- hallmark_gmt %>% lapply(function(i) i[3:length(i)])
length(hallmark)
```

```
## [1] 50
```

```r
# set element name to the name of the gene set
names(hallmark) <- hallmark_gmt %>%
    lapply(function(i) i[[1]]) %>%
    unlist(use.names = FALSE)

saveRDS(hallmark, file = glue("{out}/hallmark_genelists.Rds"))

cnmf_intersect_hallmark <- sapply(cnmf_top_genes,
                                  function(x) sapply(hallmark, function(y) length(intersect(x, y)) / length(y)))
dim(cnmf_intersect_hallmark)
```

```
## [1]  50 261
```

```r
cnmf_intersect_hallmark_long <- cnmf_intersect_hallmark %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Signature") %>%
    gather(Program, Overlap, 2:ncol(.))

cnmf_intersect_hallmark_long2 <- cnmf_intersect_hallmark_long %>%
    filter(Program %in% hm_metaprograms_filt)
```

For the MSigDB collections, our universe of genes used for sampling will be the
set of genes expressed in single-cell tumor datasets.


```r
expr_genes_per_sample <- purrr::map(sc_samples, function(i) {
    
    message("@ ", basename(i))
    
    id <- basename(i)
    load(glue("{i}/seurat.Rda"))
    
    rownames(seurat)
    
})

# save expressed genes
expr_genes <- Reduce(union, expr_genes_per_sample)
length(expr_genes)
```

```
## [1] 51413
```

```r
length(unique(expr_genes))
```

```
## [1] 51413
```

```r
saveRDS(expr_genes, file = glue("{out}/expressed_genes_tumour_RNA.Rds"))
```

Use these as the universe for computing null distributions:


```r
hallmark_p_values <- imap_dfr(hallmark,
                              ~ calc_overlap_pvalues(olaps_df = cnmf_intersect_hallmark_long2,
                                                     ref_signature = .y,
                                                     N = N,
                                                     L = length(.x),
                                                     universe = expr_genes))

cnmf_intersect_hallmark_long2 <- cnmf_intersect_hallmark_long2 %>%
    left_join(hallmark_p_values, by = c("Signature", "Program")) %T>% 
    write_tsv(glue("{out}/cNMF_program_hallmark_overlap_with_pvalues.tsv"))
```

Line plot:


```r
hallmark_sigs_signif <- cnmf_intersect_hallmark_long2 %>%
    mutate(Signif = ifelse(p_value < p_value_threshold & Overlap >= 0.1, TRUE, FALSE)) %>% 
    group_by(Signature) %>% 
    summarise(Signif_in_any_comparison = any(Signif)) %>% 
    filter(Signif_in_any_comparison) %>%
    pull(Signature)

cnmf_intersect_hallmark_long2 %>%
    filter(Signature %in% hallmark_sigs_signif) %>% 
    mutate(Overlap = ifelse(p_value < p_value_threshold & Overlap >= 0.1, Overlap, 0)) %>% 
    mutate(Program = factor(Program, levels = hm_metaprograms_filt)) %>%
    ggplot(aes(x = Program, y = Overlap, colour = Signature, group = Signature)) +
    geom_line(size = 0.8, alpha = 0.5) +
    theme_min() +
    guides(colour = guide_legend(ncol = 2, title = NULL)) +
    theme(legend.position = "bottom") +
    rotate_x() +
    ggtitle(paste0("N=", length(hallmark_sigs_signif)))
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/01/hallmark_signature_overlap-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/01hallmark_signature_overlap...*]~</span>

### KEGG gene set {.tabset}

Overlap cNMF genes with MSigDB [KEGG](http://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=CP:KEGG) gene sets:


```r
kegg_gmt <- readLines(here("data/misc/MSigDB_c2.cp.kegg.v7.4.symbols.gmt.txt")) %>%
    sapply(stringr::str_split, "\t")

# extract genes
kegg <- kegg_gmt %>% lapply(function(i) i[3:length(i)])
length(kegg)
```

```
## [1] 186
```

```r
# set element name to the name of the gene set
names(kegg) <- kegg_gmt %>%
    lapply(function(i) i[[1]]) %>%
    unlist(use.names = FALSE)

saveRDS(kegg, file = glue("{out}/kegg_genelists.Rds"))

cnmf_intersect_kegg <- sapply(cnmf_top_genes,
                              function(x) sapply(kegg, function(y) length(intersect(x, y)) / length(y)))
dim(cnmf_intersect_kegg)
```

```
## [1] 186 261
```

```r
cnmf_intersect_kegg_long <- cnmf_intersect_kegg %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Signature") %>%
    gather(Program, Overlap, 2:ncol(.))

cnmf_intersect_kegg_long2 <- cnmf_intersect_kegg_long %>%
    filter(Program %in% hm_metaprograms_filt)
```


```r
kegg_p_values <- imap_dfr(kegg,
                          ~ calc_overlap_pvalues(olaps_df = cnmf_intersect_kegg_long2,
                                                 ref_signature = .y,
                                                 N = N,
                                                 L = length(.x),
                                                 universe = expr_genes))

cnmf_intersect_kegg_long2 <- cnmf_intersect_kegg_long2 %>%
    left_join(kegg_p_values, by = c("Signature", "Program")) %T>% 
    write_tsv(glue("{out}/cNMF_program_kegg_overlap_with_pvalues.tsv"))
```

Line plot:


```r
kegg_sigs_signif <- cnmf_intersect_kegg_long2 %>%
    mutate(Signif = ifelse(p_value < p_value_threshold & Overlap >= 0.1, TRUE, FALSE)) %>% 
    group_by(Signature) %>% 
    summarise(Signif_in_any_comparison = any(Signif)) %>% 
    filter(Signif_in_any_comparison) %>%
    pull(Signature)

cnmf_intersect_kegg_long2 %>%
    filter(Signature %in% kegg_sigs_signif) %>% 
    mutate(Overlap = ifelse(p_value < p_value_threshold & Overlap >= 0.1, Overlap, 0)) %>% 
    mutate(Program = factor(Program, levels = hm_metaprograms_filt)) %>%
    ggplot(aes(x = Program, y = Overlap, colour = Signature, group = Signature)) +
    geom_line(size = 0.8, alpha = 0.5) +
    theme_min() +
    guides(colour = guide_legend(ncol = 2, title = NULL)) +
    theme(legend.position = "bottom") +
    rotate_x() +
    ggtitle(paste0("N=", length(kegg_sigs_signif)))
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/01/kegg_signature_overlap-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/01kegg_signature_overlap...*]~</span>


### PID

Overlap cNMF genes with MSigDB [PID](http://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=CP:PID) gene sets:


```r
pid_gmt <- readLines(here("data/misc/MSigDB_c2.cp.pid.v7.4.symbols.gmt.txt")) %>%
    sapply(stringr::str_split, "\t")

# extract genes
pid <- pid_gmt %>% lapply(function(i) i[3:length(i)])
length(pid)
```

```
## [1] 196
```

```r
# set element name to the name of the gene set
names(pid) <- pid_gmt %>%
    lapply(function(i) i[[1]]) %>%
    unlist(use.names = FALSE)

saveRDS(pid, file = glue("{out}/pid_genelists.Rds"))

cnmf_intersect_pid <- sapply(cnmf_top_genes,
                             function(x) sapply(pid, function(y) length(intersect(x, y)) / length(y)))
dim(cnmf_intersect_pid)
```

```
## [1] 196 261
```

```r
cnmf_intersect_pid_long <- cnmf_intersect_pid %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Signature") %>%
    gather(Program, Overlap, 2:ncol(.))

cnmf_intersect_pid_long2 <- cnmf_intersect_pid_long %>%
    filter(Program %in% hm_metaprograms_filt)
```


```r
pid_p_values <- imap_dfr(pid,
                         ~ calc_overlap_pvalues(olaps_df = cnmf_intersect_pid_long2,
                                                ref_signature = .y,
                                                N = N,
                                                L = length(.x),
                                                universe = expr_genes))

cnmf_intersect_pid_long2 <- cnmf_intersect_pid_long2 %>%
    left_join(pid_p_values, by = c("Signature", "Program")) %T>% 
    write_tsv(glue("{out}/cNMF_program_pid_overlap_with_pvalues.tsv"))
```


Line plot:


```r
pid_sigs_signif <- cnmf_intersect_pid_long2 %>%
    mutate(Signif = ifelse(p_value < p_value_threshold & Overlap >= 0.1, TRUE, FALSE)) %>% 
    group_by(Signature) %>% 
    summarise(Signif_in_any_comparison = any(Signif)) %>% 
    filter(Signif_in_any_comparison) %>%
    pull(Signature)

cnmf_intersect_pid_long2 %>%
    filter(Signature %in% pid_sigs_signif) %>% 
    mutate(Overlap = ifelse(p_value < p_value_threshold & Overlap >= 0.1, Overlap, 0)) %>% 
    mutate(Program = factor(Program, levels = hm_metaprograms_filt)) %>%
    ggplot(aes(x = Program, y = Overlap, colour = Signature, group = Signature)) +
    geom_line(size = 0.8, alpha = 0.5) +
    theme_min() +
    guides(colour = guide_legend(ncol = 2, title = NULL)) +
    theme(legend.position = "bottom") +
    rotate_x() +
    ggtitle(paste0("N=", length(pid_sigs_signif)))
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/01/pid_signature_overlap-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/01pid_signature_overlap...*]~</span>

### Assemble table


```r
TABLE_ref_cnmf_overlaps <- bind_rows(
    cnmf_intersect_atlas_long_input_filt_meta %>% tibble::add_column(Source = "Atlas", .before = 1),
    cnmf_intersect_hallmark_long2 %>% tibble::add_column(Source = "Hallmark", .before = 1),
    cnmf_intersect_kegg_long2 %>% tibble::add_column(Source = "KEGG", .before = 1),
    cnmf_intersect_pid_long2 %>% tibble::add_column(Source = "PID", .before = 1)
)

rr_write_tsv(TABLE_ref_cnmf_overlaps,
             glue("{out}/TABLE_reference_cnmf_program_overlaps.tsv"),
             "Summary of overlaps between all reference gene signatures and all tumour programs")
```

```
## ...writing description of TABLE_reference_cnmf_program_overlaps.tsv to public/R-4/output/01/TABLE_reference_cnmf_program_overlaps.desc
```

## Metapgrogram signatures

To identify the genes characterizing each module, we selected the 50 genes most
frequently associated with programs belonging to the module.


```r
cnmf_metaprograms_filt_sigs <- map(
    cnmf_metaprograms_filt_idx,
    # get the top genes for all programs in the metaprogram
    ~ cnmf_top_genes[hm_programs_filt[.x]] %>% 
        # flatten
        unlist() %>%
        # count how many programs each gene appears in
        table() %>%
        # sort from most to least common
        sort(decreasing = TRUE) %>% 
        # get the top 50
        head(50) %>%
        names())                              

# save
saveRDS(cnmf_metaprograms_filt_sigs, file = glue("{out}/cNMF_metaprogram_signatures.malignant_filt.Rds"))

# show table
(cnmf_metaprograms_filt_sigs_tbl <- enframe(map_chr(cnmf_metaprograms_filt_sigs, ~ glue_collapse(.x, ","))))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["name"],"name":[1],"type":["chr"],"align":["left"]},{"label":["value"],"name":[2],"type":["chr"],"align":["left"]}],"data":[{"1":"1","2":"RPL13,RPL15,RPS14,RPS23,EEF1A1,RPL29,RPL3,RPL32,RPLP1,RPS18,RPS2,RPS3,RPS4X,RPL10,RPL13A,RPL18,RPL19,RPL30,RPL37A,RPL7A,RPL8,RPS11,RPS12,RPS15,RPS27A,RPS6,RPS8,RPS9,RPL10A,RPL11,RPL23A,RPL35A,RPL41,RPL6,RPL7,RPS13,RPS19,RPS24,RPS25,RPS27,RPS5,FAU,RPL28,RPL34,RPL36,RPLP0,RPLP2,RPS16,RPS3A,UBA52"},{"1":"2","2":"MMP16,DSCAM,LHFPL3,NLGN1,SOX6,NOVA1,DCC,FGF14,CA10,CADM2,DNM3,FGF12,GRIA2,GRID2,KCND2,LPHN3,MAP2,NRXN1,NXPH1,PCDH15,TNR,GRIA4,SNTG1,XYLT1,CSMD1,DPP6,LRRC4C,OPCML,SEZ6L,CHST11,CNTN1,CTNNA2,ETV1,SH3D19,SMOC1,ALCAM,LINC00511,LSAMP,MAML2,PTPRZ1,SOX5,CSMD3,DLGAP1,EPN2,GALNT13,GLCCI1,MEGF11,NCAM1,SEMA5A,UNC80"},{"1":"3","2":"ASPM,AURKA,BIRC5,CASC5,CCNB1,CDCA8,CENPE,CENPF,CKAP2,DEPDC1,GTSE1,KIF23,KIF2C,KIF4A,KIFC1,KPNA2,LMNB1,MKI67,NDC80,NUF2,NUSAP1,PRC1,SMC4,TACC3,TOP2A,TPX2,TTK,ARHGAP11A,ARL6IP1,BUB1,BUB1B,CCNA2,CCNB2,CDC20,CDK1,CEP55,CKAP2L,CKS2,DLGAP5,ECT2,FAM64A,FAM83D,H2AFZ,HMGB2,HN1,KIF14,KNSTRN,NCAPG,NEK2,NUCKS1"},{"1":"4","2":"DCLK1,FMN2,KCNN3,SORBS1,DTNA,LINC00478,CTNND2,EEPD1,MGAT4C,SLC4A4,ADCY2,ANKFN1,GLIS3,PITPNC1,PLEKHA5,TNC,CD44,DCLK2,GPR98,NFIA,NRG3,RFX4,RNF219-AS1,ASTN2,CADPS,KCNQ5,LIFR,NKAIN3,AC002429.5,FAT3,GFAP,LAMA2,LRIG1,MALAT1,NEAT1,NPAS3,PDE8A,RGMA,RP11-627D16.1,SPARCL1,TNIK,AHCYL1,APOE,AQP4,GRAMD3,ITPKB,LINC01088,ROBO2,RORA,TENM4"},{"1":"5","2":"LSAMP,BCAS1,CDK18,FYN,GPR17,MBP,PKP4,ARHGAP5,BMPER,CADM2,FRMD4B,IGSF11,KIF21A,MPZL1,NFASC,PPP1R16B,SEMA4D,SH3RF3,SHROOM4,SIRT2,SLC44A1,TMEM108,TNS3,TUBB4A,CADM1,CLDN11,ENPP6,EPB41L2,FRMD4A,GPC3,LINC00511,MIR219-2,MOB3B,MYO5A,PARD3B,PLCB1,PLP1,PPFIBP1,PRICKLE1,RAPGEF1,RASGEF1B,SEMA5A,SEMA6D,TMEM132B,ADAM23,ADAMTSL1,ADRA1A,ARHGAP35,BIN1,C9orf3"},{"1":"6","2":"AKAP12,TRIO,GAP43,NDRG1,NEAT1,SAMD4A,TNC,VEGFA,EMP1,LYST,PLOD2,SPOCD1,VIM,CD44,ITGB1,LPP,PGK1,SERPINE1,GBE1,GLIS3,HILPDA,INSIG2,NRCAM,PRKCB,S100A6,VMP1,ADM,BNIP3L,CA12,CAMK2D,DGKD,ELL2,FAM162A,FRMD5,GRK5,IGFBP5,IRS2,ITGA3,NAMPT,SCG2,SLC2A1,SLC2A3,SOD2,SVIL,VCL,ABCA1,BNIP3,CA9,CEBPD,CLIC4"},{"1":"7","2":"AUTS2,MAGI1,MALAT1,CHD9,KLHDC8A,SLC24A3,AC007682.1,AC159540.1,ANKRD36,CDK6,CHD7,CSMD1,KIRREL3,MAP2,MTSS1,NOVA1,NRXN3,PAK3,PCSK2,PSD3,PTPRS,RP11-782C8.2,SEZ6L,SNTG1,TCF4,TMEM2,ANKRD36C,ARHGAP11B,ASXL3,BEST3,BICD1,CLASP2,CREB5,CTNNA2,DCX,DLX6-AS1,DPP6,DSCAM,DST,FAM65B,FLJ27365,GAD1,GLCCI1,GRIA2,GRID2,HIP1,KCNQ1OT1,MAML2,MOXD1,MSRA"},{"1":"8","2":"SPARCL1,AQP4,CLU,CST3,ITM2C,NTRK2,SLC1A3,ANXA1,CD81,GFAP,GJA1,GPM6B,GRAMD3,ID4,ITM2B,NCAN,NDRG2,RHOB,S100B,ADCY2,ATP1A2,BCAN,C1orf61,CTNND2,EZR,HIF3A,ID3,LIMCH1,LINC00844,MAOB,NLRP1,NTM,PHYHIPL,QKI,TIMP3,TSC22D4,TTYH1,AGT,AK1,ALDOC,AMER2,AP1S2,APLNR,APOC1,APOE,AQP1,ATP1B1,ATP1B2,B2M,C1orf192"},{"1":"9","2":"FIGN,LPAR1,ADCY2,BMPR1A,GLIS3,GNA14,PRRX1,ROBO2,TPST1,ZBTB20,ANK2,ARHGEF3,BMPR1B,CACNB2,DCLK2,DPP10,DTNA,ETV6,FAT3,IFI16,IMMP2L,IQGAP2,L3MBTL4,LINC00478,LINC01057,LRP1B,MAN1C1,NAMPT,NFIA,NKAIN3,NRG3,PARD3B,PBX1,QKI,SLC44A5,TANC1,TNIK,ACSS3,ACTN1,ADAMTS12,AKAP6,ANO6,ATRNL1,BAZ2B,BCL6,CACHD1,CACNA2D1,CCDC85A,CD44,CDC14B"},{"1":"10","2":"C9orf117,CCDC146,EFCAB2,LRRIQ1,RP11-356K23.1,AC005281.1,ADGB,C4orf22,C8orf34,CAPS,CCDC173,CCDC30,CCDC39,CCDC42B,DCDC1,DNAAF1,DNAH11,DNAH12,DNAH6,DNAH7,DNAH9,DNAI1,DTHD1,EFCAB6,GABRG1,NEK11,PTPRQ,RSPH1,SLC47A2,SPAG16,SPAG17,SPATA17,SPEF2,SYNE1,TMEM232,TTC29,ULK4,VWA3A,WDR49,WDR78,WDR96,ZBBX,AGBL4,AK1,ARMC3,C12orf55,C1orf173,C20orf26,C5orf49,C6orf118"},{"1":"11","2":"MT-ATP6,MT-CO1,MT-CO2,MT-CO3,MT-CYB,MT-ND1,MT-ND2,MT-ND3,MT-ND4,MT-ND4L,MT-ND5,DHFR,MT-ATP8,CLU,CST3,GFAP,MT-ND6,MTRNR2L12,S100A6,APOE,CAMK2A,FTH1,HFM1,HLA-C,MTRNR2L8,NPAS3,OLFML3,PCSK1N,RP11-146F11.1,RP5-857K21.4,SPP1,TMEM158,TMSB10,TNFRSF12A,VIM,ABTB2,AC002056.5,AC004067.5,AC005237.4,AC005488.1,AC005594.3,AC006335.13,AC007249.3,AC008278.3,AC009133.14,AC012462.2,AC079163.1,AC092580.4,AC093157.1,AC093381.2"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
write_tsv(cnmf_metaprograms_filt_sigs_tbl, glue("{out}/cNMF_metaprogram_signatures.malignant_filt.tsv"))

# get unique ones
cnmf_metaprograms_filt_sigs_uniq_flat <- unique(unlist(cnmf_metaprograms_filt_sigs))

# retrieve the same lists subsetted to unique genes
unique_genes <- c()
cnmf_metaprograms_filt_sigs_uniq <- list()

for (i in seq_along(cnmf_metaprograms_filt_sigs)) {
    
    # for each metaprogram, keep the genes which haven't been seen before
    cnmf_metaprograms_filt_sigs_uniq[[i]] <- cnmf_metaprograms_filt_sigs[[i]][
        !(cnmf_metaprograms_filt_sigs[[i]] %in% unique_genes)
        ]
    unique_genes <- c(unique_genes, cnmf_metaprograms_filt_sigs_uniq[[i]]) 
    
}
```

Finally, let's create a heatmap, showing the NMF score for each signature gene in across all programs:


```r
extract_signature_scores <- function(signatures, output_dir, program_string) {
    
    map_dfr(sc_samples, function(i) {
        
        message("@ ", basename(i))
        
        gene_scores <- get_cnmf_gene_scores(i, output_dir, program_string, genes_keep = signatures)
        
    })
    
}

cnmf_metaprogram_gene_scores <- extract_signature_scores(cnmf_metaprograms_filt_sigs_uniq_flat, "output_ngenes2000_niter100_malignant", "_cNMF_program_malignant_")
write_tsv(cnmf_metaprogram_gene_scores, glue("{out}/cNMF_metaprogram_gene_scores.malignant_filt.tsv"))
```

Heatmap:


```r
plot_metaprogram_scores <- function(...) {
    
    x <- cnmf_metaprogram_gene_scores %>% 
        tibble::column_to_rownames(var = "Program") %>% 
        t() %>% 
        # rows --> unique genes appearing in any metaprogram genes
        # cols --> all programs which belong to a metaprogram
        .[cnmf_metaprograms_filt_sigs_uniq_flat, hm_metaprograms_filt] %>% 
        set_colnames(1:ncol(.))
    
    x[is.na(x)] <- 0
    
    pheatmap(x,
             color = rdbu3, 
             border_color = NA,
             scale = "row",
             gaps_col = cumsum(map(cnmf_metaprograms_filt_idx, length)),
             gaps_row = cumsum(map(cnmf_metaprograms_filt_sigs_uniq, length)),
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             show_rownames = FALSE,
             show_colnames = TRUE,
             fontsize_col = 3,
             cellwidth = 2,
             cellheight = 0.5,
             ...)
    
}

plot_metaprogram_scores(filename = glue("{figout}/cNMF_metaprogram_scores.malignant_filt.png"))
plot_metaprogram_scores(filename = glue("{figout}/cNMF_metaprogram_scores.malignant_filt.pdf"))

knitr::include_graphics(glue("{figout}/cNMF_metaprogram_scores.malignant_filt.png"))
```

<img src="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/01/cNMF_metaprogram_scores.malignant_filt.png" width="1402" />

## Prep table of top genes / metaprogram assignment of each program

Finally, for the supplementary materials, we prepare tables with top
genes per program for each sample, and the metaprogram it was assigned to (if any):


```r
# make dataframe with metaprogram membership
cnmf_df_meta <- imap_dfr(cnmf_metaprograms_filt, ~ data.frame("Metaprogram" = .y, "Program" = .x)) %>% 
    arrange(Metaprogram)

# make dataframe with top genes
cnmf_df_top <- imap_dfr(cnmf_top_genes, ~ data.frame("Program" = .y, "Top_100_program_associated_genes" = paste0(.x, collapse = ",")))

cnmf_df_all <- qc_cnmf_correlations %>% 
    left_join(cnmf_df_meta, by = "Program") %>% 
    left_join(cnmf_program_usages, by = c("Program", "Sample")) %>% 
    mutate(Note = case_when(
        # these programs were included in the analysis, and were contained in a metaprogram
        Prop >= 0.05 & !is.na(Metaprogram) ~ "Included in analysis",
        # these programs were included in the analysis, but were part of subtrees
        # that did not pass the metaprogram definition step
        Prop >= 0.05 & is.na(Metaprogram) ~ "Included in analysis, not part of module",
        # these programs (with Prop < 0.05) were used in small # of cells, 
        # and thus were excluded
        TRUE ~ "Excluded, rare program"
    )) %>%
    rename(Prop_cells_most_active = Prop) %>% 
    left_join(cnmf_df_top, by = "Program")
    
# sanity checks
cnmf_df_all %>% filter(Note == "Included in analysis, not part of module") %>% pull(Metaprogram) %>% is.na() %>% all()
```

```
## [1] TRUE
```

```r
cnmf_df_all %>% filter(Note == "Included in analysis") %>% pull(Metaprogram) %>% is.na() %>% any()
```

```
## [1] FALSE
```

```r
cnmf_df_all %>% filter(Note == "Excluded, rare program") %>% pull(Metaprogram) %>% is.na() %>% all()
```

```
## [1] TRUE
```

```r
table(cnmf_df_all$Note)
```

```
## 
##                   Excluded, rare program 
##                                      115 
##                     Included in analysis 
##                                      124 
## Included in analysis, not part of module 
##                                       22
```

```r
nrow(cnmf_df_all)
```

```
## [1] 261
```

```r
rr_write_tsv(cnmf_df_all,
             glue("{out}/TABLE_cNMF_programs_per_sample.tsv"),
             "Overview of per-sample cNMF programs")
```

```
## ...writing description of TABLE_cNMF_programs_per_sample.tsv to public/R-4/output/01/TABLE_cNMF_programs_per_sample.desc
```

```r
nrow(cnmf_df_all)
```

```
## [1] 261
```

# Number of malignant cells

Count the total number of malignant cells -- including patients with samples profiled
by multiple technologies, and including ATAC cells:


```r
rna_samples <- c(list.files(here("data/scRNAseq/pipeline_10X/"), full.names = TRUE))
rna_samples <- rna_samples[!grepl("Makefile", rna_samples)]

multi_samples <- list.files(here("R-4/data/scMultiome/pipeline_10X_Multiome/"), full.names = TRUE)
multi_samples <- multi_samples[!grepl("Makefile", multi_samples)]

# combine
sc_samples_all <- c(rna_samples, multi_samples)

n_cells_per_sample_rna <- purrr::map_dfr(sc_samples_all, function(i) {
    
    message("@ ", basename(i))
    
    id <- basename(i)
    load(glue("{i}/seurat.Rda"))
    
    data.frame("Sample" = i,
               "N_cells" = nrow(seurat@meta.data),
               "N_cells_malignant" = seurat@meta.data %>% 
                   filter(Malignant_normal_consensus %in% c("Malignant", "Likely malignant")) %>% 
                   nrow())
    
})

atac_samples <- c(list.files(here("R-4/data/scATACseq/pipeline_10X_ATAC"), full.names = TRUE))
atac_samples <- atac_samples[!grepl("Makefile", atac_samples)]

n_cells_per_sample_atac <- purrr::map_dfr(atac_samples, function(i) {
    
    message("@ ", basename(i))
    
    id <- basename(i)
    load(glue("{i}/seurat.Rda"))
    
    # for ATACseq, cells which are not projected to immune/vascular are treated
    # as malignant
    data.frame("Sample" = i,
               "N_cells" = nrow(seurat_atac@meta.data),
               "N_cells_malignant" = seurat_atac@meta.data %>% 
                   filter(cluster_predicted.id != "Microglia/macrophages") %>% 
                   nrow())
    
})

n_cells_per_sample_all <- bind_rows(n_cells_per_sample_atac, n_cells_per_sample_rna)

write_tsv(n_cells_per_sample_all, glue("{out}/n_cells_per_sample.tsv"))

# calculate totals
sum(n_cells_per_sample_all$N_cells)
```

```
## [1] 226011
```

```r
sum(n_cells_per_sample_all$N_cells_malignant)
```

```
## [1] 181282
```



<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2022-06-29 11:03:38
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
## [1] stats     graphics  grDevices datasets  utils     methods   base     
## 
## other attached packages:
##  [1] magrittr_2.0.1     viridis_0.5.1      viridisLite_0.3.0  RColorBrewer_1.1-2
##  [5] SeuratObject_4.0.4 Seurat_4.0.0       cowplot_1.1.1      ggrastr_0.2.3     
##  [9] dendextend_1.15.1  ape_5.5            pheatmap_1.0.12    purrr_0.3.4       
## [13] tibble_3.1.6       glue_1.6.1         readr_2.1.1        ggrepel_0.9.1     
## [17] ggplot2_3.3.5      dplyr_1.0.7        tidyr_1.1.4        biomaRt_2.50.2    
## [21] here_1.0.1        
## 
## loaded via a namespace (and not attached):
##   [1] BiocFileCache_2.2.1  plyr_1.8.6           igraph_1.2.11       
##   [4] lazyeval_0.2.2       splines_4.1.2        listenv_0.8.0       
##   [7] scattermore_0.7      digest_0.6.29        htmltools_0.5.2     
##  [10] fansi_1.0.2          memoise_2.0.1        tensor_1.5          
##  [13] cluster_2.1.2        ROCR_1.0-11          tzdb_0.2.0          
##  [16] globals_0.14.0       Biostrings_2.58.0    matrixStats_0.61.0  
##  [19] vroom_1.5.7          prettyunits_1.1.1    colorspace_2.0-2    
##  [22] blob_1.2.2           rappdirs_0.3.3       xfun_0.29           
##  [25] crayon_1.4.2         jsonlite_1.7.3       spatstat_1.64-1     
##  [28] spatstat.data_2.1-2  survival_3.2-13      zoo_1.8-9           
##  [31] polyclip_1.10-0      gtable_0.3.0         zlibbioc_1.36.0     
##  [34] XVector_0.30.0       leiden_0.3.9         future.apply_1.8.1  
##  [37] BiocGenerics_0.36.1  abind_1.4-5          scales_1.1.1        
##  [40] DBI_1.1.2            miniUI_0.1.1.1       Rcpp_1.0.8          
##  [43] xtable_1.8-4         progress_1.2.2       reticulate_1.23     
##  [46] bit_4.0.4            stats4_4.1.2         htmlwidgets_1.5.4   
##  [49] httr_1.4.2           ellipsis_0.3.2       ica_1.0-2           
##  [52] pkgconfig_2.0.3      XML_3.99-0.8         uwot_0.1.11         
##  [55] deldir_1.0-6         sass_0.4.0           dbplyr_2.1.1        
##  [58] utf8_1.2.2           tidyselect_1.1.1     rlang_0.4.12        
##  [61] reshape2_1.4.4       later_1.3.0          AnnotationDbi_1.56.2
##  [64] munsell_0.5.0        tools_4.1.2          cachem_1.0.6        
##  [67] cli_3.1.1            generics_0.1.1       RSQLite_2.2.9       
##  [70] ggridges_0.5.3       evaluate_0.14        stringr_1.4.0       
##  [73] fastmap_1.1.0        yaml_2.2.1           goftest_1.2-3       
##  [76] knitr_1.37           bit64_4.0.5          fitdistrplus_1.1-6  
##  [79] RANN_2.6.1           KEGGREST_1.34.0      pbapply_1.5-0       
##  [82] future_1.23.0        nlme_3.1-153         mime_0.12           
##  [85] xml2_1.3.3           compiler_4.1.2       beeswarm_0.4.0      
##  [88] plotly_4.10.0        filelock_1.0.2       curl_4.3.2          
##  [91] png_0.1-7            spatstat.utils_2.3-0 bslib_0.3.1         
##  [94] stringi_1.7.6        lattice_0.20-45      Matrix_1.3-4        
##  [97] vctrs_0.3.8          pillar_1.6.4         lifecycle_1.0.1     
## [100] BiocManager_1.30.15  lmtest_0.9-39        jquerylib_0.1.4     
## [103] RcppAnnoy_0.0.19     data.table_1.14.2    irlba_2.3.5         
## [106] httpuv_1.6.5         patchwork_1.1.1      R6_2.5.1            
## [109] promises_1.2.0.1     renv_0.15.5          KernSmooth_2.23-20  
## [112] gridExtra_2.3        vipor_0.4.5          IRanges_2.24.1      
## [115] parallelly_1.30.0    codetools_0.2-18     MASS_7.3-54         
## [118] assertthat_0.2.1     rprojroot_2.0.2      withr_2.4.3         
## [121] sctransform_0.3.3    S4Vectors_0.28.1     mgcv_1.8-38         
## [124] parallel_4.1.2       hms_1.1.1            rpart_4.1-15        
## [127] grid_4.1.2           rmarkdown_2.11       Rtsne_0.15          
## [130] git2r_0.29.0         Biobase_2.54.0       shiny_1.7.1         
## [133] ggbeeswarm_0.6.0
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
