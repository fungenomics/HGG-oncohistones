---
title: "01B - Prepare human fetal hindbrain data"
author: "Selin Jessa [[selin.jessa@mail.mcgill.ca](mailto:selin.jessa@mail.mcgill.ca)]"
date: "29 June, 2022"
params:
  resources: "NOT SPECIFIED"
output:
  html_document:
    keep_md: yes
    code_folding: show
    theme: flatly
    css: ../include/style.css
    toc: yes
    toc_depth: 3
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
## Document index: 01B
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
## public/output/01B
```

```
## public/figures/01B
```



Setting a random seed:



```r
set.seed(100)
```



***



<!-- END OF FRONT MATTER -->


# Overview

In this analysis, we'll prepare the reference from the human fetal hindbrain data
which has been processed with our single-cell workflow, and gather the
required formats and outputs needed for comparisons to tumors:

* annotated Seurat objects, per-sample and integrated across samples
* gene signatures, per-cluster
* mean expression profiles, per-cluster

This data comes from [Bhaduri et al](http://dx.doi.org/10.1038/s41586-021-03910-8), and [Eze et al](http://dx.doi.org/10.1038/s41593-020-00794-1).

The structure of this document largely follows the analogous document for the thalamus (`01A`),
however there are no published labels used for this brain region.


# Libraries



```r
# Load libraries here
library(here)
library(tidyr)
library(dplyr)
library(readxl)
library(readr)
library(glue)
library(ggplot2)
library(magrittr)
library(purrr)
library(cowplot)
library(feather)
library(Seurat)
library(icytobox)
library(ggrepel)

source(here("include/style.R"))
source(here("code/functions/scRNAseq.R"))
ggplot2::theme_set(theme_min())
```




## Metadata



```r
info_samples <- read_tsv(here("data/metadata/metadata_human_fetal_10X.tsv")) %>% 
    filter(Region %in% c("Hindbrain", "Cerebellum"))

palette_sample <- info_samples %>% select(Sample, Sample_color) %>% tibble::deframe()
palette_alias  <- info_samples %>% select(Alias, Sample_color) %>% tibble::deframe()

ages <- info_samples$Age %>% unique %>% sort()
palette_age <- colorRampPalette(brewer.pal(8, "YlGnBu"))(n = length(ages))
names(palette_age) <- ages
```



## QC



```r
qc_stats <- map_dfr(1:nrow(info_samples), ~ read_tsv(here(file.path(info_samples[.x, ]$Path, "preprocessing/seurat_metrics.tsv"))) %>% 
                        tibble::add_column(.before = 1, Sample = info_samples[.x, ]$Sample))

#' a function to plot a scatterplot of 2 column in the aggregated qc dataframe
#' 
#' @example qc_scatter(nCount_RNA_min.postQC, nCount_RNA_mean.postQC)
qc_scatter <- function(col1, col2) {
    
    col1_quo <- enquo(col1)
    col2_quo <- enquo(col2)
    
    qc_stats %>% 
        ggplot(aes(x = !! col1_quo, y = !! col2_quo)) +
        geom_point(aes(color = Sample), size = 4, alpha = 0.5) +
        geom_text_repel(aes(label = Sample, color = Sample), size = 4) +
        scale_color_manual(values = palette_sample) +
        theme_min() +
        no_legend()
    
}

qc_scatter(N_cells_before, N_cells_after)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/01B/qc-1.png)<!-- -->

## Prep objects

For each sample:

* load and plot the Seurat object
* count the number of cells in each cluster, and give each cluster a label prefixed w/ sample name
* add the correlations predictions to the object
* save the object



```r
prep_seurat <- function(path_to_sample) {
    
    message("@ ", basename(path_to_sample))
    
    # load the object as output by the SC workflow
    load(file.path(path_to_sample, "preprocessing/seurat.Rda"))
    
    # get the labels by correlations and add them to the Seurat object
    labels_cor <- read_tsv(file.path(path_to_sample, "correlations/ref.joint_mouse_extended/correlation.predicted_labels.tsv")) %>%
        tibble::column_to_rownames(var = "cell")
    seurat <- AddMetaData(seurat,
                          metadata = labels_cor,
                          col.name = c("COR_ref.joint_mouse_extended", "COR_ref.joint_mouse_extended_score"))
    
    # get some summary stats per cluster
    seurat[["Cluster_number"]] <- Idents(object = seurat)
    
    info_samples_cluster <- seurat@meta.data %>%
        group_by(Cluster_number) %>%
        mutate(N_cells = n()) %>%
        group_by(Cluster_number, N_cells) %>%
        summarize_at(vars(nCount_RNA, nFeature_RNA, percent.mito, percent.ribo, S.Score, G2M.Score), mean) %>%
        tibble::add_column(.before = 1, Sample = basename(path_to_sample)) %>%
        # add sample info table, and prefix the cluster number by the alias, so that it's unique
        # in the form, [Alias]_[Cluster]
        left_join(select(info_samples, Sample, Age, Alias), by = "Sample") %>%
        mutate(Cluster = paste0(Alias, "_", Cluster_number)) %>%
        select(Sample, Age, Alias, Cluster_number, Cluster, everything())
    
    # get the most common projection from the correlations (excluding the information
    # of the timepoint from the mouse clusters)
    best_match <- seurat@meta.data %>% rowwise() %>% 
        mutate(Type = stringr::str_split_fixed(COR_ref.joint_mouse_extended, "_", 2) %>% getElement(2)) %>%
        group_by(Cluster_number, Type) %>%
        mutate(n = n()) %>%
        group_by(Cluster_number) %>% top_n(1, n) %>%
	  slice(1) %>%
	  select(Cluster_number, n, Type, COR_ref.joint_mouse_extended)
    
    # sanity check
    nrow(info_samples_cluster) == nrow(best_match)
    
    info_samples_cluster <- info_samples_cluster %>%
        left_join(best_match, by = "Cluster_number") %>%
        mutate(Prop_matched_mouse_atlas = round(n/N_cells, 2)) %>%
        select(-n) %>% 
        mutate(ID_20210610 = paste0(Cluster, "_", Type)) %>% 
        select(-Type)
    
    # put the unique cluster labels in the Seurat object
    seurat$Cluster <- plyr::mapvalues(seurat$Cluster_number,
                                      from = info_samples_cluster$Cluster_number,
                                      to = info_samples_cluster$Cluster)
    Idents(seurat) <- "Cluster"
    
    seurat$ID_20210610 <- plyr::mapvalues(seurat$Cluster_number,
                                          from = info_samples_cluster$Cluster_number,
                                          to = info_samples_cluster$ID_20210610)
    
    # correct the palette
    names(seurat@misc$colours) <- info_samples_cluster$Cluster
    
    # save the new seurat object
    # save(seurat, file = here(file.path(path_to_sample, "seurat.Rda")))
    
    return(info_samples_cluster)
    
}

# gather all cluster summary info
info_clusters <- map_dfr(1:nrow(info_samples), ~ prep_seurat(info_samples[.x, ]$Path))
write_tsv(info_clusters, glue("{out}/info_clusters.tsv"))
```





```r
info_clusters <- read_tsv(glue("{out}/info_clusters.tsv"))
```



## Sample exclusions

No samples will be excluded due to low QC.

## Per-cluster QC



```r
info_clusters %>% nrow()
```



```
## [1] 161
```



```r
info_clusters2 <- info_clusters %>%
    mutate(Note_sample = "OK",
           Note_cluster = case_when(
               N_cells < 30 ~ "Less than 30 cells",
               nCount_RNA < 1000 ~ "Low number of UMIs",
               percent.ribo > 50 ~ "High ribosomal content",
               TRUE ~ "OK"
           ))

info_clusters2 %>% filter(Note_sample == "OK" & Note_cluster == "OK") %>% nrow()
```



```
## [1] 141
```



```r
write_tsv(info_clusters2, glue("{out}/info_clusters2.tsv"))

# total # of cells
sum(info_clusters2$N_cells)
```



```
## [1] 79428
```



## QC tables



```r
TABLE_hindbrain_qc <- info_samples %>% 
    left_join(qc_stats, by = "Sample") %>% 
    mutate(Publication = case_when(
        grepl("CS", Sample) ~ "Eze et al, Nature Neuroscience, 2021",
        grepl("GW", Sample) ~ "Bhaduri et al, Nature, 2021"
    )) %>% 
    mutate(Note_sample = "OK") %>% 
    select(
        # Sample info
        Sample, Age, Alias, Region, Publication,
        # QC
        # N_reads = `Number of Reads`,
        # N_cells_estimated = `Estimated Number of Cells`,
        N_cells_after_filtering = N_cells_after,
        # Reads_mapped_to_genome = `Reads Mapped to Genome`,
        # Reads_mapped_to_transcriptome = `Reads Mapped Confidently to Transcriptome`,
        Mean_mitochondrial_content_after_filtering = percent.mito_mean.postQC,
        Mean_UMIs_after_filtering = nCount_RNA_mean.postQC,
        Mean_N_genes_after_filtering = nFeature_RNA_mean.postQC,
        Max_mito_threshold = percent.mito_max.threshold,
        Min_N_genes_threshold = nFeature_RNA_min.threshold,
        Max_N_genes_threshold = nFeature_RNA_max.threshold,
        Max_UMIs_threshold = nCount_RNA_max.threshold,
        Note_sample)

rr_write_tsv(TABLE_hindbrain_qc,
             glue("{out}/TABLE_hindbrain_QC.tsv"),
             "Summary of sample info and QC for hindbrain samples")
```



```
## ...writing description of TABLE_hindbrain_QC.tsv to public/output/01B/TABLE_hindbrain_QC.desc
```



# Gather gene signatures

First, load all cluster markers and filter out markers for samples/cluters that we'll exclude:



```r
clusters_keep <- info_clusters2 %>% filter(Note_cluster == "OK") %>% pull(Cluster)

cluster_markers <- map_dfr(1:nrow(info_samples), ~ read_tsv(here(file.path(info_samples[.x, ]$Path, "preprocessing/cluster_markers.tsv"))) %>%
                               tibble::add_column(Sample = info_samples[.x, ]$Sample)) %>%
    # get positive markers
    filter(avg_logFC > 0) %>%
    mutate(cluster = as.numeric(cluster)) %>%
    # add the cluster/sample info
    left_join(info_clusters2 %>% select(Sample, Cluster_number, Cluster), by = c("Sample" = "Sample", "cluster" = "Cluster_number")) %>%
    rename(Cluster_number = cluster) %>%
    # filter out excluded samples and clusters
    filter(Cluster %in% clusters_keep) %>%
    # get the diff between pct.1 and pct.2
    mutate(pct_diff = pct.1 - pct.2) %>%
    # reorder
    select(gene, pct_diff, everything())

# sanity checks
unique(cluster_markers$Cluster) %>% length()
```



```
## [1] 141
```



```r
all(cluster_markers$Cluster %in% clusters_keep)
```



```
## [1] TRUE
```



```r
head(cluster_markers)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["gene"],"name":[1],"type":["chr"],"align":["left"]},{"label":["pct_diff"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["p_val"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["avg_logFC"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["pct.1"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["pct.2"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["p_val_adj"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["Cluster_number"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["Sample"],"name":[9],"type":["chr"],"align":["left"]},{"label":["Cluster"],"name":[10],"type":["chr"],"align":["left"]}],"data":[{"1":"ARL4A","2":"0.090","3":"9.936107e-57","4":"0.3064309","5":"0.960","6":"0.870","7":"1.816618e-52","8":"0","9":"CS13_hindbrain","10":"CS13-H_0"},{"1":"DDIT4","2":"0.061","3":"1.025524e-49","4":"0.2540863","5":"0.977","6":"0.916","7":"1.874966e-45","8":"0","9":"CS13_hindbrain","10":"CS13-H_0"},{"1":"LDHA","2":"0.063","3":"7.898211e-47","4":"0.2648140","5":"0.952","6":"0.889","7":"1.444030e-42","8":"0","9":"CS13_hindbrain","10":"CS13-H_0"},{"1":"FGFBP3","2":"0.130","3":"6.428063e-39","4":"0.2520210","5":"0.865","6":"0.735","7":"1.175243e-34","8":"0","9":"CS13_hindbrain","10":"CS13-H_0"},{"1":"SLC2A1","2":"0.170","3":"2.197159e-33","4":"0.2939371","5":"0.556","6":"0.386","7":"4.017066e-29","8":"0","9":"CS13_hindbrain","10":"CS13-H_0"},{"1":"MCM5","2":"0.149","3":"4.477445e-33","4":"0.2760484","5":"0.355","6":"0.206","7":"8.186113e-29","8":"0","9":"CS13_hindbrain","10":"CS13-H_0"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# save
write_tsv(cluster_markers, glue("{out}/cluster_markers.tsv"))
```



Second, filter these down to 100-gene signatures:



```r
#' Create 100-gene signature given the markers for one cluster
filter_markers <- function(markers, sp = "hg", n_top = 100) {
    
    markers %>%
        # filter out mitochondrial and ribosomal genes
        dplyr::filter(., !grepl("RPS|RPL|MRPS|MRPL|^MT-", gene)) %>%
        group_by(Cluster) %>%
        # get positive markers
        dplyr::filter(avg_logFC > 0) %>%
        # sort by p-val
        arrange(p_val_adj) %>%
        # make sure they're unique
        distinct(gene, .keep_all = TRUE) %>%
        dplyr::slice(1:n_top) %>%
        ungroup()
    
}

# read gene annotation
anno <- read_tsv(here("data/misc/gene_annotation_with_homologous_mice_genes.txt"))
```



```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   hsapiens_ensembl_gene_id = col_character(),
##   hsapiens_external_gene_id = col_character(),
##   mmusculus_homolog_ensembl_gene = col_character(),
##   mmusculus_external_gene_id = col_character(),
##   gene_biotype = col_character(),
##   description = col_character()
## )
```



```r
hindbrain_signatures_df <- cluster_markers %>%
    left_join(select(anno,
                     hsapiens_external_gene_id,
                     hsapiens_ensembl_gene_id,
                     mmusculus_external_gene_id,
                     mmusculus_ensembl_gene_id = mmusculus_homolog_ensembl_gene),
              by = c("gene" = "hsapiens_external_gene_id")) %>%
    filter_markers()
```



## Signature QC

Compute some stats to assess how specific these are:



```r
# Sanity check
sort(table(hindbrain_signatures_df$Cluster))
```



```
## 
##    CS13-H_0    CS15-H_2    CS13-H_1   CS14-H1_0    CS13-H_7    CS15-H_0 
##           6           8          10          12          20          20 
##   CS14-H1_1    CS19-H_5    CS15-H_3   CS14-H1_8   CS14-H1_9    CS13-H_4 
##          23          23          27          33          34          35 
##    CS19-H_9   CS14-H1_3   CS14-H1_6    CS13-H_5   CS13-H_11  CS14-H1_11 
##          39          45          47          49          50          51 
##  CS22-Cb_13   CS14-H2_2   CS22-Cb_3    CS13-H_3   CS13-H_13    CS22-H_0 
##          52          54          56          63          73          73 
##    CS15-H_7   GW20-Cb_0  CS20-Cb1_5  CS20-Cb2_0   CS14-H1_5   CS14-H1_4 
##          76          76          84          84          85          86 
##    CS19-H_3  CS14-H1_15    CS22-H_7    CS13-H_2  CS20-Cb1_0  CS14-H1_12 
##          88          89          89          93          95          96 
##   CS14-H2_9   CS14-H2_0   CS19-H_13  GW20-Cb_10   CS22-Cb_1   CS13-H_10 
##          96          98          98          98          99         100 
##   CS13-H_12    CS13-H_6    CS13-H_8    CS13-H_9  CS14-H1_10  CS14-H1_13 
##         100         100         100         100         100         100 
##  CS14-H1_14  CS14-H1_16  CS14-H1_17  CS14-H1_18   CS14-H1_2   CS14-H1_7 
##         100         100         100         100         100         100 
##   CS14-H2_1  CS14-H2_10   CS14-H2_3   CS14-H2_4   CS14-H2_5   CS14-H2_6 
##         100         100         100         100         100         100 
##   CS14-H2_7   CS14-H2_8    CS15-H_1    CS15-H_4    CS15-H_5    CS15-H_6 
##         100         100         100         100         100         100 
##    CS15-H_8    CS19-H_0    CS19-H_1   CS19-H_10   CS19-H_11   CS19-H_12 
##         100         100         100         100         100         100 
##   CS19-H_14   CS19-H_15   CS19-H_16   CS19-H_17   CS19-H_18   CS19-H_19 
##         100         100         100         100         100         100 
##    CS19-H_2    CS19-H_4    CS19-H_6    CS19-H_7    CS19-H_8  CS20-Cb1_1 
##         100         100         100         100         100         100 
## CS20-Cb1_10 CS20-Cb1_11  CS20-Cb1_2  CS20-Cb1_3  CS20-Cb1_4  CS20-Cb1_6 
##         100         100         100         100         100         100 
##  CS20-Cb1_7  CS20-Cb1_8  CS20-Cb1_9  CS20-Cb2_1 CS20-Cb2_10  CS20-Cb2_2 
##         100         100         100         100         100         100 
##  CS20-Cb2_3  CS20-Cb2_4  CS20-Cb2_6  CS20-Cb2_7  CS20-Cb2_8  CS20-Cb2_9 
##         100         100         100         100         100         100 
##   CS22-Cb_0  CS22-Cb_10  CS22-Cb_11  CS22-Cb_12   CS22-Cb_2   CS22-Cb_4 
##         100         100         100         100         100         100 
##   CS22-Cb_5   CS22-Cb_6   CS22-Cb_7   CS22-Cb_8   CS22-Cb_9    CS22-H_1 
##         100         100         100         100         100         100 
##   CS22-H_10   CS22-H_11   CS22-H_12   CS22-H_13    CS22-H_2    CS22-H_3 
##         100         100         100         100         100         100 
##    CS22-H_4    CS22-H_5    CS22-H_6    CS22-H_8    CS22-H_9   GW20-Cb_1 
##         100         100         100         100         100         100 
##  GW20-Cb_11  GW20-Cb_12  GW20-Cb_13  GW20-Cb_14  GW20-Cb_15  GW20-Cb_16 
##         100         100         100         100         100         100 
##  GW20-Cb_17  GW20-Cb_18  GW20-Cb_19  GW20-Cb_20   GW20-Cb_3   GW20-Cb_4 
##         100         100         100         100         100         100 
##   GW20-Cb_7   GW20-Cb_8   GW20-Cb_9 
##         100         100         100
```



```r
hindbrain_signatures_df_qc <- hindbrain_signatures_df %>%
    group_by(Sample, Cluster) %>%
    summarise(median_logFC = median(avg_logFC),
              max_logFC = max(avg_logFC),
              n_gene_FC_above_1.5 = sum(avg_logFC > log2(1.5)),
              median_pct.1 = median(pct.1),
              median_pct.2 = median(pct.2)) %>%
    arrange(median_logFC)
```



```
## `summarise()` has grouped output by 'Sample'. You can override using the `.groups` argument.
```



```r
# save intermediate files
write_tsv(hindbrain_signatures_df, glue("{out}/hindbrain_signatures_df.tsv"))
write_tsv(hindbrain_signatures_df_qc, glue("{out}/hindbrain_signatures_df_qc.tsv"))
```



We'll also want to remove erythrocyte & melanocyte clusters, which are quite biologically distinct:



```r
erythro_markers <- c("HBB", "HBA1", "HBA2", "HBZ", "HBG1")
(erythro_count <- cluster_markers %>% filter(gene %in% erythro_markers) %>% pull(Cluster) %>% table())
```



```
## .
## CS14-H1_15   CS15-H_7 GW20-Cb_20 
##          4          4          4
```



```r
# remove the clusters with 4-5 of these markers:
erythro_clusters <- names(erythro_count)[erythro_count >= 4]

melano_markers <- c("DCT", "MITF", "PMEL")
(melano_count <- cluster_markers %>% filter(gene %in% melano_markers) %>% pull(Cluster) %>% table())
```



```
## .
##  CS13-H_12   CS13-H_8   CS13-H_9 CS14-H1_14 CS14-H2_10  CS14-H2_5  CS14-H2_7 
##          1          1          1          1          1          1          1 
## GW20-Cb_19 
##          1
```



```r
# remove the clusters with all three of these:
melano_clusters <- names(melano_count)[melano_count == 3]
```




Filter out signatures which have outlier QC stats, indicative of non-specific signatures:



```r
hindbrain_signatures_df_qc_filt <- hindbrain_signatures_df_qc %>%
    filter(median_pct.2 < 0.4) %>%
    filter(n_gene_FC_above_1.5 > 10)

signatures_keep <- hindbrain_signatures_df_qc_filt$Cluster
length(signatures_keep)
```



```
## [1] 106
```



```r
info_clusters3 <- info_clusters2 %>%
    mutate(Note_cluster = case_when(
        Cluster %in% erythro_clusters ~ "Erythrocyte cluster, excluded",
        Cluster %in% melano_clusters ~ "Melanocyte cluster, excluded",
        TRUE ~ Note_cluster
    )) %>%
    mutate(Note_signature = case_when(
        Note_cluster != "OK" ~ "Excluded, see Note_sample/Note_cluster",
        !(Cluster %in% signatures_keep) ~ "Non-specific signature",
        TRUE ~ "OK"
    ))

write_tsv(info_clusters3, glue("{out}/info_clusters3.tsv"))

clusters_drop <- c(erythro_clusters, melano_clusters)
```



Convert to a list:



```r
hindbrain_signatures_df_filt <- hindbrain_signatures_df %>%
    filter(Cluster %in% signatures_keep & !(Cluster %in% clusters_drop))

hf_hindbrain_signatures <- list('hg_sym' = split(hindbrain_signatures_df, f = hindbrain_signatures_df$Cluster) %>% map(~ pull(.x, gene) %>% .[!is.na(.)]),
                               'hg_ens' = split(hindbrain_signatures_df, f = hindbrain_signatures_df$Cluster) %>% map(~ pull(.x, hsapiens_ensembl_gene_id) %>% .[!is.na(.)]))

save(hf_hindbrain_signatures, file = glue("{out}/hindbrain_signatures.Rda"))
```



# Join data

Next, let's write out the config for integrating these samples:



```r
info_samples2 <- info_samples %>%
    select(Path, Sample = Alias, Age, Color = Sample_color) %>%
    mutate(Path = paste0(Path, "/seurat.Rda")) %>%
    write_tsv(glue("{out}/info.samples.tsv"))

info_groups <- data.frame("Covariate" = "Age",
                          "Levels" = unique(info_samples$Age),
                          "Color" = c("#7bccc4", "#ccebc5", "#fdd49e", "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000")) %T>%
    write_tsv(glue("{out}/info.groups.tsv"))
```




## Annotation

Make the annotation for the cell type projection scripts:



```r
# set cluster color as best-matching cluster in mouse atlas
palette_mouse_df <- palette_joint_mouse_extended %>% tibble::enframe()
cluster_palette_mouse_match <- info_clusters3 %>%
  left_join(palette_mouse_df, by = c("COR_ref.joint_mouse_extended" = "name")) %>% 
  select(Cluster, value) %>%
  tibble::deframe()

# write annotation
cluster_palette_mouse_match %>% 
  tibble::enframe("Cell_type", "Color") %>% 
  write_tsv(glue("{out}/hf_hindbrain_annotation.tsv"))
```



<!-- END MATTER, insert reproducibility info -->




***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:



```
## 2022-06-29 10:03:17
```



The git repository and last commit:



```
## Local:    master /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public
## Remote:   master @ origin (git@github.com:fungenomics/HGG-oncohistones.git)
## Head:     [009cdf0] 2022-06-29: Add README and update infrastructure
```



The random seed was set with `set.seed(100)`

The R session info:
<details>



```
## Error in get(genname, envir = envir) : object 'testthat_print' not found
```



```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value                           
##  version  R version 3.6.1 (2019-07-05)    
##  os       Rocky Linux 8.5 (Green Obsidian)
##  system   x86_64, linux-gnu               
##  ui       X11                             
##  language (EN)                            
##  collate  en_CA.UTF-8                     
##  ctype    en_CA.UTF-8                     
##  tz       EST5EDT                         
##  date     2022-06-29                      
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  ! package        * version    date       lib
##  P abind            1.4-5      2016-07-21 [?]
##  P assertthat       0.2.1      2019-03-21 [?]
##  P bslib            0.2.5      2021-05-12 [?]
##  P callr            3.7.0      2021-04-20 [?]
##  P cellranger       1.1.0      2016-07-27 [?]
##  P cli              2.5.0      2021-04-26 [?]
##  P cluster          2.1.0      2019-06-19 [?]
##  P codetools        0.2-16     2018-12-24 [?]
##  P colorspace       2.0-1      2021-05-04 [?]
##  P cowplot        * 1.1.1      2020-12-30 [?]
##  P crayon           1.4.1      2021-02-08 [?]
##  P data.table       1.14.0     2021-02-21 [?]
##  P DBI              1.1.1      2021-01-15 [?]
##  P deldir           0.2-10     2021-02-16 [?]
##  P desc             1.2.0      2018-05-01 [?]
##  P devtools         2.3.0      2020-04-10 [?]
##  P digest           0.6.27     2020-10-24 [?]
##  P dplyr          * 1.0.6      2021-05-05 [?]
##  P ellipsis         0.3.2      2021-04-29 [?]
##  P evaluate         0.14       2019-05-28 [?]
##  P fansi            0.4.2      2021-01-15 [?]
##  P fastmap          1.1.0      2021-01-25 [?]
##  P feather        * 0.3.5      2019-09-15 [?]
##  P fitdistrplus     1.1-3      2020-12-05 [?]
##  P fs               1.5.0      2020-07-31 [?]
##  P future           1.21.0     2020-12-10 [?]
##  P future.apply     1.7.0      2021-01-04 [?]
##  P generics         0.1.0      2020-10-31 [?]
##  P ggplot2        * 3.3.3      2020-12-30 [?]
##  P ggrepel        * 0.9.1      2021-01-15 [?]
##  P ggridges         0.5.3      2021-01-08 [?]
##  P git2r            0.27.1     2020-05-03 [?]
##  P globals          0.14.0     2020-11-22 [?]
##  P glue           * 1.4.2      2020-08-27 [?]
##  P goftest          1.2-2      2019-12-02 [?]
##  P gridExtra        2.3        2017-09-09 [?]
##  P gtable           0.3.0      2019-03-25 [?]
##  P here           * 0.1        2017-05-28 [?]
##  P hms              1.0.0      2021-01-13 [?]
##  P htmltools        0.5.1.1    2021-01-22 [?]
##  P htmlwidgets      1.5.3      2020-12-10 [?]
##  P httpuv           1.6.1      2021-05-07 [?]
##  P httr             1.4.2      2020-07-20 [?]
##  P ica              1.0-2      2018-05-24 [?]
##  P icytobox       * 1.0.1      2021-12-07 [?]
##  P igraph           1.2.6      2020-10-06 [?]
##  P irlba            2.3.3      2019-02-05 [?]
##  P jquerylib        0.1.4      2021-04-26 [?]
##  P jsonlite         1.7.2      2020-12-09 [?]
##  P KernSmooth       2.23-15    2015-06-29 [?]
##  P knitr            1.33       2021-04-24 [?]
##  P later            1.0.0      2019-10-04 [?]
##  P lattice          0.20-44    2021-05-02 [?]
##  P lazyeval         0.2.2      2019-03-15 [?]
##  P leiden           0.3.7      2021-01-26 [?]
##  P lifecycle        1.0.0      2021-02-15 [?]
##  P listenv          0.8.0      2019-12-05 [?]
##  P lmtest           0.9-38     2020-09-09 [?]
##  P magrittr       * 2.0.1      2020-11-17 [?]
##  P MASS             7.3-54     2021-05-03 [?]
##  P Matrix           1.2-18     2019-11-27 [?]
##  P matrixStats      0.58.0     2021-01-29 [?]
##  P memoise          1.1.0      2017-04-21 [?]
##  P mgcv             1.8-35     2021-04-18 [?]
##  P mime             0.10       2021-02-13 [?]
##  P miniUI           0.1.1.1    2018-05-18 [?]
##  P munsell          0.5.0      2018-06-12 [?]
##  P nlme             3.1-152    2021-02-04 [?]
##  P parallelly       1.25.0     2021-04-30 [?]
##  P patchwork        1.1.1      2020-12-17 [?]
##  P pbapply          1.4-3      2020-08-18 [?]
##  P pillar           1.6.0      2021-04-13 [?]
##  P pkgbuild         1.0.8      2020-05-07 [?]
##  P pkgconfig        2.0.3      2019-09-22 [?]
##  P pkgload          1.0.2      2018-10-29 [?]
##  P plotly           4.9.3      2021-01-10 [?]
##  P plyr             1.8.6      2020-03-03 [?]
##  P png              0.1-7      2013-12-03 [?]
##  P polyclip         1.10-0     2019-03-14 [?]
##  P prettyunits      1.1.1      2020-01-24 [?]
##  P processx         3.5.2      2021-04-30 [?]
##  P promises         1.1.0      2019-10-04 [?]
##  P ps               1.6.0      2021-02-28 [?]
##  P purrr          * 0.3.4      2020-04-17 [?]
##  P R6               2.5.0      2020-10-28 [?]
##  P RANN             2.6.1      2019-01-08 [?]
##  P RColorBrewer   * 1.1-2      2014-12-07 [?]
##  P Rcpp             1.0.6      2021-01-15 [?]
##  P RcppAnnoy        0.0.18     2020-12-15 [?]
##  P readr          * 1.4.0      2020-10-05 [?]
##  P readxl         * 1.3.1      2019-03-13 [?]
##  P remotes          2.1.1      2020-02-15 [?]
##    renv             0.14.0     2021-07-21 [1]
##  P reshape2         1.4.4      2020-04-09 [?]
##  P reticulate       1.20       2021-05-03 [?]
##  P rlang            0.4.11     2021-04-30 [?]
##  P rmarkdown        2.8        2021-05-07 [?]
##  P ROCR             1.0-11     2020-05-02 [?]
##  P rpart            4.1-15     2019-04-12 [?]
##  P rprojroot        2.0.2      2020-11-15 [?]
##  P rstudioapi       0.13       2020-11-12 [?]
##  P rsvd             1.0.3      2020-02-17 [?]
##  P Rtsne            0.15       2018-11-10 [?]
##  P sass             0.4.0      2021-05-12 [?]
##  P scales           1.1.1      2020-05-11 [?]
##  P sctransform      0.3.2      2020-12-16 [?]
##  P sessioninfo      1.1.1      2018-11-05 [?]
##  P Seurat         * 3.2.1      2020-09-07 [?]
##  P shiny            1.6.0      2021-01-25 [?]
##  P spatstat         1.64-1     2020-05-12 [?]
##  P spatstat.data    2.1-0      2021-03-21 [?]
##  P spatstat.utils   2.1-0      2021-03-15 [?]
##  P stringi          1.6.1      2021-05-10 [?]
##  P stringr          1.4.0      2019-02-10 [?]
##  P survival         3.2-11     2021-04-26 [?]
##  P tensor           1.5        2012-05-05 [?]
##  P testrmd          0.0.1.9000 2021-12-06 [?]
##  P testthat         2.3.2      2020-03-02 [?]
##  P tibble           3.1.1      2021-04-18 [?]
##  P tidyr          * 1.1.3      2021-03-03 [?]
##  P tidyselect       1.1.1      2021-04-30 [?]
##  P usethis          1.6.1      2020-04-29 [?]
##  P utf8             1.2.1      2021-03-12 [?]
##  P uwot             0.1.10     2020-12-15 [?]
##  P vctrs            0.3.8      2021-04-29 [?]
##  P viridis        * 0.5.1      2018-03-29 [?]
##  P viridisLite    * 0.4.0      2021-04-13 [?]
##  P withr            2.4.2      2021-04-18 [?]
##  P xfun             0.22       2021-03-11 [?]
##  P xtable           1.8-4      2019-04-21 [?]
##  P yaml             2.2.1      2020-02-01 [?]
##  P zoo              1.8-9      2021-03-09 [?]
##  source                               
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  Github (fungenomics/icytobox@730e8b8)
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  Github (rmflight/testrmd@0735c20)    
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
## 
## [1] /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/renv/library/R-3.6/x86_64-pc-linux-gnu
## [2] /tmp/Rtmpt7dKUV/renv-system-library
## 
##  P ── Loaded and on-disk path mismatch.
```



</details>

The resources requested when this document was last rendered:



```
## #SBATCH --time=01:00:00
## #SBATCH --cpus-per-task=1
## #SBATCH --mem=60G
```




***



<!-- END OF END MATTER -->
