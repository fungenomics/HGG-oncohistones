---
title: "01A - Prepare human fetal thalamus data"
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
## Document index: 01A
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
## public/output/01A
```

```
## public/figures/01A
```



Setting a random seed:



```r
set.seed(100)
```



***



<!-- END OF FRONT MATTER -->


# Overview

In this analysis, we'll prepare the reference from the human fetal thalamus data
which has been processed with our single-cell workflow, and gather the
required formats and outputs needed for comparisons to tumors:

* labels from the original publications, as well as manual additions/corrections based on known markers
* Seurat objects, per-sample and integrated across samples
* gene signatures, per-cluster
* mean expression profiles, per-cluster

This data comes from [Bhaduri et al](http://dx.doi.org/10.1038/s41586-021-03910-8), and [Eze et al](http://dx.doi.org/10.1038/s41593-020-00794-1).


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
library(pvclust)
library(dendextend)
library(feather)
library(Seurat)
library(icytobox)
library(ggrepel)

source(here("include/style.R"))
source(here("code/functions/scRNAseq.R"))
ggplot2::theme_set(theme_min())
```




# Overview of individual samples

Each sample has a folder in the pipeline at:

```
data/scRNAseq/pipeline_human_fetal/
|-- sample1/
--- |-- preprocessing # output of the SC pipeline
----|-- correlations  # output of projections

```

## Metadata

Load per sample metadata:



```r
info_samples <- read_tsv(here("data/metadata/metadata_human_fetal_10X.tsv")) %>% 
    filter(Region == "Thalamus") %>% 
    filter(grepl("GW", Age))

palette_sample <- info_samples %>% select(Sample, Sample_color) %>% tibble::deframe()
palette_alias  <- info_samples %>% select(Alias, Sample_color) %>% tibble::deframe()

ages <- info_samples$Age %>% unique %>% sort()
palette_age <- colorRampPalette(brewer.pal(8, "YlGnBu"))(n = length(ages))
names(palette_age) <- ages
```



Load the labels from the authors of the original studies:



```r
metadata_pub_gw <- read_excel(here("data/scRNAseq/references/Bhaduri2021_fetal_brain_GW/GW_sample_labels_Nature_2021_Supp_Table_1.xlsx")) %>% 
    filter(structure == "thalamus")

# get sample names in the same format
samples_pub <- metadata_pub_gw %>% filter(structure == "thalamus") %>%
    pull(cell.name) %>%
    {gsub('.{17}$', '', .)} %>%
    unique()

all(info_samples$Sample_public %in% samples_pub)
```



```
## [1] TRUE
```



## QC

Verify QC of each sample, to determine whether to include all:



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

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/01A/qc-1.png)<!-- -->

## Prep objects

For each sample:

* load the Seurat object
* count the number of cells in each cluster, and give each cluster a label prefixed w/ sample name
* add the labels from the published metadata



```r
prep_seurat <- function(path_to_sample, sample_public_id) {
    
    message("@ ", basename(path_to_sample))
    
    # load the object as output by the SC workflow
    load(file.path(path_to_sample, "preprocessing/seurat.Rda"))
    
    # get the labels from the published metadata and add them to the Seurat object
    public_labels <- seurat@meta.data %>%
        tibble::rownames_to_column(var = "cell.name") %>%
        rowwise() %>%
        mutate(barcode_public = paste0(sample_public_id, "_", gsub("-1", "", cell.name, fixed = TRUE))) %>%
        left_join(metadata_pub_gw %>%
                      select(barcode_public = cell.name, Label_public = `cell type`), by = c("barcode_public")) %>% 
        select(cell.name, Label_public) %>%
        tibble::column_to_rownames(var = "cell.name")
    
    seurat <- AddMetaData(seurat,
                          metadata = public_labels,
                          col.name = c("Label_published"))
    
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
    
    # get the most common projection from the public data
    cluster_labels_published <- seurat@meta.data %>% rowwise() %>% 
        group_by(Cluster_number) %>% 
        mutate(N_cells_in_cluster = n()) %>% 
        filter(!is.na(Label_published)) %>% 
        group_by(Cluster_number, N_cells_in_cluster, Label_published) %>%
        summarize(N_celltype_in_cluster = n()) %>% 
        mutate(Prop_celltype_in_cluster = round(N_celltype_in_cluster/N_cells_in_cluster, 2)) %>%
        group_by(Cluster_number) %>% top_n(1, Prop_celltype_in_cluster) %>%
        slice(1) %>%
        filter(Prop_celltype_in_cluster >= 0.25) %>% 
        select(Cluster_number, N_celltype_in_cluster, Prop_celltype_in_cluster, Label_published)
    
    info_samples_cluster <- info_samples_cluster %>%
        left_join(cluster_labels_published, by = "Cluster_number") %>%
        mutate(ID_20220321 = paste0(Cluster, "_", Label_published)) %>% 
        dplyr::rename(Prop_cluster_with_published_label = Prop_celltype_in_cluster) %>% 
        select(-N_celltype_in_cluster)
    
    return(info_samples_cluster)
    
}

# gather all cluster summary info
info_clusters <- map_dfr(1:nrow(info_samples), ~ prep_seurat(info_samples[.x, ]$Path, info_samples[.x, ]$Sample_public))
write_tsv(info_clusters, glue("{out}/info_clusters.tsv"))
```





```r
info_clusters <- read_tsv(glue("{out}/info_clusters.tsv"))
```



## Sample exclusions

We'll exclude the following samples:

* GW19_thalamus, due to low # of cells
* GW14_thalamus, due to low structure
* GW18_2_ventral_thalamus, due to high mito content
* GW20_34_dorsalthalamus, due to low coverage

## Per-cluster QC



```r
samples_drop <- c("GW19_thalamus", "GW14_thalamus", "GW18_2_ventral_thalamus", "GW20_34_dorsalthalamus")

info_clusters %>% nrow()
```



```
## [1] 233
```



```r
info_clusters2 <- info_clusters %>% 
    mutate(Note_sample = case_when(
        Sample == "GW19_thalamus" ~ "Sample excluded, low number of cells",
        Sample == "GW14_thalamus" ~ "Sample excluded, low structure",
        Sample == "GW18_2_ventral_thalamus" ~ "Sample excluded, high mitochondrial content",
        Sample == "GW20_34_dorsalthalamus" ~ "Sample excluded, low coverage",
        TRUE ~ "OK"),
        Note_cluster = case_when(
            Note_sample != "OK" ~ "Sample excluded",
            N_cells < 30 ~ "Less than 30 cells",
            nCount_RNA < 1000 ~ "Low number of UMIs",
            percent.ribo > 50 ~ "High ribosomal content",
            TRUE ~ "OK"
        ))

info_clusters2 %>% filter(Note_sample == "OK" & Note_cluster == "OK") %>% nrow()
```



```
## [1] 167
```



```r
write_tsv(info_clusters2, glue("{out}/info_clusters2.tsv"))

# total # of cells
sum(info_clusters2$N_cells)
```



```
## [1] 90226
```



## QC tables

Produce a QC table for the supplementary materials:



```r
TABLE_thalamus_qc <- info_samples %>% 
    left_join(qc_stats, by = "Sample") %>% 
    mutate(Note_sample = case_when(
        Sample == "GW19_thalamus" ~ "Sample excluded, low number of cells",
        Sample == "GW14_thalamus" ~ "Sample excluded, low structure",
        Sample == "GW18_2_ventral_thalamus" ~ "Sample excluded, high mitochondrial content",
        Sample == "GW20_34_dorsalthalamus" ~ "Sample excluded, low coverage",
        TRUE ~ "OK")) %>% 
    mutate(Publication = case_when(
        grepl("CS", Sample) ~ "Eze et al, Nature Neuroscience, 2021",
        grepl("GW", Sample) ~ "Bhaduri et al, Nature, 2021"
    )) %>% 
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

rr_write_tsv(TABLE_thalamus_qc,
             glue("{out}/TABLE_thalamus_QC.tsv"),
             "Summary of sample info and QC for thalamus samples")
```



```
## ...writing description of TABLE_thalamus_QC.tsv to public/output/01A/TABLE_thalamus_QC.desc
```



# Gather gene signatures

First, load all cluster markers and filter out markers for samples/clusters that we'll exclude:



```r
clusters_keep <- info_clusters2 %>% filter(Note_cluster == "OK") %>% pull(Cluster)

cluster_markers <- map_dfr(1:nrow(info_samples), ~ read_tsv(here(file.path(info_samples[.x, ]$Path, "preprocessing/cluster_markers.tsv"))) %>% 
                               tibble::add_column(Sample = info_samples[.x, ]$Sample)) %>% 
    # get positive markers
    filter(avg_logFC > 0) %>%
    mutate(cluster = as.numeric(cluster)) %>% 
    # add the cluster/sample info
    left_join(info_clusters2 %>% select(Sample, Cluster_number, Cluster) %>% mutate(Cluster_number = as.numeric(as.character(Cluster_number))),
              by = c("Sample" = "Sample", "cluster" = "Cluster_number")) %>% 
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
## [1] 167
```



```r
any(cluster_markers$Sample %in% samples_drop)
```



```
## [1] FALSE
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
{"columns":[{"label":["gene"],"name":[1],"type":["chr"],"align":["left"]},{"label":["pct_diff"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["p_val"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["avg_logFC"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["pct.1"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["pct.2"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["p_val_adj"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["Cluster_number"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["Sample"],"name":[9],"type":["chr"],"align":["left"]},{"label":["Cluster"],"name":[10],"type":["chr"],"align":["left"]}],"data":[{"1":"GADD45B","2":"0.537","3":"5.648400e-227","4":"1.2249812","5":"0.864","6":"0.327","7":"1.111605e-222","8":"0","9":"GW18_thalamus","10":"GW18-T_0"},{"1":"JUNB","2":"0.483","3":"6.641893e-192","4":"0.9973744","5":"0.910","6":"0.427","7":"1.307125e-187","8":"0","9":"GW18_thalamus","10":"GW18-T_0"},{"1":"C1orf61","2":"0.059","3":"1.497058e-190","4":"0.7234438","5":"1.000","6":"0.941","7":"2.946210e-186","8":"0","9":"GW18_thalamus","10":"GW18-T_0"},{"1":"IER2","2":"0.400","3":"1.196228e-184","4":"0.9341488","5":"0.963","6":"0.563","7":"2.354177e-180","8":"0","9":"GW18_thalamus","10":"GW18-T_0"},{"1":"PTPRZ1","2":"0.162","3":"4.318492e-169","4":"0.7589503","5":"1.000","6":"0.838","7":"8.498793e-165","8":"0","9":"GW18_thalamus","10":"GW18-T_0"},{"1":"FABP7","2":"0.124","3":"1.561269e-157","4":"0.7030516","5":"0.999","6":"0.875","7":"3.072578e-153","8":"0","9":"GW18_thalamus","10":"GW18-T_0"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# save
write_tsv(cluster_markers, glue("{out}/cluster_markers.tsv"))
```



We'll also want to remove erythrocyte & melanocyte clusters, which are quite biologically distinct:



```r
erythro_markers <- c("HBB", "HBA1", "HBA2", "HBZ", "HBG1")
(erythro_count <- cluster_markers %>% filter(gene %in% erythro_markers) %>% pull(Cluster) %>% table())
```



```
## .
##  GW19-T2_20 GW20-DT2_13  GW22-T1_14  GW22-T2_16   GW25-T_18    GW25-T_8 
##           5           1           4           4           5           1
```



```r
# remove the clusters with 4-5 of these markers:
erythro_clusters <- names(erythro_count)[erythro_count >= 4]

melano_markers <- c("DCT", "MITF", "PMEL")
(melano_count <- cluster_markers %>% filter(gene %in% melano_markers) %>% pull(Cluster) %>% table())
```



```
## .
## GW22-T1_15 GW22-T2_17  GW25-T_15 
##          3          3          1
```



```r
# remove the clusters with all three of these:
melano_clusters <- names(melano_count)[melano_count == 3]
```



# Complete and refine cluster labels

We will manually update the cluster labels (provided by the author) in order to ensure all populations 
of interest are labelled. We will use the following markers:

* Astrocytes: FABP7, S100B, CLU, AQP4
* Neurons: STMN2, TUBB, SYT1
* Microglia: C1QC, LY86
* Ependymal: FOXJ1,  DNAH10, KIF19, DNAH12, TEK1
* OPC: PDGFRA, OLIG1, OLIG2
* Glial progenitors: SOX2, OLIG1/2

After using these markers to update the labels outside of R, we load in the updated labels:



```r
labels_manual <- read_tsv(here("data/scRNAseq/references/Bhaduri2021_fetal_brain_GW/2022-03-21-human_fetal_thalamus_labels_manual_SJ.tsv")) %>% 
    select(Cluster, Label_manual)

info_clusters3 <- info_clusters2 %>% 
    left_join(labels_manual) %>% 
    mutate(Note_cluster = case_when(
        grepl("Outlier", Label_published) ~ "Exclude cluster, outlier",
        Label_manual == "Erythrocytes" ~ "Exclude cluster, erythrocytes",
        Label_manual == "Melanocytes" ~ "Exclude cluster, melanocytes",
        TRUE ~ Note_cluster
    )) %>% 
    mutate(ID_20220321 = case_when(
        Note_sample != "OK" | Note_cluster != "OK" ~ paste0(Cluster, "_EXCLUDE"),
        is.na(Label_manual) ~ ID_20220321,
        !is.na(Label_manual) ~ paste0(Cluster, "_", Label_manual)
    )) %>% 
    select(-Label_manual)

write_tsv(info_clusters3, glue("{out}/info_clusters3.tsv"))

clusters_keep <- info_clusters3 %>% filter(!grepl("EXCLUDE", ID_20220321)) %>% pull(Cluster)
length(clusters_keep)    
```



```
## [1] 160
```




# Make reference for cell type projection

In this section, we prepare the reference matrix for cell type projection,
containing the mean expression of each gene in each cluster.

## Load data



```r
load(here("data/scRNAseq/integrations/human_fetal_thalamus/output/seurat_joint.Rda"))
```




## Mean expression

Calculate mean expression per cluster:



```r
Idents(seurat_joint) <- "Cluster"
hf_thalamus_meanexp <- AverageExpression(seurat_joint, return.seurat = TRUE, assays = "RNA", verbose = FALSE)
hf_thalamus_meanexp <- GetAssayData(hf_thalamus_meanexp)

# peek at the corner
dim(hf_thalamus_meanexp)
```



```
## [1] 23955   217
```



```r
hf_thalamus_meanexp[1:5, 1:5]
```



```
##                  CS19-T_1    CS19-T_0    CS19-T_4   CS19-T_6   CS19-T_2
## FO538757.2    0.642409234 0.671287431 0.659322505 0.72569247 0.67059287
## AP006222.2    0.379104823 0.457272354 0.424646100 0.53114856 0.58883959
## RP4-669L17.10 0.008780798 0.006391665 0.008229348 0.01997191 0.00430778
## RP11-206L10.9 0.116477049 0.102566086 0.039760352 0.15006909 0.03784705
## LINC00115     0.018969886 0.016771483 0.030140318 0.00000000 0.03769330
```



```r
hf_thalamus_meanexp <- hf_thalamus_meanexp[, clusters_keep]

# sanity checks
dim(hf_thalamus_meanexp)
```



```
## [1] 23955   160
```



```r
length(clusters_keep)
```



```
## [1] 160
```



```r
save(hf_thalamus_meanexp, file = glue("{out}/hf_thalamus_mean_expression.Rda"))
saveRDS(hf_thalamus_meanexp, glue("{out}/hf_thalamus_mean_expression.Rds"))
```



## Annotation

Make the annotation for the cell type projection scripts:



```r
info_clusters3 %>%
    filter(Cluster %in% clusters_keep) %>% 
    select(Cell_type = Cluster, ID_20220321) %>% 
    mutate(Color = case_when(
        grepl("Proliferating OPC", ID_20220321) ~ palette_type[["Proliferating OPC"]],
        grepl("OPC", ID_20220321) ~ palette_type[["OPC"]],
        grepl("Oligo", ID_20220321) ~ palette_type[["Oligodendrocytes"]],
        grepl("[Nn]euron", ID_20220321) ~ palette_type[["Neurons"]],
        grepl("Astro", ID_20220321) ~ palette_type[["Astrocytes"]],
        grepl("Ependymal", ID_20220321) ~ palette_type[["Ependymal"]],
        grepl("Micro", ID_20220321) ~ palette_type[["Immune"]],
        grepl("Endo", ID_20220321) ~ palette_type[["Vascular & other"]],
        grepl("Dividing|RG", ID_20220321) ~ palette_type[["RGC"]],
    )) %>% 
    select(-ID_20220321) %>%
    write_tsv(glue("{out}/hf_thalamus_annotation.tsv"))
```



<!-- END MATTER, insert reproducibility info -->




***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:



```
## 2022-09-08 14:54:28
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
##  P dendextend     * 1.14.0     2020-08-26 [?]
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
##  P pvclust        * 2.2-0      2019-11-19 [?]
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
## [2] /tmp/RtmpRjPYXN/renv-system-library
## 
##  P ── Loaded and on-disk path mismatch.
```



</details>

The resources requested when this document was last rendered:



```
## #SBATCH --time=01:00:00
## #SBATCH --cpus-per-task=1
## #SBATCH --mem=50G
```




***



<!-- END OF END MATTER -->
