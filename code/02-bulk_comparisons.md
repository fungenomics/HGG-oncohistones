---
title: "02 - Bulk RNAseq/ChIPseq comparisons"
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
## Document index: 02
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
## public/output/02
```

```
## public/figures/02
```



Setting a random seed:



```r
set.seed(100)
```



***



<!-- END OF FRONT MATTER -->


# Overview

Here we interrogate differences in bulk RNAseq as well as H3K27ac and H3K27me3
ChIPseq data between tumor types (H3.1 vs H3.3K27M pontine gliomas, and pontine vs
thalamic H3.3K27M gliomas). This document generates the scatter plots integrating
H3K27ac/H3K27me3/RNAseq, as shown in Figures 2 and 4.

# Libraries



```r
# Load libraries here
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(glue)
library(tibble)
library(ggplot2)
library(ggrepel)
library(stringr)
library(purrr)
library(cowplot)
library(feather)
library(icytobox)

source(here("include/style.R")) # contains palettes & plotting utils
source(here("code/functions/RNAseq.R"))
ggplot2::theme_set(theme_min())
```






# Load metadata

Load the sample metadata for the project:



```r
meta      <- read_tsv(here("data/metadata/metadata_patient_samples_NGS.tsv"))
```






```r
meta_chip <- read_tsv(here("data/metadata/metadata_chip_all.tsv")) %>% 
    right_join(meta, by = c("BioID", "ID_paper", "Material")) %>% 
    mutate(Factor = gsub("-M", "", Factor)) %>%
    filter(grepl("Cell-of-origin", Analyses))
```



# Differential expression (pons vs. thalamus, and H3.1 vs H3.3)


Load the differential expression analyses between H3.1 vs H3.3 and pons vs. thalamic tumors:



```r
# H3.1 has log2FoldChange > 0
info_h31_vs_H33 <- read_tsv(here("data/RNAseq/pipeline_l3/DGE/HGG-H3.1.2K27M-Pons_vs_HGG-H3.3K27M-Pons/info.samples.tsv"))
dge_h31_vs_h33 <- read_tsv(here("data/RNAseq/pipeline_l3/DGE/HGG-H3.1.2K27M-Pons_vs_HGG-H3.3K27M-Pons/diff/Ensembl.ensGene.exon/HGG-H3.3K27M-PonsvsHGG-H3.1.2K27M-Pons.tsv")) %>% 
    separate(ID, into = c("ENSID", "symbol"), sep = ":")
table(info_h31_vs_H33$Group)
```



```
## 
## HGG-H3.1.2K27M-Pons   HGG-H3.3K27M-Pons 
##                   9                  11
```



```r
# Thalamus has log2FoldChange > 0
info_thalamus_vs_pons <- read_tsv(here("data/RNAseq/pipeline_l3/DGE/HGG-H3.3K27M-Thal._vs_HGG-H3.3K27M-Pons/info.samples.tsv"))
dge_thalamus_vs_pons <- read_tsv(here("data/RNAseq/pipeline_l3/DGE/HGG-H3.3K27M-Thal._vs_HGG-H3.3K27M-Pons/diff/Ensembl.ensGene.exon/HGG-H3.3K27M-PonsvsHGG-H3.3K27M-Thal..tsv")) %>% 
    separate(ID, into = c("ENSID", "symbol"), sep = ":")
table(info_thalamus_vs_pons$Group)
```



```
## 
##  HGG-H3.3K27M-Pons HGG-H3.3K27M-Thal. 
##                 11                  6
```



## Prepare tables

We'll assemble 3 tables:

1. raw counts
2. DGE H3.1 vs H3.3
3. DGE thalamus vs pons



```r
# COUNTS -----------------------------------------------------------------------

# get sample order as Table 1
sample_order <- data.table::fread(here("data/metadata/bulk_info.samples.tsv"), data.table = FALSE) %>% .$Nickname

counts_RNA <- read.table(here("data/RNAseq/pipeline_l3/all/counts/Ensembl.ensGene.exon.raw.tsv.gz"),
                         header = T, sep = "\t", check.names = FALSE) %>%
    tibble::rownames_to_column(var = "ID") %>% 
    separate(ID, into = c("ENSID", "symbol"), sep = ":") %>% 
    arrange(ENSID) %>%
    .[, c("ENSID", "symbol", sample_order)] %T>%
    rr_write_tsv(glue("{out}/TABLE_bulk_counts.tsv"),
                 "Raw counts for bulk RNAseq data")
```



```
## ...writing description of TABLE_bulk_counts.tsv to public/output/02/TABLE_bulk_counts.desc
```



```r
# sanity checks
length(colnames(counts_RNA)) - 2
```



```
## [1] 64
```



```r
dim(counts_RNA)
```



```
## [1] 60234    66
```



```r
# DGE --------------------------------------------------------------------------
dge_h31_vs_h33 %>% arrange(ENSID) %>%
    rr_write_tsv(glue("{out}/TABLE_dge_H3.1_vs_H3.3.tsv"),
                 "Differential gene expression analysis for H3.1K27M pons vs H3.3K27M pons (LFC > 0 indicates up in H3.1K27M)")
```



```
## ...writing description of TABLE_dge_H3.1_vs_H3.3.tsv to public/output/02/TABLE_dge_H3.1_vs_H3.3.desc
```



```r
dge_thalamus_vs_pons %>% arrange(ENSID) %>%
    rr_write_tsv(glue("{out}/TABLE_dge_thal_vs_pons.tsv"),
                 "Differential gene expression analysis for H3.3K27M thalamic vs H3.3K27M pons (LFC > 0 indicates up in thalamic)")
```



```
## ...writing description of TABLE_dge_thal_vs_pons.tsv to public/output/02/TABLE_dge_thal_vs_pons.desc
```



```r
# sanity checks
dim(dge_h31_vs_h33)
```



```
## [1] 60234     8
```



```r
dim(dge_thalamus_vs_pons)
```



```
## [1] 60234     8
```



# ChIPseq

## Load data



```r
# slow
data_chip <- read_xlsx(here("data/ChIPseq/quantifications_@mhulswit/20210521_H3.3H3.3EZHIP_H3K27ac_H3K27me3_RefSeq_Promoters_SJ.xlsx"), skip = 1)
```



Let's separate out the data into the summary stats (group medians and z-scores)
and per-sample quantifications, and convert to tidy (long) format. Starting
with the per-sample quantifications:



```r
data_chip_tidy <- data_chip %>%
    dplyr::select(1:100, matches("Median_"),
                  # only keep relevant Z-scores
                  Z_H3.3K27M_Thalamus_H3.3K27M_Pons_H3K27ac,
                  Z_H3.3K27M_Pons_H3.1K27M_Pons_H3K27ac,
                  Z_H3.3K27M_Thalamus_H3.3K27M_Pons_H3K27me3,
                  Z_H3.3K27M_Pons_H3.1K27M_Pons_H3K27me3) %>%
    gather(dataset, value, 8:ncol(.)) %>%
    mutate(Mark = case_when(
        grepl("K27ac", dataset) ~ "H3K27ac",
        grepl("K27me3", dataset) ~ "H3K27me3"
    )) %>%
    mutate(Sample = gsub("_Count|_RPKM", "", dataset)) %>% 
    left_join(meta_chip %>% mutate(bw = basename(Path_bw_internal)) %>% dplyr::select(ID_paper, bw, Replicate),
              # column name in ChIP quantification matches bw filename
              by = c("Sample" = "bw"))  %>% 
    filter(!grepl("NormalPons", Sample)) %>% 
    # only keep the relevant Z
    filter(!is.na(ID_paper) | grepl("^Median|^Z", Sample)) %>% 
    # tidy
    mutate(Sample = case_when(
        grepl("Median_PFA_PF", Sample) ~ gsub("PFA_PF", "EZHIP_PFA", Sample),
        grepl("Median_H3.1K27M_PF", Sample) ~ gsub("PF", "PFA", Sample),
        TRUE ~ Sample))

# sanity check: ensure that no samples present in this dataset are *NOT* included
# in the ChIPseq metadata file
base::setdiff(unique(data_chip_tidy$ID_paper), unique(meta_chip$ID_paper))
```



```
## [1] NA
```



```r
save(data_chip_tidy, file = glue("{out}/data_chip_tidy.Rda"))

data_chip_tidy %>%
    distinct(Mark, Sample, ID_paper) %>%
    left_join(meta %>% select(ID_paper, Group)) %>%
    count(Mark, Group)
```



```
## Joining, by = "ID_paper"
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Mark"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Group"],"name":[2],"type":["chr"],"align":["left"]},{"label":["n"],"name":[3],"type":["int"],"align":["right"]}],"data":[{"1":"H3K27ac","2":"HGG-H3.1/2K27M-Pons","3":"16"},{"1":"H3K27ac","2":"HGG-H3.3K27M-Pons","3":"16"},{"1":"H3K27ac","2":"HGG-H3.3K27M-Thal.","3":"5"},{"1":"H3K27ac","2":"HGG-H3WT-Cortex","3":"3"},{"1":"H3K27ac","2":"HGG-H3WT-NA","3":"2"},{"1":"H3K27ac","2":"PFA-EZHIP-PF","3":"4"},{"1":"H3K27ac","2":"PFA-H3.1K27M-PF","3":"1"},{"1":"H3K27ac","2":"NA","3":"9"},{"1":"H3K27me3","2":"HGG-H3.1/2K27M-Pons","3":"3"},{"1":"H3K27me3","2":"HGG-H3.3K27M-Pons","3":"6"},{"1":"H3K27me3","2":"HGG-H3.3K27M-Thal.","3":"4"},{"1":"H3K27me3","2":"HGG-H3WT-Cortex","3":"4"},{"1":"H3K27me3","2":"PFA-EZHIP-PF","3":"4"},{"1":"H3K27me3","2":"PFA-H3.1K27M-PF","3":"1"},{"1":"H3K27me3","2":"NA","3":"9"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# make table
data_chip_tidy %>%
    mutate(ID_paper = case_when(
        grepl("^Z|^Median", Sample) ~ Sample,
        TRUE ~ paste0(ID_paper, "-", Replicate, "__", Mark)
    )) %>% 
    dplyr::select(ID, gene_symbol = name, chr, start, end, strand, ID_paper, value) %>%
    pivot_wider(names_from = ID_paper, values_from = value) %>% 
    arrange(chr, start) %>%
    rr_write_tsv(glue("{out}/TABLE_promoter_H3K27ac_H3K27me3_per_sample.tsv"),
                 "Promoter H3K27ac and H3K27me3 enrichment (RPKM) in each sample, and summary stats")
```



```
## ...writing description of TABLE_promoter_H3K27ac_H3K27me3_per_sample.tsv to public/output/02/TABLE_promoter_H3K27ac_H3K27me3_per_sample.desc
```



Tidy summary stats:



```r
data_chip_stats <- data_chip_tidy %>% 
    filter(grepl("Z_", Sample)) %>% 
    rename(statistic = Sample) %>% 
    select(-ID_paper, -Replicate, -Mark, -dataset) %>% 
    group_by(chr, strand, name, statistic) %>% 
    top_n(1, `Exon Length`) %>% 
    ungroup() %>% 
    pivot_wider(names_from = statistic, values_from = value)

save(data_chip_stats, file = glue("{out}/data_chip_stats.Rda"))
```




## Scatterplot for H3.1 vs H3.3 pons

Here, we'll combine the DGE data with each of different source of ChIP-seq analysis,
following code used to do a similar analysis for H3.3G34 gliomas here: https://github.com/fungenomics/G34-gliomas/blob/master/bulk_transcriptome_epigenome/analysis/03-ChIPseq.Rmd#L151 (from [Chen et al, 2020](https://www.cell.com/cell/fulltext/S0092-8674(20)31529-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867420315294%3Fshowall%3Dtrue)).

The input for the scatter plot is essentially a data frame where, for each gene,
we have:

* A score for differential expression (e.g. log2FoldChange from DESeq2 differential expression analysis)
* A score for differential H3K27me3 (e.g. Z-score for promoter H3K27me3)
* A score for differential H3K27ac (e.g. Z-score for promoter H3K27ac)



```r
dge_h31_vs_h33_and_chip <- data_chip_stats %>%
    select(name, H3K27ac = Z_H3.3K27M_Pons_H3.1K27M_Pons_H3K27ac, H3K27me3 = Z_H3.3K27M_Pons_H3.1K27M_Pons_H3K27me3) %>%
    distinct() %>% 
    mutate(H3K27ac = -as.numeric(H3K27ac),
           H3K27me3 = -as.numeric(H3K27me3)) %>% 
    left_join(dge_h31_vs_h33, by = c("name" = "symbol")) %>% 
    filter(baseMean > 100) %>% 
    # remove mitochondrial and Y-RNA
    filter(!grepl("Y_RNA|^MT_", name))

dim(dge_h31_vs_h33_and_chip)
```



```
## [1] 14954    10
```






```r
# get limits of z-scores for x/y axes
xlim <- range(dge_h31_vs_h33_and_chip$H3K27ac)
ylim <- range(dge_h31_vs_h33_and_chip$H3K27me3)

# all genes
to_label <- dge_h31_vs_h33_and_chip %>%
    filter((H3K27ac > 1 | H3K27ac < -0.5) |
               (H3K27me3 > 2.5 | H3K27me3 < -3))
nrow(to_label)
```



```
## [1] 241
```



```r
dge_h31_vs_h33_and_chip %>% 
    arrange(abs(log2FoldChange)) %>%
    mutate(name = factor(name, levels = unique(.$name))) %>%
    rr_ggplot(aes(x = H3K27ac, y = H3K27me3), plot_num = 1) +
    geom_hline(yintercept = c(-0.5, 0.5), size = 0.5, colour = "gray90") +
    geom_vline(xintercept = c(-0.5, 0.5), size = 0.5, colour = "gray90") +
    geom_point(aes(colour = log2FoldChange, size = -log10(padj)), alpha = 0.8) +
    scale_colour_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0)  +
    geom_text_repel(aes(label = name), to_label,
                    inherit.aes = TRUE, size = 2, segment.color = "gray70") +
    xlab("Z score H3K27ac (H3.1K27M pons / H3.3K27M pons)") +
    ylab("Z score H3K27me3 (H3.1K27M pons / H3.3K27M pons)")
```



```
## ...writing source data of ggplot to public/figures/02/scatterplot_h31_vs_h33-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/02/scatterplot_h31_vs_h33-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/02/scatterplot_h31_vs_h33...*]~</span>


## Scatterplot for H3.3 pons vs H3.3 thalamus



```r
dge_pons_vs_thal_and_chip <- data_chip_stats %>%
    select(name,
           H3K27ac = Z_H3.3K27M_Thalamus_H3.3K27M_Pons_H3K27ac,
           H3K27me3 = Z_H3.3K27M_Thalamus_H3.3K27M_Pons_H3K27me3) %>%
    distinct() %>% 
    mutate(H3K27ac = as.numeric(H3K27ac),
           H3K27me3 = as.numeric(H3K27me3)) %>% 
    left_join(dge_thalamus_vs_pons, by = c("name" = "symbol")) %>% 
    filter(baseMean > 100) %>% 
    # remove mitochondrial and Y-RNA
    filter(!grepl("Y_RNA|^MT_", name))

dim(dge_pons_vs_thal_and_chip)
```



```
## [1] 14809    10
```






```r
# get limits of z-scores for x/y axes
xlim <- range(dge_pons_vs_thal_and_chip$H3K27ac)
ylim <- range(dge_pons_vs_thal_and_chip$H3K27me3)

to_label <- dge_pons_vs_thal_and_chip %>%
    filter((H3K27ac > 0.5 | H3K27ac < -0.2) &
               (H3K27me3 > 1 | H3K27me3 < -0.5) & abs(log2FoldChange) > 2) %>% 
    filter(-log10(padj) > 4 | log2FoldChange > 0)
nrow(to_label)
```



```
## [1] 19
```



```r
dge_pons_vs_thal_and_chip %>% 
    arrange(abs(log2FoldChange)) %>%
    mutate(name = factor(name, levels = unique(.$name))) %>%
    rr_ggplot(aes(x = H3K27ac, y = H3K27me3), plot_num = 1) +
    geom_hline(yintercept = c(-0.5, 0.5), size = 0.5, colour = "gray90") +
    geom_vline(xintercept = c(-0.5, 0.5), size = 0.5, colour = "gray90") +
    geom_point(aes(colour = log2FoldChange, size = -log10(padj)), alpha = 0.8) +
    scale_colour_gradient2(low = "navy", mid = "white", high = "red", midpoint = 0)  +
    geom_text_repel(aes(label = name), to_label,
                    inherit.aes = TRUE, size = 2, segment.color = "gray70") +
    xlab("Z score H3K27ac (H3.3K27M thalamus / H3.3K27M pons)") +
    ylab("Z score H3K27me3 (H3.3K27M thalamus / H3.3K27M pons)")
```



```
## ...writing source data of ggplot to public/figures/02/scatterplot_thal_vs_pons-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/02/scatterplot_thal_vs_pons-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/02/scatterplot_thal_vs_pons...*]~</span>


<!-- END MATTER, insert reproducibility info -->




***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:



```
## 2022-06-29 10:10:41
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
##  P Seurat           3.2.1      2020-09-07 [?]
##  P shiny            1.6.0      2021-01-25 [?]
##  P spatstat         1.64-1     2020-05-12 [?]
##  P spatstat.data    2.1-0      2021-03-21 [?]
##  P spatstat.utils   2.1-0      2021-03-15 [?]
##  P stringi          1.6.1      2021-05-10 [?]
##  P stringr        * 1.4.0      2019-02-10 [?]
##  P survival         3.2-11     2021-04-26 [?]
##  P tensor           1.5        2012-05-05 [?]
##  P testrmd          0.0.1.9000 2021-12-06 [?]
##  P testthat         2.3.2      2020-03-02 [?]
##  P tibble         * 3.1.1      2021-04-18 [?]
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
## [2] /tmp/Rtmp21MFFn/renv-system-library
## 
##  P ── Loaded and on-disk path mismatch.
```



</details>

The resources requested when this document was last rendered:



```
## #SBATCH --time=00:20:00
## #SBATCH --cpus-per-task=1
## #SBATCH --mem=10G
```




***



<!-- END OF END MATTER -->
