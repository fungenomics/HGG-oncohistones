---
title: "05 - Summary figures for extended mouse atlas"
author: "Selin Jessa [[selin.jessa@mail.mcgill.ca](mailto:selin.jessa@mail.mcgill.ca)]"
date: "12 September, 2022"
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
## Document index: 05
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
## public/output/05
```

```
## public/figures/05
```



Setting a random seed:



```r
set.seed(100)
```



***



<!-- END OF FRONT MATTER -->


# Overview

This document summarizes the extended mouse scRNAseq developmental atlas and generates
the dendrograms over single-cell clusters for each brain region, as shown in 
Extended Data Figure 1.

# Libraries



```r
# Load libraries here
library(here)
library(tidyr)
library(dplyr)
library(readr)
library(readxl)
library(glue)
library(purrr)
library(pvclust)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(dendextend)
library(ggrastr)
library(Seurat)

source(here("include/style.R"))
source(here("code/functions/scRNAseq.R"))
ggplot2::theme_set(theme_min())

# set seed to reproduce earlier analysis
set.seed(23)
```



# Load data



```r
atlas_path <- read_lines(here("data/scRNAseq/references/mouse_atlas_extended/atlas_path_hydra.tsv"))

mean_expression_profile <- readRDS(file.path(atlas_path, "data/joint_mouse_extended/mean_expression_per_cluster.Rds"))

cluster_signatures <- readRDS(file.path(atlas_path, "data/joint_mouse_extended/joint_mouse_extended.signatures_ID_20210710.Rds"))

mouse_signatures <- cluster_signatures$mm_sym
```





# Compute dendrograms

Compute dendrograms based on correlation between mean transcriptome of each cluster. We use [pvclust](https://cran.r-project.org/web/packages/pvclust/index.html) for constructing dendrograms.



```r
# get the unique signature genes from mouse only
atlas_unique_genes_mm <- mouse_signatures %>% unlist %>% unique
f_signatures <- mouse_signatures[grepl("F-", names(mouse_signatures))] %>% unlist %>% unique
p_signatures <- mouse_signatures[grepl("P-.+_", names(mouse_signatures))] %>% unlist %>% unique

meanexp_mm_uniq <- mean_expression_profile[atlas_unique_genes_mm] 
rownames(meanexp_mm_uniq) <- mean_expression_profile$ID_20210710

# subset the gene signatures to the clusters used in the dendrogram
meanexp_f_uniq <- meanexp_mm_uniq %>% filter(grepl("F-", rownames(meanexp_mm_uniq))) %>% .[f_signatures]
meanexp_p_uniq <- meanexp_mm_uniq %>% filter(grepl("P-.+_", rownames(meanexp_mm_uniq))) %>% .[p_signatures]

# transpose
meanexp_f_uniq <- as.data.frame(t(meanexp_f_uniq))
meanexp_p_uniq <- as.data.frame(t(meanexp_p_uniq))

# take correlations
cor_f <- cor(meanexp_f_uniq, meanexp_f_uniq, method = "spearman", use = "complete.obs")
cor_p <- cor(meanexp_p_uniq, meanexp_p_uniq, method = "spearman", use = "complete.obs")

# remove NAs
meanexp_f_uniq_no_NA <- meanexp_f_uniq
meanexp_f_uniq_no_NA[is.na(meanexp_f_uniq_no_NA)] <- 0

meanexp_p_uniq_no_NA <- meanexp_p_uniq
meanexp_p_uniq_no_NA[is.na(meanexp_p_uniq_no_NA)] <- 0

result_f <- pvclust(meanexp_f_uniq_no_NA, method.dist = spearman,
                    method.hclust = "complete", nboot = 100)
```



```
## Bootstrap (r = 0.5)... Done.
## Bootstrap (r = 0.6)... Done.
## Bootstrap (r = 0.7)... Done.
## Bootstrap (r = 0.8)... Done.
## Bootstrap (r = 0.9)... Done.
## Bootstrap (r = 1.0)... Done.
## Bootstrap (r = 1.1)... Done.
## Bootstrap (r = 1.2)... Done.
## Bootstrap (r = 1.3)... Done.
## Bootstrap (r = 1.4)... Done.
```



```r
result_p <- pvclust(meanexp_p_uniq_no_NA, method.dist = spearman,
                    method.hclust = "complete", nboot = 100)
```



```
## Bootstrap (r = 0.5)... Done.
## Bootstrap (r = 0.6)... Done.
## Bootstrap (r = 0.7)... Done.
## Bootstrap (r = 0.8)... Done.
## Bootstrap (r = 0.9)... Done.
## Bootstrap (r = 1.0)... Done.
## Bootstrap (r = 1.1)... Done.
## Bootstrap (r = 1.2)... Done.
## Bootstrap (r = 1.3)... Done.
## Bootstrap (r = 1.4)... Done.
```






```r
dend_f <- as.dendrogram(result_f)
dend_order_f <- colnames(cor_f)[order.dendrogram(dend_f)]

dend_p <- as.dendrogram(result_p)
dend_order_p <- colnames(cor_p)[order.dendrogram(dend_p)]

saveRDS(dend_order_f,
        file = glue("{out}/dendrogram_order_joint_extended_forebrain.Rds"))
saveRDS(dend_order_p,
        file = glue("{out}/dendrogram_order_joint_extended_pons.Rds"))

f_palette <- palette_joint_mouse_extended_full[dend_order_f]
p_palette <- palette_joint_mouse_extended_full[dend_order_p]

dend_f <- dend_f %>% set("labels_col", f_palette) %>% set("labels_cex", 0.3)
dend_p <- dend_p %>% set("labels_col", p_palette) %>% set("labels_cex", 0.3)

par(mar = c(20,2,2,2))
plot(dend_f)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/05/plot_dend-1.png)<!-- -->

```r
plot(dend_p)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/05/plot_dend-2.png)<!-- -->

```r
dev.off()
```



```
## null device 
##           1
```






# Generate tables

Per-sample data:



```r
mouse_new_sample_qc <- read_xlsx(here("data/scRNAseq/references/mouse_atlas_extended/new_timepoints_stats.xlsx"))

mouse_sample_summary <- read_xlsx(here("data/metadata/2021-10-20_Omega_table-singlecell.xlsx")) %>% 
    filter(Diagnosis_1 == "Normal Mouse Brain") %>% 
    filter(!(Sample %in% c("PT_CT_0F", "fresh_mouse_brain_E18"))) %>% 
    left_join(mouse_new_sample_qc, by = c("Aliases" = "Alias")) %>% 
    rename(Alias = Aliases)
```



```
## New names:
## * `library yield (nM)` -> `library yield (nM)...33`
## * `copies/ul` -> `copies/ul...34`
## * `copies/ul` -> `copies/ul...35`
## * `library yield (nM)` -> `library yield (nM)...36`
```



```r
# get QC stats after filtering
mouse_sample_qc_after_filt <- map_dfr(mouse_sample_summary$Alias,
         ~ data.table::fread(glue("{atlas_path}/data/{.x}/{.x}.metadata.tsv"), data.table = FALSE) %>%
             group_by(alias) %>% 
             summarise_at(.vars = c("nCount_RNA", "nFeature_RNA"), .funs = mean)
         )

mouse_sample_summary_with_qc <- mouse_sample_summary %>%
    left_join(mouse_sample_qc_after_filt, by = c("Alias" = "alias")) %>% 
    select(Path, Original_sample = Sample, Alias, Protocol, Publication,
           Kit = `kit version`, Species, Age, Location,
           Starting_material = `Starting material`,
           N_reads = `Number of reads`,
           Reads_mapped_to_genome = `Reads mapped to genome`,
           Reads_mapped_to_transcriptome = `Reads mapped to transcriptome`,
           N_cells_estimated = `Estimated number of cells`,
           N_cells_after_filtering = `number of cells post-filtering`,
           Mean_UMIs_after_filtering = nCount_RNA,
           Mean_N_genes_after_filtering = nFeature_RNA,
           Min_cells_threshold = min_cells,
           Min_N_genes_threshold = min_nGene,
           Max_N_genes_threshold = max_nGene,
           Min_UMIs_threshold = min_nUMI,
           Max_UMIs_threshold = max_nUMI,
           Number_of_PCs = `number of PCs`,
           Clustering_resolution = `clustering resolution`) %>% 
    mutate(Publication = ifelse(is.na(Publication), "This study", Publication)) %>% 
    arrange(Publication, Age, Location)

mouse_sample_summary_with_qc %>%
    rr_write_tsv(glue("{out}/TABLE_mouse_sample_info.tsv"),
                 "Sample info and QC stats for mouse brain samples")
```



```
## ...writing description of TABLE_mouse_sample_info.tsv to public/output/05/TABLE_mouse_sample_info.desc
```



Per-cluster data:



```r
mouse_info_clusters <- data.table::fread(file.path(atlas_path, "data/metadata_extended/metadata_20210710_with_qc.tsv"), data.table = FALSE) %>% 
    left_join(mouse_sample_summary %>% select(Original_sample = Sample, Alias), by = "Alias") %>% 
    select(Sample, Original_sample, Alias, everything())

# sanity check that all the signatures are in the table
all(names(cluster_signatures$mm_sym) %in% mouse_info_clusters$Label)
```



```
## [1] TRUE
```



```r
# now, add the signatures for each cluster to the table
cluster_signatures_df <- data.frame("Signature_mm_symbols" = map_chr(cluster_signatures$mm_sym, ~ glue_collapse(.x, sep = ",")),
                                    "Signature_mm_ensembl" = map_chr(cluster_signatures$mm_ens, ~ glue_collapse(.x, sep = ",")),
                                    "Signature_hg_symbols" = map_chr(cluster_signatures$hg_sym, ~ glue_collapse(.x, sep = ",")),
                                    "Signature_hg_ensembl" = map_chr(cluster_signatures$hg_ens, ~ glue_collapse(.x, sep = ","))) %>% 
    tibble::rownames_to_column(var = "Label")

mouse_info_clusters_with_sigs <- mouse_info_clusters %>% 
    left_join(cluster_signatures_df, by = "Label")

# add cell class
# mouse_info_clusters_with_sigs_class <- mouse_info_clusters_with_sigs %>% 
#     summarize_cell_types("Level4_short") %>% 
#     rename(Cell_ontological_class = Type) %>% 
#     relocate(Cell_ontological_class, .before = "Level1_type")

# add dend order
mouse_info_clusters_with_sigs_class_dend <- mouse_info_clusters_with_sigs %>% 
    left_join(data.frame("Label" = dend_order_f) %>%
                  tibble::rowid_to_column(var = "Dendrogram_order_forebrain")) %>% 
    left_join(data.frame("Label" = dend_order_p) %>%
                  tibble::rowid_to_column(var = "Dendrogram_order_pons")) %>% 
    relocate(Dendrogram_order_forebrain, Dendrogram_order_pons, .after = Exclude)
```



```
## Joining, by = "Label"
## Joining, by = "Label"
```



```r
mouse_info_clusters_with_sigs_class_dend %>%
    rr_write_tsv(glue("{out}/TABLE_mouse_cluster_info.tsv"),
                 "Cluster info, QC stats, signatures, and dendrogram order for each mouse brain single-cell cluster")
```



```
## ...writing description of TABLE_mouse_cluster_info.tsv to public/output/05/TABLE_mouse_cluster_info.desc
```




# Supplementary figures

Generate annotations for the dendrograms:



```r
anno_dend_f <- mouse_info_clusters_with_sigs_class_dend %>% 
    select(Label, Dendrogram_order_forebrain, Timepoint) %>% 
    filter(!is.na(Dendrogram_order_forebrain)) %>% 
    arrange(Dendrogram_order_forebrain) %>% 
    mutate(Label = factor(Label, levels = .$Label))

anno_dend_f %>% 
    ggplot(aes(x = Label, y = 1)) +
    geom_tile(aes(fill = Timepoint), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = palette_timepoint) +
    theme_min() +
    theme(panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank()) +
    rotate_x()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/05/dend_f_anno-1.png)<!-- -->



```r
anno_dend_p <- mouse_info_clusters_with_sigs_class_dend %>% 
    select(Label, Dendrogram_order_pons, Timepoint) %>% 
    filter(!is.na(Dendrogram_order_pons)) %>% 
    arrange(Dendrogram_order_pons) %>% 
    mutate(Label = factor(Label, levels = .$Label))

anno_dend_p %>% 
    ggplot(aes(x = Label, y = 1)) +
    geom_tile(aes(fill = Timepoint), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = palette_timepoint) +
    theme_min() +
    theme(panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank()) +
    rotate_x()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/05/dend_p_anno-1.png)<!-- -->

Plot the # cells per sample:



```r
joint_mouse_meta <- readRDS(file.path(atlas_path,
                                      "2020-10_analyze_extended_dataset/integration_joint_mouse_extended/output/metadata.Rds"))

joint_mouse_meta %>% 
    separate(Sample, into = c("Location2", "Timepoint"), sep = " ") %>% 
    group_by(Timepoint, Location) %>% 
    count() %>% 
    ggplot(aes(x = Timepoint, y = n)) +
    geom_col(aes(fill = Timepoint), width = 0.5) +
    facet_wrap(~ Location) +
    geom_text(aes(label = n),  angle = 90, size = 4, hjust = -0.05) +
    scale_fill_manual(values = palette_timepoint) +
    rotate_x() +
    ylim(c(0, 13000))
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/05/barplot_n_cells-1.png)<!-- -->


<!-- END MATTER, insert reproducibility info -->




***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:



```
## 2022-09-12 15:22:22
```



The git repository and last commit:



```
## Local:    master /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public
## Remote:   master @ origin (git@github.com:fungenomics/HGG-oncohistones.git)
## Head:     [1a06382] 2022-09-08: Update comments, documentation, etc, based on lab feedback
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
##  date     2022-09-12                      
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  ! package        * version    date       lib source                           
##  P abind            1.4-5      2016-07-21 [?] CRAN (R 3.6.1)                   
##  P assertthat       0.2.1      2019-03-21 [?] CRAN (R 3.6.1)                   
##  P beeswarm         0.3.1      2021-03-07 [?] CRAN (R 3.6.1)                   
##  P bslib            0.2.5      2021-05-12 [?] CRAN (R 3.6.1)                   
##  P callr            3.7.0      2021-04-20 [?] CRAN (R 3.6.1)                   
##  P cellranger       1.1.0      2016-07-27 [?] CRAN (R 3.6.1)                   
##  P cli              2.5.0      2021-04-26 [?] CRAN (R 3.6.1)                   
##  P cluster          2.1.0      2019-06-19 [?] CRAN (R 3.6.1)                   
##  P codetools        0.2-16     2018-12-24 [?] CRAN (R 3.6.1)                   
##  P colorspace       2.0-1      2021-05-04 [?] CRAN (R 3.6.1)                   
##  P cowplot        * 1.1.1      2020-12-30 [?] CRAN (R 3.6.1)                   
##  P crayon           1.4.1      2021-02-08 [?] CRAN (R 3.6.1)                   
##  P data.table       1.14.0     2021-02-21 [?] CRAN (R 3.6.1)                   
##  P DBI              1.1.1      2021-01-15 [?] CRAN (R 3.6.1)                   
##  P deldir           0.2-10     2021-02-16 [?] CRAN (R 3.6.1)                   
##  P dendextend     * 1.14.0     2020-08-26 [?] CRAN (R 3.6.1)                   
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
##  P hms              1.0.0      2021-01-13 [?] CRAN (R 3.6.1)                   
##  P htmltools        0.5.1.1    2021-01-22 [?] CRAN (R 3.6.1)                   
##  P htmlwidgets      1.5.3      2020-12-10 [?] CRAN (R 3.6.1)                   
##  P httpuv           1.6.1      2021-05-07 [?] CRAN (R 3.6.1)                   
##  P httr             1.4.2      2020-07-20 [?] CRAN (R 3.6.1)                   
##  P ica              1.0-2      2018-05-24 [?] CRAN (R 3.6.1)                   
##  P igraph           1.2.6      2020-10-06 [?] CRAN (R 3.6.1)                   
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
##  P parallelly       1.25.0     2021-04-30 [?] CRAN (R 3.6.1)                   
##  P patchwork        1.1.1      2020-12-17 [?] CRAN (R 3.6.1)                   
##  P pbapply          1.4-3      2020-08-18 [?] CRAN (R 3.6.1)                   
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
##  P promises         1.1.0      2019-10-04 [?] CRAN (R 3.6.1)                   
##  P ps               1.6.0      2021-02-28 [?] CRAN (R 3.6.1)                   
##  P purrr          * 0.3.4      2020-04-17 [?] CRAN (R 3.6.1)                   
##  P pvclust        * 2.2-0      2019-11-19 [?] CRAN (R 3.6.1)                   
##  P R6               2.5.0      2020-10-28 [?] CRAN (R 3.6.1)                   
##  P RANN             2.6.1      2019-01-08 [?] CRAN (R 3.6.1)                   
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
##  P rsvd             1.0.3      2020-02-17 [?] CRAN (R 3.6.1)                   
##  P Rtsne            0.15       2018-11-10 [?] CRAN (R 3.6.1)                   
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
##  P tibble           3.1.1      2021-04-18 [?] CRAN (R 3.6.1)                   
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
##  P xtable           1.8-4      2019-04-21 [?] CRAN (R 3.6.1)                   
##  P yaml             2.2.1      2020-02-01 [?] CRAN (R 3.6.1)                   
##  P zoo              1.8-9      2021-03-09 [?] CRAN (R 3.6.1)                   
## 
## [1] /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/renv/library/R-3.6/x86_64-pc-linux-gnu
## [2] /tmp/RtmppM1df5/renv-system-library
## 
##  P ── Loaded and on-disk path mismatch.
```



</details>

The resources requested when this document was last rendered:



```
## #SBATCH --time=01:00:00
## #SBATCH --cpus-per-task=1
## #SBATCH --mem=20G
```




***



<!-- END OF END MATTER -->
