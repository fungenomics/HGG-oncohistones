---
title: "06 - Effect of ACVR1"
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
## Document index: 06
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
## public/output/06
```

```
## public/figures/06
```



Setting a random seed:



```r
set.seed(100)
```



***



<!-- END OF FRONT MATTER -->


# Overview

This document contains analyses to understand the dependency of ACVR1 in K27M
and especially H3.1K27M gliomas, and to interrogate the changes that result from
knocking-out ACVR1 using CRISPR in ACVR1 mutant HGG cell lines, as shown in Figure 6.

# Libraries



```r
# Load libraries here
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
library(stringr)
library(cowplot)
library(fgsea)
library(icytobox)

source(here("include/style.R")) # contains palettes and custom style elements
source(here("code/functions/RNAseq.R"))
source(here("code/functions/ssGSEA.R"))
ggplot2::theme_set(theme_min())
```





# ID genes ddPCR

ddPCR assay of ID gene expression.



```r
# load data
ddpcr <- suppressMessages(read_tsv(here("data/experimental/2022-03-10-ddPCR_ID_genes_@mhulswit.txt")))

# barplot
ddpcr_tidy <- ddpcr %>% 
  # remove xenografts
  filter(!grepl("[Xx]eno", Sample)) %>%
  gather(Stat, Value, matches("ID")) %>% 
  separate(Stat, into = c("Gene", "Stat"), sep = "_") %>% 
  spread(Stat, Value) %>% 
  separate(Sample, into = c("Cell_line", "CRISPR"), sep = "_") %>% 
  mutate(Condition = factor(Condition, levels = c("ACVR1 mutant", "ACVR1 KO")))

ddpcr_tidy %>% 
  rr_ggplot(aes(x = Condition, y = fold), plot_num = 1) +
  geom_bar(aes(fill = Condition), stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = fold - sd, ymax = fold + sd), width = 0.25) +
  scale_fill_manual(values = palette_acvr1) +
  facet_grid(Cell_line ~ Gene, scales = "free_x", space = "free", drop = TRUE) + 
  ylim(c(0, 1.2)) +
  rotate_x()
```



```
## ...writing source data of ggplot to public/figures/06/ddpcr_id_genes-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/06/ddpcr_id_genes-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/06/ddpcr_id_genes...*]~</span>


# Clonogenesis

Data for clone formation assay for ACVR1 mutant and KO cell lines.



```r
# load data
clono <- suppressMessages(read_tsv(here("data/experimental/clonogenic_assay_DIPG36_@am.txt")))

# define standard error to quantify uncertainty around the mean
se <- function(x) sqrt(var(x)/length(x))

# tidy & plot
clono_tidy <- clono %>% 
    set_colnames(c("ACVR1 mutant", "ACVR1 KO")) %>% 
    gather(Condition, N_clones) %>% 
    mutate(Condition = factor(Condition, levels = names(palette_acvr1))) %>% 
    filter(!is.na(N_clones))

clono_tidy %>% 
    group_by(Condition) %>% 
    summarize(N_clones_mean = mean(N_clones),
              N_clones_se   = se(N_clones)) %>% 
    rr_ggplot(aes(x = Condition, y = N_clones_mean), plot_num = 1) +
    geom_bar(aes(fill = Condition), stat = "identity", width = 0.5) +
    geom_errorbar(aes(ymin = N_clones_mean - N_clones_se, ymax = N_clones_mean + N_clones_se), width = 0.25) +
    geom_jitter(data = clono_tidy, aes(x = Condition, y = N_clones), width = 0.2, size = 1) +
    scale_fill_manual(values = palette_acvr1) +
    ylim(c(0, 150))
```



```
## ...writing source data of ggplot to public/figures/06/clonogenesis-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/06/clonogenesis-1.png)<!-- -->

```r
t.test(clono$Parental, clono$`ACVR1 KO`, alternative = "two.sided", var.equal = TRUE)
```



```
## 
## 	Two Sample t-test
## 
## data:  clono$Parental and clono$`ACVR1 KO`
## t = 2.7602, df = 7, p-value = 0.02809
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##   6.855717 88.810950
## sample estimates:
## mean of x mean of y 
## 116.50000  68.66667
```

<br><span style="color:#0d00ff">~[figure @ *public/figures/06/clonogenesis...*]~</span>


# Doubling times

Doubling times of H3.1 and H3.3K27M cell lines.



```r
# load data
doubling_H3.1_vs_H3.3 <- suppressMessages(read_tsv(here("data/experimental/doubling_time_H3.3_vs_H3.1_@am.txt")))

# tidy & plot
doubling_H3.1_vs_H3.3_tidy <- doubling_H3.1_vs_H3.3 %>% 
    set_colnames(c("DIPGXIII", "HSJ-019", "DIPG36", "DIPGIV")) %>% 
    gather(Cell_line, Doubling_time) %>% 
    filter(!is.na(Doubling_time))

# use standard deviation to assesss variability of the measurements
doubling_H3.1_vs_H3.3_tidy %>% 
    group_by(Cell_line) %>% 
    summarize(Doubling_mean = mean(Doubling_time),
              Doubling_sd   = sd(Doubling_time)) %>% 
    rr_ggplot(aes(x = Cell_line, y = Doubling_mean), plot_num = 1) +
    geom_bar(aes(fill = Cell_line), stat = "identity", width = 0.5) +
    geom_errorbar(aes(ymin = Doubling_mean - Doubling_sd, ymax = Doubling_mean + Doubling_sd), width = 0.25) +
    geom_jitter(data = doubling_H3.1_vs_H3.3_tidy, aes(x = Cell_line, y = Doubling_time), width = 0.1, size = 1) +
    scale_fill_manual(values = c("DIPG36" = "orange", "DIPGIV" = "orange", "DIPGXIII" = "red3", "HSJ-019" = "red3")) + 
    no_legend()
```



```
## ...writing source data of ggplot to public/figures/06/doubling_H3.1_vs_H3.3-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/06/doubling_H3.1_vs_H3.3-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/06/doubling_H3.1_vs_H3.3...*]~</span>



```r
# load data
doubling_ACVR1 <- suppressMessages(read_tsv(here("data/experimental/doubling_time_DIPGIV_ACVR1mut_vs_KO_@am.txt")))

# tidy & plot
doubling_ACVR1_tidy <- doubling_ACVR1 %>% 
    set_colnames(c("ACVR1 mutant", "ACVR1 KO")) %>% 
    gather(Condition, Doubling_time) %>% 
    mutate(Condition = factor(Condition, levels = names(palette_acvr1))) %>% 
    filter(!is.na(Condition))

doubling_ACVR1_tidy %>% 
    group_by(Condition) %>% 
    summarize(Doubling_mean = mean(Doubling_time),
              Doubling_sd   = sd(Doubling_time)) %>% 
    rr_ggplot(aes(x = Condition, y = Doubling_mean), plot_num = 1) +
    geom_bar(aes(fill = Condition), stat = "identity", width = 0.5) +
    geom_errorbar(aes(ymin = Doubling_mean - Doubling_sd, ymax = Doubling_mean + Doubling_sd), width = 0.25) +
    geom_jitter(data = doubling_ACVR1_tidy, aes(x = Condition, y = Doubling_time), width = 0.1, size = 1) +
    scale_fill_manual(values = palette_acvr1) + 
    no_legend() +
    ylim(c(0, 50))
```



```
## ...writing source data of ggplot to public/figures/06/doubling_ACVR1-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/06/doubling_ACVR1-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/06/doubling_ACVR1...*]~</span>


<!-- END MATTER, insert reproducibility info -->




***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:



```
## 2022-09-12 15:23:43
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
##  ! package        * version    date       lib
##  P abind            1.4-5      2016-07-21 [?]
##  P assertthat       0.2.1      2019-03-21 [?]
##  P BiocParallel     1.20.1     2019-12-21 [?]
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
##  P farver           2.1.0      2021-02-28 [?]
##  P fastmap          1.1.0      2021-01-25 [?]
##  P fastmatch        1.1-0      2017-01-28 [?]
##  P fgsea          * 1.12.0     2019-10-29 [?]
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
##  P highr            0.9        2021-04-16 [?]
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
##  P labeling         0.4.2      2020-10-20 [?]
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
##  P pbapply        * 1.4-3      2020-08-18 [?]
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
##  P Rcpp           * 1.0.6      2021-01-15 [?]
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
##  Bioconductor                         
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  Bioconductor                         
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
##  CRAN (R 3.6.1)                       
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
## [2] /tmp/RtmpYBfkw4/renv-system-library
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
