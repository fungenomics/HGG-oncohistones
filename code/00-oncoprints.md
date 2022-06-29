---
title: "00 - Generate oncoprints"
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
## Document index: 00
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
## public/output/00
```

```
## public/figures/00
```



Setting a random seed:



```r
set.seed(100)
```



***



<!-- END OF FRONT MATTER -->


# Overview

This document plots summaries of the samples included in this analysis as oncoprints as shown in Figure 1.

# Libraries



```r
# Load libraries here
library(here)
library(tidyr)
library(dplyr)
library(readr)
library(glue)
library(tibble)
library(ggplot2)
library(cowplot)

source(here("include/style.R")) # contains palettes & plotting utils
ggplot2::theme_set(theme_min())
```



# Oncoprint

We'll generate an oncoprint for unique patient tumor samples, and patient-derived
cell lines.

## Tumors



```r
oncoprint_input_tumors <- read_tsv(here("data/metadata/oncoprint_input_tumors.tsv")) %>% 
  mutate(Molecular = factor(Molecular, levels = names(palette_molecular)),
           Location = factor(Location, levels = names(palette_location))) %>%
  arrange(desc(Type), Molecular, Location, GrowthFactorReceptor,
            desc(scRNAseq), desc(scATACseq), desc(RNAseq), desc(H3K27ac), desc(H3K27me3), desc(H3K27me2)) %>%
    mutate(ID_patient = factor(ID_patient, levels = .$ID_patient))

p0 <- oncoprint_input_tumors %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = Material), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Tumor" = "azure2", "Cell line" = "gray80", "Both" = "gray50"),
                      guide = guide_legend(ncol = 2)) +
    theme_min() + # the theme_row() theme is defined in include/style.R
    theme_row()

p1 <- oncoprint_input_tumors %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = Molecular), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = palette_molecular, guide = guide_legend(ncol = 2)) +
    theme_min() +
    theme_row()

p2 <- oncoprint_input_tumors %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = GrowthFactorReceptor), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = palette_gfr, na.value = "gray90", guide = guide_legend(ncol = 2)) +
    theme_min() +
    theme_row()

p3 <- oncoprint_input_tumors %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = Location), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = palette_location, na.value = "gray90", guide = guide_legend(ncol = 2)) +
    theme_min() +
    theme_row()

p4 <- oncoprint_input_tumors %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = Age), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"), na.value = "gray90") +
    theme_min() +
    theme_row()

p5 <- oncoprint_input_tumors %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = scRNAseq), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Y (new)" = "#2a9135", "Y" = "#6fc978", "-" = "gray90")) +
    theme_min() +
    theme_row()

p6 <- oncoprint_input_tumors %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = scATACseq), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Y (new)" = "#0ad1b7", "-" = "gray90")) +
    theme_min() +
    theme_row()

p7 <- oncoprint_input_tumors %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = RNAseq), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Y (new)" = "#2064f7", "Y" = "#6c88c4", "-" = "gray90")) +
    theme_min() +
    theme_row()

p8 <- oncoprint_input_tumors %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = H3K27ac), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Y (new)" = "#f72076", "Y" = "#c95583", "-" = "gray90")) +
    theme_min() +
    theme_row()

p9 <- oncoprint_input_tumors %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = H3K27me3), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Y (new)" = "#650ac7", "Y" = "#9363c7", "-" = "gray90")) +
    theme_min() +
    theme_row()

p10 <- oncoprint_input_tumors %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = H3K27me2), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Y (new)" = "#e04902", "Y" = "#d18460", "-" = "gray90")) +
    theme(panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank()) +
    rotate_x()

cowplot::plot_grid(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
                   align = "v", axis = "r", ncol = 1,
                   rel_heights = c(rep(0.08, 10), 0.2))
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/00/oncoprint_tumors-1.png)<!-- -->

## Cell lines



```r
oncoprint_input_cl <- read_tsv(here("data/metadata/oncoprint_input_cl.tsv")) %>% 
  mutate(Molecular = factor(Molecular, levels = names(palette_molecular)),
           Location = factor(Location, levels = names(palette_location))) %>%
  arrange(desc(Type), Molecular, Location, GrowthFactorReceptor,
            desc(scRNAseq), desc(scATACseq), desc(RNAseq), desc(H3K27ac), desc(H3K27me3), desc(H3K27me2)) %>%
    mutate(ID_patient = factor(ID_patient, levels = .$ID_patient))

p0 <- oncoprint_input_cl %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = Material), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Tumor" = "azure2", "Cell line" = "gray80", "Both" = "gray50"),
                      guide = guide_legend(ncol = 2)) +
    theme_min() +
    theme_row()

p1 <- oncoprint_input_cl %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = Molecular), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = palette_molecular, guide = guide_legend(ncol = 2)) +
    theme_min() +
    theme_row()

p2 <- oncoprint_input_cl %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = GrowthFactorReceptor), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = palette_gfr, na.value = "gray90", guide = guide_legend(ncol = 2)) +
    theme_min() +
    theme_row()

p3 <- oncoprint_input_cl %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = Location), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = palette_location, na.value = "gray90", guide = guide_legend(ncol = 2)) +
    theme_min() +
    theme_row()

p4 <- oncoprint_input_cl %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = Age), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"), na.value = "gray90", limits = c(0, 20)) +
    theme_min() +
    theme_row()

p5 <- oncoprint_input_cl %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = scRNAseq), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Y (new)" = "#2a9135", "Y" = "#6fc978", "-" = "gray90")) +
    theme_min() +
    theme_row()

p6 <- oncoprint_input_cl %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = scATACseq), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Y (new)" = "#0ad1b7", "-" = "gray90")) +
    theme_min() +
    theme_row()

p7 <- oncoprint_input_cl %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = RNAseq), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Y (new)" = "#2064f7", "Y" = "#6c88c4", "-" = "gray90")) +
    theme_min() +
    theme_row()

p8 <- oncoprint_input_cl %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = H3K27ac), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Y (new)" = "#f72076", "Y" = "#c95583", "-" = "gray90")) +
    theme_min() +
    theme_row()

p9 <- oncoprint_input_cl %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = H3K27me3), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Y (new)" = "#650ac7", "Y" = "#9363c7", "-" = "gray90")) +
    theme_min() +
    theme_row()

p10 <- oncoprint_input_cl %>%
    ggplot(aes(x = ID_patient, y = 1)) +
    geom_tile(aes(fill = H3K27me2), colour = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("Y (new)" = "#e04902", "Y" = "#d18460", "-" = "gray90")) +
    theme(panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank()) +
    rotate_x()

cowplot::plot_grid(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
                   align = "v", axis = "r", ncol = 1,
                   rel_heights = c(rep(0.08, 10), 0.2))
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/00/oncoprint_cl-1.png)<!-- -->



<!-- END MATTER, insert reproducibility info -->




***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:



```
## 2022-06-29 09:45:56
```



The git repository and last commit:



```
## Local:    master /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public
## Remote:   master @ origin (git@github.com:fungenomics/HGG-oncohistones.git)
## Head:     [056f679] 2022-06-14: Add R-4 renv lockfile
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
##  ! package      * version    date       lib source                           
##  P assertthat     0.2.1      2019-03-21 [?] CRAN (R 3.6.1)                   
##  P bslib          0.2.5      2021-05-12 [?] CRAN (R 3.6.1)                   
##  P callr          3.7.0      2021-04-20 [?] CRAN (R 3.6.1)                   
##  P cli            2.5.0      2021-04-26 [?] CRAN (R 3.6.1)                   
##  P colorspace     2.0-1      2021-05-04 [?] CRAN (R 3.6.1)                   
##  P cowplot      * 1.1.1      2020-12-30 [?] CRAN (R 3.6.1)                   
##  P crayon         1.4.1      2021-02-08 [?] CRAN (R 3.6.1)                   
##  P DBI            1.1.1      2021-01-15 [?] CRAN (R 3.6.1)                   
##  P desc           1.2.0      2018-05-01 [?] CRAN (R 3.6.1)                   
##  P devtools       2.3.0      2020-04-10 [?] CRAN (R 3.6.1)                   
##  P digest         0.6.27     2020-10-24 [?] CRAN (R 3.6.1)                   
##  P dplyr        * 1.0.6      2021-05-05 [?] CRAN (R 3.6.1)                   
##  P ellipsis       0.3.2      2021-04-29 [?] CRAN (R 3.6.1)                   
##  P evaluate       0.14       2019-05-28 [?] CRAN (R 3.6.1)                   
##  P fansi          0.4.2      2021-01-15 [?] CRAN (R 3.6.1)                   
##  P farver         2.1.0      2021-02-28 [?] CRAN (R 3.6.1)                   
##  P fs             1.5.0      2020-07-31 [?] CRAN (R 3.6.1)                   
##  P generics       0.1.0      2020-10-31 [?] CRAN (R 3.6.1)                   
##  P ggplot2      * 3.3.3      2020-12-30 [?] CRAN (R 3.6.1)                   
##  P git2r          0.27.1     2020-05-03 [?] CRAN (R 3.6.1)                   
##  P glue         * 1.4.2      2020-08-27 [?] CRAN (R 3.6.1)                   
##  P gridExtra      2.3        2017-09-09 [?] CRAN (R 3.6.1)                   
##  P gtable         0.3.0      2019-03-25 [?] CRAN (R 3.6.1)                   
##  P here         * 0.1        2017-05-28 [?] CRAN (R 3.6.1)                   
##  P highr          0.9        2021-04-16 [?] CRAN (R 3.6.1)                   
##  P hms            1.0.0      2021-01-13 [?] CRAN (R 3.6.1)                   
##  P htmltools      0.5.1.1    2021-01-22 [?] CRAN (R 3.6.1)                   
##  P jquerylib      0.1.4      2021-04-26 [?] CRAN (R 3.6.1)                   
##  P jsonlite       1.7.2      2020-12-09 [?] CRAN (R 3.6.1)                   
##  P knitr          1.33       2021-04-24 [?] CRAN (R 3.6.1)                   
##  P labeling       0.4.2      2020-10-20 [?] CRAN (R 3.6.1)                   
##  P lifecycle      1.0.0      2021-02-15 [?] CRAN (R 3.6.1)                   
##  P magrittr     * 2.0.1      2020-11-17 [?] CRAN (R 3.6.1)                   
##  P memoise        1.1.0      2017-04-21 [?] CRAN (R 3.6.1)                   
##  P munsell        0.5.0      2018-06-12 [?] CRAN (R 3.6.1)                   
##  P pillar         1.6.0      2021-04-13 [?] CRAN (R 3.6.1)                   
##  P pkgbuild       1.0.8      2020-05-07 [?] CRAN (R 3.6.1)                   
##  P pkgconfig      2.0.3      2019-09-22 [?] CRAN (R 3.6.1)                   
##  P pkgload        1.0.2      2018-10-29 [?] CRAN (R 3.6.1)                   
##  P prettyunits    1.1.1      2020-01-24 [?] CRAN (R 3.6.1)                   
##  P processx       3.5.2      2021-04-30 [?] CRAN (R 3.6.1)                   
##  P ps             1.6.0      2021-02-28 [?] CRAN (R 3.6.1)                   
##  P purrr          0.3.4      2020-04-17 [?] CRAN (R 3.6.1)                   
##  P R6             2.5.0      2020-10-28 [?] CRAN (R 3.6.1)                   
##  P RColorBrewer * 1.1-2      2014-12-07 [?] CRAN (R 3.6.1)                   
##  P readr        * 1.4.0      2020-10-05 [?] CRAN (R 3.6.1)                   
##  P remotes        2.1.1      2020-02-15 [?] CRAN (R 3.6.1)                   
##    renv           0.14.0     2021-07-21 [1] CRAN (R 3.6.1)                   
##  P rlang          0.4.11     2021-04-30 [?] CRAN (R 3.6.1)                   
##  P rmarkdown      2.8        2021-05-07 [?] CRAN (R 3.6.1)                   
##  P rprojroot      2.0.2      2020-11-15 [?] CRAN (R 3.6.1)                   
##  P rstudioapi     0.13       2020-11-12 [?] CRAN (R 3.6.1)                   
##  P sass           0.4.0      2021-05-12 [?] CRAN (R 3.6.1)                   
##  P scales         1.1.1      2020-05-11 [?] CRAN (R 3.6.1)                   
##  P sessioninfo    1.1.1      2018-11-05 [?] CRAN (R 3.6.1)                   
##  P stringi        1.6.1      2021-05-10 [?] CRAN (R 3.6.1)                   
##  P stringr        1.4.0      2019-02-10 [?] CRAN (R 3.6.1)                   
##  P testrmd        0.0.1.9000 2021-12-06 [?] Github (rmflight/testrmd@0735c20)
##  P testthat       2.3.2      2020-03-02 [?] CRAN (R 3.6.1)                   
##  P tibble       * 3.1.1      2021-04-18 [?] CRAN (R 3.6.1)                   
##  P tidyr        * 1.1.3      2021-03-03 [?] CRAN (R 3.6.1)                   
##  P tidyselect     1.1.1      2021-04-30 [?] CRAN (R 3.6.1)                   
##  P usethis        1.6.1      2020-04-29 [?] CRAN (R 3.6.1)                   
##  P utf8           1.2.1      2021-03-12 [?] CRAN (R 3.6.1)                   
##  P vctrs          0.3.8      2021-04-29 [?] CRAN (R 3.6.1)                   
##  P viridis      * 0.5.1      2018-03-29 [?] CRAN (R 3.6.1)                   
##  P viridisLite  * 0.4.0      2021-04-13 [?] CRAN (R 3.6.1)                   
##  P withr          2.4.2      2021-04-18 [?] CRAN (R 3.6.1)                   
##  P xfun           0.22       2021-03-11 [?] CRAN (R 3.6.1)                   
##  P yaml           2.2.1      2020-02-01 [?] CRAN (R 3.6.1)                   
## 
## [1] /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/renv/library/R-3.6/x86_64-pc-linux-gnu
## [2] /tmp/RtmpZIxQ4e/renv-system-library
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
