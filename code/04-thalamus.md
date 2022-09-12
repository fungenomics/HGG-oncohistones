---
title: "04 - Analysis of thalamic H3K27M tumors"
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
## Document index: 04
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
## public/output/04
```

```
## public/figures/04
```



Setting a random seed:



```r
set.seed(100)
```



***



<!-- END OF FRONT MATTER -->


# Overview

In this document, we investigate the activation of diencephalon patterning genes
in thalamic gliomas as shown in Figure 3.

# Libraries



```r
# Load libraries here
library(here)
library(tidyr)
library(dplyr)

library(GenomicRanges)
library(BentoBox)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)

library(ggrepel)
library(readr)
library(glue)
library(tibble)
library(ggplot2)
library(purrr)
library(cowplot)

source(here("include/style.R"))
source(here("code/functions/BentoBox_helpers.R"))
source(here("code/functions/RNAseq.R"))
source(here("code/functions/scRNAseq.R"))
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




# Activation of prosomere markers

To map thalamic tumors to the prosomeres of the diencephalon, we'll assess the activation
of various prosomere 2/3 markers in the tumors based on RNAseq, H3K27ac and H3K27me3:



```r
# set genes to examine
prosomere_markers <- list(
    GRanges("chr3",  147093346:147150046, name = "ZIC1/4"),
    GRanges("chr2",  63273543:63289360,   name = "OTX1"),
    GRanges("chr2",  237071880:237079029, name = "GBX2"),
    GRanges("chr16", 54315217:54322699,   name = "IRX3"),
    GRanges("chr7",  2544962:2574432,     name = "LFNG"),
    GRanges("chr10", 119295203:119315811, name = "EMX2"),
    GRanges("chr2",  45163986:45175965,   name = "SIX3"),
    GRanges("chr7",  121940057:121946104, name = "FEZF1"))

# make params
params_prosomere <- map(prosomere_markers, ~ bb_params(chrom = as.character(seqnames(.x)),
                                                       chromstart = start(.x),
                                                       chromend = end(.x),
                                                       assembly = "hg19"))
names(params_prosomere) <- map(prosomere_markers, ~ .x$name)
```






```r
# load in track configuration
prosomere_config <- prep_track_input(here("data/BentoBox_config/prosomere_thalamus.config.tsv"))
```



Make a summary figure of the tracks at these regions. In pseudocode:

```

for each region:
-- for each sample:
----- plot signal & annotation
----- plot genome label

```

Generate figure:



```r
x_positions <- seq(0.5, 15, by = 1.75)
y_positions <- seq(0.5, by = 0.4, length.out = nrow(prosomere_config))

bb_pageCreate(width = 15, height = 7, default.units = "inches")

for (i in seq_along(params_prosomere)) {
    
    x_i <- x_positions[i]
    params_i <- params_prosomere[[i]]

    # iterating through samples,
    # plot H3K27ac and H3K27me3
    pwalk(list(prosomere_config$bw, prosomere_config$ID_paper, prosomere_config$Data, prosomere_config$Ymax, y_positions),
          ~ bb_placeSignalAndLabel(data = ..1,
                                   annotation = ..2,
                                   color = palette_tracks[..3],
                                   range = c(0, ..4),
                                   y = ..5, x = x_i,
                                   params = params_i,
                                   width = 1.5,
                                   height = 0.35,
                                   fontsize = 4))
    
    # place 0.03in below last plot
    bb_plotGenomeLabel(params = params_i, scale = "Kb", x = x_i, y = "0.03b", length = 1.5, fontsize = 8)
    
    bb_plotGenes(params = params_i, x = x_i, y = "0.05b", height = 1.5, width = 1.5,
                 strandcolors = c("navy", "black"), fontcolors = c("navy", "black"),
                 stroke = 0.05,
                 just = c("left", "top"), default.units = "inches")
    
}

bb_pageGuideHide()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/04/prosomere_tracks-1.png)<!-- -->




<!-- END MATTER, insert reproducibility info -->




***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:



```
## 2022-09-12 15:22:30
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
##  ! package                           * version    date       lib
##  P AnnotationDbi                     * 1.48.0     2019-10-29 [?]
##  P askpass                             1.1        2019-01-13 [?]
##  P assertthat                          0.2.1      2019-03-21 [?]
##  P BentoBox                          * 0.1.0      2022-06-14 [?]
##  P Biobase                           * 2.46.0     2019-10-29 [?]
##  P BiocFileCache                       1.10.2     2019-11-08 [?]
##  P BiocGenerics                      * 0.32.0     2019-10-29 [?]
##  P BiocManager                         1.30.15    2021-05-11 [?]
##  P BiocParallel                        1.20.1     2019-12-21 [?]
##  P biomaRt                             2.42.1     2020-03-26 [?]
##  P Biostrings                        * 2.54.0     2019-10-29 [?]
##  P bit                                 4.0.4      2020-08-04 [?]
##  P bit64                               4.0.5      2020-08-30 [?]
##  P bitops                              1.0-7      2021-04-24 [?]
##  P blob                                1.2.1      2020-01-20 [?]
##  P BSgenome                          * 1.54.0     2019-10-29 [?]
##  P BSgenome.Hsapiens.UCSC.hg19       * 1.4.0      2021-12-06 [?]
##  P bslib                               0.2.5      2021-05-12 [?]
##  P callr                               3.7.0      2021-04-20 [?]
##  P cli                                 2.5.0      2021-04-26 [?]
##  P codetools                           0.2-16     2018-12-24 [?]
##  P colorspace                          2.0-1      2021-05-04 [?]
##  P cowplot                           * 1.1.1      2020-12-30 [?]
##  P crayon                              1.4.1      2021-02-08 [?]
##  P curl                                4.3.1      2021-04-30 [?]
##  P data.table                        * 1.14.0     2021-02-21 [?]
##  P DBI                                 1.1.1      2021-01-15 [?]
##  P dbplyr                              2.1.1      2021-04-06 [?]
##  P DelayedArray                        0.12.3     2020-04-09 [?]
##  P desc                                1.2.0      2018-05-01 [?]
##  P devtools                            2.3.0      2020-04-10 [?]
##  P digest                              0.6.27     2020-10-24 [?]
##  P dplyr                             * 1.0.6      2021-05-05 [?]
##  P ellipsis                            0.3.2      2021-04-29 [?]
##  P evaluate                            0.14       2019-05-28 [?]
##  P fansi                               0.4.2      2021-01-15 [?]
##  P fs                                  1.5.0      2020-07-31 [?]
##  P generics                            0.1.0      2020-10-31 [?]
##  P GenomeInfoDb                      * 1.22.1     2020-03-27 [?]
##  P GenomeInfoDbData                    1.2.2      2021-12-06 [?]
##  P GenomicAlignments                   1.22.1     2019-11-12 [?]
##  P GenomicFeatures                   * 1.38.2     2020-02-15 [?]
##  P GenomicRanges                     * 1.38.0     2019-10-29 [?]
##  P ggplot2                           * 3.3.3      2020-12-30 [?]
##  P ggplotify                         * 0.0.7      2021-05-11 [?]
##  P ggrepel                           * 0.9.1      2021-01-15 [?]
##  P git2r                               0.27.1     2020-05-03 [?]
##  P glue                              * 1.4.2      2020-08-27 [?]
##  P gmp                               * 0.6-2      2021-01-07 [?]
##  P gridExtra                           2.3        2017-09-09 [?]
##  P gridGraphics                        0.5-1      2020-12-13 [?]
##  P gtable                              0.3.0      2019-03-25 [?]
##  P here                              * 0.1        2017-05-28 [?]
##  P highr                               0.9        2021-04-16 [?]
##  P hms                                 1.0.0      2021-01-13 [?]
##  P htmltools                           0.5.1.1    2021-01-22 [?]
##  P httr                                1.4.2      2020-07-20 [?]
##  P IRanges                           * 2.20.2     2020-01-13 [?]
##  P jquerylib                           0.1.4      2021-04-26 [?]
##  P jsonlite                            1.7.2      2020-12-09 [?]
##  P knitr                               1.33       2021-04-24 [?]
##  P lattice                             0.20-44    2021-05-02 [?]
##  P lifecycle                           1.0.0      2021-02-15 [?]
##  P magrittr                          * 2.0.1      2020-11-17 [?]
##  P Matrix                              1.2-18     2019-11-27 [?]
##  P matrixStats                         0.58.0     2021-01-29 [?]
##  P memoise                             1.1.0      2017-04-21 [?]
##  P munsell                             0.5.0      2018-06-12 [?]
##  P openssl                             1.4.4      2021-04-30 [?]
##  P org.Hs.eg.db                      * 3.10.0     2021-12-06 [?]
##  P pillar                              1.6.0      2021-04-13 [?]
##  P pkgbuild                            1.0.8      2020-05-07 [?]
##  P pkgconfig                           2.0.3      2019-09-22 [?]
##  P pkgload                             1.0.2      2018-10-29 [?]
##  P prettyunits                         1.1.1      2020-01-24 [?]
##  P processx                            3.5.2      2021-04-30 [?]
##  P progress                            1.2.2      2019-05-16 [?]
##  P ps                                  1.6.0      2021-02-28 [?]
##  P purrr                             * 0.3.4      2020-04-17 [?]
##  P R6                                  2.5.0      2020-10-28 [?]
##  P rappdirs                            0.3.3      2021-01-31 [?]
##  P RColorBrewer                      * 1.1-2      2014-12-07 [?]
##  P Rcpp                                1.0.6      2021-01-15 [?]
##  P RCurl                               1.98-1.3   2021-03-16 [?]
##  P readr                             * 1.4.0      2020-10-05 [?]
##  P remotes                             2.1.1      2020-02-15 [?]
##    renv                                0.14.0     2021-07-21 [1]
##  P rlang                               0.4.11     2021-04-30 [?]
##  P rmarkdown                           2.8        2021-05-07 [?]
##  P Rmpfr                             * 0.8-4      2021-04-11 [?]
##  P rprojroot                           2.0.2      2020-11-15 [?]
##  P Rsamtools                           2.2.3      2020-02-23 [?]
##  P RSQLite                             2.2.1      2020-09-30 [?]
##  P rstudioapi                          0.13       2020-11-12 [?]
##  P rtracklayer                       * 1.46.0     2019-10-29 [?]
##  P rvcheck                             0.1.8      2020-03-01 [?]
##  P S4Vectors                         * 0.24.4     2020-04-09 [?]
##  P sass                                0.4.0      2021-05-12 [?]
##  P scales                              1.1.1      2020-05-11 [?]
##  P sessioninfo                         1.1.1      2018-11-05 [?]
##  P stringi                             1.6.1      2021-05-10 [?]
##  P stringr                             1.4.0      2019-02-10 [?]
##  P SummarizedExperiment                1.16.1     2019-12-19 [?]
##  P testrmd                             0.0.1.9000 2021-12-06 [?]
##  P testthat                            2.3.2      2020-03-02 [?]
##  P tibble                            * 3.1.1      2021-04-18 [?]
##  P tidyr                             * 1.1.3      2021-03-03 [?]
##  P tidyselect                          1.1.1      2021-04-30 [?]
##  P TxDb.Hsapiens.UCSC.hg19.knownGene * 3.2.2      2021-12-06 [?]
##  P usethis                             1.6.1      2020-04-29 [?]
##  P utf8                                1.2.1      2021-03-12 [?]
##  P vctrs                               0.3.8      2021-04-29 [?]
##  P viridis                           * 0.5.1      2018-03-29 [?]
##  P viridisLite                       * 0.4.0      2021-04-13 [?]
##  P withr                               2.4.2      2021-04-18 [?]
##  P xfun                                0.22       2021-03-11 [?]
##  P XML                                 3.99-0.3   2020-01-20 [?]
##  P XVector                           * 0.26.0     2019-10-29 [?]
##  P yaml                                2.2.1      2020-02-01 [?]
##  P zlibbioc                            1.32.0     2019-10-29 [?]
##  source                                
##  Bioconductor                          
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  Github (PhanstielLab/BentoBox@72701d1)
##  Bioconductor                          
##  Bioconductor                          
##  Bioconductor                          
##  CRAN (R 3.6.1)                        
##  Bioconductor                          
##  Bioconductor                          
##  Bioconductor                          
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  Bioconductor                          
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
##  Bioconductor                          
##  Bioconductor                          
##  Bioconductor                          
##  Bioconductor                          
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
##  Bioconductor                          
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  Bioconductor                          
##  CRAN (R 3.6.1)                        
##  Bioconductor                          
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  Bioconductor                          
##  Github (rmflight/testrmd@0735c20)     
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
##  Bioconductor                          
##  CRAN (R 3.6.1)                        
##  Bioconductor                          
## 
## [1] /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/renv/library/R-3.6/x86_64-pc-linux-gnu
## [2] /tmp/RtmpBnm3KZ/renv-system-library
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
