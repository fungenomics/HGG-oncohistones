---
title: "07 - Analysis in K27M and K27M-KO cell lines"
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
## Document index: 07A
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
## public/output/07A
```

```
## public/figures/07A
```



Setting a random seed:



```r
set.seed(100)
```



***



<!-- END OF FRONT MATTER -->


# Overview

Here we will investigate the histone modifications H3K27me2 and H3K27me3 genome-wide
in H3K27M HGG and PFA-EP tumor and cell lines, and in isogenic K27M and K27M-KO contexts,
as shown in Figures 7 and 8.


# Libraries



```r
# Load libraries here
library(here)
library(tidyr)
library(dplyr)
library(readr)
library(glue)
library(purrr)
library(readxl)
library(ggplot2)
library(ggrastr)
library(ggExtra)
library(ggrepel)
library(cowplot)
library(ensurer)
library(plotly)
library(rtracklayer)
library(GenomicRanges)

source(here("include/style.R")) # contains palettes & plotting utils
source(here("code/functions/RNAseq.R"))
source(here("code/functions/testing.R")) # contains ensurer contracts
ggplot2::theme_set(theme_min())
```




# Load metadata

Load the sample metadata for the project:




```r
meta_chip <- read_tsv(here("data/metadata/metadata_chip_all.tsv"))
```



```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   .default = col_character(),
##   Replicate = col_double()
## )
## ℹ Use `spec()` for the full column specifications.
```




# Profile of H3K27me2/me3 enrichment

## Heatmaps

Prepare data for heatmaps, which are then generated int he scripts `07D` and `07E` using [deeptools](https://deeptools.readthedocs.io/en/develop/).



```r
bams_path <- read_lines(here("data/misc/bams_path.tsv"))

meta_chip_parental_hm_input <- meta_chip %>%
    # filter out G477/pcGBM2 here
    filter(!(ID_paper %in% c("G477", "pcGBM2"))) %>% 
    # filter to samples included in the K27M effects analysis
    filter(grepl("K27M effects", Analyses)) %>%
    # filter to H3K27me2/3 marks
    filter(Factor %in% c("H3K27me3", "H3K27me2")) %>% 
    # only keep parental cell lines
    filter(CRISPR == "Parental" | Material == "Tumor") %>% 
    dplyr::select(-matches("Filename"), -Path_bw) %>% 
    dplyr::rename(Path_bw = Path_bw_internal) %>% 
    arrange(Factor, Group, BioID) %>% 
    mutate(Path_bam_with_index = file.path(
        bams_path,
        basename(Path_bam))) %T>%
    write_tsv(glue("{out}/ChIP_metadata_parental.heatmap_input.tsv"))
```



preserve5be4c4dad4d3481a


```r
ensure_present(meta_chip_parental_hm_input$Path_bam)
ensure_present(meta_chip_parental_hm_input$Path_bam_with_index)
```


  </div>
</div>

# Histone mass-spectrometry

Load and plot the MS data:



```r
data_ms <- read_excel(here("data/experimental/2022-06-20-mass_spec_data_points_@ah.xlsx"), sheet = 1)
cols_ms <- read_excel(here("data/experimental/2022-06-20-mass_spec_data_points_@ah.xlsx"), sheet = 2)

data_ms_tidy <- data_ms %>% 
    gather(Sample, Ratio, 2:ncol(.)) %>% 
    left_join(cols_ms, by = "Sample") %>% 
    filter(Peptide %in% c("H3_27_40 K27ac", "H3.1-2-H3K27me2-total", "H3.1-2-H3K27me3-total")) %>% 
    filter(Group %in% c("H3.1K27M", "H3.3K27M", "EZHIP", "H3WT", "ST")) %>% 
    mutate(Group = factor(Group, levels = c("H3.1K27M", "H3.3K27M", "EZHIP", "H3WT", "ST")),
           Peptide = recode(Peptide,
                            "H3_27_40 K27ac" = "H3K27ac",
                            "H3.1-2-H3K27me2-total" = "H3K27me2",
                            "H3.1-2-H3K27me3-total" = "H3K27me3"),
           Peptide = factor(Peptide, levels = c("H3K27me3", "H3K27me2", "H3K27ac")),
           Percentage = Ratio * 100)

data_ms_tidy %>% distinct(Sample, Group) %>% group_by(Group) %>% count()
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Group"],"name":[1],"type":["fct"],"align":["left"]},{"label":["n"],"name":[2],"type":["int"],"align":["right"]}],"data":[{"1":"H3.1K27M","2":"10"},{"1":"H3.3K27M","2":"20"},{"1":"EZHIP","2":"6"},{"1":"H3WT","2":"10"},{"1":"ST","2":"4"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Plot levels per group:



```r
data_ms_tidy %>% 
    group_by(Peptide, Group) %>% 
    summarize(Percentage_mean = mean(Percentage),
              Percentage_sd = sd(Percentage)) %>% 
    ggplot(aes(x = Group, y = Percentage_mean)) +
    geom_bar(aes(fill = Group), stat = "identity", width = 0.5) +
    geom_errorbar(aes(ymin = Percentage_mean - Percentage_sd,
                      ymax = Percentage_mean + Percentage_sd), width = 0.25)  +
    geom_jitter(data = data_ms_tidy, aes(x = Group, y = Percentage), width = 0.2, size = 1) +
    scale_fill_manual(values = c(palette_molecular, "ST" = "gray50")) +
    facet_wrap(~ Peptide, scales = "free_y") +
    rotate_x() +
    ylab("%") +
    no_legend()
```



```
## `summarise()` has grouped output by 'Peptide'. You can override using the `.groups` argument.
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/07A/ms-1.png)<!-- -->

Statistical tests:



```r
stats_ms <- map_dfr(c("H3K27me3", "H3K27me2", "H3K27ac"), function(peptide){
    
    h3.1 <- data_ms_tidy %>% filter(Peptide == peptide & Group == "H3.1K27M") %>% pull(Percentage)
    h3.3 <- data_ms_tidy %>% filter(Peptide == peptide & Group == "H3.3K27M") %>% pull(Percentage)
    pfa  <- data_ms_tidy %>% filter(Peptide == peptide & Group == "EZHIP") %>% pull(Percentage)
    wt   <- data_ms_tidy %>% filter(Peptide == peptide & Group == "H3WT") %>% pull(Percentage)
    st   <- data_ms_tidy %>% filter(Peptide == peptide & Group == "ST") %>% pull(Percentage)
    
    tribble(
        ~ "Comp", ~ "Test", ~ "pvalue",
        "H3.1 vs H3.3", t.test(h3.1, h3.3)$method, t.test(h3.1, h3.3)$p.value,
        "H3.1 vs PFA",  t.test(h3.1, pfa)$method,  t.test(h3.1, pfa)$p.value,
        "H3.3 vs PFA",  t.test(h3.3, pfa)$method,  t.test(h3.3, pfa)$p.value,
        "H3.1 vs H3WT", t.test(h3.1, wt)$method,   t.test(h3.1, wt)$p.value,
        "H3.3 vs H3WT", t.test(h3.3, wt)$method,   t.test(h3.3, wt)$p.value,
        "PFA vs ST",    t.test(pfa, st)$method,    t.test(pfa, st)$p.value,
    ) %>% 
        tibble::add_column("Peptide" = peptide, .before = 1) %>% 
        mutate(Signif = ifelse(pvalue < 0.05, "*", NA)) %>% 
        mutate(pvalue = format(pvalue, digits = 3))
    
}) %T>% 
    write_tsv(glue("{out}/ms_stat_testing.tsv"))

stats_ms
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Peptide"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Comp"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Test"],"name":[3],"type":["chr"],"align":["left"]},{"label":["pvalue"],"name":[4],"type":["chr"],"align":["left"]},{"label":["Signif"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"H3K27me3","2":"H3.1 vs H3.3","3":"Welch Two Sample t-test","4":"5.28e-08","5":"*"},{"1":"H3K27me3","2":"H3.1 vs PFA","3":"Welch Two Sample t-test","4":"1.56e-05","5":"*"},{"1":"H3K27me3","2":"H3.3 vs PFA","3":"Welch Two Sample t-test","4":"2.62e-05","5":"*"},{"1":"H3K27me3","2":"H3.1 vs H3WT","3":"Welch Two Sample t-test","4":"2.14e-07","5":"*"},{"1":"H3K27me3","2":"H3.3 vs H3WT","3":"Welch Two Sample t-test","4":"3.27e-07","5":"*"},{"1":"H3K27me3","2":"PFA vs ST","3":"Welch Two Sample t-test","4":"1.07e-06","5":"*"},{"1":"H3K27me2","2":"H3.1 vs H3.3","3":"Welch Two Sample t-test","4":"8.44e-12","5":"*"},{"1":"H3K27me2","2":"H3.1 vs PFA","3":"Welch Two Sample t-test","4":"1.04e-07","5":"*"},{"1":"H3K27me2","2":"H3.3 vs PFA","3":"Welch Two Sample t-test","4":"1.81e-08","5":"*"},{"1":"H3K27me2","2":"H3.1 vs H3WT","3":"Welch Two Sample t-test","4":"4.07e-11","5":"*"},{"1":"H3K27me2","2":"H3.3 vs H3WT","3":"Welch Two Sample t-test","4":"7.97e-12","5":"*"},{"1":"H3K27me2","2":"PFA vs ST","3":"Welch Two Sample t-test","4":"1.57e-06","5":"*"},{"1":"H3K27ac","2":"H3.1 vs H3.3","3":"Welch Two Sample t-test","4":"9.82e-01","5":"NA"},{"1":"H3K27ac","2":"H3.1 vs PFA","3":"Welch Two Sample t-test","4":"4.40e-03","5":"*"},{"1":"H3K27ac","2":"H3.3 vs PFA","3":"Welch Two Sample t-test","4":"1.75e-05","5":"*"},{"1":"H3K27ac","2":"H3.1 vs H3WT","3":"Welch Two Sample t-test","4":"7.42e-03","5":"*"},{"1":"H3K27ac","2":"H3.3 vs H3WT","3":"Welch Two Sample t-test","4":"5.32e-05","5":"*"},{"1":"H3K27ac","2":"PFA vs ST","3":"Welch Two Sample t-test","4":"2.59e-01","5":"NA"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>




# Pervasive H3K27ac

Load and plot the quantifications of genome-wide H3K27ac:



```r
iso_k27ac <- read_tsv(here("data/ChIPseq/quantifications_@svaradharajan/ISOGENIC_BigwigSignals_H31andH33.txt")) %>% 
    filter(category == "H31")
```






```r
iso_k27ac %>% 
    mutate(RPKM = log10(RPKM + 1)) %>% 
    rr_ggplot(aes(x = variable, y = RPKM), plot_num = 1) +
    geom_boxplot(aes(fill = cat), outlier.shape = NA) +
    facet_wrap(~ sample) +
    scale_fill_manual(values = c("H3K27M" = "orange", "KO" = "blue")) +
    ylim(c(0, 2)) +
    no_legend() +
    ylab("log10(RPKM + 1)") + xlab("Group")
```



```
## ...writing source data of ggplot to public/figures/07A/h3k27ac_pervasive-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/07A/h3k27ac_pervasive-1.png)<!-- -->


# Genomic context of H3K27me3

Here, we'll investigate where in the genome H3K27em3 is retained in K27M gliomas.

[Harutyunyan et al, 2019](https://www.nature.com/articles/s41467-019-09140-x), have 
previously shown that in K27M mutants, H3K27me3 is retained at unmethylated CGIs.
Here, we investigate whether this is consistent for H3.1K27M gliomas by quantifying
H3K27me3 across all tumor types in different genomic regions (CGIs, SUZ12 peaks, inter-
and intra-genic regions, and gene promoters).

The code used in this section was adapted from work by Nicolas De Jay in
[Khazaei et al, Cancer Discovery, 2020](https://cancerdiscovery.aacrjournals.org/content/10/12/1968.abstract).

## Prepare inputs

For this analysis, we'll use the same samples as were included in the heatmaps
visualizations, with the addtions of H3WT cell lines with H3K27me3 data from
Harutyunyan et al, 2019.



```r
meta_chip_parental_deeptools_input <- meta_chip_parental_hm_input %>%
    tibble::add_row(BioID = "G477-WT1",   ID_paper = "G477-WT1",   Factor = "H3K27me3", Group = "HGG-H3WT-Cortex",
                    Path_bam_with_index = here("data/ChIPseq/bam/cell_parental/H3K27me3/G477-WT1_H3K27me3.bam")) %>% 
    tibble::add_row(BioID = "pcGBM2-WT1", ID_paper = "pcGBM2-WT1", Factor = "H3K27me3", Group = "HGG-H3WT-Cortex",
                    Path_bam_with_index = here("data/ChIPseq/bam/cell_parental/H3K27me3/pcGBM2-WT1_H3K27me3.bam")) %T>%
    write_tsv(glue("{out}/ChIP_metadata_parental.deeptools_input.tsv"))
```




## Quantify H3K27me2/3 ChIPseq enrichment

Quantification of H3K27me3 and H3K27me2 in tumors/cell lines is done genome-wide
in 10kb bins using deeptools, by the script `code/07B-deeptools_counts.sh`. The
total number of mapped reads for each mark in each sample is calculated by the script
`code/07C-count_reads.sh`. The outputs from both are loaded below.

## Load counts produced by deeptools

Here, we'll load the counts produced by deeptools:.



```r
# config parameters
input.counts.path = here("output/07B/blacklist_bs10000_raw_MAPQ0_noDup/counts.tab")
input.normalize   = "rpkm"              # either RPKM or RPM
input.bs          = 10000               # bin size
df.chrsToKeep     = paste0("chr", 1:29) # chromosomes to keep
df.signalSums     = NULL                # if ChIP-Rx were being used

# metadata
sample_ids <- read_lines(here("output/07B/samples.txt"))
df.meta <- as.data.frame(meta_chip_parental_deeptools_input)
rownames(df.meta) <- sample_ids
df.meta$sample_id <- sample_ids
df.meta$n_mapped_reads <- as.numeric(read_lines(here("output/07C/n_mapped_reads.txt"), skip = 1))

# load counts
df.counts <- data.table::fread(input.counts.path, header = TRUE, data.table = FALSE) %>% distinct()
colnames(df.counts)[[1]] = "chr"
colnames(df.counts)   %<>% gsub("'", "", .)

# peek
df.counts[1:5, 1:5]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["chr"],"name":[1],"type":["chr"],"align":["left"]},{"label":["start"],"name":[2],"type":["int"],"align":["right"]},{"label":["end"],"name":[3],"type":["int"],"align":["right"]},{"label":["DIPG21__P28__H3K27me2"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["DIPG36__P20__H3K27me2"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"chr3","2":"73161131","3":"73171131","4":"162","5":"166","_rn_":"1"},{"1":"chr3","2":"73171131","3":"73181131","4":"168","5":"256","_rn_":"2"},{"1":"chr3","2":"73181131","3":"73191131","4":"157","5":"391","_rn_":"3"},{"1":"chr3","2":"73191131","3":"73201131","4":"137","5":"436","_rn_":"4"},{"1":"chr3","2":"73201131","3":"73211131","4":"107","5":"118","_rn_":"5"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
dim(df.counts)
```



```
## [1] 312982     33
```



```r
# extract genomic locations and signal counts into different data frames
# (df.loc and df.signal, respectively)
df.loc    <- df.counts[, 1:3]
df.loc %<>% mutate(location = paste0(chr, ":", start, "-", end))
df.signal <- df.counts[, 4:ncol(df.counts), drop = F]

# save in one data structure
ldc <- list(location = df.loc,    # genomic locations
            signal   = df.signal) # per-sample signal for each location
```



preservec9192003cb3371c5


```r
# intersect metadatas and check the same samples exist:
ensure_that(colnames(df.signal), all(. == rownames(df.meta)))
```


  </div>
</div>

Next we normalize using RPKM:



```r
# let's ignore all but the autosomal chromosomes.
df.filter            <- (df.loc[, "chr"] %in% df.chrsToKeep) & (df.signal %>% rowSums(na.rm = TRUE) > 0)

# calculate normalization factor
df.signalSums        <- df.meta$n_mapped_reads / 1e6
names(df.signalSums) <- rownames(df.meta)

# adjust normalization factor using RPKM
df.signalSums        <- df.signalSums * (input.bs / 1000)

# gather results
ldc[["raw"]]         <- ldc[["signal"]]
ldc[["signal"]]      <- t(((t(df.signal) / df.signalSums) )) %>% as.data.frame()
ldc[["filter"]]      <- df.filter
ldc[["signalSums"]]  <- df.signalSums

save(ldc, file = glue("{out}/deeptools_counts.Rda"))
```



## Annotate genomic bins

Next, we aggregate genomic annotations, to annotate each bin based on its genomic context.



```r
# load genomic annotations
# references copied from NDJ @ /project/kleinman/nicolas.dejay/from_hydra/2020/201115-H31vsH33/data/04-input_beds
genes <- here("data/ChIPseq/references/Ensembl.ensGene.whole.hg19.bed")
cgis  <- here("data/ChIPseq/references/IGV.annotations.cpgIslands.hg19.bed")

bed.genic.gr <- rtracklayer::import(genes)

# gene promoters
promoter_distance <- 2500
bed.promoter.gr <- GRanges(
    seqnames = seqnames(bed.genic.gr),
    ranges   = IRanges(start = ifelse(strand(bed.genic.gr) == "+",
                                      start(bed.genic.gr) - promoter_distance, end(bed.genic.gr) - promoter_distance),
                       end   = ifelse(strand(bed.genic.gr) == "+",
                                      start(bed.genic.gr) + promoter_distance, end(bed.genic.gr) + promoter_distance)),
    strand   = strand(bed.genic.gr))
elementMetadata(bed.promoter.gr) = elementMetadata(bed.genic.gr)

# convert dataframe with bin locations to GRanges object where each bin is an interval
# with() evaluates the expression within environment of that dataframe
loc.gr <- with(df.loc, GRanges(chr, IRanges(start = start, end = end)))

# CGI
bed.cgi.gr <- rtracklayer::import(cgis)

# SUZ12 -- taken from Harutyunyan et al, 2019
suz12 <- list(here("data/ChIPseq/peaks/cell_parental/SUZ12/BT245_SUZ12.narrowPeak.bed"),
              here("data/ChIPseq/peaks/cell_parental/SUZ12/DIPGXIII_SUZ12.narrowPeak.bed"),
              here("data/ChIPseq/peaks/cell_parental/SUZ12/G477_SUZ12.narrowPeak.bed"),
              here("data/ChIPseq/peaks/cell_parental/SUZ12/pcGBM2_SUZ12.narrowPeak.bed"))

extra.cols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")

bed.suz12.K27M.gr <- map(suz12[1:2], ~ rtracklayer::import(.x, extraCols = extra.cols)) %>% 
    GRangesList() %>% # covert to GenomicRangesList
    unlist()          # collapse into one set of granges

bed.suz12.WT.gr   <- map(suz12[3:4], ~ rtracklayer::import(.x, extraCols = extra.cols)) %>% 
    GRangesList() %>% # covert to GenomicRangesList
    unlist()          # collapse into one set of granges

# collate
annot.gr <- list(
    "genic"      = bed.genic.gr,
    "promoter"   = bed.promoter.gr,
    "cgi"        = bed.cgi.gr,
    "suz12_K27M" = bed.suz12.K27M.gr,
    "suz12_WT"   = bed.suz12.WT.gr
)

# write out
rtracklayer::export(annot.gr[["genic"]],    glue("{out}/genic.bed"))
rtracklayer::export(annot.gr[["promoter"]], glue("{out}/promoter.bed"))
rtracklayer::export(annot.gr[["cgi"]],      glue("{out}/cgi.bed"))
```



Now, we need to associate each genomic bin (in which ChIPseq signal was counted)
with these annotations in a dataframe that captures the possible many-to-many
relationship between genomic bins and annotation.



```r
# SLOW step
annot.map <- annot.gr %>% sapply(function(gr) {
    # output data frame gives bin indices (column 1) and genomic annotation indices(column 2)
    overlaps = findOverlaps(loc.gr, gr) 
    # construct a dataframe with actual locations of bins / annotations
    data.frame("location" = overlaps %>% queryHits %>% df.loc$location[.],            # which bin
               "gene"     = overlaps %>% subjectHits %>% elementMetadata(gr)$name[.]) # which genomic annotation
}, simplify = FALSE)

# clean up the gene/promoter names
annot.map$promoter %<>%
    dplyr::mutate(gene_short = gene %>%
                      sapply(function (i) strsplit(i %>% as.character, ":") %>%
                                 sapply(`[`, 2:3) %>%
                                 paste(collapse = ":")))

annot.map$genic %<>%
    dplyr::mutate(gene_short = gene %>%
                      sapply(function (i) strsplit(i %>% as.character, ":") %>%
                                 sapply(`[`, 2:3) %>%
                                 paste(collapse = ":")))
```



Now, let's collapse this relationship to a one-to-one relationship between genomic bins
and a list of annotations (for internal readibility).



```r
# SLOW step
annot.squash <- annot.map %>% sapply(function (map) {
    map %>%
        group_by(location) %>%
        summarise("genes" = paste(gene, collapse = ";"),
                  "genes_short" = paste(gene, collapse = ";"))
}, simplify = FALSE)
```



Now, let's also produce a simple classification scheme of bins into the following categories:
genic, promoter, intergenic.



```r
# set defaults
annot.class <- df.loc[, "location", drop = F] %>%
    dplyr::mutate("class"      = "intergenic",
                  "CGI"        = "non-CpG island",
                  "SUZ12_K27M" = "non-SUZ12 peak",
                  "SUZ12_WT"   = "non-SUZ12 peak")

annot.class[df.loc$location %in% unique(annot.squash$genic$location),    "class"] <- "genic"
annot.class[df.loc$location %in% unique(annot.squash$promoter$location), "class"] <- "promoter"
annot.class %<>% dplyr::mutate("class" = factor(class, levels = c("promoter", "genic", "intergenic")))

annot.class[df.loc$location %in% unique(annot.squash$cgi$location), "CGI"] <- "CpG island"
annot.class %<>% dplyr::mutate("CGI" = factor(CGI, levels = c("CpG island", "non-CpG island")))

annot.class[df.loc$location %in% unique(annot.squash$suz12_K27M$location), "SUZ12_K27M"] <- "SUZ12 peak"
annot.class %<>% dplyr::mutate("SUZ12_K27M" = factor(SUZ12_K27M, levels = c("SUZ12 peak", "non-SUZ12 peak")))

annot.class[df.loc$location %in% unique(annot.squash$suz12_WT$location), "SUZ12_WT"] <- "SUZ12 peak"
annot.class %<>% dplyr::mutate("SUZ12_WT" = factor(SUZ12_WT, levels = c("SUZ12 peak", "non-SUZ12 peak")))

# peek
head(annot.class)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["location"],"name":[1],"type":["chr"],"align":["left"]},{"label":["class"],"name":[2],"type":["fct"],"align":["left"]},{"label":["CGI"],"name":[3],"type":["fct"],"align":["left"]},{"label":["SUZ12_K27M"],"name":[4],"type":["fct"],"align":["left"]},{"label":["SUZ12_WT"],"name":[5],"type":["fct"],"align":["left"]}],"data":[{"1":"chr3:73161131-73171131","2":"promoter","3":"non-CpG island","4":"non-SUZ12 peak","5":"non-SUZ12 peak","_rn_":"1"},{"1":"chr3:73171131-73181131","2":"intergenic","3":"non-CpG island","4":"non-SUZ12 peak","5":"non-SUZ12 peak","_rn_":"2"},{"1":"chr3:73181131-73191131","2":"intergenic","3":"non-CpG island","4":"non-SUZ12 peak","5":"non-SUZ12 peak","_rn_":"3"},{"1":"chr3:73191131-73201131","2":"intergenic","3":"non-CpG island","4":"non-SUZ12 peak","5":"non-SUZ12 peak","_rn_":"4"},{"1":"chr3:73201131-73211131","2":"intergenic","3":"non-CpG island","4":"non-SUZ12 peak","5":"non-SUZ12 peak","_rn_":"5"},{"1":"chr3:73211131-73221131","2":"intergenic","3":"non-CpG island","4":"non-SUZ12 peak","5":"non-SUZ12 peak","_rn_":"6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
save(annot.gr, annot.map, annot.squash, annot.class,
     file = glue("{out}/annot.Rda"))
```



## Visualize enrichment of marks in genomic context

Now, we can join the signal, the location, and the genomic annotations.



```r
df.signal.annot <- data.frame(df.loc) %>%
    bind_cols(ldc$signal) %>% 
    left_join(annot.class, by = "location") %>% 
    dplyr::relocate(class, CGI, SUZ12_K27M, SUZ12_WT, .after = location)

# peek
df.signal.annot[1:5, 1:7]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["chr"],"name":[1],"type":["chr"],"align":["left"]},{"label":["start"],"name":[2],"type":["int"],"align":["right"]},{"label":["end"],"name":[3],"type":["int"],"align":["right"]},{"label":["location"],"name":[4],"type":["chr"],"align":["left"]},{"label":["class"],"name":[5],"type":["fct"],"align":["left"]},{"label":["CGI"],"name":[6],"type":["fct"],"align":["left"]},{"label":["SUZ12_K27M"],"name":[7],"type":["fct"],"align":["left"]}],"data":[{"1":"chr3","2":"73161131","3":"73171131","4":"chr3:73161131-73171131","5":"promoter","6":"non-CpG island","7":"non-SUZ12 peak","_rn_":"1"},{"1":"chr3","2":"73171131","3":"73181131","4":"chr3:73171131-73181131","5":"intergenic","6":"non-CpG island","7":"non-SUZ12 peak","_rn_":"2"},{"1":"chr3","2":"73181131","3":"73191131","4":"chr3:73181131-73191131","5":"intergenic","6":"non-CpG island","7":"non-SUZ12 peak","_rn_":"3"},{"1":"chr3","2":"73191131","3":"73201131","4":"chr3:73191131-73201131","5":"intergenic","6":"non-CpG island","7":"non-SUZ12 peak","_rn_":"4"},{"1":"chr3","2":"73201131","3":"73211131","4":"chr3:73201131-73211131","5":"intergenic","6":"non-CpG island","7":"non-SUZ12 peak","_rn_":"5"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

What type of genomic regions retain H3K27me3?

We can look at the top 1% of bins, for each sample, and ask what the distribution
of genomic regions is for each sample/each group.



```r
features <- c("location", "class", "CGI", "SUZ12_K27M", "SUZ12_WT")
samples <- colnames(df.signal.annot)[9:ncol(df.signal.annot)]
samples_k27me3 <- samples[grepl("K27me3", samples)]

df.top.bins.k27me3 <- map_dfr(samples_k27me3, function(i) {
    
    df.signal.annot.sample <- df.signal.annot[, c(features, i)]
    colnames(df.signal.annot.sample)[6] <- "RPKM"
    df.signal.annot.sample %>% filter(RPKM > quantile(.$RPKM, 0.99)) %>% mutate(Sample = i)
    
})

palette_genomic_annotation <- c("CpG island"     = "red3",
                                "non-CpG island" = "gray70",
                                "SUZ12 peak"     = "red3",
                                "non-SUZ12 peak" = "gray70",
                                "genic"          = "darksalmon",
                                "promoter"       = "red3",
                                "intergenic"     = "gray70"
)

df.top.bins.k27me3_long <- df.top.bins.k27me3 %>%
    gather(annotation, status, class, CGI, SUZ12_K27M, SUZ12_WT) %>% 
    dplyr::left_join(df.meta, by = c("Sample" = "sample_id")) %>% 
    mutate(Group = factor(Group, levels = names(palette_groups))) %>% 
    arrange(Group) %>% 
    mutate(Sample = factor(Sample, levels = unique(.$Sample))) %>% 
    mutate(status = factor(status, levels = names(palette_genomic_annotation)))

axis_colours <- map_chr(df.meta %>% filter(Factor == "H3K27me3") %>% pull(Group), ~ palette_groups[.x])
```




Plot by group:



```r
df.top.bins.k27me3_long_prop <- df.top.bins.k27me3_long %>% 
    group_by(Sample) %>% 
    mutate(N_bins_per_sample = length(unique(location))) %>%
    group_by(Group, annotation, Sample, status, N_bins_per_sample) %>% 
    count() %>%
    mutate(Prop = 100 * n / N_bins_per_sample) %>% 
    filter(status %in% c("CpG island", "SUZ12 peak") &
               annotation %in% c("CGI", "SUZ12_K27M")) %>% 
    mutate(
        Group2 = case_when(
            Group == "HGG-H3.1/2K27M-Pons" ~ "H3.1K27M",
            grepl("H3.3", Group) ~ "H3.3K27M",
            grepl("PFA", Group) ~ "PFA",
            Group == "HGG-H3WT-Cortex" ~ "H3WT"),
        Group2 = factor(Group2, levels = c("H3.1K27M", "H3.3K27M", "PFA", "H3WT")),
        Genotype = case_when(
            Group == "HGG-H3.1/2K27M-Pons" ~ "H3.1K27M",
            grepl("H3.3", Group) ~ "H3.3K27M",
            grepl("EZHIP", Group) ~ "EZHIP-PFA",
            Group == "PFA-H3.1K27M-PF" ~ "H3.1K27M-PFA",
            Group == "HGG-H3WT-Cortex" ~ "H3WT"
        ))

df.top.bins.k27me3_long_prop %>% 
    rr_ggplot(aes(x = Group2, y = Prop), plot_num = 1) +
    geom_jitter(aes(fill = Genotype, shape = Genotype), size = 4, alpha = 0.8, width = 0.1) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5) +
    scale_fill_manual(values = palette_genotype) +
    scale_shape_manual(values = c("EZHIP-PFA" = 21, "H3.1K27M-PFA" = 24,
                                  "H3.1K27M" = 21, "H3.3K27M" = 21, "H3WT" = 21)) +
    facet_wrap(~ annotation, scales = "free_y") +
    rotate_x() +
    no_legend() +
    ylim(c(0, 100)) +
    ylab("% of top bins") +
    ggtitle("Restriction of H3K27me3 to CGIs & PRC2 landing sites")
```



```
## ...writing source data of ggplot to public/figures/07A/dotplot_genomic_regions-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/07A/dotplot_genomic_regions-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/07A/dotplot_genomic_regions...*]~</span>

Compute p-values:



```r
# CGIs
h3.1 <- df.top.bins.k27me3_long_prop %>% filter(Genotype == "H3.1K27M"  & annotation == "CGI") %>% pull(n)
h3.3 <- df.top.bins.k27me3_long_prop %>% filter(Genotype == "H3.3K27M"  & annotation == "CGI") %>% pull(n)
pfa  <- df.top.bins.k27me3_long_prop %>% filter(Genotype == "EZHIP-PFA" & annotation == "CGI") %>% pull(n)
wt   <- df.top.bins.k27me3_long_prop %>% filter(Genotype == "H3WT"      & annotation == "CGI") %>% pull(n)

tribble(
    ~ "Comp", ~ "pvalue",
    "H3.1 vs H3WT", t.test(h3.1, wt)$p.value,
    "H3.3 vs H3WT", t.test(h3.3, wt)$p.value,
    "PFA  vs H3WT", t.test(pfa,  wt)$p.value
)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Comp"],"name":[1],"type":["chr"],"align":["left"]},{"label":["pvalue"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"H3.1 vs H3WT","2":"0.0241672794"},{"1":"H3.3 vs H3WT","2":"0.0257190223"},{"1":"PFA  vs H3WT","2":"0.0009043772"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# SUZ12
h3.1 <- df.top.bins.k27me3_long_prop %>% filter(Genotype == "H3.1K27M"  & annotation == "SUZ12_K27M") %>% pull(n)
h3.3 <- df.top.bins.k27me3_long_prop %>% filter(Genotype == "H3.3K27M"  & annotation == "SUZ12_K27M") %>% pull(n)
pfa  <- df.top.bins.k27me3_long_prop %>% filter(Genotype == "EZHIP-PFA" & annotation == "SUZ12_K27M") %>% pull(n)
wt   <- df.top.bins.k27me3_long_prop %>% filter(Genotype == "H3WT"      & annotation == "SUZ12_K27M") %>% pull(n)

tribble(
    ~ "Comp", ~ "pvalue",
    "H3.1 vs H3WT", t.test(h3.1, wt)$p.value,
    "H3.3 vs H3WT", t.test(h3.3, wt)$p.value,
    "PFA  vs H3WT", t.test(pfa,  wt)$p.value
)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Comp"],"name":[1],"type":["chr"],"align":["left"]},{"label":["pvalue"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"H3.1 vs H3WT","2":"0.062988723"},{"1":"H3.3 vs H3WT","2":"0.022642899"},{"1":"PFA  vs H3WT","2":"0.001617436"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

# Quantification of H3K27me3/2 genome-wide

## Load data



```r
# load the quantification in tumors & cell lines
data_k27me3_cgi <- data.table::fread(here("data/ChIPseq/quantifications_@aharutyunyan/20220128-GBM_PFA_K27me3-CGI-v3.txt"),
                                     data.table = FALSE) %>% 
    set_colnames(data.table::fread(here("data/ChIPseq/quantifications_@aharutyunyan/20220128-GBM_PFA_K27me3-CGI-v3_columns.txt"),
                                   data.table = FALSE, header = FALSE)[[2]]) %>% 
    mutate(H3.1_PFA_K27me3 = `P-2077_S-2077_K27me3`)

rr_write_tsv(data_k27me3_cgi,
             glue("{out}/TABLE_K27me3_CGIs.tsv"),
             "Quantification of H3K27me3 ChIPseq signal at CGIs, genome-wide, for each sample and per group")
```



```
## ...writing description of TABLE_K27me3_CGIs.tsv to public/output/07A/TABLE_K27me3_CGIs.desc
```



```r
# load the quantification in ACVR1 mutant and KO lines
data_k27me3_cgi_acvr1 <- data.table::fread(here("data/ChIPseq/quantifications_@aharutyunyan/20220206-DIPGIV-ACVR1-KO_K27me3_CGI_calculated.txt"),
                                           data.table = FALSE)
colnames(data_k27me3_cgi_acvr1) <- gsub("-input-CGI", "", colnames(data_k27me3_cgi_acvr1))

# load the quantifications for K27M OE
data_k27me3_cgi_oe <- data.table::fread(here("data/ChIPseq/quantifications_@aharutyunyan/20220206-K27M_OE_K27me3_CGI_calculated.txt"),
                                        data.table = FALSE)
colnames(data_k27me3_cgi_oe) <- gsub("-input-CGI", "", colnames(data_k27me3_cgi_oe))
colnames(data_k27me3_cgi_oe) <- gsub("-", "_", colnames(data_k27me3_cgi_oe))

# sanity check that features are the same
all(data_k27me3_cgi$name == data_k27me3_cgi_acvr1$name)
```



```
## [1] TRUE
```



```r
all(data_k27me3_cgi$name == data_k27me3_cgi_oe$name)
```



```
## [1] TRUE
```



```r
data_k27me3_cgi_exp <- bind_cols(data_k27me3_cgi_acvr1 %>%
                                     dplyr::rename(DIPGIV_ACVR1_mut_K27me3 = "DIPGIV-K27me3",
                                                   DIPGIV_ACVR1_KO_K27me3  = "DIPGIV-ACVR1-KO-K27me3"),
                                 data_k27me3_cgi_oe %>% select(7:ncol(.)))
```





```r
data_k27me2_100kb <- data.table::fread(here("data/ChIPseq/quantifications_@aharutyunyan/20220128-GBM_PFA_K27me2-100kb-v3.txt"),
                                       data.table = FALSE) %>% 
    set_colnames(data.table::fread(here("data/ChIPseq/quantifications_@aharutyunyan/20220128-GBM_PFA_K27me2-100kb-v3_columns.txt"),
                                   data.table = FALSE, header = FALSE)[[2]]) %>% 
    mutate(H3.1_PFA_K27me2 = `P-2077_S-2077_K27me2`)

rr_write_tsv(data_k27me2_100kb,
             glue("{out}/TABLE_K27me2_100kb_bins.tsv"),
             "Quantification of H3K27me2 ChIPseq signal in 100kb bins, genome-wide, for each sample and per group")
```



```
## ...writing description of TABLE_K27me2_100kb_bins.tsv to public/output/07A/TABLE_K27me2_100kb_bins.desc
```



```r
dim(data_k27me2_100kb)
```



```
## [1] 31419    26
```



```r
save(data_k27me3_cgi, data_k27me2_100kb, file = glue("{out}/data_k27me3_k27me2_quant.Rda"))

# load the quantification in ACVR1 mutant and KO lines
data_k27me2_100kb_acvr1 <- data.table::fread(here("data/ChIPseq/quantifications_@aharutyunyan/20220206-DIPGIV-ACVR1-KO_K27me2-3_K36me2_100kb_calculated.txt"),
                                             data.table = FALSE)
colnames(data_k27me2_100kb_acvr1) <- gsub("-input", "", colnames(data_k27me2_100kb_acvr1))

# load the quantifications for K27M OE
data_k27me2_100kb_oe <- data.table::fread(here("data/ChIPseq/quantifications_@aharutyunyan/20220206-K27M_OE_K27me2_100kb_calculated.txt"),
                                          data.table = FALSE)
colnames(data_k27me2_100kb_oe) <- gsub("-input", "", colnames(data_k27me2_100kb_oe))
colnames(data_k27me2_100kb_oe) <- gsub("-", "_", colnames(data_k27me2_100kb_oe))

data_k27me2_100kb_exp <- bind_cols(data_k27me2_100kb_acvr1 %>%
                                       dplyr::rename(DIPGIV_ACVR1_mut_K27me2 = "DIPGIV-K27me2",
                                                     DIPGIV_ACVR1_KO_K27me2  = "DIPGIV-ACVR1-KO-K27me2") %>%
                                       select(1:6),
                                   data_k27me2_100kb_oe %>% select(5:ncol(.)))
```





## Pairwise comparisons between tumor groups

Here, we'll compare the distribution of marks over genomic regions (CGIs for H3K27me3,
or 100kb-bins for H3K27me2) in different groups.

### H3K27me3



```r
#' Plot a scatterplot of H3K27me3 in CGIs in one group vs. another
#'
#' The data is plot on a log scale, but the labels are in the original
#' un-transformed data space.
#'
#' @param x Unquoted name of column containing H3K27me3 values for x-axis
#' @param y Unquoted name of column containing H3K27me3 values for y-axis
#' @param colour Unquoted column name for binary variable indicating whether
#' a region is marked or not
#' @param title Character, plot title
#' @param x_colour Character, x-axis label colour
#' @param y_colour Character, y-axis label colour
#' @param density Logical, whether to draw 2D density curves on top of points
#' (if TRUE, then \code{rasterize} must be FALSE)
#' @param rasterize Logical, whether to rasterize points layer
#' @param marginal Logical, whether to draw marginal distributions
plot_scatter_k27me3 <- function(x, y, colour, title, x_colour = NULL, y_colour = NULL,
                                density = TRUE,
                                rasterize = TRUE,
                                marginal = FALSE) {
    
    x <- enquo(x)
    y <- enquo(y)
    colour <- enquo(colour)
    
    gg <- data_k27me3_cgi %>%
        ggplot(aes(x = log2(!!x), y = log2(!!y))) +
        geom_abline(slope = 1, colour = "red") +
        geom_hline(yintercept = 0, colour = "gray90") +
        geom_vline(xintercept = 0, colour = "gray90")
    
    if (rasterize) gg <- gg +
        rasterize(geom_point(size = 0.1, alpha = 0.2, aes(colour = !!colour)), dpi = 500)
    else gg <- gg + geom_point(size = 0.1, alpha = 0.2, aes(colour = !!colour))
    
    if (density) gg <- gg +
        geom_density_2d(data = data_k27me3_cgi %>% filter(!!colour), aes(colour = !!colour))
    
    breaks <- c(0.03, 1, 32, 1024)
    
    gg <- gg +
        scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "gray70")) +
        # re-label axes in original data space
        scale_x_continuous(breaks = log2(breaks),
                           labels = breaks,
                           limits = log2(c(0.01, 1024))) +
        scale_y_continuous(breaks = log2(breaks),
                           labels = breaks,
                           limits = log2(c(0.01, 1024))) +
        ggtitle(title) +
        theme(legend.position = "none",
              axis.title.x = element_text(colour = x_colour),
              axis.title.y = element_text(colour = y_colour))
    
    data_marginal <- data_k27me3_cgi %>% filter(!!colour)
    
    if (marginal & !rasterize) ggMarginal(gg, groupColour = TRUE, groupFill = TRUE)
    else if (marginal & rasterize) stop("Cannot produce marginal distributions w/ point rasterization")
    else return(gg)
    
}

data_k27me3_cgi <- data_k27me3_cgi %>%
    mutate(Shown = ifelse(H3.3_HGG_K27me3 < 1024 & H3.1_HGG_K27me3 < 1024, TRUE, FALSE),
           Marked_H3.1_H3.3 = ifelse(H3.1_HGG_K27me3 > 1 | H3.3_HGG_K27me3 > 1, TRUE, FALSE),
           Marked_H3.1_PFA  = ifelse(H3.1_HGG_K27me3 > 1 | EZHIP_PFA_K27me3 > 1, TRUE, FALSE),
           Marked_H3.3_PFA  = ifelse(H3.3_HGG_K27me3 > 1 | EZHIP_PFA_K27me3 > 1, TRUE, FALSE),
           Marked_PFA_PFA   = ifelse(H3.1_PFA_K27me3 > 1 | EZHIP_PFA_K27me3 > 1, TRUE, FALSE))

p1 <- plot_scatter_k27me3(x = H3.3_HGG_K27me3, y = H3.1_HGG_K27me3, colour = Marked_H3.1_H3.3,
                          title = "H3.1 vs. H3.3",
                          "red4", "orange",
                          rasterize = FALSE, marginal = TRUE)


p2 <- plot_scatter_k27me3(x = EZHIP_PFA_K27me3, y = H3.1_HGG_K27me3, colour = Marked_H3.1_PFA,
                          title = "H3.1 vs. PFA",
                          "darkorchid4", "orange",
                          rasterize = FALSE, marginal = TRUE)

p3 <- plot_scatter_k27me3(x = EZHIP_PFA_K27me3, y = H3.3_HGG_K27me3, colour = Marked_H3.3_PFA,
                          title = "H3.3 vs. PFA",
                          "darkorchid4", "red4",
                          rasterize = FALSE, marginal = TRUE)

p4 <- plot_scatter_k27me3(x = EZHIP_PFA_K27me3, y = H3.1_PFA_K27me3, colour = Marked_PFA_PFA,
                          title = "H3.1-PFA vs. EZHIP-PFA",
                          "darkorchid4", "darkorchid4",
                          rasterize = FALSE, marginal = TRUE)

plot_grid(p1, p2, p3, p4, nrow = 1)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/07A/H3K27me3_scatter-1.png)<!-- -->

### H3K27me2



```r
data_k27me2_100kb <- data_k27me2_100kb %>%
    mutate(Shown = ifelse(H3.3_HGG_K27me2 < 1024 & H3.1_HGG_K27me2 < 1024, TRUE, FALSE),
           Marked_H3.1_H3.3 = ifelse(H3.1_HGG_K27me2 > 1 | H3.3_HGG_K27me2 > 1, TRUE, FALSE),
           Marked_H3.1_PFA  = ifelse(H3.1_HGG_K27me2 > 1 | EZHIP_PFA_K27me2 > 1, TRUE, FALSE),
           Marked_H3.3_PFA  = ifelse(H3.3_HGG_K27me2 > 1 | EZHIP_PFA_K27me2 > 1, TRUE, FALSE),
           Marked_PFA_PFA   = ifelse(H3.1_PFA_K27me2 > 1 | EZHIP_PFA_K27me2 > 1, TRUE, FALSE))

#' Plot a scatterplot of H3K27me2 in 100kb-bins in one group vs. another
#' This function follows exactly plot_scatter_k27me3() above
plot_scatter_k27me2 <- function(x, y, colour, title, x_colour = NULL, y_colour = NULL,
                                density = TRUE,
                                rasterize = TRUE,
                                marginal = FALSE) {
    
    x <- enquo(x)
    y <- enquo(y)
    colour <- enquo(colour)
    
    gg <- data_k27me2_100kb %>%
        ggplot(aes(x = log2(!!x), y = log2(!!y))) +
        geom_abline(slope = 1, colour = "red") +
        geom_hline(yintercept = 0, colour = "gray90") +
        geom_vline(xintercept = 0, colour = "gray90")
    
    if (rasterize) gg <- gg +
        rasterize(geom_point(size = 0.1, alpha = 0.2, aes(colour = !!colour)), dpi = 500)
    else gg <- gg + geom_point(size = 0.1, alpha = 0.2, aes(colour = !!colour))
    
    if (density) gg <- gg +
        geom_density_2d(data = data_k27me2_100kb %>% filter(!!colour), aes(colour = !!colour))
    
    breaks <- c(0.2, 0.5, 1, 2, 4, 8)
    
    gg <- gg +
        scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "gray70")) +
        # re-label axes in original data space
        scale_x_continuous(breaks = log2(breaks),
                           labels = breaks,
                           limits = log2(c(0.2, 8))) +
        scale_y_continuous(breaks = log2(breaks),
                           labels = breaks,
                           limits = log2(c(0.2, 8))) +
        ggtitle(title) +
        theme(legend.position = "none", # c(0.1, 0.8),
              axis.title.x = element_text(colour = x_colour),
              axis.title.y = element_text(colour = y_colour))
    
    data_marginal <- data_k27me2_100kb %>% filter(!!colour)
    
    if (marginal & !rasterize) ggMarginal(gg, groupColour = TRUE, groupFill = TRUE)
    else if (marginal & rasterize) stop("Cannot produce marginal distributions w/ point rasterization")
    else return(gg)
    
}

p1 <- plot_scatter_k27me2(x = H3.3_HGG_K27me2, y = H3.1_HGG_K27me2, colour = Marked_H3.1_H3.3,
                          title = "H3.1 vs. H3.3",
                          "red4", "orange",
                          rasterize = FALSE, marginal = TRUE)


p2 <- plot_scatter_k27me2(x = EZHIP_PFA_K27me2, y = H3.1_HGG_K27me2, colour = Marked_H3.1_PFA,
                          title = "H3.1 vs. PFA",
                          "darkorchid4", "orange",
                          rasterize = FALSE, marginal = TRUE)

p3 <- plot_scatter_k27me2(x = EZHIP_PFA_K27me2, y = H3.3_HGG_K27me2, colour = Marked_H3.3_PFA,
                          title = "H3.3 vs. PFA",
                          "darkorchid4", "red4",
                          rasterize = FALSE, marginal = TRUE)

p4 <- plot_scatter_k27me2(x = EZHIP_PFA_K27me2, y = H3.1_PFA_K27me2, colour = Marked_PFA_PFA,
                          title = "H3.1-PFA vs. EZHIP-PFA",
                          "darkorchid4", "darkorchid4",
                          rasterize = FALSE, marginal = TRUE)

plot_grid(p1, p2, p3, p4, nrow = 1)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/07A/H3K27me2_scatter-1.png)<!-- -->

## Numbers of marked regions

We can also binarize the signal of each mark at each region by taking regions with RPKM > $x$.
For H3K27me3, we can compute for each sample, how many CGIs genome-wide have RPKM > $x$.
For H3K27me2, we can compute for each sample, how many 100kb bins genome-wide have RPKM > $x$.

In this case, $x$ is the threshold on RPKM used to consider a region as "marked".
We can sweep $x$ across different values to investigate how the result changes for
different definitions of a region being "marked".

### H3K27me3

For cell lines/tumors:



```r
# convert to long format, with columns CGI, Sample, RPKM, Group
data_k27me3_cgi_long <- data_k27me3_cgi %>% 
    select(ID, matches("K27me3"), -H3.1_HGG_K27me3, -H3.3_HGG_K27me3, -EZHIP_PFA_K27me3, -H3.1_PFA_K27me3, -PFA_K27me3) %>%
    gather(Sample, RPKM, 2:ncol(.)) %>% 
    filter(!is.na(RPKM)) %>%
    mutate(Group_granular = case_when(
        Sample == "P-2077_S-2077_K27me3" ~ "PFA-H3.1K27M",
        Sample %in% c("P-5425_S-6886_K27me3", "P-5426_S-6887_K27me3",
                      "P-5428_S-6889_K27me3", "P-5430_S-6891_K27me3",
                      "PFA4_K27me3", "PFA9_K27me3") ~ "PFA-EZHIP",
        Sample %in% c("DIPGIV_K27me3", "DIPG21_K27me3", "DIPG36_K27me3") ~ "HGG-H3.1/2K27M",
        Sample %in% c("BT869_K27me3", "DIPG007_K27me3",
                      "HSJ019_K27me3", "BT245_K27me3", "DIPGXIII_K27me3",
                      "HSJ-051_K27me3") ~ "HGG-H3.3K27M"
    ), Group = case_when(
        Sample %in% c("P-5425_S-6886_K27me3", "P-5426_S-6887_K27me3",
                      "P-5428_S-6889_K27me3", "P-5430_S-6891_K27me3",
                      "P-2077_S-2077_K27me3", "PFA4_K27me3", "PFA9_K27me3") ~ "PFA",
        Sample %in% c("DIPGIV_K27me3", "DIPG21_K27me3", "DIPG36_K27me3") ~ "HGG-H3.1/2K27M",
        Sample %in% c("BT869_K27me3", "DIPG007_K27me3",
                      "HSJ019_K27me3", "BT245_K27me3", "DIPGXIII_K27me3",
                      "HSJ-051_K27me3") ~ "HGG-H3.3K27M"
    ))

summary(data_k27me3_cgi_long$RPKM)
```



```
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
##  0.000e+00  0.000e+00  1.000e+00 4.716e+306  3.000e+00 1.798e+308
```



```r
n_cgi_marked <- map_dfr(1:10, ~ data_k27me3_cgi_long %>% 
                            # for each sample (in each group), count how many CGIs have RPKM > threshold
                            group_by(Sample, Group, Group_granular) %>% 
                            summarize(N_marked = sum(RPKM > .x), .groups = "drop") %>% 
                            mutate(Threshold = .x))

palette_shape <- c("BT245_K27me3" = 21,
                   "BT869_K27me3" = 21,
                   "DIPG007_K27me3" = 21,
                   "DIPG21_K27me3" = 23,
                   "DIPG36_K27me3" = 25,
                   "DIPGIV_K27me3" = 21,
                   "DIPGXIII_K27me3" = 21,
                   "HSJ-051_K27me3" = 21,
                   "HSJ019_K27me3" = 21,
                   "P-2077_S-2077_K27me3" = 24,
                   "P-5425_S-6886_K27me3" = 21,
                   "P-5426_S-6887_K27me3" = 21,
                   "P-5428_S-6889_K27me3" = 21,
                   "P-5430_S-6891_K27me3" = 21,
                   "PFA4_K27me3" = 21,
                   "PFA9_K27me3" = 21)

n_cgi_marked_threshold_3 <- n_cgi_marked %>% 
    filter(Threshold == 3)

h3.1 <- n_cgi_marked_threshold_3 %>% filter(Group == "HGG-H3.1/2K27M") %>% pull(N_marked)
h3.3 <- n_cgi_marked_threshold_3 %>% filter(Group == "HGG-H3.3K27M") %>% pull(N_marked)
pfa  <- n_cgi_marked_threshold_3 %>% filter(Group == "PFA") %>% pull(N_marked)

tribble(
    ~ "Comp", ~ "Test", ~ "pvalue",
    "H3.1 vs H3.3", t.test(h3.1, h3.3)$method, t.test(h3.1, h3.3)$p.value,
    "H3.1 vs PFA",  t.test(h3.1, pfa)$method, t.test(h3.1, pfa)$p.value,
    "H3.3 vs PFA",  t.test(h3.3, pfa)$method, t.test(h3.3, pfa)$p.value
)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Comp"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Test"],"name":[2],"type":["chr"],"align":["left"]},{"label":["pvalue"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"H3.1 vs H3.3","2":"Welch Two Sample t-test","3":"0.54658352"},{"1":"H3.1 vs PFA","2":"Welch Two Sample t-test","3":"0.02164242"},{"1":"H3.3 vs PFA","2":"Welch Two Sample t-test","3":"0.01003385"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div><br><span style="color:#0d00ff">~[figure @ *public/figures/07A/dotplots_cgi_H3K27me3...*]~</span>




```r
# dotplot per sample, grouped by molecular group
n_cgi_marked %>%
    filter(Threshold == 3) %>% 
    ggplot(aes(x = Group, y = N_marked, label = Sample)) +
    geom_jitter(aes(fill = Group_granular, shape = Sample), width = 0.1, alpha = 0.8, size = 4) +
    scale_fill_manual(values = palette_groups_simple, guide = FALSE) +
    scale_shape_manual(values = palette_shape) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5) +
    rotate_x() +
    scale_y_continuous(labels = scales::comma) +
    ggtitle("Number of CGIs genome-wide marked by H3K27me3 \n(RPKM > 3)") +
    ylab("# CGIs") +
    ylim(c(0, 13000)) +
    theme(legend.position = "bottom")
```



```
## Scale for 'y' is already present. Adding another scale for 'y', which will
## replace the existing scale.
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/07A/dotplots_cgi_H3K27me3_exp-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/07A/dotplots_cgi_H3K27me3_exp...*]~</span>



### H3K27me2



```r
data_k27me2_100kb_long <- data_k27me2_100kb %>% 
    mutate(ID = paste0(chr, ":", start, "-", end)) %>% 
    select(ID, matches("K27me2"), -H3.1_HGG_K27me2, -H3.3_HGG_K27me2,
           -EZHIP_PFA_K27me2, -H3.1_PFA_K27me2, -PFA_K27me2) %>%
    gather(Sample, RPKM, 2:ncol(.)) %>% 
    filter(!is.na(RPKM)) %>%
    mutate(Group_granular = case_when(
        Sample == "P-2077_S-2077_K27me2" ~ "PFA-H3.1K27M",
        Sample %in% c("P-5425_S-6886_K27me2", "P-5426_S-6887_K27me2",
                      "P-5428_S-6889_K27me2", "P-5430_S-6891_K27me2",
                      "PFA2_K27me2", "PFA4_K27me2", "PFA5_K27me2") ~ "PFA-EZHIP",
        Sample %in% c("DIPGIV_K27me2", "DIPG21_K27me2", "DIPG36_K27me2") ~ "HGG-H3.1/2K27M",
        Sample %in% c("BT869_K27me2", "DIPG007_K27me2", "HSJ019_K27me2", "BT245_K27me2",
                      "DIPGXIII_K27me2", "HSJ-051_K27me2") ~ "HGG-H3.3K27M"
    ), Group = case_when(
        Sample %in% c("P-5425_S-6886_K27me2", "P-5426_S-6887_K27me2",
                      "P-5428_S-6889_K27me2", "P-5430_S-6891_K27me2",
                      "P-2077_S-2077_K27me2", "PFA2_K27me2",
                      "PFA4_K27me2", "PFA5_K27me2") ~ "PFA",
        Sample %in% c("DIPGIV_K27me2", "DIPG21_K27me2", "DIPG36_K27me2") ~ "HGG-H3.1/2K27M",
        Sample %in% c("BT869_K27me2", "DIPG007_K27me2",
                      "HSJ019_K27me2", "BT245_K27me2", "DIPGXIII_K27me2",
                      "HSJ-051_K27me2") ~ "HGG-H3.3K27M"
    ))

n_100kb_marked <- map_dfr(seq(0.6, 1.5, by = 0.1), ~ data_k27me2_100kb_long %>% 
                              # for each sample (in each group), count how many CGIs have RPKM > threshold
                              group_by(Sample, Group, Group_granular) %>% 
                              summarize(N_marked = sum(RPKM > .x), .groups = "drop") %>% 
                              mutate(Threshold = .x))

palette_shape <- c("BT245_K27me2" = 21,
                   "BT869_K27me2" = 21,
                   "DIPG007_K27me2" = 21,
                   "DIPG21_K27me2" = 23,
                   "DIPG36_K27me2" = 25,
                   "DIPGIV_K27me2" = 21,
                   "DIPGXIII_K27me2" = 21,
                   "HSJ-051_K27me2" = 21,
                   "HSJ019_K27me2" = 21,
                   "P-2077_S-2077_K27me2" = 24,
                   "P-5425_S-6886_K27me2" = 21,
                   "P-5426_S-6887_K27me2" = 21,
                   "P-5428_S-6889_K27me2" = 21,
                   "P-5430_S-6891_K27me2" = 21,
                   "PFA4_K27me2" = 21,
                   "PFA9_K27me2" = 21,
                   "PFA2_K27me2" = 21,
                   "PFA5_K27me2" = 21)

n_100kb_marked_threshold_1 <- n_100kb_marked %>% 
    filter(Threshold == 1)

h3.1 <- n_100kb_marked_threshold_1 %>% filter(Group == "HGG-H3.1/2K27M") %>% pull(N_marked)
h3.3 <- n_100kb_marked_threshold_1 %>% filter(Group == "HGG-H3.3K27M") %>% pull(N_marked)
pfa <- n_100kb_marked_threshold_1 %>% filter(Group == "PFA") %>% pull(N_marked)

tribble(
    ~ "Comp", ~ "Test", ~ "pvalue",
    "H3.1 vs H3.3", t.test(h3.1, h3.3)$method, t.test(h3.1, h3.3)$p.value,
    "H3.1 vs PFA",  t.test(h3.1, pfa)$method, t.test(h3.1, pfa)$p.value,
    "H3.3 vs PFA",  t.test(h3.3, pfa)$method, t.test(h3.3, pfa)$p.value
)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Comp"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Test"],"name":[2],"type":["chr"],"align":["left"]},{"label":["pvalue"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"H3.1 vs H3.3","2":"Welch Two Sample t-test","3":"1.709237e-01"},{"1":"H3.1 vs PFA","2":"Welch Two Sample t-test","3":"1.353987e-02"},{"1":"H3.3 vs PFA","2":"Welch Two Sample t-test","3":"1.807836e-05"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div><br><span style="color:#0d00ff">~[figure @ *public/figures/07A/dotplots_100kb_H3K27me2...*]~</span>




```r
# dotplot per sample, grouped by molecular group
n_100kb_marked %>%
    filter(Threshold == 1) %>% 
    ggplot(aes(x = Group, y = N_marked, label = Sample)) +
    geom_jitter(aes(fill = Group_granular, shape = Sample), width = 0.1, alpha = 0.8, size = 4) +
    scale_fill_manual(values = palette_groups_simple) +
    scale_shape_manual(values = palette_shape) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5) +
    rotate_x() +
    scale_y_continuous(labels = scales::comma) +
    ggtitle("Number of 100kb-bins genome-wide marked by H3K27me2 \n(RPKM > 1)") +
    ylab("# bins") +
    no_legend() +
    ylim(c(0, 20000)) +
    theme(legend.position = "bottom")
```



```
## Scale for 'y' is already present. Adding another scale for 'y', which will
## replace the existing scale.
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/07A/dotplots_100kb_H3K27me2_exp-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/07A/dotplots_100kb_H3K27me2_exp...*]~</span>


## H3K27me2 domains

Plot the total length of the genome covered by H3K27me2 domains:



```r
data_k27me2_domains <- read_xlsx(here("data/ChIPseq/quantifications_@mhulswit/20211006-H3K27me2_domains.xlsx")) %>% 
    dplyr::rename(Total_coverage = `Total genome pasted with K27me2 domains`) %>% 
    mutate(Genotype = case_when(
        Category == "H3.1K27M_KO"  ~ "H3.1K27M-KO",
        Category == "H3.1K27M_NKO" ~ "H3.1K27M",
        Category == "H3.3K27M_KO"  ~ "H3.3K27M-KO",
        Category == "H3.3K27M_NKO" ~ "H3.3K27M",
        TRUE ~ Category
    )) %>% 
    separate(Sample, into = c("Cell_line", "Extra"), sep = "_") %>% 
    # DIPG21 excluded from all analysies in the paper due to high background
    filter(Cell_line != "DIPG21") %>% 
    mutate(Genotype = factor(Genotype, levels =
                                 c("H3.1K27M", "H3.1K27M-KO", "H3.3K27M", "H3.3K27M-KO", "EZHIP-PFA", "H3.1K27M-PFA")))

palette_cell_line <- c(palette_cell_line,
                       c(21, 21, 21) %>% setNames(setdiff(unique(data_k27me2_domains$Cell_line), names(palette_cell_line))))

palette_genotype <- c("H3.1K27M"     = "orange",
                      "H3.3K27M"     = "red3",
                      "KO"           = "blue",
                      "H3.1K27M-KO"  = "orange",
                      "H3.3K27M-KO"  = "red3",
                      "EZHIP-PFA"    = "darkorchid4",
                      "H3.1K27M-PFA" = "mediumorchid")

data_k27me2_domains %>% 
    rr_ggplot(aes(x = Genotype, y = Total_coverage), plot_num = 1) +
    geom_jitter(aes(fill = Genotype, shape = Cell_line), size = 4, alpha = 0.8, width = 0.1) +
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.5) +
    scale_fill_manual(values = palette_genotype) +
    scale_shape_manual(values = palette_cell_line) +
    scale_y_continuous(labels = scales::comma) +
    ggtitle("Total length of genome marked \nby H3K27me2 domains") +
    rotate_x()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/07A/H3K27me2_domains-1.png)<!-- -->

Run T-tests to compare total H3K27me2 domain coverage between groups:



```r
# H3.1
h3.1 <- data_k27me2_domains %>% filter(Genotype == "H3.1K27M") %>% pull(Total_coverage)
h3.1_ko <- data_k27me2_domains %>% filter(Genotype == "H3.1K27M-KO") %>% pull(Total_coverage)
t.test(h3.1, h3.1_ko)
```



```
## 
## 	Welch Two Sample t-test
## 
## data:  h3.1 and h3.1_ko
## t = -9.1735, df = 4.1651, p-value = 0.0006471
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -1643975182  -889144985
## sample estimates:
##  mean of x  mean of y 
##  756343667 2022903750
```



```r
# H3.3
h3.3 <- data_k27me2_domains %>% filter(Genotype == "H3.3K27M") %>% pull(Total_coverage)
h3.3_ko <- data_k27me2_domains %>% filter(Genotype == "H3.3K27M-KO") %>% pull(Total_coverage)
t.test(h3.3, h3.3_ko)
```



```
## 
## 	Welch Two Sample t-test
## 
## data:  h3.3 and h3.3_ko
## t = -6.5637, df = 6.2096, p-value = 0.0005205
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -955905885 -439814781
## sample estimates:
##  mean of x  mean of y 
## 1141829833 1839690167
```



Total domain length:



```r
meta_chip_iso <- read_xlsx(here("data/metadata/2022-02-15-Isogenic_ChIP_metadata.xlsx"))
rx_chip_iso <- map_dfr(c("H3K27me3", "H3K27me2"), ~ read_tsv(here(glue("data/ChIPseq/2022-02-24-Rx_values_{.x}.tsv"))))

data_k27me2_domains_oe <- read_tsv(here("data/ChIPseq/quantifications_@mhulswit/20220228-TotalLengthMe2Domains_OE.txt")) %>% 
    dplyr::rename(Total_coverage = Total_Length_H3K27me2_domains) %>% 
    left_join(meta_chip_iso, by = c("Cell line" = "SampleID", "Condition")) %>% 
    mutate(Condition = gsub("-", " ", Condition)) %>% 
    mutate(Condition = factor(Condition, levels = names(palette_condition)),
           `Cell line` = factor(`Cell line`, levels = c("G477", "DIPG36", "BT245", "DIPGXIII", "DIPGIV"))) %>% 
    arrange(`Cell line`, Condition) %>% 
    mutate(Sample = factor(Sample, levels = unique(.$Sample)))
```



Plot distribution of domains:



```r
domains_paths <- list.files(here("data/ChIPseq/domains_@mhulswit/"), full.names = TRUE, pattern = "*.bed")
domains <- map(domains_paths, ~ rtracklayer::import(.x))
names(domains) <- basename(domains_paths)

domains_lengths <- imap_dfr(domains, ~ data.frame("Sample" = .y,
                                                  "Widths" = GenomicRanges::width(.x))) %>%
    left_join(data_k27me2_domains_oe, by = "Sample") %T>% 
    write_tsv(glue("{figout}/H3K27me2_domains_OE_distribution.source_data.tsv"))

plot_dist <- function(cl, ymax) {
    
    p1 <- domains_lengths %>% 
        filter(`Cell line` == cl) %>% 
        # swap the order because we're rotating coordinates
        mutate(Condition = factor(Condition, levels = rev(levels(.$Condition)))) %>% 
        ggplot(aes(x = Condition, y = Widths)) +
        geom_boxplot(aes(fill = Condition), alpha = 0.7, width = 0.5, outlier.colour = NA) +
        scale_fill_manual(values = palette_condition) +
        coord_flip(ylim = c(0, ymax)) +
        scale_y_continuous(labels = scales::comma,
                           name = "Distribution of H3K27me2 domain length") +
        no_legend() +
        rotate_x() +
        ggtitle(cl)
    
    show(p1)
    ggsave(plot = p1, glue("{figout}/H3K27me2_domains_OE_distribution_{cl}.pdf"), width = 5, height = 3)
    
}

plot_dist("BT245",    2750000)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/07A/domain_length_dist-1.png)<!-- -->

```r
plot_dist("DIPGIV",   350000)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/07A/domain_length_dist-2.png)<!-- -->

```r
plot_dist("DIPGXIII", 3000000)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/07A/domain_length_dist-3.png)<!-- -->

```r
# print N
domains_lengths %>%
    filter(`Cell line` %in% c("BT245", "DIPGIV", "DIPGXIII")) %>%
    group_by(`Cell line`, Condition) %>%
    count()
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Cell line"],"name":[1],"type":["fct"],"align":["left"]},{"label":["Condition"],"name":[2],"type":["fct"],"align":["left"]},{"label":["n"],"name":[3],"type":["int"],"align":["right"]}],"data":[{"1":"BT245","2":"H3.3K27M","3":"8436"},{"1":"BT245","2":"H3.3K27M KO + H3.1WT OE","3":"3782"},{"1":"BT245","2":"H3.3K27M KO + H3.1K27M OE","3":"8200"},{"1":"DIPGXIII","2":"H3.3K27M","3":"16630"},{"1":"DIPGXIII","2":"H3.3K27M KO + H3.1WT OE","3":"3388"},{"1":"DIPGXIII","2":"H3.3K27M KO + H3.1K27M OE","3":"11568"},{"1":"DIPGIV","2":"H3.1K27M + ACVR1 mut","3":"11736"},{"1":"DIPGIV","2":"H3.1K27M + ACVR1 KO","3":"10614"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


<!-- END MATTER, insert reproducibility info -->




***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:



```
## 2022-09-12 15:34:31
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
##  ! package              * version    date       lib
##  P assertthat             0.2.1      2019-03-21 [?]
##  P beeswarm               0.3.1      2021-03-07 [?]
##  P Biobase                2.46.0     2019-10-29 [?]
##  P BiocGenerics         * 0.32.0     2019-10-29 [?]
##  P BiocParallel           1.20.1     2019-12-21 [?]
##  P Biostrings             2.54.0     2019-10-29 [?]
##  P bitops                 1.0-7      2021-04-24 [?]
##  P bslib                  0.2.5      2021-05-12 [?]
##  P callr                  3.7.0      2021-04-20 [?]
##  P cellranger             1.1.0      2016-07-27 [?]
##  P cli                    2.5.0      2021-04-26 [?]
##  P codetools              0.2-16     2018-12-24 [?]
##  P colorspace             2.0-1      2021-05-04 [?]
##  P cowplot              * 1.1.1      2020-12-30 [?]
##  P crayon                 1.4.1      2021-02-08 [?]
##  P data.table             1.14.0     2021-02-21 [?]
##  P DBI                    1.1.1      2021-01-15 [?]
##  P DelayedArray           0.12.3     2020-04-09 [?]
##  P desc                   1.2.0      2018-05-01 [?]
##  P devtools               2.3.0      2020-04-10 [?]
##  P digest                 0.6.27     2020-10-24 [?]
##  P dplyr                * 1.0.6      2021-05-05 [?]
##  P ellipsis               0.3.2      2021-04-29 [?]
##  P ensurer              * 1.0        2021-12-06 [?]
##  P evaluate               0.14       2019-05-28 [?]
##  P fansi                  0.4.2      2021-01-15 [?]
##  P farver                 2.1.0      2021-02-28 [?]
##  P fastmap                1.1.0      2021-01-25 [?]
##  P fs                     1.5.0      2020-07-31 [?]
##  P generics               0.1.0      2020-10-31 [?]
##  P GenomeInfoDb         * 1.22.1     2020-03-27 [?]
##  P GenomeInfoDbData       1.2.2      2021-12-06 [?]
##  P GenomicAlignments      1.22.1     2019-11-12 [?]
##  P GenomicRanges        * 1.38.0     2019-10-29 [?]
##  P ggbeeswarm             0.6.0      2017-08-07 [?]
##  P ggExtra              * 0.10.0     2022-03-23 [?]
##  P ggplot2              * 3.3.3      2020-12-30 [?]
##  P ggrastr              * 0.2.3      2021-03-01 [?]
##  P ggrepel              * 0.9.1      2021-01-15 [?]
##  P git2r                  0.27.1     2020-05-03 [?]
##  P glue                 * 1.4.2      2020-08-27 [?]
##  P gridExtra              2.3        2017-09-09 [?]
##  P gtable                 0.3.0      2019-03-25 [?]
##  P here                 * 0.1        2017-05-28 [?]
##  P highr                  0.9        2021-04-16 [?]
##  P hms                    1.0.0      2021-01-13 [?]
##  P htmltools              0.5.1.1    2021-01-22 [?]
##  P htmlwidgets            1.5.3      2020-12-10 [?]
##  P httpuv                 1.6.1      2021-05-07 [?]
##  P httr                   1.4.2      2020-07-20 [?]
##  P IRanges              * 2.20.2     2020-01-13 [?]
##  P jquerylib              0.1.4      2021-04-26 [?]
##  P jsonlite               1.7.2      2020-12-09 [?]
##  P knitr                  1.33       2021-04-24 [?]
##  P labeling               0.4.2      2020-10-20 [?]
##  P later                  1.0.0      2019-10-04 [?]
##  P lattice                0.20-44    2021-05-02 [?]
##  P lazyeval               0.2.2      2019-03-15 [?]
##  P lifecycle              1.0.0      2021-02-15 [?]
##  P magrittr             * 2.0.1      2020-11-17 [?]
##  P Matrix                 1.2-18     2019-11-27 [?]
##  P matrixStats            0.58.0     2021-01-29 [?]
##  P memoise                1.1.0      2017-04-21 [?]
##  P mime                   0.10       2021-02-13 [?]
##  P miniUI                 0.1.1.1    2018-05-18 [?]
##  P munsell                0.5.0      2018-06-12 [?]
##  P pillar                 1.6.0      2021-04-13 [?]
##  P pkgbuild               1.0.8      2020-05-07 [?]
##  P pkgconfig              2.0.3      2019-09-22 [?]
##  P pkgload                1.0.2      2018-10-29 [?]
##  P plotly               * 4.9.3      2021-01-10 [?]
##  P prettyunits            1.1.1      2020-01-24 [?]
##  P processx               3.5.2      2021-04-30 [?]
##  P promises               1.1.0      2019-10-04 [?]
##  P ps                     1.6.0      2021-02-28 [?]
##  P purrr                * 0.3.4      2020-04-17 [?]
##  P R6                     2.5.0      2020-10-28 [?]
##  P RColorBrewer         * 1.1-2      2014-12-07 [?]
##  P Rcpp                   1.0.6      2021-01-15 [?]
##  P RCurl                  1.98-1.3   2021-03-16 [?]
##  P readr                * 1.4.0      2020-10-05 [?]
##  P readxl               * 1.3.1      2019-03-13 [?]
##  P remotes                2.1.1      2020-02-15 [?]
##    renv                   0.14.0     2021-07-21 [1]
##  P rlang                  0.4.11     2021-04-30 [?]
##  P rmarkdown              2.8        2021-05-07 [?]
##  P rprojroot              2.0.2      2020-11-15 [?]
##  P Rsamtools              2.2.3      2020-02-23 [?]
##  P rstudioapi             0.13       2020-11-12 [?]
##  P rtracklayer          * 1.46.0     2019-10-29 [?]
##  P S4Vectors            * 0.24.4     2020-04-09 [?]
##  P sass                   0.4.0      2021-05-12 [?]
##  P scales                 1.1.1      2020-05-11 [?]
##  P sessioninfo            1.1.1      2018-11-05 [?]
##  P shiny                  1.6.0      2021-01-25 [?]
##  P stringi                1.6.1      2021-05-10 [?]
##  P stringr                1.4.0      2019-02-10 [?]
##  P SummarizedExperiment   1.16.1     2019-12-19 [?]
##  P testrmd                0.0.1.9000 2021-12-06 [?]
##  P testthat               2.3.2      2020-03-02 [?]
##  P tibble                 3.1.1      2021-04-18 [?]
##  P tidyr                * 1.1.3      2021-03-03 [?]
##  P tidyselect             1.1.1      2021-04-30 [?]
##  P usethis                1.6.1      2020-04-29 [?]
##  P utf8                   1.2.1      2021-03-12 [?]
##  P vctrs                  0.3.8      2021-04-29 [?]
##  P vipor                  0.4.5      2017-03-22 [?]
##  P viridis              * 0.5.1      2018-03-29 [?]
##  P viridisLite          * 0.4.0      2021-04-13 [?]
##  P whisker                0.4        2019-08-28 [?]
##  P withr                  2.4.2      2021-04-18 [?]
##  P xfun                   0.22       2021-03-11 [?]
##  P XML                    3.99-0.3   2020-01-20 [?]
##  P xtable                 1.8-4      2019-04-21 [?]
##  P XVector                0.26.0     2019-10-29 [?]
##  P yaml                   2.2.1      2020-02-01 [?]
##  P zlibbioc               1.32.0     2019-10-29 [?]
##  source                           
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
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
##  Bioconductor                     
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
##  Github (smbache/ensurer@feb1def) 
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
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
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
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
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
##  Bioconductor                     
##  Bioconductor                     
##  CRAN (R 3.6.1)                   
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
##  CRAN (R 3.6.1)                   
##  CRAN (R 3.6.1)                   
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
##  Bioconductor                     
## 
## [1] /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/renv/library/R-3.6/x86_64-pc-linux-gnu
## [2] /tmp/RtmpsrGsxB/renv-system-library
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
