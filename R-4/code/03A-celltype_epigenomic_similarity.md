---
title: "03A - Epigenomic similarity between tumors & COO"
author: "Selin Jessa [[selin.jessa@mail.mcgill.ca](mailto:selin.jessa@mail.mcgill.ca)]"
date: "16 June, 2022"
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
## Document index: 03A
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
## public/R-4/output/03A
```

```
## public/R-4/figures/03A
```



Setting a random seed:

```r
set.seed(100)
```

***

<!-- END OF FRONT MATTER -->


# Overview

Here we'll examine the similarity between HGGs and PFAs to their respective
cells-of-origin based on epigenomic data across genomic features.

For this analysis, we obtained Paired-Tag data for different histone modifications
and corresponding cell type specific epigenome prifles from [Zhu et al, Nature Methods, 2021](https://doi.org/10.1038/s41592-021-01060-3).

# Libraries


```r
library(here)

# For genome analysis / references
library(GenomicRanges)
library(plotgardener)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(rtracklayer)

# general purpose
library(tidyr)
library(dplyr)
library(ggrepel)
library(ggrastr)
library(readr)
library(glue)
library(tibble)
library(ggplot2)
library(purrr)
library(pheatmap)
library(stringr)
library(fgsea)
library(icytobox)
library(cowplot)

source(here("include/style.R"))
source(here("R-4/code/functions/plotgardener_helpers.R"))
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



```r
meta_chip <- read_tsv(here("data/metadata/metadata_chip_all.tsv"))
```

```
## Rows: 172 Columns: 23
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (22): BioID, ID_paper, Location, Group, Factor, Material, CRISPR, Clone_...
## dbl  (1): Replicate
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```


# Load data

Import ChromHMM state calls from Zhu et al using `rtracklayer::import())`:


```r
# OPC is cluster 15 and ependymal is cluster 22
opc_chromhmm <- import(
    here("data/Paired-Tag/Zhu_2021/ChromHMM/15_8_segments.bed"))
epen_chromhmm <- import(
    here("data/Paired-Tag/Zhu_2021/ChromHMM/22_8_segments.bed"))

# subset to autosomes
chrs_keep <- paste0("chr", 1:29)
opc_chromhmm <- opc_chromhmm[seqnames(opc_chromhmm) %in% chrs_keep]
epen_chromhmm <- epen_chromhmm[seqnames(epen_chromhmm) %in% chrs_keep]

# load segment info
read_tsv(here("data/Paired-Tag/Zhu_2021/ChromHMM/README_segment.txt"))
```

```
## Rows: 7 Columns: 2
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (2): 2-Promoter-Weak, E1
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["2-Promoter-Weak"],"name":[1],"type":["chr"],"align":["left"]},{"label":["E1"],"name":[2],"type":["chr"],"align":["left"]}],"data":[{"1":"1-Promoter-Active","2":"E2"},{"1":"3-Enhancer-Active","2":"E3"},{"1":"4-Enhancer-Poised","2":"E4"},{"1":"7-No signal","2":"E5"},{"1":"8-No signal-Gene Desert","2":"E6"},{"1":"6-Heterochromatin-H3K9me3","2":"E7"},{"1":"5-Heterochromatin-H3K27me3","2":"E8"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

# Get cell-type specific genomic regions

First, subset to active/inactive genomic regions for each cell type:


```r
active_states <- c("E1", "E2", "E3")
inactive_states <- c("E7", "E8")

opc_active   <- opc_chromhmm[(elementMetadata(opc_chromhmm)[, "name"] %in% active_states)]
opc_inactive <- opc_chromhmm[(elementMetadata(opc_chromhmm)[, "name"] %in% inactive_states)]

epen_active   <- epen_chromhmm[(elementMetadata(epen_chromhmm)[, "name"] %in% active_states)]
epen_inactive <- epen_chromhmm[(elementMetadata(epen_chromhmm)[, "name"] %in% inactive_states)]
```

Next, look for overlaps between active regions in one cell type and inactive in the other,
requiring 50% overlap


```r
# regions active in OPC and inactive in epen
hits <- findOverlaps(opc_active, epen_inactive)
overlaps <- pintersect(opc_active[queryHits(hits)], epen_inactive[subjectHits(hits)])
length(overlaps)
```

```
## [1] 418
```

```r
percentOverlap <- width(overlaps) / width(opc_active[subjectHits(hits)])
hits <- hits[percentOverlap > 0.1]
opc_active_specific <- opc_active[queryHits(hits)]
length(opc_active_specific)
```

```
## [1] 416
```

```r
# regions active in epen and inactive in OPC
hits <- findOverlaps(epen_active, opc_inactive)
overlaps <- pintersect(epen_active[queryHits(hits)], opc_inactive[subjectHits(hits)])
length(overlaps)
```

```
## [1] 291
```

```r
percentOverlap <- width(overlaps) / width(epen_active[subjectHits(hits)])
hits <- hits[percentOverlap > 0.1]
epen_active_specific <- epen_active[queryHits(hits)]
length(epen_active_specific)
```

```
## [1] 291
```

That gives us 200-300 features which are specific to one cell type or another.
The union over the two will be our genomic features:


```r
elementMetadata(opc_active_specific)$celltype <- "OPC"
elementMetadata(epen_active_specific)$celltype <- "Ependymal"
```

Load BED file of mm10 genes, and look for overlaps:


```r
mm10_genes_bed <- rtracklayer::import(here("data/ChIPseq/references/Ensembl.ensGene.mm10.collapsed.bed"))
mm10_genes_bed <- mm10_genes_bed[seqnames(mm10_genes_bed) %in% chrs_keep]

# add 5kb on either side
start(mm10_genes_bed) <- start(mm10_genes_bed) - 2500
# end(mm10_genes_bed)   <- end(mm10_genes_bed) + 5000

# find nearest gene for each OPC element
opc_nearest_genes_idx <- nearest(opc_active_specific, mm10_genes_bed)
opc_nearest_genes <- mm10_genes_bed[opc_nearest_genes_idx]
opc_nearest_genes_hg <- opc_nearest_genes$name %>% sapply(str_split, ":") %>% 
    sapply(getElement, 2) %>% 
    mm2hg() %>% 
    unlist() %>%
    unique()

# how many?
length(unique(opc_nearest_genes$name))
```

```
## [1] 300
```

```r
# same thing for ependymal elements
epen_nearest_genes_idx <- nearest(epen_active_specific, mm10_genes_bed)
epen_nearest_genes <- mm10_genes_bed[epen_nearest_genes_idx]
length(unique(epen_nearest_genes$name))
```

```
## [1] 197
```

```r
epen_nearest_genes_hg <- epen_nearest_genes$name %>% sapply(str_split, ":") %>% 
    sapply(getElement, 2) %>% 
    mm2hg() %>% 
    unlist() %>%
    unique()


# mm
epen_nearest_genes_mm <- epen_nearest_genes$name %>% sapply(str_split, ":") %>%
    sapply(getElement, 2) %>% unname() %>% unique()
opc_nearest_genes_mm <- opc_nearest_genes$name %>% sapply(str_split, ":") %>%
    sapply(getElement, 2) %>% unname() %>% unique()

(common_mm <- base::intersect(opc_nearest_genes_mm, epen_nearest_genes_mm))
```

```
## [1] "Cbx8"          "Rbfox1"        "6030443J06Rik"
```

```r
opc_nearest_genes_mm <- base::setdiff(opc_nearest_genes_mm, common_mm)
epen_nearest_genes_mm <- base::setdiff(epen_nearest_genes_mm, common_mm)

# hg
(common <- base::intersect(opc_nearest_genes_hg, epen_nearest_genes_hg))
```

```
## [1] "RBFOX1" "CBX8"   ""
```

```r
opc_nearest_genes_hg <- base::setdiff(opc_nearest_genes_hg, common)
epen_nearest_genes_hg <- base::setdiff(epen_nearest_genes_hg, common)

nearest_genes_hg_uniq <- c(opc_nearest_genes_hg,
                           epen_nearest_genes_hg)

length(nearest_genes_hg_uniq)
```

```
## [1] 413
```

```r
# make a dataframe to annotate
feature_anno <- bind_rows(
    data.frame("Feature"   =  opc_nearest_genes_hg,
               "Cell_type" = "OPC"),
    data.frame("Feature"   =  epen_nearest_genes_hg,
               "Cell_type" = "Ependymal"))

# save
save(opc_nearest_genes_mm, epen_nearest_genes_mm,
     opc_nearest_genes_hg, epen_nearest_genes_hg, nearest_genes_hg_uniq,
     feature_anno,
     file = glue("{out}/cell_type_specific_features.Rda"))
```

Next, we need to obtain the corresponding human features, quantify H3K27ac/H3K27me3/scATAC
there, and then see first if the information in these features is enough to separate
tumor types.

# Evaluate samples at cell-type features

## Format inputs

Loading the quantifications of ChIPseq for H3K27ac/me3 at gene promoters for tumors,
we can use the cell-type specific genes as a reduced feature space, and perform PCA.


```r
data_chip_wide <- read_tsv(here("output/02/TABLE_promoter_H3K27ac_H3K27me3_per_sample.tsv"))
```

```
## Rows: 66970 Columns: 93
```

```
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr  (4): ID, gene_symbol, chr, strand
## dbl (88): start, end, P-1411_S-1411-1__H3K27ac, P-1425_S-1425-2__H3K27ac, P-...
## lgl  (1): Median_H3.1K27M_Thalamus_H3K27me3
```

```
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
# reformat promoter H3K27ac to a matrix
promoter_mat <- data_chip_wide %>% 
    select(-matches("Median|Z")) %>% 
    select(gene_symbol, 7:ncol(.)) %>% 
    filter(gene_symbol %in% nearest_genes_hg_uniq) %>% 
    distinct(gene_symbol, .keep_all = TRUE) %>% 
    tibble::column_to_rownames(var = "gene_symbol") %>% 
    t()

# sanity checks
dim(promoter_mat)
```

```
## [1]  69 400
```

```r
length(unique(rownames(promoter_mat)))
```

```
## [1] 69
```

```r
nearest_genes_hg_uniq[!(nearest_genes_hg_uniq %in% data_chip_wide$gene_symbol)]
```

```
##  [1] "AC058822.1"   "AL592490.1"   "AC092329.3"   "ZNF66"        "AC115220.1"  
##  [6] "AL359922.1"   "AC010615.4"   "AC005154.6"   "FP236240.1"   "CU639417.1"  
## [11] "AL117339.5"   "FAM47E-STBD1" "MDFIC2"
```

```r
feature_anno <- feature_anno %>% 
    filter(Feature %in% colnames(promoter_mat))

save(promoter_mat,
     feature_anno,
     file = glue("{out}/promoter_H3K27ac_H3K27me3_mat.Rda"))

dim(feature_anno)
```

```
## [1] 400   2
```

## H3K27ac @ cell-type features

Next, we can verify that these features, which were selected
by ChromHMM, do indeed correspond to high H3K27ac in the requisite cell type,
and low in the other.

Load quantifications of promoter H3K27ac in OPC and ependymal cells using deeptools:


```r
mm10_promoters <- rtracklayer::import(here("data/ChIPseq/references/Ensembl.ensGene.mm10.collapsed.promoter5kb.bounded.bed")) %>% 
    as.data.frame()
celltype_k27ac <- data.table::fread(here("R-4/output/03B/MBS_Ensembl.ensGene.mm10.collapsed.promoter5kb.tab"), data.table = FALSE)
colnames(celltype_k27ac) <- c("chr", "start", "end", "OPC", "Ependymal")

celltype_k27ac_anno <- celltype_k27ac %>% 
    # don't join by start, because import does start+1
    left_join(mm10_promoters, by = c("chr" = "seqnames", "end" = "end")) %>% 
    separate(name, into = c("ENSID", "symbol"), sep = ":")
```



```r
celltype_k27ac_anno_df <- celltype_k27ac_anno %>%
    filter(symbol %in% c(opc_nearest_genes_mm, epen_nearest_genes_mm)) %>% 
    mutate(Feature_type = case_when(
        symbol %in% opc_nearest_genes_mm ~ "OPC",
        symbol %in% epen_nearest_genes_mm ~ "Ependymal"
    )) %>% 
    select(symbol, Feature_type, OPC, Ependymal) %>% 
    distinct(symbol, .keep_all = TRUE)

p1 <- celltype_k27ac_anno_df %>% 
    arrange(OPC) %>% 
    mutate(symbol = factor(symbol, levels = unique(.$symbol))) %>% 
    ggplot(aes(x = symbol, y = OPC)) +
    geom_point(stat = "identity", aes(colour = Feature_type), alpha = 0.5) +
    geom_text_repel(data = celltype_k27ac_anno_df %>% top_n(30, OPC),
                    aes(label = symbol, colour = Feature_type), size = 3, max.overlaps = 30) +
    scale_colour_manual(values = palette_type) +
    rotate_x() +
    theme(axis.text.x = element_blank()) +
    ylab("Promoter H3K27ac of all cell-type specific features in OPCs") +
    no_legend()

p2 <- celltype_k27ac_anno_df %>% 
    arrange(Ependymal) %>% 
    mutate(symbol = factor(symbol, levels = unique(.$symbol))) %>% 
    ggplot(aes(x = symbol, y = Ependymal)) +
    geom_point(stat = "identity", aes(colour = Feature_type), alpha = 0.5) +
    geom_text_repel(data = celltype_k27ac_anno_df %>% top_n(30, Ependymal),
                    aes(label = symbol, colour = Feature_type), size = 3, max.overlaps = 30) +
    scale_colour_manual(values = palette_type) +
    rotate_x() +
    theme(axis.text.x = element_blank()) +
    ylab("Promoter H3K27ac of all cell-type specific features in Ependymal cells") +
    no_legend()

plot_grid(p1, p2, nrow = 1)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/03A/opc_epen_k27ac_ranked_plot-1.png)<!-- -->

```r
save(celltype_k27ac_anno_df, file = glue("{out}/celltype_k27ac_anno_df.Rda"))
```

This indicates that not _all_ the features we selected are activated by H3K27ac in the
appropriate cell type. We can select the top 15 in each:


```r
discrim_features_mm <- c(celltype_k27ac_anno_df %>% top_n(20, OPC) %>% pull(symbol),
                         celltype_k27ac_anno_df %>% top_n(20, Ependymal) %>% pull(symbol))

discrim_features_hg <- discrim_features_mm %>% mm2hg() %>% 
    unlist() %>%
    unique() %>% 
    discard(. == "") %>% 
    # make sure these features are present in tumors
    keep(. %in% feature_anno$Feature)
```


## Tumors and cell lines (HGG and PFA)

In this analysis, we'll include all HGG and PFA tumors and cell lines:


```r
samples_1 <- meta_chip %>%
    filter(Group %in% c("HGG-H3.1/2K27M-Pons",
                        "HGG-H3.3K27M-Pons",
                        "HGG-H3.3K27M-Thalamus",
                        "PFA-EZHIP-PF",
                        "PFA-H3.1K27M-PF") &
               (Material == "Tumor" | CRISPR == "Parental") &
               !grepl("Nagaraja", Source)) %>%
    filter(Factor == "H3K27ac") %>% 
    mutate(bw = gsub(".bw", "", basename(Path_bw))) %>% 
    select(ID_paper, Group, bw) %>% 
    filter(bw %in% rownames(promoter_mat))

promoter_H3K27ac_mat_1 <- promoter_mat[samples_1$bw, ]
dim(promoter_H3K27ac_mat_1)
```

```
## [1]  16 400
```



Representing these as a heatmap:


```r
hm_fun <- purrr::partial(
    pheatmap,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA,
    treeheight_row = 20,
    treeheight_col = 20,
    scale = "none",
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = samples_1 %>%
        tibble::column_to_rownames(var = "bw") %>%
        select(-ID_paper),
    annotation_row = feature_anno %>% tibble::column_to_rownames(var = "Feature"),
    annotation_colors = list("Group" = palette_groups,
                             "Cell_type" = palette_type),
    color = rdbu2,
    breaks = seq(-6, 6, length.out = 101),
    fontsize_row = 6,
    cellwidth = 10,
    cellheight = 6)
```


```r
# without H3.3K27Ms
samples_2 <- samples_1 %>% filter(Group %in% c("HGG-H3.1/2K27M-Pons",
                                               "PFA-EZHIP-PF",
                                               "PFA-H3.1K27M-PF")) %>% 
    filter(bw %in% rownames(promoter_mat))

hm_fun(mat = t(log2(promoter_H3K27ac_mat_1)[samples_2$bw, discrim_features_hg]),
       cutree_cols = 2,
       filename = glue("{figout}/celltype_H3K27ac_discrim_heatmap2.png"),
       title = "Cell-type features discriminative of H3.1K27M HGG vs EZHIP PFA")
hm_fun(mat = t(log2(promoter_H3K27ac_mat_1)[samples_2$bw, discrim_features_hg]),
       cutree_cols = 2,
       filename = glue("{figout}/celltype_H3K27ac_discrim_heatmap2.pdf"),
       title = "Cell-type features discriminative of H3.1K27M HGG vs EZHIP PFA")

knitr::include_graphics(glue("{figout}/celltype_H3K27ac_discrim_heatmap2.png"))
```

<img src="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/03A/celltype_H3K27ac_discrim_heatmap2.png" width="1620" />



# Visualization of chromatin state @ COO features

## Mouse reference data

First, set the features in the mm10 genome:


```r
opc_epen_features_mm10 <- list(
    GRanges("chr5",  37777685:37864531,   name = "Msx1"),
    GRanges("chr9",  91336804:91408070,   name = "Zic1/4"),
    GRanges("chr14", 122442135:122513576, name = "Zic2/5"),
    GRanges("chr16", 91197862:91254035,   name = "Olig2"),
    GRanges("chr2",  147137444:147241609, name = "Nkx2-2"),
    GRanges("chr15", 79149288:79189457,   name = "Sox10"))

params_mm <- map(opc_epen_features_mm10, ~ pgParams(chromstart = start(.x),
                                                    chromend = end(.x),
                                                    assembly = "mm10"))
names(params_mm) <- map(opc_epen_features_mm10, ~ .x$name)
```

Build the config file row-wise:


```r
mm10_config <- tribble(
    ~Data,      ~ID,               ~Ymax,  ~bw,
    "RNAseq",   "OPC scRNA",       151,    "data/Paired-Tag/Zhu_2021/bigwig/Paired-Tag_RNA_OPC.bw",
    "H3K27ac",  "OPC scH3K27ac",   112,    "data/Paired-Tag/Zhu_2021/bigwig/Paired-Tag_H3K27ac_OPC.bw",
    "H3K27me3", "OPC H3K27me3",    70,     "data/ChIPseq/references/Bartosovic_2021/Bulk/GSM4980053_P13556_1026_dedup.bw",
    "RNAseq",   "Epen. scRNA",     151,    "data/Paired-Tag/Zhu_2021/bigwig/Paired-Tag_RNA_Ependymal.bw",
    "H3K27ac",  "Epen. scH3K27ac", 112,    "data/Paired-Tag/Zhu_2021/bigwig/Paired-Tag_H3K27ac_Ependymal.bw"
) %>% 
    mutate(bw = here(bw))

# sanity check
all(file.exists(mm10_config$bw))
```

```
## [1] TRUE
```

```r
# import chromHMM
chromhmm_epen <- rtracklayer::import(here("data/Paired-Tag/Zhu_2021/ChromHMM/cluster22_dense.bed"))
chromhmm_opc  <- rtracklayer::import(here("data/Paired-Tag/Zhu_2021/ChromHMM/cluster15_dense.bed"))
(palette_chromhmm <- as.data.frame(chromhmm_opc) %>% select(name, itemRgb) %>% distinct() %>% deframe())
```

```
##         8         4         3         1         7         2         6         5 
## "#C8C8C8" "#00688A" "#04B1EE" "#B3EE39" "#727171" "#008B45" "#C30D19" "#956134"
```

Construct the figure:


```r
x_positions <- seq(0.5, 9, by = 1.25) + 1.5
y_positions <- seq(0.5, by = 0.4, length.out = nrow(mm10_config))
width <- 1

pageCreate(width = 11.5, height = 4, default.units = "inches")

# in this loop, for each region, we will
# 1. plot OPC bw
# 2. plot OPC ChromHMM
# 3. plot Epen. bw
# 4. plot Epen ChromHMM
# 5. plot genome labels
for (i in seq_along(params_mm)) {
    
    x_i <- x_positions[i]
    params_i <- params_mm[[i]]
    chrom <- as.character(seqnames(opc_epen_features_mm10[[i]]))
    # for the RNA & H3K27me3 bw files, chromsome names don't have the "chr" prefix...
    chroms_i <- unlist(map(mm10_config$Data,
                           ~ ifelse(.x %in% c("RNAseq", "H3K27me3"),
                                    as.numeric(stringr::str_extract(chrom, "[0-9]+")),
                                    chrom)))
    
    pwalk(list(mm10_config$bw[1:3], mm10_config$Data[1:3], mm10_config$Ymax[1:3], y_positions[1:3], chroms_i[1:3]),
          ~ pg_placeSignalAndLabel(data = ..1,
                                   color = palette_tracks[..2],
                                   range = c(0, ..3),
                                   y = ..4,
                                   x = x_i,
                                   chr = ..5,
                                   params = params_i,
                                   width = width,
                                   height = 0.35))
    
    plotRanges(chromhmm_opc, fill = chromhmm_opc$itemRgb, params = params_i, chrom = chrom, collapse = TRUE,
               x = x_i, y = "0.05b", height = 0.1, width = width)
    
    pwalk(list(mm10_config$bw[4:5], mm10_config$Data[4:5], mm10_config$Ymax[4:5], y_positions[4:5], chroms_i[4:5]),
          ~ pg_placeSignalAndLabel(data = ..1,
                                   color = palette_tracks[..2],
                                   range = c(0, ..3),
                                   y = ..4 + 0.3, # add a bit to account for chromHMM track
                                   x = x_i,
                                   chr = ..5,
                                   params = params_i,
                                   width = width,
                                   height = 0.35))
    
    plotRanges(chromhmm_epen, fill = chromhmm_epen$itemRgb, params = params_i, chrom = chrom, collapse = TRUE,
               x = x_i, y = "0.05b", height = 0.1, width = width)
    
    # place 0.1in below last plot
    plotGenomeLabel(params = params_i, chrom = chrom, scale = "Kb", x = x_i, y = "0.1b", length = width, fontsize = 8)
    
    plotGenes(params = params_i, chrom = chrom, x = x_i, y = "0.25b", height = 0.5, width = width,
              fontcolor = c("navy", "black"), fill = c("navy", "black"),
              stroke = 0.05,
              just = c("left", "top"), default.units = "inches")
    
}

# add data labels at left
pmap(list(c(mm10_config$ID[1:3], "OPC ChromHMM", mm10_config$ID[4:5], "Epen. ChromHMM"),
          c(y_positions[1:3], 1.5, y_positions[4:5] + 0.3, 2.7)),
     ~ plotText(label = ..1, fonsize = 3, fontcolor = "black",
                x = 0.2, y = ..2 + 0.2, just = c("left", "top"), default.units = "inches"))
```

```
## [[1]]
## text[text1] 
## 
## [[2]]
## text[text2] 
## 
## [[3]]
## text[text2] 
## 
## [[4]]
## text[text2] 
## 
## [[5]]
## text[text2] 
## 
## [[6]]
## text[text2] 
## 
## [[7]]
## text[text2]
```

```r
pageGuideHide()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/03A/mm10_tracks-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/03Amm10_tracks...*]~</span>

## Human tumors / cell lines

Define the regions in the human genome:


```r
opc_epen_features_hg19 <- list(
    GRanges("chr4",  4836193:4890861,     name = "MSX1"),
    GRanges("chr3",  147066046:147179311, name = "ZIC1/4"),
    GRanges("chr13", 100595157:100661032, name = "ZIC2/5"),
    GRanges("chr21", 34372339:34427397,   name = "OLIG2"),
    GRanges("chr20", 21475155:21511204,   name = "NKX2-2"),
    GRanges("chr22", 38344095:38403785,   name = "SOX10"))

params_hg <- map(opc_epen_features_hg19, ~ pgParams(
    chrom = as.character(seqnames(.x)),
    chromstart = start(.x),
    chromend = end(.x),
    assembly = "hg19"))
names(params_hg) <- map(opc_epen_features_hg19, ~ .x$name)
```


Define the input samples:


```r
hg19_config <- meta_chip %>%
    filter(Group %in% c("HGG-H3.1/2K27M-Pons", "PFA-EZHIP-PF", "PFA-H3.1K27M-PF") &
               Factor == "H3K27ac" &
               grepl("Cell-of-origin", Analyses) &
               Material == "Tumor") %>%
    arrange(Group) %>%
    distinct(ID_paper, .keep_all = TRUE) %>% 
    select(Data = Factor, ID_paper, Group, bw = Path_bw_internal) %>% 
    mutate(Ymax = case_when(
        Group == "HGG-H3.1/2K27M-Pons" ~ 10,
        grepl("PFA", Group) ~ 25))

# sanity check
all(file.exists(hg19_config$bw))
```

```
## [1] TRUE
```

Construct the figure:


```r
x_positions <- seq(0.5, 9, by = 1.25) + 1.5
y_positions <- seq(0.5, by = 0.25, length.out = nrow(hg19_config))
width <- 1

pageCreate(width = 11.5, height = 6.5, default.units = "inches")

for (i in seq_along(params_hg)) {
    
    x_i <- x_positions[i]
    params_i <- params_hg[[i]]
    
    pwalk(list(hg19_config$bw, hg19_config$Group, hg19_config$Ymax, y_positions),
          ~ pg_placeSignalAndLabel2(data = ..1,
                                    color = palette_groups[..2],
                                    range = c(0, ..3),
                                    y = ..4,
                                    x = x_i,
                                    params = params_i,
                                    width = width,
                                    height = 0.2))
    
    plotGenomeLabel(params = params_i, scale = "Kb", x = x_i, y = "0.1b", length = width, fontsize = 8)
    
    plotGenes(params = params_i, x = x_i, y = "0.25b", height = 0.75, width = width,
              fontcolor = c("navy", "black"), fill = c("navy", "black"),
              stroke = 0.05,
              just = c("left", "top"), default.units = "inches")
    
}

# add data labels at left
pmap(list(hg19_config$ID, y_positions),
     ~ plotText(label = ..1, fonsize = 3, fontcolor = "black",
                x = 0.2, y = ..2 + 0.1, just = c("left", "top"), default.units = "inches"))
```

```
## list()
```

```r
pageGuideHide()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/R-4/figures/03A/hg19_tracks-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/R-4/figures/03Ahg19_tracks...*]~</span>



<!-- END MATTER, insert reproducibility info -->


***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:

```
## 2022-06-16 14:45:55
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
## [1] parallel  stats4    stats     graphics  grDevices datasets  utils    
## [8] methods   base     
## 
## other attached packages:
##  [1] magrittr_2.0.1                           
##  [2] viridis_0.5.1                            
##  [3] viridisLite_0.3.0                        
##  [4] RColorBrewer_1.1-2                       
##  [5] cowplot_1.1.1                            
##  [6] icytobox_1.0.1                           
##  [7] fgsea_1.20.0                             
##  [8] stringr_1.4.0                            
##  [9] pheatmap_1.0.12                          
## [10] purrr_0.3.4                              
## [11] tibble_3.1.6                             
## [12] glue_1.6.1                               
## [13] readr_2.1.1                              
## [14] ggrastr_0.2.3                            
## [15] ggrepel_0.9.1                            
## [16] ggplot2_3.3.5                            
## [17] dplyr_1.0.7                              
## [18] tidyr_1.1.4                              
## [19] rtracklayer_1.54.0                       
## [20] org.Mm.eg.db_3.14.0                      
## [21] TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
## [22] org.Hs.eg.db_3.14.0                      
## [23] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2  
## [24] GenomicFeatures_1.46.4                   
## [25] AnnotationDbi_1.56.2                     
## [26] Biobase_2.54.0                           
## [27] plotgardener_1.0.14                      
## [28] GenomicRanges_1.42.0                     
## [29] GenomeInfoDb_1.26.7                      
## [30] IRanges_2.24.1                           
## [31] S4Vectors_0.28.1                         
## [32] BiocGenerics_0.36.1                      
## [33] here_1.0.1                               
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2                  reticulate_1.23            
##   [3] tidyselect_1.1.1            RSQLite_2.2.9              
##   [5] htmlwidgets_1.5.4           grid_4.1.2                 
##   [7] BiocParallel_1.24.1         Rtsne_0.15                 
##   [9] strawr_0.0.9                munsell_0.5.0              
##  [11] codetools_0.2-18            ica_1.0-2                  
##  [13] future_1.23.0               miniUI_0.1.1.1             
##  [15] withr_2.4.3                 colorspace_2.0-2           
##  [17] filelock_1.0.2              knitr_1.37                 
##  [19] Seurat_4.0.0                ROCR_1.0-11                
##  [21] tensor_1.5                  listenv_0.8.0              
##  [23] MatrixGenerics_1.6.0        git2r_0.29.0               
##  [25] GenomeInfoDbData_1.2.4      polyclip_1.10-0            
##  [27] bit64_4.0.5                 rprojroot_2.0.2            
##  [29] parallelly_1.30.0           vctrs_0.3.8                
##  [31] generics_0.1.1              xfun_0.29                  
##  [33] BiocFileCache_2.2.1         R6_2.5.1                   
##  [35] ggbeeswarm_0.6.0            spatstat.utils_2.3-0       
##  [37] bitops_1.0-7                cachem_1.0.6               
##  [39] gridGraphics_0.5-1          DelayedArray_0.16.1        
##  [41] assertthat_0.2.1            vroom_1.5.7                
##  [43] promises_1.2.0.1            BiocIO_1.4.0               
##  [45] scales_1.1.1                beeswarm_0.4.0             
##  [47] gtable_0.3.0                globals_0.14.0             
##  [49] goftest_1.2-3               rlang_0.4.12               
##  [51] splines_4.1.2               lazyeval_0.2.2             
##  [53] plyranges_1.14.0            abind_1.4-5                
##  [55] BiocManager_1.30.15         yaml_2.2.1                 
##  [57] reshape2_1.4.4              httpuv_1.6.5               
##  [59] tools_4.1.2                 ggplotify_0.1.0            
##  [61] ellipsis_0.3.2              jquerylib_0.1.4            
##  [63] ggridges_0.5.3              Rcpp_1.0.8                 
##  [65] plyr_1.8.6                  progress_1.2.2             
##  [67] zlibbioc_1.36.0             RCurl_1.98-1.5             
##  [69] prettyunits_1.1.1           deldir_1.0-6               
##  [71] rpart_4.1-15                pbapply_1.5-0              
##  [73] zoo_1.8-9                   SeuratObject_4.0.4         
##  [75] SummarizedExperiment_1.20.0 cluster_2.1.2              
##  [77] data.table_1.14.2           scattermore_0.7            
##  [79] lmtest_0.9-39               RANN_2.6.1                 
##  [81] fitdistrplus_1.1-6          matrixStats_0.61.0         
##  [83] hms_1.1.1                   patchwork_1.1.1            
##  [85] mime_0.12                   evaluate_0.14              
##  [87] xtable_1.8-4                XML_3.99-0.8               
##  [89] gridExtra_2.3               compiler_4.1.2             
##  [91] biomaRt_2.50.2              KernSmooth_2.23-20         
##  [93] crayon_1.4.2                htmltools_0.5.2            
##  [95] mgcv_1.8-38                 later_1.3.0                
##  [97] tzdb_0.2.0                  DBI_1.1.2                  
##  [99] dbplyr_2.1.1                MASS_7.3-54                
## [101] rappdirs_0.3.3              Matrix_1.3-4               
## [103] cli_3.1.1                   igraph_1.2.11              
## [105] pkgconfig_2.0.3             GenomicAlignments_1.26.0   
## [107] plotly_4.10.0               xml2_1.3.3                 
## [109] vipor_0.4.5                 bslib_0.3.1                
## [111] XVector_0.30.0              yulab.utils_0.0.4          
## [113] digest_0.6.29               sctransform_0.3.3          
## [115] RcppAnnoy_0.0.19            spatstat.data_2.1-2        
## [117] Biostrings_2.58.0           rmarkdown_2.11             
## [119] leiden_0.3.9                fastmatch_1.1-3            
## [121] uwot_0.1.11                 restfulr_0.0.13            
## [123] curl_4.3.2                  shiny_1.7.1                
## [125] Rsamtools_2.6.0             rjson_0.2.21               
## [127] nlme_3.1-153                lifecycle_1.0.1            
## [129] jsonlite_1.7.3              fansi_1.0.2                
## [131] pillar_1.6.4                lattice_0.20-45            
## [133] KEGGREST_1.34.0             fastmap_1.1.0              
## [135] httr_1.4.2                  survival_3.2-13            
## [137] spatstat_1.64-1             png_0.1-7                  
## [139] bit_4.0.4                   stringi_1.7.6              
## [141] sass_0.4.0                  blob_1.2.2                 
## [143] memoise_2.0.1               renv_0.15.5                
## [145] irlba_2.3.5                 future.apply_1.8.1
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
