---
title: "03A - Systematic analysis of HOX activation"
author: "Selin Jessa [[selin.jessa@mail.mcgill.ca](mailto:selin.jessa@mail.mcgill.ca)]"
date: "16 June, 2022"
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
## Document index: 03A
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
## public/output/03A
```

```
## public/figures/03A
```



Setting a random seed:



```r
set.seed(100)
```



***



<!-- END OF FRONT MATTER -->


# Overview

This document systematically analyzes HOX expression/activation across tumor types and data types.

* In the section _Assemble HOX coordinates_, we describe how to get the coordinates
of the HOX genes/transcripts/promoters
* In the section _Load data_, we load RNAseq and ChIPseq quantifications
at the HOX genes, and generate the associated supplementary tables
* In the following two sections, we visualize bulk RNAseq/ChIPseq quantifications
at the HOX loci
* In the section _Chromatin state_, we examine RNAseq, ChIPseq, and scATAC
data together at the HOXD locus
* In the section _3D chromatin conformation_, we examine HiC data as well
as RNAseq/ChIPseq data at the HOXD locus


# Libraries



```r
# Load libraries here
library(here)

# For genome analysis / references
library(GenomicRanges)
library(BentoBox)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
library(org.Mm.eg.db)
library(rtracklayer)

# General purpose
library(tidyr)
library(dplyr)
library(ggrepel)
library(readr)
library(readxl)
library(glue)
library(tibble)
library(ggplot2)
library(purrr)
library(pheatmap)
library(cowplot)

# colour palettes for the HOX genes are defined in the style file
source(here("include/style.R")) 
source(here("code/functions/BentoBox_helpers.R"))
ggplot2::theme_set(theme_min())
```





# Load metadata

Load the sample metadata for the project:




```r
meta      <- read_tsv(here("data/metadata/metadata_patient_samples_NGS.tsv"))
```



```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   .default = col_character(),
##   Exclude_entirely = col_logical(),
##   Age = col_double(),
##   Smartseq2 = col_logical(),
##   Smartseq2_path = col_logical(),
##   Smartseq2_ID2 = col_logical(),
##   RING1B = col_double()
## )
## ℹ Use `spec()` for the full column specifications.
```






```r
meta_chip <- read_tsv(here("data/metadata/metadata_chip_all.tsv")) %>% 
    right_join(meta, by = c("BioID", "ID_paper", "Material")) %>% 
    mutate(Factor = gsub("-M", "", Factor)) %>%
    filter(grepl("Cell-of-origin", Analyses))
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




Load the RNAseq metadata:



```r
info_samples <- data.table::fread(here("data/RNAseq/pipeline_l3/HGG-H3K27M_and_PFA.with_HOX/info.samples.tsv"), data.table = FALSE) %>% 
    mutate(Group = factor(Group, levels = names(palette_groups)))
info_groups  <- data.table::fread(here("data/RNAseq/pipeline_l3/HGG-H3K27M_and_PFA.with_HOX/info.groups.tsv"), data.table = FALSE)

table(info_samples$Group) %>% discard(. == 0)
```



```
## 
##        PFA-EZHIP-PF HGG-H3.1/2K27M-Pons   HGG-H3.3K27M-Pons  HGG-H3.3K27M-Thal. 
##                  13                   9                  11                   6
```




# Assemble HOX coordinates

This analysis involves 3 types of features:

* genes (for quantifying expression)
* transcripts (for quantifying tumor subgroup-specific expression of different HOX transcripts)
* promoters (for quantifying enrichment of H3K27ac and H3K27me3 histone marks)

In some cases, assembly of these features was done by othe scripts; where that's the case, I provide the path to those scripts from within the project repository/directory.

## Transcripts

To get the coordinates for the HOX **transcripts**, I extracted all HOX exons from the Ensembl
_exon_ annotation used in the in-house bulk RNAseq pipeline. This is performed in the script `code/scripts/make_HOX_gtf.R`.

The GTF contains exons, so we compile the full HOX transcripts from the GTF file by
taking, for each transcript, the longest interval spanning all exons associated with
that transcipt:



```r
# load Ensembl exon GTF subsetted to HOX genes
hox_exons <- rtracklayer::import(here("data/RNAseq/references/Ensembl.ensGene.exon.HOX.gtf"))

hox_exons_df <- as.data.frame(hox_exons)

# for each transcript, defined by Ensembl transcript ID,
# get the longest interval spanning all exons
hox_transcripts_df <- hox_exons_df %>%
    # for each transcript...
    group_by(ensembl_transcript) %>%
    # get longest interval
    mutate(start = min(start)) %>%
    mutate(end = max(end)) %>%
    # collapse by transcript
    group_by(seqnames, start, end, ensembl_transcript, gene_symbol, strand) %>%
    summarise() %>% 
    ungroup() %>% 
    # fix order of the transcripts based on HOX cluster (HOXA, HOXB, HOXC, HOXD)
    mutate(gene_symbol = factor(gene_symbol, levels = c(hoxa, hoxb, hoxc, hoxd))) %>% 
    filter(!is.na(gene_symbol)) %>% 
    arrange(gene_symbol, seqnames, start) %>% 
    # make a unique ID for each transcript of ENSEMBL_TX_ID:GENE_SYMBOL
    mutate(id = paste0(ensembl_transcript, ":", gene_symbol)) %>% 
    mutate(id = factor(id, levels = unique(.$id))) %>% 
    mutate(name = id) %>% 
    # get the length of the transcript
    mutate(length = abs(start - end)) %>% 
    mutate(Hox_cluster = case_when(
        grepl("HOXA", gene_symbol) ~ "HOXA",
        grepl("HOXB", gene_symbol) ~ "HOXB",
        grepl("HOXC", gene_symbol) ~ "HOXC",
        grepl("HOXD", gene_symbol) ~ "HOXD"
    ))

rr_write_tsv(hox_transcripts_df,
             glue("{out}/HOX_transcripts_coordinates.tsv"),
             desc = "Genomic coordinates of all HOX transcripts (longest interval spanning all exons) used in the analysis.")
```

<br><span style="color:#0d00ff">~[output @ *public/output/03A/HOX_transcripts_coordinates.tsv*]~</span>

## Genes

Next, to get the coordinates for the HOX **genes**, I extracted them from the Ensembl
_whole gene_ annotation used in the in-house bulk RNAseq pipeline. This is done in the following scripts:

1. `code/scripts/ENSEMBL_00-gtf_to_bed.sh` (convert whole-gene GTF to BED file)
2. `code/scripts/ENSEMBL_01-collapse_genes.R` (for each gene, defined by ENSEMBL ID, get the coordinates of the gene i.e. the longest interval spanning all exons)
3. `code/scripts/ENSEMBL_02-keep_strand.sh` (retrieve the strand information for each gene)

Here, we load the per-gene BED file produced by those scripts, and filter to HOX genes.



```r
hox_genes_bed <- rtracklayer::import(here("data/ChIPseq/references/Ensembl.ensGene.hg19.collapsed.bed")) %>% 
    as.data.frame() %>% 
    separate(name, into = c("ENSG", "gene_symbol"), sep = ":") %>% 
    filter(gene_symbol %in% hox_transcripts_df$gene_symbol)

rtracklayer::export(hox_genes_bed, glue("{out}/HOX_genes_coordinates.bed"))
```

<br><span style="color:#0d00ff">~[output @ *public/output/03A/HOX_genes_coordinates.bed*]~</span>

## Promoters

Finally, we'll define **promoters around each HOX transcript** as TSS +/- 2.5kbp,
which will be used as input for quantifications of ChIPseq data for H3K27ac and H3K27me3.
For this, I use the `GenomicRanges::promoters()` function to add flanking regions
to the TSS.



```r
# get promoters
hox_transcripts <- GRanges(hox_transcripts_df %>% dplyr::select(-length))
hox_transcripts_promoters <- GenomicRanges::promoters(hox_transcripts, upstream = 2500, downstream = 2500)

# sanity check
all(width(hox_transcripts_promoters) == 5000)
```



```
## [1] TRUE
```



```r
rtracklayer::export(hox_transcripts_promoters, glue("{out}/HOX_transcripts_promoters_5kb.bed"))
```

<br><span style="color:#0d00ff">~[output @ *public/output/03A/HOX_transcripts_promoters_5kb.bed*]~</span>


# Load data

_**NOTE:**_ RNAseq and ChIPseq quantifications were done outside this script. RNAseq
quantification was done using the bulk RNAseq pipeline, using the HOX transcripts as
features for custom feature counting. ChIPseq quantifications were performed by the Jabado
lab.

Load the RNAseq quantifications per transcript. The quantification was performed
with the in-house bulk RNAseq pipeline, configured here:
`/lustre03/project/6004736/pipeline/v0/levels1-2/sjessa/2021-04-21-HGG_HOX_transcript_counts`,
using the JSON file `.HOX_plus_Ensembl.json` in that directory



```r
load(here("data/RNAseq/pipeline_l3/HGG-H3K27M_and_PFA.with_HOX/.counts.RData"))
# get an idea of the shape of the data (transcript x sample)
counts$norm$plus_HOX_tx.exon[1:5, 1:5]
```



```
##                             EPT0508     EPT0509   EPT0562     EPT0578
## ENSG00000118473:SGIP1     114.01181    70.75683  272.8021   310.77025
## ENSG00000162426:SLC45A1    78.00808    99.89199  156.9547    77.69256
## ENSG00000157191:NECAP2   2361.24463  1846.61449  855.7766  1533.40582
## ENSG00000169504:CLIC4   28657.46896 13504.84255 8810.0145 19267.75526
## ENSG00000142920:ADC        99.01026   169.26143  151.3491   284.19121
##                             EPT0805
## ENSG00000118473:SGIP1     107.35236
## ENSG00000162426:SLC45A1    80.08827
## ENSG00000157191:NECAP2   1318.90042
## ENSG00000169504:CLIC4   10936.30869
## ENSG00000142920:ADC       294.79299
```



```r
# convert to long format, with columns Sample, Transcript, Expression
# using DESeq2-normalized expression data
counts_rna_long <- counts$norm$plus_HOX_tx.exon %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "id") %>% 
    filter(grepl(":HOX", id)) %>% 
    gather("Sample", "Norm_expression", 2:ncol(.)) %>% 
    separate(id, into = c("ensembl_transcript", "gene_symbol"), sep = ":") %>% 
    # join with sample metadata to get tumor group
    left_join(info_samples, by = c("Sample" = "Nickname")) %>% 
    # join with HOX transcripts annotation to get HOX cluster
    right_join(hox_transcripts_df)
```



```
## Joining, by = c("ensembl_transcript", "gene_symbol")
```



```r
# calculate median expression per group
counts_rna_long_medians <- counts_rna_long %>% 
    group_by(Group, gene_symbol, id, Hox_cluster) %>% 
    summarize(Median = median(Norm_expression)) %>% 
    ungroup()
```



```
## `summarise()` has grouped output by 'Group', 'gene_symbol', 'id'. You can override using the `.groups` argument.
```



```r
# make supplementary table
counts_rna_long %>%
    dplyr::select(ensembl_transcript, gene_symbol, chr = seqnames, start, end, strand, length, Sample, Norm_expression, Hox_cluster) %>%
    pivot_wider(names_from = Sample, values_from = Norm_expression) %>% 
    arrange(Hox_cluster, start) %>%
    rr_write_tsv(glue("{out}/TABLE_HOX_expression_per_transcript.tsv"),
                 "Normalized (DESeq2) expression of each HOX transcript in each sample")
```



```
## ...writing description of TABLE_HOX_expression_per_transcript.tsv to public/output/03A/TABLE_HOX_expression_per_transcript.desc
```




Load the ChIPseq quantifications per transcript. ChIPseq data for histone marks
H3K27ac and H3K27me3 was quantified at the HOX promoters in a per-transcript manner, using
the TSS +/- 2.5kbp.



```r
# the ChIPseq quantifications has two parts:
# 1. the counts (saved into the variable counts_chip)
# 2. the column annotations for the counts (saved into the variable meta_chipq)
counts_chip <- read_xlsx(here("data/ChIPseq/quantifications_@mhulswit/20210511_H3.3H3.1EZHIP_H3K27ac_H3K27me3_Ensembl.ensGene.exon.Hox.promoters_5kb_SJ.xlsx"), sheet = 1)
meta_chipq  <- read_xlsx(here("data/ChIPseq/quantifications_@mhulswit/20210511_H3.3H3.1EZHIP_H3K27ac_H3K27me3_Ensembl.ensGene.exon.Hox.promoters_5kb_SJ.xlsx"), sheet = 2) %>% 
    mutate(Sample = gsub("_Count|_RPKM", "", Sample)) %>% 
    distinct(Sample, Group, Source) %>% 
    # make the group annotation consistent with the rest of the analysis/code
    mutate(Group = case_when(
        Group == "EPN-PFA_Posterior fossa" ~ "PFA-EZHIP-PF",
        Group == "H3.1K27M_Pons"           ~ "HGG-H3.1/2K27M-Pons",
        Group == "H3.3K27M_Pons"           ~ "HGG-H3.3K27M-Pons",
        Group == "H3.3K27M_Thalamus"       ~ "HGG-H3.3K27M-Thal.",
        TRUE ~ "Drop")) %>% 
    # join with the project ChIPseq metadata table to get data paths
    left_join(meta_chip %>%
                  mutate(bw = basename(Path_bw_internal)) %>%
                  dplyr::select(ID_paper, bw, Replicate),
              by = c("Sample" = "bw"))

# convert to long format
counts_chip_long <- counts_chip %>%
    dplyr::select(-chr, -strand, -name) %>%
    dplyr::rename(id = ID) %>% 
    gather(Stat, Value, 4:ncol(.)) %>% 
    mutate(
        Statistic = case_when(
            grepl("Count", Stat) ~ "Count",
            grepl("RPKM",  Stat) ~ "RPKM" ),
        Sample = gsub("_Count|_RPKM", "", Stat),
        Mark = ifelse(grepl("K27me3", Sample), "H3K27me3", "H3K27ac")
    ) %>% 
    left_join(meta_chipq, by = "Sample") %>% 
    # exlcude tumor groups not included in the study
    filter(Group != "Drop") %>% 
    right_join(hox_transcripts_df %>% dplyr::select(-start, -end, -length), by = "id")

counts_chip_long %>% distinct(Sample, Mark, Group) %>% count(Mark, Group)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Mark"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Group"],"name":[2],"type":["chr"],"align":["left"]},{"label":["n"],"name":[3],"type":["int"],"align":["right"]}],"data":[{"1":"H3K27ac","2":"HGG-H3.1/2K27M-Pons","3":"16"},{"1":"H3K27ac","2":"HGG-H3.3K27M-Pons","3":"16"},{"1":"H3K27ac","2":"HGG-H3.3K27M-Thal.","3":"5"},{"1":"H3K27ac","2":"PFA-EZHIP-PF","3":"4"},{"1":"H3K27me3","2":"HGG-H3.1/2K27M-Pons","3":"3"},{"1":"H3K27me3","2":"HGG-H3.3K27M-Pons","3":"6"},{"1":"H3K27me3","2":"HGG-H3.3K27M-Thal.","3":"4"},{"1":"H3K27me3","2":"PFA-EZHIP-PF","3":"2"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# calculate median promoter ChIP RPKM per group
counts_chip_long_medians <- counts_chip_long %>% 
    filter(Statistic == "RPKM") %>% 
    group_by(Mark, Group, id, Hox_cluster) %>% 
    summarize(Median_RPKM = median(Value)) %>% 
    ungroup()
```



```
## `summarise()` has grouped output by 'Mark', 'Group', 'id'. You can override using the `.groups` argument.
```



```r
# make supplementary table
counts_chip_long %>%
    filter(Statistic == "RPKM") %>% 
    mutate(Sample = paste0(ID_paper, "__", Replicate, "__", Mark)) %>% 
    mutate(length = end - start) %>% 
    dplyr::select(id, chr = seqnames, start, end, length, strand, Sample, Value, Hox_cluster) %>%
    separate(id, into = c("ensembl_transcript", "gene_symbol"), sep = ":") %>% 
    pivot_wider(names_from = Sample, values_from = Value) %>% 
    arrange(Hox_cluster, start) %>%
    rr_write_tsv(glue("{out}/TABLE_HOX_H3K27ac_H3K27me3_per_transcript.tsv"),
                 "Promoter (TSS +/- 5kbp) H3K27ac and H3K27me3 enrichment (RPKM) of each HOX transcript in each sample")
```



```
## ...writing description of TABLE_HOX_H3K27ac_H3K27me3_per_transcript.tsv to public/output/03A/TABLE_HOX_H3K27ac_H3K27me3_per_transcript.desc
```




# Gene expression by bulk RNAseq


For a per-gene analysis, we'll select, for each gene, **the most highly-expressed 
transcript, followed by the longest if there are ties**.



```r
# this variable will store the most highly expressed transcript for each HOX gene,
# which will be used for any downstream per-gene analyses
keep_transcripts <- counts_rna_long %>%
    group_by(id, gene_symbol, length) %>%
    summarize(mean_norm_expression = mean(Norm_expression)) %>% 
    ungroup() %>% 
    group_by(gene_symbol) %>% 
    # sort by descending expression
    dplyr::arrange(desc(mean_norm_expression), length) %>% 
    # get the most highly expressed transcript
    dplyr::slice(1) %>% 
    pull(id)
```



```
## `summarise()` has grouped output by 'id', 'gene_symbol'. You can override using the `.groups` argument.
```



```r
# sanity check: we should have one selected transcript per gene, 39 in total (as many as there are HOX genes)
length(keep_transcripts)
```



```
## [1] 39
```



```r
length(keep_transcripts) == length(unique(hox_transcripts_df$gene_symbol))
```



```
## [1] TRUE
```



```r
# make a similar plot as above, but per gene
counts_rna_long %>% 
    filter(id %in% keep_transcripts) %>%
    mutate(gene_symbol = factor(gene_symbol, levels = names(palette_hox))) %>% 
    dplyr::select(id, Norm_expression, gene_symbol, Group, Hox_cluster) %>% 
    rr_ggplot(plot_num = 1, aes(x = gene_symbol, y = Norm_expression)) +
    geom_boxplot(aes(fill = gene_symbol), alpha = 0.5, width = 0.5, outlier.color = NA) +
    geom_jitter(aes(color = gene_symbol), alpha = 0.5, width = 0.1, size = 0.8, stroke = 0.5) +
    facet_grid(Group ~ Hox_cluster, scales = "free_x") +
    scale_fill_manual(values = palette_hox) +
    scale_color_manual(values = palette_hox) +
    theme(strip.text.y = element_text(size = rel(0.8)),
          strip.text.x = element_blank()) +
    rotate_x() +
    no_legend() 
```



```
## ...writing source data of ggplot to public/figures/03A/hox_expression_boxplots_per_gene-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03A/hox_expression_boxplots_per_gene-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/03A/hox_expression_boxplots_per_gene...*]~</span>


## Highlight HOX genes in DGE genes

We highlight the HOX genes among the volcano plot of differentially-expressed
genes between H3.3K27M thalamus and pons, to demonstrate that they are among
the most significantly upregulated in pons tumors. The differential expression
analysis was done with the in-house bulk RNAseq pipeline using the expression-tools
module.



```r
# load differentially-expressed genes from DESeq2 analysis from in-house pipeline
dge_thalamus_vs_pons <- read_tsv(here("data/RNAseq/pipeline_l3/DGE/HGG-H3.3K27M-Thal._vs_HGG-H3.3K27M-Pons/diff/Ensembl.ensGene.exon/HGG-H3.3K27M-PonsvsHGG-H3.3K27M-Thal..tsv")) %>% 
    separate(ID, into = c("ENSID", "symbol"), sep = ":") %>% 
    # note the HOX genes
    mutate(HOX = ifelse(symbol %in% hox_genes_bed$gene_symbol, TRUE, FALSE)) %>% 
    # put HOX genes (points) on top
    arrange(HOX) %>% 
    # filter to expressed genes
    filter(baseMean > 100)
```



```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   ID = col_character(),
##   baseMean = col_double(),
##   log2FoldChange = col_double(),
##   lfcSE = col_double(),
##   stat = col_double(),
##   pvalue = col_double(),
##   padj = col_double()
## )
```



```r
# make volcano plot
dge_thalamus_vs_pons %>% 
    rr_ggplot(aes(x = log2FoldChange, y = -log10(padj)), plot_num = 1) +
    geom_point(alpha = 0.5, aes(colour = HOX, size = HOX)) + 
    scale_colour_manual(values = c("TRUE" = "darkorchid4", "FALSE" = "gray70")) +
    scale_size_manual(values = c("TRUE" = 2, "FALSE" = 1)) +
    geom_text_repel(data = dge_thalamus_vs_pons %>%
                        filter(HOX & baseMean > 100 & abs(log2FoldChange) > 2),
                    inherit.aes = TRUE, aes(label = symbol), colour = "darkorchid4") +
    ggtitle("H3.3K27M thalamus vs. H3.3K27M pons")
```



```
## ...writing source data of ggplot to public/figures/03A/dge_volcano_hox-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03A/dge_volcano_hox-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/03A/dge_volcano_hox...*]~</span>

# Heatmap quantification {.tabset}

To complement the RNAseq, we'll also examine the quantification of HOX expression,
again per gene and per transcript, in both ChIPseq data for H3K27ac and H3K27me3.

To put all data types in the same range, for each data type, **we scale the values
for all HOX genes across samples to [0, 1]**. This will mean that **any values for the
same data type are comparable across HOX clusters/genes and across samples**. However,
values between data types are not comparable.

For the per-gene analysis, we first filter to the most highly expressed transcript
per gene, and then do the scaling. For the per-transcript analysis, we scale on
the per-transcript values.



```r
# per gene
rna_per_gene_scaled <- counts_rna_long_medians %>% 
    # filter to most highly expressed transcript first
    filter(id %in% keep_transcripts) %>% 
    mutate(Median_scaled = scales::rescale(Median))

k27ac_per_gene_scaled <- counts_chip_long_medians %>%
    filter(Mark == "H3K27ac") %>% 
    filter(id %in% keep_transcripts) %>% 
    mutate(Median_scaled = scales::rescale(Median_RPKM))

k27me3_per_gene_scaled <- counts_chip_long_medians %>%
    filter(Mark == "H3K27me3") %>% 
    filter(id %in% keep_transcripts) %>% 
    mutate(Median_scaled = scales::rescale(Median_RPKM))

# per tx
rna_per_tx_scaled <- counts_rna_long_medians %>% 
    mutate(Median_scaled = scales::rescale(Median))

k27ac_per_tx_scaled <- counts_chip_long_medians %>%
    filter(Mark == "H3K27ac") %>% 
    mutate(Median_scaled = scales::rescale(Median_RPKM))

k27me3_per_tx_scaled <- counts_chip_long_medians %>%
    filter(Mark == "H3K27me3") %>% 
    mutate(Median_scaled = scales::rescale(Median_RPKM))
```



Function for generating the heatmaps across data types:



```r
# for each unique group, add a suffix for each data type
group_order <- paste0(
    rep(names(palette_groups)[names(palette_groups) %in% unique(rna_per_gene_scaled$Group)], each = 3),
    c(":RNA", ":H3K27ac", ":H3K27me3")
)

# define order of genes and transcripts
gene_order <- hox_transcripts_df %>% filter(id %in% keep_transcripts) %>% pull(id) %>% 
    as.character()

tx_order <- levels(hox_transcripts_df$id)

#' Plot heatmap of HOX activation based on RNA/H3K27ac/H3K27me3
#'
#' @param hox_cluster Character, one of "HOXA", "HOXB", "HOXC", "HOXD"
#' @param feature_order Character, vector of genes/transcripts in the order they
#' should be plot along the columns of the heatmap
#' @param type Character, one of "per_gene" or "per_transcript", indicating whether
#' to generate heatmaps by genes or by transcripts
#'
#' @examples
#' plot_hox_heatmap("HOXA", gene_order, "per_gene")
plot_hox_heatmap <- function(hox_cluster, feature_order, type) {
    
    # extract the data for the genes in this cluster,
    # for either the per_gene quantification or per_transcript quantification
    if (type == "per_gene") {
        
        d1 <- rna_per_gene_scaled    %>% filter(Hox_cluster == hox_cluster) %>% 
            dplyr::select(id, Group, RNA      = Median_scaled)
        d2 <- k27ac_per_gene_scaled  %>% filter(Hox_cluster == hox_cluster) %>%
            dplyr::select(id, Group, H3K27ac  = Median_scaled)
        d3 <- k27me3_per_gene_scaled %>% filter(Hox_cluster == hox_cluster) %>%
            dplyr::select(id, Group, H3K27me3 = Median_scaled)
        
    } else if (type == "per_tx") {
        
        d1 <- rna_per_tx_scaled    %>% filter(Hox_cluster == hox_cluster) %>% 
            dplyr::select(id, Group, RNA      = Median_scaled)
        d2 <- k27ac_per_tx_scaled  %>% filter(Hox_cluster == hox_cluster) %>%
            dplyr::select(id, Group, H3K27ac  = Median_scaled)
        d3 <- k27me3_per_tx_scaled %>% filter(Hox_cluster == hox_cluster) %>%
            dplyr::select(id, Group, H3K27me3 = Median_scaled)
        
    }
    
    # merge the RNA and ChIPseq data
    d4 <- d1 %>% left_join(d2) %>% left_join(d3) %>%
        ungroup() %>% 
        gather(stat, value, 3:5) %>% 
        mutate(stat = paste0(Group, ":", stat)) %>% 
        dplyr::select(-Group) %>% 
        spread(stat, value) %>% 
        tibble::column_to_rownames(var = "id")
    
    # plot in the order provided
    input <- d4 %>% 
        t() %>% 
        .[group_order, feature_order[feature_order %in% colnames(.)] ]
    
    # function to generate the heatmap -- these parameters will be kept constant
    hm_fun <- purrr::partial(pheatmap,
                             input,
                             scale = "none",
                             cluster_rows = FALSE,
                             cluster_cols = FALSE,
                             border_color = "black",
                             color = purples,
                             cellwidth = 50, cellheight = 25,
                             breaks = seq(0, 1, length.out = 101))
    
    # make both file types
    hm_fun(filename = glue("{figout}/heatmap_{hox_cluster}_{type}.pdf"))
    hm_fun(filename = glue("{figout}/heatmap_{hox_cluster}_{type}.png"))
    
    # show the PNG in the HTML
    knitr::include_graphics(glue("{figout}/heatmap_{hox_cluster}_{type}.png"))
    
}
```




Now make a heatmap for each cluster per gene. The rows of the heatmaps are re-ordered
using Illustrator for the final figure.



```r
plot_hox_heatmap("HOXA", gene_order, "per_gene")
```

<img src="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03A/heatmap_HOXA_per_gene.png" width="3185" />

```r
plot_hox_heatmap("HOXB", gene_order, "per_gene")
```

<img src="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03A/heatmap_HOXB_per_gene.png" width="2977" />

```r
plot_hox_heatmap("HOXC", gene_order, "per_gene")
```

<img src="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03A/heatmap_HOXC_per_gene.png" width="2769" />

```r
plot_hox_heatmap("HOXD", gene_order, "per_gene")
```

<img src="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03A/heatmap_HOXD_per_gene.png" width="2769" />


# Chromatin state {.tabset}

Next, we show that the tumors exhibit bipartite functional domains at the HOX
clusters based on chromatin state data. We'll assemble representative data tracks
for each tumor group and for normal cells at the HOXD cluster.

For this section, I heavily use helper functions defined in `code/functions/BentoBox_helpers.R`.
To configure the plots, I define the samples to include in configuration TSV
files at `data/BentoBox_config/hoxd_*.config.tsv`

First, define the region and the BentoBox genomic plotting parameters:



```r
# define region
hoxd <- GRanges("chr2", 176927000:177075000)
params_hoxd <- bb_params(chrom = as.character(seqnames(hoxd)),
                         chromstart = start(hoxd),
                         chromend = end(hoxd),
                         assembly = "hg19")

# set up gene track & genomic ranges object
hoxd_genes <- hox_transcripts_df %>% filter(Hox_cluster == "HOXD") %>%
  pull(gene_symbol) %>%
  unique() %>%
  as.character()
hoxd_genes_gr <- hox_genes_bed %>% filter(gene_symbol %in% hoxd_genes) %>%
    dplyr::rename(name = gene_symbol) %>%
    mutate(name = factor(name, levels = hoxd_genes)) %>%
    arrange(name) %>%
    GRanges()

hoxd_genes_gr
```



```
## GRanges object with 9 ranges and 3 metadata columns:
##       seqnames              ranges strand |            ENSG     name     score
##          <Rle>           <IRanges>  <Rle> |     <character> <factor> <numeric>
##   [1]     chr2 177053307-177055688      + | ENSG00000128645    HOXD1         0
##   [2]     chr2 177001340-177037830      + | ENSG00000128652    HOXD3         0
##   [3]     chr2 177015950-177017954      + | ENSG00000170166    HOXD4         0
##   [4]     chr2 176994422-176997423      + | ENSG00000175879    HOXD8         0
##   [5]     chr2 176987088-176989853      + | ENSG00000128709    HOXD9         0
##   [6]     chr2 176973518-176984670      + | ENSG00000128710   HOXD10         0
##   [7]     chr2 176968944-176974722      + | ENSG00000128713   HOXD11         0
##   [8]     chr2 176964458-176966408      + | ENSG00000170178   HOXD12         0
##   [9]     chr2 176957619-176960666      + | ENSG00000128714   HOXD13         0
##   -------
##   seqinfo: 68 sequences from an unspecified genome; no seqlengths
```




## H3.3K27M pons



```r
# define tracks (function defined in code/functions/BentoBox_helpers.R)
hoxd_config <- prep_track_input(here("data/BentoBox_config/hoxd_H33.config.tsv"))

# make figure
x <- 0.5
y_positions <- seq(0.5, by = 0.5, length.out = nrow(hoxd_config))

bb_pageCreate(width = 4, height = 5, default.units = "inches")

# add bw tracks
pwalk(list(hoxd_config$bw, hoxd_config$ID_paper, hoxd_config$Data, y_positions),
      # helper function defined in code/functions/BentoBox_helpers.R
      ~ bb_placeSignalAndLabel(data = ..1,
                               annotation = ..2,
                               color = palette_tracks[..3],
                               y = ..4, x = x,
                               params = params_hoxd,
                               ymax = 1,
                               fontsize = 4))

bb_plotGenomeLabel(params = params_hoxd, scale = "Kb", x = 0.5, y = "0.03b", length = 3)

# plot ticks for the genes in colour
# bb_plotBed(params = params_hoxd, data = hoxd_genes_gr, fill = palette_hox[hoxd_genes], colorby = colorby("name"),
#            width = 3, height = 0.6, y = 3.5, x = 0.5,
#            baseline.lwd = 2, boxHeight = unit(3, "mm"), collapse = FALSE)

bb_plotGenes(params = params_hoxd, x = 0.5, y = "0.1b", height = 1, width = 3,
             strandcolors = c("navy", "black"), fontcolors = c("navy", "black"),
             stroke = 0.05,
             just = c("left", "top"), default.units = "inches")

bb_pageGuideHide()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03A/hoxd_tracks_H3.3-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/03A/hoxd_tracks_H3.3...*]~</span>



## H3.1K27M pons



```r
# define tracks
hoxd_config <- prep_track_input(here("data/BentoBox_config/hoxd_H31.config.tsv"))

# make figure
x <- 0.5
y_positions <- seq(0.5, by = 0.5, length.out = nrow(hoxd_config))

bb_pageCreate(width = 4, height = 5.5, default.units = "inches")

# add bw tracks
pwalk(list(hoxd_config$bw, hoxd_config$ID_paper, hoxd_config$Data, y_positions),
      ~ bb_placeSignalAndLabel(data = ..1,
                               annotation = ..2,
                               color = palette_tracks[..3],
                               y = ..4, x = x,
                               params = params_hoxd,
                               ymax = 1,
                               fontsize = 4))

bb_plotGenomeLabel(params = params_hoxd, scale = "Kb", x = 0.5, y = 3.53, length = 3)

bb_plotGenes(params = params_hoxd, x = 0.5, y = "0.1b", height = 1, width = 3,
             strandcolors = c("navy", "black"), fontcolors = c("navy", "black"),
             stroke = 0.05,
             just = c("left", "top"), default.units = "inches")

bb_pageGuideHide()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03A/hoxd_tracks_H3.1-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/03A/hoxd_tracks_H3.1...*]~</span>


## PFA



```r
# define tracks
hoxd_config <- prep_track_input(here("data/BentoBox_config/hoxd_PFA.config.tsv"))

# make figure
x <- 0.5
y_positions <- seq(0.5, by = 0.5, length.out = nrow(hoxd_config))

bb_pageCreate(width = 4, height = 4, default.units = "inches")

# add bw tracks
pwalk(list(hoxd_config$bw, hoxd_config$ID_paper, hoxd_config$Data, y_positions),
      ~ bb_placeSignalAndLabel(data = ..1,
                               annotation = ..2,
                               color = palette_tracks[..3],
                               y = ..4, x = x,
                               params = params_hoxd,
                               ymax = 1,
                               fontsize = 4))

bb_plotGenomeLabel(params = params_hoxd, scale = "Kb", x = 0.5, y = "0.03b", length = 3)

bb_plotGenes(params = params_hoxd, x = 0.5, y = "0.1b", height = 1, width = 3,
             strandcolors = c("navy", "black"), fontcolors = c("navy", "black"),
             stroke = 0.05,
             just = c("left", "top"), default.units = "inches")

bb_pageGuideHide()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03A/hoxd_tracks_PFA-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/03A/hoxd_tracks_PFA...*]~</span>



## Motor neurons & ESCs

This is data from [Narendra et al, Science, 2015](https://pubmed.ncbi.nlm.nih.gov/25722416/),
downloaded from GEO ID [GSE60232](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60232).



```r
# set region in mm9
hoxd_mm9 <- GRanges("chr2", 74480000:74610000)
params_hoxd_mm9 <- bb_params(chrom = as.character(seqnames(hoxd_mm9)),
                             chromstart = start(hoxd_mm9),
                             chromend = end(hoxd_mm9),
                             assembly = "mm9")

# build the config data frame row-wise
hoxd_config <- tribble(
    ~Data,       ~ID,              ~Ymax, ~bw,
    "H3K27me3",  "H3K27me3 (ESC)", 8,     "GSM1468396_H3K27me3.ESC.WT.bw",
    "Pol2",      "Pol2 (MN)",      20,    "GSM1468406_Pol2.MN.WT.bw",
    "H3K4me3",   "H3K4me3 (MN)",   150,   "GSM1468401_H3K4me3.MN.WT.bw",
    "H3K27me3",  "H3K27me3 (MN)",  120,   "GSM1468398_H3K27me3.MN.WT.bw",
    "EZH2",      "EZH2 (MN)",      150,   "GSM1468404_EZH2.MN.WT.bw",
    "H2AK119Ub", "H2AK119Ub (MN)", 20,    "GSM1468408_H2AK119Ub.MN.WT.bw",
    "CTCF",      "CTCF (MN)",      30,    "GSM1468394_CTCF.MN.WT.bw"
) %>%
    mutate(bw = here(file.path("data/ChIPseq/bigwig/Narendra_et_al_2015/", bw)))

# make figure
x <- 0.5
y_positions <- seq(0.5, by = 0.5, length.out = nrow(hoxd_config))

bb_pageCreate(width = 4, height = 6, default.units = "inches")

# add bw tracks
pwalk(list(hoxd_config$bw, hoxd_config$ID, hoxd_config$Data, y_positions, hoxd_config$Ymax),
      ~ bb_placeSignalAndLabel(data = ..1,
                               annotation = ..2,
                               color = palette_tracks[..3],
                               y = ..4, x = x,
                               range = c(0, ..5),
                               params = params_hoxd_mm9,
                               fontsize = 4))

bb_plotGenomeLabel(params = params_hoxd_mm9, scale = "Kb", x = 0.5, y = "0.03b", length = 3)

bb_plotGenes(params = params_hoxd_mm9, x = 0.5, y = "0.1b", height = 1, width = 3,
             strandcolors = c("navy", "black"), fontcolors = c("navy", "black"),
             stroke = 0.05,
             just = c("left", "top"), default.units = "inches")

bb_pageGuideHide()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03A/hoxd_tracks_MN_Narendra-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/03A/hoxd_tracks_MN_Narendra...*]~</span>




# 3D chromatin conformation

Finally, we interrogate 3D chromatin conformation at the HOX locus using HiC data,
and plot this together with other genomic data. 

First, we'll assemble the plotting parameters:



```r
# define gene palette and region
palette_hoxd <- palette_hox[grepl("HOXD", names(palette_hox))]
hoxd_genes_gr <- GRanges(hox_genes_bed[grepl("HOXD", hox_genes_bed$gene_symbol),])
region <- GRanges("chr2", 176725000:177320000)

params_hoxd_hic <- bb_params(chrom = "chr2",
                         chromstart = start(region),
                         chromend = end(region),
                         assembly = "hg19",
                         zrange = c(-4, 4),
                         palette = colorRampPalette(colors = c("blue", "white", "red")))

# get x/y-positions; heatmaps can be placed 2 inches apart starting at y = 0.5
y_positions <- 0.5
x_positions <- c(2, 5.5, 9)
x_starts    <- c(0.5, 4, 7.5)
```



Next, assemble the inputs:



```r
# define input data ------------------------------------------------------------
hic_files <- file.path(here("data/HiC/hic/"), c("tumor/E823_hg19.hic",
                                                "cell/HSJ-031_hg19.hic",
                                                "cell/BT245_hg19.hic"))

# load normalized HiC counts for this region, for all samples
hic_counts <- map(hic_files, function(hic) {
    
    out <- bb_readHic(file = hic,
                      chrom = "chr2",
                      chromstart = start(region),
                      chromend = end(region),
                      assembly = "hg19",
                      norm = "KR",
                      matrix = "oe",
                      resolution = 10000)
    
    # get log2(O/E)
    out$counts <- log2(out$counts)
    
    return(out)
    
})
```



```
## Read in hic file with KR normalization at 10000 BP resolution.
## Read in hic file with KR normalization at 10000 BP resolution.
## Read in hic file with KR normalization at 10000 BP resolution.
```



```r
# HiC
inputs <- list(data = hic_counts,
               annotation = c("E823 (PFA)", "HSJ-031 \n(H3.3 pons)", "BT245 \n(H3.3 thal)"),
               x = x_positions)

# add ChIP bw tracks
inputs_ctcf <- list("ctcf" = c(here("data/ChIPseq/bigwig/cell_parental/CTCF/DIPGXIII-1__CTCF.bw"),
                               here("data/ChIPseq/bigwig/cell_parental/CTCF/BT245-1__CTCF.bw")),
                    "annotation" = c("DIPGXIII", "BT245"),
                    x = x_starts[2:3])

inputs_suz12 <- list("suz12" = c(here("data/ChIPseq/bigwig/cell_parental/SUZ12/DIPGXIII-1__SUZ12.bw"),
                                 here("data/ChIPseq/bigwig/cell_parental/SUZ12/BT245-1__SUZ12.bw")),
                     "annotation" = c("DIPGXIII", "BT245"),
                     x = x_starts[2:3])

# helper function to get the internal pipleine path for a sample without exposing ID
get_rna_path <- function(id, material = "Tumor") {
    
    RNAseq_path <- meta %>% filter(ID_paper == id & Material == material) %>% pull(RNAseq_path)
    file.path("/lustre06/project/6004736/pipeline/v0/levels1-2", RNAseq_path,
              "star", paste0(basename(RNAseq_path), ".sorted.bw"))
    
}

inputs_rna <- list("rna" = c(get_rna_path("P-5426_S-6887"),
                             get_rna_path("P-3407_S-3447"),
                             get_rna_path("P-4111_S-4496")),
                   "annotation" = c("P-5426_S-6887", "P-3407_S-3447", "P-4111_S-4496"),
                   x = x_starts)

inputs_k27ac <- list("k27ac" = c(here("data/ChIPseq/bigwig/tumor/H3K27ac/P-5430_S-6891-1__H3K27ac.bw"),
                                 here("data/ChIPseq/bigwig/cell_parental/H3K27ac/DIPGXIII-1__H3K27ac.bw"),
                                 here("data/ChIPseq/bigwig/cell_parental/H3K27ac/BT245-1__H3K27ac.bw")),
                     "annotation" = c("P-5430_S-6891", "DIPGXIII", "BT245"),
                     x = x_starts)

inputs_k27me3 <- list("k27me3" = c(here("data/ChIPseq/bigwig/tumor/H3K27me3/P-5430_S-6891-1__H3K27me3.bw"),
                                   here("data/ChIPseq/bigwig/cell_parental/H3K27me3/DIPGXIII-1__H3K27me3.bw"),
                                   here("data/ChIPseq/bigwig/cell_parental/H3K27me3/BT245-1__H3K27me3.bw")),
                      "annotation" = c("P-5430_S-6891", "DIPGXIII", "BT245"),
                      x = x_starts)
```



Assemble the plot:



```r
# make figure ------------------------------------------------------------------

# create BB page
bb_pageCreate(width = 11.5, height = 6, default.units = "inches")

# place the heatmap & annotations for each sample, 1 x 3
pwalk(inputs,
      # helper function defined in code/functions/BentoBox_helpers.R
      ~ bb_placeHeatmapAndLegend(counts = ..1,
                                 annotation = ..2,
                                 y = 0.5, x = ..3,
                                 params = params_hoxd_hic))

# add ChIP using helper function defined in code/functions/BentoBox_helpers.R
pwalk(inputs_rna,    ~ bb_placeSignalAndLabel(..1, y = 2.05, x = ..3, annotation = ..2,
                                              params = params_hoxd_hic, ymax = 1, color = "black"))
pwalk(inputs_k27ac,  ~ bb_placeSignalAndLabel(..1, y = 2.55, x = ..3, annotation = ..2,
                                              params = params_hoxd_hic, ymax = 1, color = "blue"))
pwalk(inputs_k27me3, ~ bb_placeSignalAndLabel(..1, y = 3.05, x = ..3, annotation = ..2,
                                              params = params_hoxd_hic, ymax = 1, color = "red"))
pwalk(inputs_suz12,  ~ bb_placeSignalAndLabel(..1, y = 3.55, x = ..3, annotation = ..2,
                                              params = params_hoxd_hic, ymax = 1, color = "red"))
pwalk(inputs_ctcf,   ~ bb_placeSignalAndLabel(..1, y = 4.05, x = ..3, annotation = ..2,
                                              params = params_hoxd_hic, ymax = 1, color = "black"))

# annotate data types
bb_plotText("RNA",      x = 10.55, y = 2.25, just = c("left", "top"), fontcolor = "black")
bb_plotText("H3K27ac",  x = 10.55, y = 2.75, just = c("left", "top"), fontcolor = "blue")
bb_plotText("H3K27me3", x = 10.55, y = 3.25, just = c("left", "top"), fontcolor = "red")
bb_plotText("SUZ12",    x = 10.55, y = 3.75, just = c("left", "top"), fontcolor = "red")
bb_plotText("CTCF",     x = 10.55, y = 4.25, just = c("left", "top"), fontcolor = "black")

# annotate x-axis genome label
walk(x_starts, ~ bb_plotGenomeLabel(params = params_hoxd_hic,
                                    scale = "Kb", x = .x, y = 4.55, length = 3))

# plot ticks for the genes in colour
walk(x_starts, ~ bb_plotBed(params = params_hoxd_hic, data = hoxd_genes_gr,
                            fill = palette_hoxd, colorby = colorby("gene_symbol"),
                            width = 3, height = 0.6, y = 4.75, x = .x,
                            baseline.lwd = 2, boxHeight = unit(3, "mm"), collapse = FALSE))

# plot the gene names in colour
x_starts_ticks <- seq(1.6, by = 0.1, length.out = 9)
x_starts_ticks <- list(x_starts_ticks, x_starts_ticks + 3.5, x_starts_ticks + 7)
walk(x_starts_ticks, ~ bb_plotText(label = rev(names(palette_hoxd)),
                                   fontcolor = rev(unname(palette_hoxd)),
                                   rot = 90, fontsize = 6,
                                   x = .x, y = 5.5,
                                   just = c("right"), default.units = "inches"))

# done!
bb_pageGuideHide()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03A/hoxd_hic-1.png)<!-- --><br><span style="color:#0d00ff">~[figure @ *public/figures/03A/hoxd_hic...*]~</span>


<!-- END MATTER, insert reproducibility info -->




***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:



```
## 2022-06-16 14:57:36
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
##  date     2022-06-16                      
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
##  P cellranger                          1.1.0      2016-07-27 [?]
##  P cli                                 2.5.0      2021-04-26 [?]
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
##  P org.Mm.eg.db                      * 3.10.0     2021-12-06 [?]
##  P pheatmap                          * 1.0.12     2019-01-04 [?]
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
##  P readxl                            * 1.3.1      2019-03-13 [?]
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
##  P TxDb.Mmusculus.UCSC.mm9.knownGene * 3.2.2      2021-12-06 [?]
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
## [2] /tmp/Rtmp2jElRW/renv-system-library
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
