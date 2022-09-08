---
title: "03B - Analysis of D-V patterning in the hindbrain"
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
## Document index: 03B
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
## public/output/03B
```

```
## public/figures/03B
```



Setting a random seed:



```r
set.seed(100)
```



***



<!-- END OF FRONT MATTER -->


# Overview

Here, we investigate dorsal-ventral (D-V) patterning within the hindbrain, and the 
presence of D-V signals within the tumors. We focus mainly on opposing Nkx6-1 and Pax3
signals in the H3.1 vs H3.3K27M pons tumors.


# Libraries



```r
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
library(readxl)
library(glue)
library(tibble)
library(ggplot2)
library(purrr)
library(pheatmap)
library(gprofiler2)
library(Seurat)
library(cowplot)
library(icytobox)
library(plotly)

source(here("include/style.R")) # contains palettes & plotting utils
source(here("code/functions/scRNAseq.R"))
source(here("code/functions/ssGSEA.R"))
source(here("code/functions/RNAseq.R"))
source(here("code/functions/BentoBox_helpers.R"))
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




# Mouse atlas, extended

## Load data

Examining the expression patterns of our key markers, NKX6-1 and PAX3 in the mouse atlas,
focusing on the samples from the pons/hindbrain:



```r
metadata_all <- read_xlsx(here("data/scRNAseq/references/mouse_atlas_extended/20201028-sample_summary_mouse_extended.xlsx"))

metadata <- metadata_all %>% 
    filter(Location == "Hindbrain/pons")

metadata_per_cluster <- data.table::fread(here("data/scRNAseq/references/mouse_atlas_extended/metadata_20210710_with_qc.tsv"), data.table = FALSE)

atlas_data_path <- read_lines(here("data/scRNAseq/references/mouse_atlas_extended/atlas_path_hydra.tsv"))

titles <- metadata$Timepoint
```



Load objects:



```r
samples_pons <- list()

for (i in metadata$Alias) {
    
    message("@ ", i)
    
    x <- get(load(glue("{atlas_data_path}/{i}/{i}.seurat.Rda")))
    x <- UpdateSeuratObject(x)
    
    # run UMAP if it hasn't been
    if (!("umap" %in% names(names(x@reductions)))) {
        
        x <- RunUMAP(x, dims = 1:20, seed.use = 42, verbose = FALSE)
        
    }
    
    samples_pons[[i]] <- x
    
}

names(samples_pons) <- metadata$Alias

save(samples_pons, file = glue("{out}/mouse_atlas_pons_extended.Rda"))
```



Load our saved data:



```r
load(glue("{out}/mouse_atlas_pons_extended.Rda"))
```





## Extract Nkx6-1/Pax3 regulons

Regulons for the normal mouse pons scRNAseq samples were inferred using [SCENIC](https://github.com/aertslab/SCENIC) as described in Jessa et al, Nature Genetics, 2019. We load the inferred regulons here.



```r
regulons_pons <- map_dfr(metadata$Alias, ~ readRDS(glue("{atlas_data_path}/data/{.x}/{.x}.regulon_target_info.Rds")) %>% 
                             filter(grepl("Nkx6-1|Pax3|Nkx6-2|Pax7", TF)) %>% 
                             filter(TF != gene) %>% 
                             as.data.frame() %>% 
                             mutate(Alias = .x)) %>% 
    mutate(Weight = ifelse(!is.na(Genie3Weight), Genie3Weight, Genie3Weight.weight)) %>% 
    select(Alias, TF, Target = gene, NES, Weight, highConfAnnot) %>% 
    left_join(metadata, by = "Alias")

# check if the # of unique targets is tractable...
length(unique(regulons_pons$Target))
```



```
## [1] 280
```



```r
saveRDS(regulons_pons, file = glue("{out}/SCENIC_regulons_pons.Rds"))
```



How many targets are shared between Nkx6-1 and Pax3 regulons?
And between Nkx6-1/Nkx6-2 regulons?



```r
regulons_pons_per_tf <- split(regulons_pons, regulons_pons$TF)

# nkx6-1 and pax3
(shared_targets_nkx_pax <- base::intersect(regulons_pons_per_tf$`Nkx6-1`$Target, regulons_pons_per_tf$Pax3$Target))
```



```
## [1] "E130114P18Rik" "Sp8"           "Bcl7a"         "Hes5"         
## [5] "Chd7"          "Neurog2"       "Nfia"          "Dhrs3"        
## [9] "Cyp26b1"
```



```r
length(unique(regulons_pons_per_tf$`Nkx6-1`$Target))
```



```
## [1] 118
```



```r
length(unique(regulons_pons_per_tf$Pax3$Target))
```



```
## [1] 100
```



```r
length(shared_targets_nkx_pax)
```



```
## [1] 9
```




This indicates there are very few target genes shared within the two networks (9/200).

Let's try compact representation, where for a given TF, we plot its targets across
development like a "clock" using polar coordinates, with the height being
the weight of the target-TF association. Each target is only plot once, at the first
timepoint at which it is identified as a target of the TF.

We can then colour genes by an external annotation:



```r
#' Generate a "clock plot" for targets of a given TF across timepoints
polar_plot2 <- function(regulons, tf, plot_num, genes_highlight, genes_color) {
    
    network_input <- regulons %>% 
        filter(TF == tf) %>% 
        select(TF, Target, Weight, Timepoint) %>% 
        mutate(Timepoint = factor(Timepoint, levels = names(palette_timepoint)))
    
    # for each target, find the earliest timepoint at which it appears
    earliest_appearance <- network_input %>%
        group_by(TF, Target) %>% 
        top_n(-1, Timepoint) %>%
        rename(Earliest_timepoint = Timepoint) %>% 
        select(-Weight)
    
    network_input2 <- network_input %>% 
        left_join(earliest_appearance, by = c("TF", "Target")) %>% 
        arrange(Earliest_timepoint, desc(Weight)) %>% 
        mutate(Target = factor(Target, levels = unique(.$Target)))
    
    # for each target, only shown the bar for the highest weight across timepoints
    network_input3 <- network_input2 %>% 
        group_by(TF, Target) %>% 
        top_n(1, Weight)
    
    n_targets <- length(unique(network_input3$Target))
    
    # get x-label angle to follow the coordinates, thanks to
    # https://stackoverflow.com/a/36109189
    angles <- 360/(2*pi)*rev( pi/2 + seq( pi/n_targets, 2*pi-pi/n_targets, len = n_targets ))
    # rotate by 180degrees the labels on the left-hand side of the plot
    angles[ ceiling(n_targets/2) : n_targets] <- angles[ ceiling(n_targets/2) : n_targets] + 180
    
    # make palette
    # palette_genes <- rep("gray80", n_targets)
    # names(palette_genes) <- unique(network_input3$Target)
    # palette_genes[genes_highlight] <- highlight_color
    
    other_genes <- base::setdiff(unique(network_input3$Target), genes_highlight)
    palette_genes_other <- rep("gray80", length(other_genes))
    names(palette_genes_other) <- other_genes
    palette_genes <- c(genes_color, palette_genes_other)
    
    network_input3 %>% 
        rr_ggplot(aes(x = Target, y = log10(Weight)), plot_num = plot_num) +
        geom_bar(stat = "identity", aes(fill = Target),
                 position = position_dodge2(width = 2, preserve = "total"),
                 alpha = 0.8) +
        scale_fill_manual(values = palette_genes) +
        coord_polar() +
        theme(panel.border = element_blank(),
              axis.text.x = element_text(
                  # color = palette_genes, 
                  angle = angles)) +
        ggtitle(tf) +
        no_legend()
    
}
```



### BMP/SHH genes

Highlighting genes of interest/involved in SHH or BMP pathways:



```r
# define genes to highlight / palette
shh_notch_genes <- c("Sulf1", "Sulf2", "Dll3", "Hes6", "Ptch1", "Hes5", "Notch1", "Gas1")
bmp_wnt_genes <- c(c("Tcf12", "Tcf4", "Wnt7b", "Id2", "Ctnnb1", "Tcf3", "Dcn"))
genes_highlight <- c(shh_notch_genes, bmp_wnt_genes)
palette_genes_highlight <- c(rep("orange", length(shh_notch_genes)),
                             rep("red3", length(bmp_wnt_genes)))
names(palette_genes_highlight) <- c(shh_notch_genes, bmp_wnt_genes)

# generate the plots
nkx6_polar_plot2 <- polar_plot2(regulons_pons, "Nkx6-1", plot_num = 1,
                                genes_highlight = palette_genes_highlight,
                                genes_color = palette_genes_highlight)
```



```
## ...writing source data of ggplot to public/figures/03B/polar_plot2-1.source_data.tsv
```



```r
pax3_polar_plot2 <- polar_plot2(regulons_pons, "Pax3", plot_num = 2,
                                genes_highlight = palette_genes_highlight,
                                genes_color = palette_genes_highlight)
```



```
## ...writing source data of ggplot to public/figures/03B/polar_plot2-2.source_data.tsv
```



```r
plot_grid(nkx6_polar_plot2, pax3_polar_plot2, nrow = 1)
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/polar_plot2-1.png)<!-- -->

### Ependymal genes 

Let's compute the cell type specificity of the Nkx6-1 targets. This function is defined in `code/functions/scRNAseq.R`



```r
# prep
pct1_feather <- here("data/scRNAseq/references/mouse_atlas_extended/pct_per_cluster.feather")
metadata_mouse <- data.table::fread(here("data/scRNAseq/references/mouse_atlas_extended/metadata_20210710_with_qc.tsv"), data.table = FALSE)
n_cells_df <- metadata_mouse %>% dplyr::select(Cluster = Label, N_cells)

# compute for NKX6-1 and PAX3 targets
pct1_nkx6.1 <- prep_celltype_specificity(regulons_pons_per_tf$`Nkx6-1`$Target)

specificity_nkx6.1 <- compute_celltype_specificity(
    pct1_nkx6.1, n_cells_df, max_across_timepoints = TRUE) %T>%
    rr_write_tsv(glue("{out}/pons_nkx6.1_targets.specificity.tsv"),
                 "Cell-type specificity scores of Nkx6-1 targets in the mouse brain")
```



```
## ...writing description of pons_nkx6.1_targets.specificity.tsv to public/output/03B/pons_nkx6.1_targets.specificity.desc
```



```r
pct1_pax3 <- prep_celltype_specificity(regulons_pons_per_tf$Pax3$Target)

specificity_pax3 <- compute_celltype_specificity(
    pct1_pax3, n_cells_df, max_across_timepoints = TRUE) %T>%
    rr_write_tsv(glue("{out}/pons_pax3_targets.specificity.tsv"),
                 "Cell-type specificity scores of Pax3 targets in the mouse brain")
```



```
## ...writing description of pons_pax3_targets.specificity.tsv to public/output/03B/pons_pax3_targets.specificity.desc
```



Visualize the results, filtering to genes with a specifcity score > 0.5:



```r
specificity_nkx6.1 %>% 
    arrange(desc(Score)) %>%
    filter(Score > 0.5) %>% 
    mutate(Gene = factor(Gene, levels = unique(.$Gene))) %>% 
    ggplot(aes(x = Gene, y = Score)) +
    geom_col(aes(fill = Type)) +
    geom_text(aes(label = Top, color = Type),
              angle = 80, size = 2, fontface = "italic", hjust = -0.05) +
    scale_fill_manual(values = palette_type) +
    scale_color_manual(values = palette_type) +
    rotate_x() +
    ggtitle("Cell-type specificity score of Nkx6-1 targets in normal pons") +
    theme(axis.text.x = element_text(size = rel(0.6)),
          legend.position = "bottom") +
    ylim(c(0, 1.2))
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/visualize_celltype_specificity-1.png)<!-- -->

```r
specificity_pax3 %>% 
    arrange(desc(Score)) %>% 
    filter(Score > 0.5) %>% 
    mutate(Gene = factor(Gene, levels = unique(.$Gene))) %>% 
    ggplot(aes(x = Gene, y = Score)) +
    geom_col(aes(fill = Type)) +
    geom_text(aes(label = Top, color = Type),
              angle = 80, size = 2, fontface = "italic", hjust = -0.05) +
    scale_fill_manual(values = palette_type) +
    scale_color_manual(values = palette_type) +
    rotate_x() +
    ggtitle("Cell-type specificity score of Pax3 targets in normal pons") +
    theme(axis.text.x = element_text(size = rel(0.6)),
          legend.position = "bottom") +
    ylim(c(0, 1.2))
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/visualize_celltype_specificity-2.png)<!-- -->

Generate bubble plot:



```r
dendrogram_order_pons <- readRDS(here("output/05/dendrogram_order_joint_extended_pons.Rds"))

plot_bubble(genes = c("Sulf1", "Metrnl", "Foxj1", "Riiad1", "Cmtm8",
                      "Ntn1", "Sulf2", "Tppp3", "Slit2", "Tnnt3"),
            cluster_col = "ID_20210710",
            meanexp_feather = here("data/scRNAseq/references/mouse_atlas_extended/mean_expression_per_cluster.feather"),
            pct_feather = pct1_feather,
            mean_exp = mouse_atlas_meanexp,
            dendrogram_order = dendrogram_order_pons)
```



```
## All clusters present: TRUE
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/epen_genes_bubbleplot-1.png)<!-- -->

Clean up to save memory:



```r
rm(samples_pons)
```



## Expression heatmaps, by class

In these heatmaps, we will show the expression of Nkx6-1 and Pax3 (rows), in each cell belonging to a given cell class (columns).



```r
load(here("data/scRNAseq/references/mouse_atlas_extended/joint_pons_extended.seurat.Rda"))

# get cell type classes
seurat_joint@meta.data$Cell_type <- seurat_joint@meta.data %>%
    summarize_cell_types("ID_20201028_with_exclude") %>%
    pull(Type)
```



Subset to Nkx6-1+ or Pax3+ cells:



```r
# subset by expression
seurat_pons_positive <- subset(seurat_joint, subset = `Nkx6-1` > 0 | Pax3 > 0)
rm(seurat_joint)

# split by cell type
seurat_pons_positive_split <- SplitObject(seurat_pons_positive, split.by = "Cell_type")
rm(seurat_pons_positive)

# rescale so we can make a heatmap
seurat_pons_positive_split <- map(seurat_pons_positive_split,
                                  ~ ScaleData(.x, vars.to.regress = c("nCount_RNA", "percent_mito")))

# drop non-neuroectodermal classes
seurat_pons_positive_split$`Vascular & other` <- NULL
seurat_pons_positive_split$Immune <- NULL
```



Make heatmaps:



```r
mouse_heatmaps_by_class <- imap(
    seurat_pons_positive_split,
    function(seurat, title) {
        
        df <- FetchData(seurat, c("Nkx6-1", "Pax3")) %>% 
            arrange(desc(`Nkx6-1`), desc(`Pax3`))
        
        n_double_pos <- df %>% filter(`Nkx6-1` > 0 & `Pax3` > 0) %>% nrow()
        
        seurat$id <- "cell"
        DoHeatmap(seurat, c("Nkx6-1", "Pax3"), cells = rownames(df),
                  group.by = "id", label = FALSE, group.bar = FALSE,
                  size = 3, angle = 90, disp.max = 2.5, disp.min = -2.5, raster = FALSE) +
            scale_fill_gradientn(colors = c("#6161B0", "white", "#FF0000"),
                                 na.value = "white", limits = c(-2.5, 2.5)) +
            ggtitle(glue("{title}"),
                    subtitle = glue("{nrow(df)} Nkx6-1+ or Pax3+ \n{n_double_pos} Nkx6-1+ and Pax3+"))
        
    })

plot_grid(plotlist = mouse_heatmaps_by_class, ncol = 2, align = "hv", axis = "rb")
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/pons_heatmaps_by_celltype-1.png)<!-- -->

Clean up to save memory:



```r
rm(seurat_pons_positive_split)
```




# Allen Brain Atlas

To quantify expression patterns observed in E13 and P56 mouse brains from the Allen Brain Atlas, we plot the quantification here. These were downloaded from the indicated links.



```r
# http://mouse.brain-map.org/experiment/show/81657661
p56 <- read_delim(here("data/misc/ABA_NKX6-1_P56.txt"), delim = " ",
                  col_names = c("Expression", "1", "2", "3", "4", "Region"))
```



```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   Expression = col_double(),
##   `1` = col_double(),
##   `2` = col_double(),
##   `3` = col_character(),
##   `4` = col_double(),
##   Region = col_character()
## )
```



```r
# https://developingmouse.brain-map.org/experiment/show/113105268
e13 <- read_delim(here("data/misc/ABA_NKX6-1_E13.5.txt"), delim = "\t",
                  col_names = c("Region", "Expression"))
```



```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   Region = col_character(),
##   Expression = col_double()
## )
```



```r
p1 <- e13 %>%
  mutate(Region = factor(Region, levels = .$Region)) %>% 
  ggplot(aes(x = Region, y = Expression)) +
  geom_bar(stat = "identity", fill = "black") +
  rotate_x() +
  ylim(c(0, max(e13$Expression)))

p2 <- p56 %>%
  mutate(Region = factor(Region, levels = .$Region)) %>% 
  ggplot(aes(x = Region, y = Expression)) +
  geom_bar(stat = "identity", fill = "black") +
  rotate_x() +
  ylim(c(0, max(e13$Expression)))

plot_grid(p1, p2, align = "h", axis = "tb")
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/aba_barplots-1.png)<!-- -->

# Human fetal hindbrain

Load metadata and palette:



```r
info_samples_hindbrain <- data.table::fread(here("data/metadata/metadata_human_fetal_10X.tsv"), data.table = FALSE) %>% 
    filter(Region %in% c("Hindbrain", "Cerebellum"))

palette_hindbrain <- data.table::fread(here("output/01B/hf_hindbrain_annotation.tsv"), data.table = FALSE) %>% 
    tibble::deframe()
```




## Expression heatmaps, by class

Subset to Nkx6-1+ or Pax3+ cells:



```r
load(here("data/scRNAseq/integrations/human_fetal_hindbrain/output/seurat_joint.Rda"))

seurat_joint@meta.data$Cell_type <- seurat_joint@meta.data %>%
    summarize_cell_types("ID_20210610") %>%
    pull(Type)

# subset by expression
seurat_hb_positive <- subset(seurat_joint, subset = `NKX6-1` > 0 | PAX3 > 0)
rm(seurat_joint)

# split by cell type
table(seurat_hb_positive$Cell_type)
```



```
## 
##                  RGC    Glial progenitors    Proliferating OPC 
##                 2635                  338                    0 
##                  OPC     Oligodendrocytes           Astrocytes 
##                    6                    0                    0 
##            Ependymal Neuronal progenitors              Neurons 
##                   22                 2330                 1155 
##               Immune     Vascular & other               Normal 
##                   13                 3653                    0
```



```r
seurat_hb_positive_split <- SplitObject(seurat_hb_positive, split.by = "Cell_type")
rm(seurat_hb_positive)

# rescale so we can make a heatmap
seurat_hb_positive_split <- map(seurat_hb_positive_split,
                                ~ ScaleData(.x, vars.to.regress = c("nCount_RNA", "percent.mito"), features = c("NKX6-1", "PAX3")))

# drop non-neuroectodermal classes
seurat_hb_positive_split$`Vascular & other` <- NULL
seurat_hb_positive_split$Immune <- NULL
```



Make heatmaps:



```r
hb_heatmaps_by_class <- imap(
    seurat_hb_positive_split,
    function(seurat, title) {
        
        df <- FetchData(seurat, c("NKX6-1", "PAX3")) %>% 
            arrange(desc(`NKX6-1`), desc(`PAX3`))
        
        n_double_pos <- df %>% filter(`NKX6-1` > 0 & `PAX3` > 0) %>% nrow()
        
        seurat$id <- "cell"
        DoHeatmap(seurat, c("NKX6-1", "PAX3"), cells = rownames(df),
                  group.by = "id", label = FALSE, group.bar = FALSE,
                  size = 3, angle = 90, disp.max = 2.5, disp.min = -2.5, raster = FALSE) +
            scale_fill_gradientn(colors = c("#6161B0", "white", "#FF0000"),
                                 na.value = "white", limits = c(-2.5, 2.5)) +
            ggtitle(glue("{title}"),
                    subtitle = glue("{nrow(df)} NKX6-1+ or PAX3+ \n{n_double_pos} NKX6-1+ and PAX3+"))
        
    })
```



```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```



```r
plot_grid(plotlist = hb_heatmaps_by_class, ncol = 2, align = "hv", axis = "rb")
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/hb_heatmaps_by_celltype-1.png)<!-- -->



# Tumor RNAseq

## H3.1 vs H3.3




```r
K27M_pons_rna <- "data/RNAseq/pipeline_l3/HGG-H3K27M-pons_with_DFCI"

K27M_pons_info <- data.table::fread(here(file.path(K27M_pons_rna, "info.samples.tsv")), data.table = FALSE) %>% 
    left_join(meta, by = c("Nickname" = "ID_paper", "Group"))

table(K27M_pons_info$Group)
```



```
## 
## HGG-H3.1/2K27M-Pons   HGG-H3.3K27M-Pons 
##                  19                  14
```



```r
K27M_pons_counts <- extract_pipeline_counts(path = here(file.path(K27M_pons_rna, "counts/Ensembl.ensGene.exon.norm.tsv.gz")),
                                            goi = c("NKX6-1", "PAX3")) %>% 
    left_join(K27M_pons_info, by = c("sample" = "Nickname")) %>% 
    mutate(ACVR1_simple = case_when(
        grepl("ACVR1-", ACVR1) ~ "Mutant",
        TRUE ~ "WT"
    ))
```



Visualize expression of NKX6-1/PAX3 genes, split by K27M status, and encoding ACVR1 status in the shape:



```r
p1 <- K27M_pons_counts %>% 
    ggplot(aes(x = Group, y = gene_expression, label = sample)) +
    geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.4) +
    geom_jitter(alpha = 0.7, size = 2, width = 0.2, aes(shape = ACVR1_simple)) +
    # Open circles: WT, closed circles: mutant
    scale_shape_manual(values = palette_acvr1_simple) +
    scale_fill_manual(values = palette_groups) +
    facet_wrap(~ gene_symbol, nrow = 1, scales = "free_y") +
    no_legend() +
    rotate_x() +
    ylab("Expression (norm)") +
    ggtitle("PAX3/NKX6-1 expression in bulk RNAseq \nfor H3.1/2 and H3.3K27M Pons")

p2 <- K27M_pons_counts %>% 
    select(sample, Group, GrowthFactorReceptor, gene_symbol, gene_expression, ACVR1_simple) %>% 
    pivot_wider(names_from = gene_symbol, values_from = gene_expression) %>% 
    ggplot(aes(x = `NKX6-1`, y = PAX3)) +
    geom_point(aes(fill = Group, color = Group, shape = ACVR1_simple), size = 4, alpha = 0.7) +
    scale_shape_manual(values = palette_acvr1_simple) +
    scale_fill_manual(values = palette_groups) +
    scale_color_manual(values = palette_groups) +
    xlab("NKX6-1 expression (Norm)") + ylab("PAX3 expression (Norm)") +
    facet_wrap(~ Group) +
    theme(legend.position = "bottom")

plot_grid(p1, p2, ncol = 2, align = "h", axis = "tb")
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/h3.1_vs_h3.3_bulk_expression_by_acvr1-1.png)<!-- -->



# Tumor ChIP/scATAC

In this section, we visualize ChIPseq and scATAC data at regions relevant
to NKX6-1/PAX3.

For this section, I heavily use helper functions defined in `code/functions/BentoBox_helpers.R`.
To configure the plots, I define the samples to include in configuration TSV
files at `data/BentoBox_config/hoxd_*.config.tsv`

## NKX6-1/PAX3 promoters {.tabset}

First, we consider the regions surrounding each gene.

### H3.1K27M

Get config for the tracks:



```r
# load config for H3.1 samples
nkx6_config <- prep_track_input(here("data/BentoBox_config/NKX6-1_H31.config.tsv"))
nkx6_region <- GRanges("chr4", 85412437:85421387, name = "NKX6-1")
params_nkx6 <- bb_params(chrom = as.character(seqnames(nkx6_region)),
                         chromstart = start(nkx6_region),
                         chromend = end(nkx6_region),
                         assembly = "hg19")

pax3_region <- GRanges("chr2", 223040199:223188125, name = "PAX3")
params_pax3 <- bb_params(chrom = as.character(seqnames(pax3_region)),
                         chromstart = start(pax3_region),
                         chromend = end(pax3_region),
                         assembly = "hg19")

# check ACVR1 status for each
nkx6_config %>%
    separate(ID_paper, into = c("ID_paper", "extra"), sep = " ") %>%
    left_join(meta %>% dplyr::select(ID_paper, GrowthFactorReceptor), by = "ID_paper") %>% 
    dplyr::select(ID_paper, Data, GrowthFactorReceptor)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ID_paper"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Data"],"name":[2],"type":["chr"],"align":["left"]},{"label":["GrowthFactorReceptor"],"name":[3],"type":["chr"],"align":["left"]}],"data":[{"1":"P-1741_S-2756","2":"H3K27ac","3":"ACVR1-G328V"},{"1":"P-1703_S-2705","2":"H3K27ac","3":"None"},{"1":"DIPGIV","2":"H3K27ac","3":"ACVR1-G328V"},{"1":"DIPG21","2":"H3K27ac","3":"ACVR1-G328W"},{"1":"DIPG36","2":"H3K27ac","3":"ACVR1-G328E"},{"1":"DIPGIV","2":"H3K27me3","3":"ACVR1-G328V"},{"1":"DIPG21","2":"H3K27me3","3":"ACVR1-G328W"},{"1":"DIPG36","2":"H3K27me3","3":"ACVR1-G328E"},{"1":"P-6253_S-8498","2":"scATACseq","3":"ACVR1-G328E"},{"1":"P-2687_S-2688","2":"scATACseq","3":"ACVR1-G328V"},{"1":"P-1780_S-1780","2":"scMultiome_ATAC","3":"None"},{"1":"P-1780_S-1780","2":"scMultiome_ATAC","3":"None"},{"1":"P-1780_S-1780","2":"scMultiome_ATAC","3":"None"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Make figure:



```r
# NKX6-1 ------
x <- 0.5
y_positions <- seq(0.5, by = 0.35, length.out = nrow(nkx6_config))

bb_pageCreate(width = 8, height = 7, default.units = "inches")

# add bw tracks
pwalk(list(nkx6_config$bw, nkx6_config$ID_paper, nkx6_config$Data, nkx6_config$Ymax, y_positions),
      ~ bb_placeSignalAndLabel(data = ..1,
                               annotation = ..2,
                               color = palette_tracks[..3],
                               range = c(0, ..4),
                               y = ..5, x = x,
                               params = params_nkx6,
                               height = 0.3,
                               fontsize = 4))

bb_plotGenomeLabel(params = params_nkx6, scale = "Kb",
                   x = x, y = "0.03b", length = 3)

bb_plotGenes(params = params_nkx6, x = x, y = "0.1b", height = 1, width = 3,
             strandcolors = c("navy", "black"), fontcolors = c("navy", "black"),
             stroke = 0.05,
             just = c("left", "top"), default.units = "inches")

# PAX3 ------
x <- 4
y_positions <- seq(0.5, by = 0.35, length.out = nrow(nkx6_config))

# add bw tracks
pwalk(list(nkx6_config$bw, nkx6_config$ID_paper, nkx6_config$Data, nkx6_config$Ymax, y_positions),
      ~ bb_placeSignalAndLabel(data = ..1,
                               annotation = ..2,
                               color = palette_tracks[..3],
                               range = c(0, ..4),
                               y = ..5, x = x,
                               params = params_pax3,
                               height = 0.3,
                               fontsize = 4))

bb_plotGenomeLabel(params = params_pax3, scale = "Kb",
                   x = x, y = "0.03b", length = 3)

bb_plotGenes(params = params_pax3, x = x, y = "0.1b", height = 1, width = 3,
             strandcolors = c("navy", "black"), fontcolors = c("navy", "black"),
             stroke = 0.05,
             just = c("left", "top"), default.units = "inches")

bb_pageGuideHide()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/nkx6_pax3_tracks-1.png)<!-- -->

### H3.3K27M

Get config for the tracks:



```r
# load config for H3.3 samples
nkx6_config2 <- prep_track_input(here("data/BentoBox_config/NKX6-1_H33.config.tsv"))

# check ACVR1 status for each
nkx6_config2 %>%
    separate(ID_paper, into = c("ID_paper", "extra"), sep = " ") %>%
    left_join(meta %>% dplyr::select(ID_paper, GrowthFactorReceptor), by = "ID_paper") %>% 
    dplyr::select(ID_paper, Data, GrowthFactorReceptor)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ID_paper"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Data"],"name":[2],"type":["chr"],"align":["left"]},{"label":["GrowthFactorReceptor"],"name":[3],"type":["chr"],"align":["left"]}],"data":[{"1":"P-1752_SX-2748","2":"H3K27ac","3":"None"},{"1":"P-3147_S-3146","2":"H3K27ac","3":"None"},{"1":"DIPGVI","2":"H3K27ac","3":"None"},{"1":"DIPG007","2":"H3K27ac","3":"ACVR1-R206H"},{"1":"DIPGXIII","2":"H3K27ac","3":"None"},{"1":"P-1752_SX-2748","2":"H3K27me3","3":"None"},{"1":"P-3147_S-3146","2":"H3K27me3","3":"None"},{"1":"DIPGVI","2":"H3K27me3","3":"None"},{"1":"DIPG007","2":"H3K27me3","3":"ACVR1-R206H"},{"1":"DIPGXIII","2":"H3K27me3","3":"None"},{"1":"P-1764_S-1766","2":"scATACseq","3":"PDGFRA mut&amp"},{"1":"P-1694_S-1694","2":"scMultiome_ATAC","3":"None"},{"1":"P-1701_S-1701","2":"scMultiome_ATAC","3":"ACVR1-R206H"},{"1":"P-1779_S-1781","2":"scMultiome_ATAC","3":"None"},{"1":"P-1764_S-1766","2":"scMultiome_ATAC","3":"PDGFRA mut&amp"},{"1":"P-6774_S-10146","2":"scMultiome_ATAC","3":"None"},{"1":"P-1764_S-1766","2":"scMultiome_ATAC","3":"PDGFRA mut&amp"},{"1":"P-6337_S-8821","2":"scMultiome_ATAC","3":"None"},{"1":"P-6774_S-10146","2":"scMultiome_ATAC","3":"None"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Make figure:



```r
# NKX6-1 ------
x <- 0.5
y_positions <- seq(0.5, by = 0.35, length.out = nrow(nkx6_config2))

bb_pageCreate(width = 8, height = 8.5, default.units = "inches")

# add bw tracks
pwalk(list(nkx6_config2$bw, nkx6_config2$ID_paper, nkx6_config2$Data, nkx6_config2$Ymax, y_positions),
      ~ bb_placeSignalAndLabel(data = ..1,
                               annotation = ..2,
                               color = palette_tracks[..3],
                               range = c(0, ..4),
                               y = ..5, x = x,
                               params = params_nkx6,
                               height = 0.3,
                               fontsize = 4))

bb_plotGenomeLabel(params = params_nkx6, scale = "Kb",
                   x = x, y = "0.03b", length = 3)

bb_plotGenes(params = params_nkx6, x = x, y = "0.1b", height = 1, width = 3,
             strandcolors = c("navy", "black"), fontcolors = c("navy", "black"),
             stroke = 0.05,
             just = c("left", "top"), default.units = "inches")

# PAX3 ------
x <- 4
y_positions <- seq(0.5, by = 0.35, length.out = nrow(nkx6_config2))

# add bw tracks
pwalk(list(nkx6_config2$bw, nkx6_config2$ID_paper, nkx6_config2$Data, nkx6_config2$Ymax, y_positions),
      ~ bb_placeSignalAndLabel(data = ..1,
                               annotation = ..2,
                               color = palette_tracks[..3],
                               range = c(0, ..4),
                               y = ..5, x = x,
                               params = params_pax3,
                               height = 0.3,
                               fontsize = 4))

bb_plotGenomeLabel(params = params_pax3, scale = "Kb",
                   x = x, y = "0.03b", length = 3)

bb_plotGenes(params = params_pax3, x = x, y = "0.1b", height = 1, width = 3,
             strandcolors = c("navy", "black"), fontcolors = c("navy", "black"),
             stroke = 0.05,
             just = c("left", "top"), default.units = "inches")

bb_pageGuideHide()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/nkx6_pax3_tracks2-1.png)<!-- -->

## NKX6-1 cis-regulatory elements

Second, we consider a region neighbouring NKX6-1 containing cis-regulatory elements
in H3.1K27M tumors.

We'll also show scATAC data from (the multiome samples) around NKX6-1:



```r
# define tracks (function defined in code/functions/BentoBox_helpers.R)
nkx6_cis_config <- prep_track_input(here("data/BentoBox_config/NKX6-1_ciselements_H31.config.tsv"))
nkx6_cis_region <- GRanges("chr4", 85230984:85514706, name = "NKX6-1_ce")
params_nkx6_cis <- bb_params(chrom = as.character(seqnames(nkx6_cis_region)),
                             chromstart = start(nkx6_cis_region),
                             chromend = end(nkx6_cis_region),
                             assembly = "hg19")

# acvr1 status
nkx6_cis_config %>%
    separate(ID_paper, into = c("ID_paper", "extra"), sep = " ") %>%
    left_join(meta %>% dplyr::select(ID_paper, GrowthFactorReceptor), by = "ID_paper") %>% 
    dplyr::select(ID_paper, Data, GrowthFactorReceptor)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ID_paper"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Data"],"name":[2],"type":["chr"],"align":["left"]},{"label":["GrowthFactorReceptor"],"name":[3],"type":["chr"],"align":["left"]}],"data":[{"1":"P-6253_S-8498","2":"scATACseq","3":"ACVR1-G328E"},{"1":"P-2687_S-2688","2":"scATACseq","3":"ACVR1-G328V"},{"1":"P-1780_S-1780","2":"scMultiome_ATAC","3":"None"},{"1":"P-6253_S-8498","2":"scMultiome_ATAC","3":"ACVR1-G328E"},{"1":"P-6253_S-8498","2":"scATACseq","3":"ACVR1-G328E"},{"1":"P-2687_S-2688","2":"scATACseq","3":"ACVR1-G328V"},{"1":"P-6253_S-8498","2":"scMultiome_ATAC","3":"ACVR1-G328E"},{"1":"P-1780_S-1780","2":"scMultiome_ATAC","3":"None"},{"1":"P-1780_S-1780","2":"scMultiome_ATAC","3":"None"},{"1":"P-1780_S-1780","2":"scMultiome_ATAC","3":"None"},{"1":"P-6253_S-8498","2":"scMultiome_RNA","3":"ACVR1-G328E"},{"1":"P-1780_S-1780","2":"scMultiome_RNA","3":"None"},{"1":"P-1780_S-1780","2":"scMultiome_RNA","3":"None"},{"1":"P-1780_S-1780","2":"scMultiome_RNA","3":"None"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Make figure:



```r
# NKX6-1 ------
x <- 0.5
y_positions <- seq(0.5, by = 0.4, length.out = nrow(nkx6_cis_config))

bb_pageCreate(width = 8, height = 8, default.units = "inches")

# add bw tracks
pwalk(list(nkx6_cis_config$bw, nkx6_cis_config$ID_paper, nkx6_cis_config$Color, nkx6_cis_config$Ymax, y_positions),
      # helper function defined in code/functions/BentoBox_helpers.R
      ~ bb_placeSignalAndLabel(data = ..1,
                               annotation = ..2,
                               color = ..3,
                               range = c(0, ..4),
                               y = ..5, x = x,
                               width = 7,
                               height = 0.35,
                               params = params_nkx6_cis,
                               fontsize = 4))

bb_plotGenomeLabel(params = params_nkx6_cis, scale = "Kb",
                   x = x, y = "0.03b", length = 7)

# plot BED file of VISTA enhancers
bb_plotBed(data = here("data/scATACseq/references/vista_enhancers_hg19.bed"),
           width = 7,
           height = 0.1, fill = "navy",
           params = params_nkx6_cis,
           y = "0.2b", x = x, collapse = TRUE)

# plot liftover of region containing NKX6-1 CRM
bb_plotBed(data = GRanges("chr4", 85272335:85297195),
           width = 7,
           height = 0.1, fill = "navy",
           params = params_nkx6_cis,
           y = "0.2b", x = x, collapse = TRUE)

bb_plotGenes(params = params_nkx6_cis, x = x, y = "0.1b", height = 1, width = 7,
             strandcolors = c("navy", "black"), fontcolors = c("navy", "black"),
             stroke = 0.05,
             just = c("left", "top"), default.units = "inches")

bb_pageGuideHide()
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/nkx6_cis_tracks-1.png)<!-- -->



## CRC/enhancer landscape

Here, we visualize results from differential enhancer analysis and core regulatory
circuitry analysis which was performed by the Mack lab.



```r
enhancers_H3.1_vs_H3.3 <- data.table::fread(here(
    "data/ChIPseq/CRC_SEs_@svaradharajan/20220202-EnhancersComparisons_H31_Pons_vs_H33_Pons.txt"),
    data.table = FALSE) %>% 
    rename(TF = TFs,
           H3.1K27M_pons_score = H31_P,
           H3.3K27M_pons_score = H33_P,
           Foldchange = FC,
           Log_foldchange = LFC)

to_label <- enhancers_H3.1_vs_H3.3 %>% 
    filter(!between(dense_rank(Log_foldchange), 10 + 1, n() - 10))

enhancers_H3.1_vs_H3.3 %>% 
    rr_ggplot(aes(x = rank, y = Log_foldchange), plot_num = 1) +
    geom_hline(yintercept = 0, color = "gray90") +
    geom_point(aes(colour = Log_foldchange), size = 3, shape = 21) +
    scale_colour_gradientn(colors = c("red3", "white", "orange"), limits = c(-9.5, 9.5)) +
    geom_text_repel(data = to_label %>% filter(Log_foldchange > 0), aes(label = TF), size = 2, max.overlaps = 100, nudge_x = 30) +
    geom_text_repel(data = to_label %>% filter(Log_foldchange < 0), aes(label = TF), size = 2, max.overlaps = 100, nudge_x = -30)
```



```
## ...writing source data of ggplot to public/figures/03B/enhancer_diff-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/enhancer_diff-1.png)<!-- -->



Visualizing the difference in in/out degree of core TFs (i.e. differences
between numbers of regulators and targets) between tumor type:



```r
crc_delta <- data.table::fread(here(
    "data/ChIPseq/CRC_SEs_@svaradharajan/20220202-CRC_h31_h33_table.txt"),
    data.table = FALSE)

to_label <- crc_delta %>% 
    filter(abs(log2FoldChange) > 1 | abs(H3K27ac) > 1) 

crc_delta %>% 
    arrange(abs(H3K27ac)) %>% 
    filter(!is.na(H3K27ac)) %>% 
    rr_ggplot(aes(x = DeltaIn, y = DeltaOut), plot_num = 1) +
    geom_hline(yintercept = 0, colour = "gray90") +
    geom_vline(xintercept = 0, colour = "gray90") +
    geom_point(aes(colour = H3K27ac, size = abs(log2FoldChange)), alpha = 0.8) +
    geom_text_repel(data = to_label, aes(label = TFs), size = 3) +
    scale_colour_gradientn(colours = c("red3", "#E37171", "white", "#FFCD71", "orange"), limits = c(-2, 2))
```



```
## ...writing source data of ggplot to public/figures/03B/crc_delta_in_out-1.source_data.tsv
```

![](/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/figures/03B/crc_delta_in_out-1.png)<!-- -->


<!-- END MATTER, insert reproducibility info -->




***

<!-- Create reproducibility receipt e.g. https://github.com/benmarwick/rrtools/blob/master/inst/templates/paper.Rmd -->

# Reproducibility

This document was last rendered on:



```
## 2022-09-08 15:19:11
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
##  ! package                           * version    date       lib
##  P abind                               1.4-5      2016-07-21 [?]
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
##  P cluster                             2.1.0      2019-06-19 [?]
##  P codetools                           0.2-16     2018-12-24 [?]
##  P colorspace                          2.0-1      2021-05-04 [?]
##  P cowplot                           * 1.1.1      2020-12-30 [?]
##  P crayon                              1.4.1      2021-02-08 [?]
##  P curl                                4.3.1      2021-04-30 [?]
##  P data.table                        * 1.14.0     2021-02-21 [?]
##  P DBI                                 1.1.1      2021-01-15 [?]
##  P dbplyr                              2.1.1      2021-04-06 [?]
##  P DelayedArray                        0.12.3     2020-04-09 [?]
##  P deldir                              0.2-10     2021-02-16 [?]
##  P desc                                1.2.0      2018-05-01 [?]
##  P devtools                            2.3.0      2020-04-10 [?]
##  P digest                              0.6.27     2020-10-24 [?]
##  P dplyr                             * 1.0.6      2021-05-05 [?]
##  P ellipsis                            0.3.2      2021-04-29 [?]
##  P evaluate                            0.14       2019-05-28 [?]
##  P fansi                               0.4.2      2021-01-15 [?]
##  P fastmap                             1.1.0      2021-01-25 [?]
##  P fitdistrplus                        1.1-3      2020-12-05 [?]
##  P fs                                  1.5.0      2020-07-31 [?]
##  P future                              1.21.0     2020-12-10 [?]
##  P future.apply                        1.7.0      2021-01-04 [?]
##  P generics                            0.1.0      2020-10-31 [?]
##  P GenomeInfoDb                      * 1.22.1     2020-03-27 [?]
##  P GenomeInfoDbData                    1.2.2      2021-12-06 [?]
##  P GenomicAlignments                   1.22.1     2019-11-12 [?]
##  P GenomicFeatures                   * 1.38.2     2020-02-15 [?]
##  P GenomicRanges                     * 1.38.0     2019-10-29 [?]
##  P ggplot2                           * 3.3.3      2020-12-30 [?]
##  P ggplotify                         * 0.0.7      2021-05-11 [?]
##  P ggrepel                           * 0.9.1      2021-01-15 [?]
##  P ggridges                            0.5.3      2021-01-08 [?]
##  P git2r                               0.27.1     2020-05-03 [?]
##  P globals                             0.14.0     2020-11-22 [?]
##  P glue                              * 1.4.2      2020-08-27 [?]
##  P gmp                               * 0.6-2      2021-01-07 [?]
##  P goftest                             1.2-2      2019-12-02 [?]
##  P gprofiler2                        * 0.2.0      2020-08-27 [?]
##  P gridExtra                           2.3        2017-09-09 [?]
##  P gridGraphics                        0.5-1      2020-12-13 [?]
##  P gtable                              0.3.0      2019-03-25 [?]
##  P here                              * 0.1        2017-05-28 [?]
##  P hms                                 1.0.0      2021-01-13 [?]
##  P htmltools                           0.5.1.1    2021-01-22 [?]
##  P htmlwidgets                         1.5.3      2020-12-10 [?]
##  P httpuv                              1.6.1      2021-05-07 [?]
##  P httr                                1.4.2      2020-07-20 [?]
##  P ica                                 1.0-2      2018-05-24 [?]
##  P icytobox                          * 1.0.1      2021-12-07 [?]
##  P igraph                              1.2.6      2020-10-06 [?]
##  P IRanges                           * 2.20.2     2020-01-13 [?]
##  P irlba                               2.3.3      2019-02-05 [?]
##  P jquerylib                           0.1.4      2021-04-26 [?]
##  P jsonlite                            1.7.2      2020-12-09 [?]
##  P KernSmooth                          2.23-15    2015-06-29 [?]
##  P knitr                               1.33       2021-04-24 [?]
##  P later                               1.0.0      2019-10-04 [?]
##  P lattice                             0.20-44    2021-05-02 [?]
##  P lazyeval                            0.2.2      2019-03-15 [?]
##  P leiden                              0.3.7      2021-01-26 [?]
##  P lifecycle                           1.0.0      2021-02-15 [?]
##  P listenv                             0.8.0      2019-12-05 [?]
##  P lmtest                              0.9-38     2020-09-09 [?]
##  P magrittr                          * 2.0.1      2020-11-17 [?]
##  P MASS                                7.3-54     2021-05-03 [?]
##  P Matrix                              1.2-18     2019-11-27 [?]
##  P matrixStats                         0.58.0     2021-01-29 [?]
##  P memoise                             1.1.0      2017-04-21 [?]
##  P mgcv                                1.8-35     2021-04-18 [?]
##  P mime                                0.10       2021-02-13 [?]
##  P miniUI                              0.1.1.1    2018-05-18 [?]
##  P munsell                             0.5.0      2018-06-12 [?]
##  P nlme                                3.1-152    2021-02-04 [?]
##  P openssl                             1.4.4      2021-04-30 [?]
##  P org.Hs.eg.db                      * 3.10.0     2021-12-06 [?]
##  P parallelly                          1.25.0     2021-04-30 [?]
##  P patchwork                           1.1.1      2020-12-17 [?]
##  P pbapply                           * 1.4-3      2020-08-18 [?]
##  P pheatmap                          * 1.0.12     2019-01-04 [?]
##  P pillar                              1.6.0      2021-04-13 [?]
##  P pkgbuild                            1.0.8      2020-05-07 [?]
##  P pkgconfig                           2.0.3      2019-09-22 [?]
##  P pkgload                             1.0.2      2018-10-29 [?]
##  P plotly                            * 4.9.3      2021-01-10 [?]
##  P plyr                                1.8.6      2020-03-03 [?]
##  P png                                 0.1-7      2013-12-03 [?]
##  P polyclip                            1.10-0     2019-03-14 [?]
##  P prettyunits                         1.1.1      2020-01-24 [?]
##  P processx                            3.5.2      2021-04-30 [?]
##  P progress                            1.2.2      2019-05-16 [?]
##  P promises                            1.1.0      2019-10-04 [?]
##  P ps                                  1.6.0      2021-02-28 [?]
##  P purrr                             * 0.3.4      2020-04-17 [?]
##  P R6                                  2.5.0      2020-10-28 [?]
##  P RANN                                2.6.1      2019-01-08 [?]
##  P rappdirs                            0.3.3      2021-01-31 [?]
##  P RColorBrewer                      * 1.1-2      2014-12-07 [?]
##  P Rcpp                                1.0.6      2021-01-15 [?]
##  P RcppAnnoy                           0.0.18     2020-12-15 [?]
##  P RCurl                               1.98-1.3   2021-03-16 [?]
##  P readr                             * 1.4.0      2020-10-05 [?]
##  P readxl                            * 1.3.1      2019-03-13 [?]
##  P remotes                             2.1.1      2020-02-15 [?]
##    renv                                0.14.0     2021-07-21 [1]
##  P reshape2                            1.4.4      2020-04-09 [?]
##  P reticulate                          1.20       2021-05-03 [?]
##  P rlang                               0.4.11     2021-04-30 [?]
##  P rmarkdown                           2.8        2021-05-07 [?]
##  P Rmpfr                             * 0.8-4      2021-04-11 [?]
##  P ROCR                                1.0-11     2020-05-02 [?]
##  P rpart                               4.1-15     2019-04-12 [?]
##  P rprojroot                           2.0.2      2020-11-15 [?]
##  P Rsamtools                           2.2.3      2020-02-23 [?]
##  P RSQLite                             2.2.1      2020-09-30 [?]
##  P rstudioapi                          0.13       2020-11-12 [?]
##  P rsvd                                1.0.3      2020-02-17 [?]
##  P rtracklayer                       * 1.46.0     2019-10-29 [?]
##  P Rtsne                               0.15       2018-11-10 [?]
##  P rvcheck                             0.1.8      2020-03-01 [?]
##  P S4Vectors                         * 0.24.4     2020-04-09 [?]
##  P sass                                0.4.0      2021-05-12 [?]
##  P scales                              1.1.1      2020-05-11 [?]
##  P sctransform                         0.3.2      2020-12-16 [?]
##  P sessioninfo                         1.1.1      2018-11-05 [?]
##  P Seurat                            * 3.2.1      2020-09-07 [?]
##  P shiny                               1.6.0      2021-01-25 [?]
##  P spatstat                            1.64-1     2020-05-12 [?]
##  P spatstat.data                       2.1-0      2021-03-21 [?]
##  P spatstat.utils                      2.1-0      2021-03-15 [?]
##  P stringi                             1.6.1      2021-05-10 [?]
##  P stringr                             1.4.0      2019-02-10 [?]
##  P SummarizedExperiment                1.16.1     2019-12-19 [?]
##  P survival                            3.2-11     2021-04-26 [?]
##  P tensor                              1.5        2012-05-05 [?]
##  P testrmd                             0.0.1.9000 2021-12-06 [?]
##  P testthat                            2.3.2      2020-03-02 [?]
##  P tibble                            * 3.1.1      2021-04-18 [?]
##  P tidyr                             * 1.1.3      2021-03-03 [?]
##  P tidyselect                          1.1.1      2021-04-30 [?]
##  P TxDb.Hsapiens.UCSC.hg19.knownGene * 3.2.2      2021-12-06 [?]
##  P usethis                             1.6.1      2020-04-29 [?]
##  P utf8                                1.2.1      2021-03-12 [?]
##  P uwot                                0.1.10     2020-12-15 [?]
##  P vctrs                               0.3.8      2021-04-29 [?]
##  P viridis                           * 0.5.1      2018-03-29 [?]
##  P viridisLite                       * 0.4.0      2021-04-13 [?]
##  P withr                               2.4.2      2021-04-18 [?]
##  P xfun                                0.22       2021-03-11 [?]
##  P XML                                 3.99-0.3   2020-01-20 [?]
##  P xtable                              1.8-4      2019-04-21 [?]
##  P XVector                           * 0.26.0     2019-10-29 [?]
##  P yaml                                2.2.1      2020-02-01 [?]
##  P zlibbioc                            1.32.0     2019-10-29 [?]
##  P zoo                                 1.8-9      2021-03-09 [?]
##  source                                
##  CRAN (R 3.6.1)                        
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
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  Github (fungenomics/icytobox@730e8b8) 
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
##  CRAN (R 3.6.1)                        
##  Bioconductor                          
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  Bioconductor                          
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
##  CRAN (R 3.6.1)                        
##  CRAN (R 3.6.1)                        
##  Bioconductor                          
##  CRAN (R 3.6.1)                        
##  Bioconductor                          
##  CRAN (R 3.6.1)                        
## 
## [1] /lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/public/renv/library/R-3.6/x86_64-pc-linux-gnu
## [2] /tmp/RtmpPoDSn2/renv-system-library
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
