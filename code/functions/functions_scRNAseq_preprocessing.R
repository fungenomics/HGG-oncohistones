# Functions used for code/scripts/scRNAseq_preprocessing.Rmd


#
# General purpose ----
#

# Gradient colour palettes, yellow-to-red and red-to-blue
ylrd <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "OrRd"))(n = 100)
rdbu <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(n = 100))



get_cells_to_filter <- function(seurat,
                                min_features,
                                max_features,
                                min_umi,
                                max_umi,
                                min_mito,
                                max_mito) {
  
  seurat@meta.data[seurat@meta.data$nFeature_RNA >= min_features &
                     seurat@meta.data$nFeature_RNA <= max_features &
                     seurat@meta.data$nCount_RNA >= min_umi &
                     seurat@meta.data$nCount_RNA <= max_umi &
                     seurat@meta.data$percent.mito >= min_mito &
                     seurat@meta.data$percent.mito <= max_mito, ] %>%
    rownames()
  
}




set_qual_pal <- function(seurat) {
  
  pal_qual <- c("#646BA8",
                "#ef893b",
                "#e2445e",
                "#5E7A41",
                "#FFFC63",
                "#5DAD3B",
                "#6E3688",
                "#A2ACD3",
                "#f2c7d8",
                "#62babd",
                "#519674",
                "#f0c992",
                "#BEBEBE",
                "#2e3082",
                "#61cfe8")
  
  n_clust <- length(levels(seurat@active.ident))
  
  if (n_clust <= length(pal_qual)) { pal_qual_ramped <- head(pal_qual, n_clust)}
  else {
    
    pal_qual_ramped <- colorRampPalette(pal_qual)(n = n_clust)
    
  }
  
  pal_seurat <- pal_qual_ramped %>% setNames(levels(seurat@active.ident))
  seurat@misc$colours <- pal_seurat
  return(seurat)
  
}


#
# For single-cell RNAseq ----
#


#' getVarianceExplained
#'
#' Compute variance explained by PCA, given a Seurat object for which the PCA
#' dim. reduction has been calculated
#'
#' @param seurat Seurat object
#' @param n Numeric, number of PCs for which variance should be reported.
#' Default: 10
#'
#' @return List with two vectors, "percent.var.explained" and "cumulative.var.explained",
#' reported for the first \code{n} PCs
#'
#' @examples
#' getVarianceExplained(pbmc, n = 5)
#'
#' @author Adapted from Yang Yang
get_variance_explained <- function(seurat, n = 10) {
  
  sdev <- seurat@reductions$pca@stdev
  variance <- sdev^2
  sum.variance <- sum(variance)
  proportion.variance <- variance/sum.variance * 100
  acc_prop_var <- cumsum(proportion.variance)
  
  return(list(percent.var.explained = head(proportion.variance, n),
              cum.var.explained     = head(acc_prop_var, n)))
  
}



#' Adapted from Alexis-Blanchet Cohen, by Marie Coutelier
compute_cell_cycle_whitfield <- function(seurat,
                                         species = "m_musculus",
                                         facets = TRUE,
                                         legend = FALSE,
                                         return_scores = FALSE) {
  
  cell.cycle.genes <- cell.cycle.genes.whitfield.2002
  if (species == "m_musculus") {
    cell.cycle.genes$gene.symbol <- cell.cycle.genes.whitfield.2002$mmusculus.gene.symbol
  } else {
    cell.cycle.genes$gene.symbol <- cell.cycle.genes.whitfield.2002$hsapiens.gene.symbol
  }
  
  g1.s.genes <- filter(cell.cycle.genes, phase == "G1/S") %>% .$gene.symbol
  g2.m.genes <- filter(cell.cycle.genes, phase == "G2/M") %>% .$gene.symbol
  
  expression.data <- as.data.frame(as.matrix(GetAssayData(object = seurat)))
  
  expression.data.g1.s.genes <- filter(expression.data, rownames(expression.data) %in% g1.s.genes)
  expression.data.g2.m.genes <- filter(expression.data, rownames(expression.data) %in% g2.m.genes)
  
  expression.data.g1.s.scores <- colMeans(expression.data.g1.s.genes)
  expression.data.g2.m.scores <- colMeans(expression.data.g2.m.genes)
  
  cell.cycle.scores <- as.data.frame(rbind(expression.data.g1.s.scores, expression.data.g2.m.scores))
  
  rownames(cell.cycle.scores) <- gsub("expression.data.", "", rownames(cell.cycle.scores))
  
  cell.cycle.scores.tidy <- as.data.frame(t(cell.cycle.scores))
  cell.cycle.scores.tidy <- tibble::rownames_to_column(cell.cycle.scores.tidy, "cell")
  cell.cycle.scores.tidy <- tibble::add_column(cell.cycle.scores.tidy, cluster=Idents(object = seurat), .after="cell")
  
  if (return_scores) return(cell.cycle.scores.tidy)
  
  # Plots
  p <- ggplot(cell.cycle.scores.tidy, aes(x = g1.s.scores, y = g2.m.scores)) +
    geom_point(aes(color = cluster)) + xlab("G1/S score") + ylab("G2/M score")
  
  if (facets) {
    p <- p + facet_grid(~cluster) }
  
  if (legend == FALSE) {
    p <- p + theme(legend.position="none")
  }
  
  return(p)
}