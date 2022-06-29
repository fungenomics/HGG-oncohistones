
# BentoBox helpers -------------------------------------------------------------

#' Given the sample IDs and data types, find the path to the .bw files
#' within the project directory
prep_track_input <- function(config) {
    
    if (is.character(config)) config <- suppressWarnings(data.table::fread(config, data.table = FALSE))
    
    config <- config %>% tibble::rowid_to_column(var = "Order")
    
    config_rna <- config %>%
        filter(Data == "RNAseq") %>%
        dplyr::select(-Replicate) %>%
        dplyr::left_join(meta, by = c("ID_paper", "Material")) %>%
        distinct(Data, ID_paper, .keep_all = TRUE) %>%
        mutate(bw = file.path("/lustre06/project/6004736/pipeline/v0/levels1-2", RNAseq_path,
                              "star", paste0(basename(RNAseq_path), ".sorted.bw")))
    
    config_chip <- config %>%
        filter(Data %in% c("H3K27ac", "H3K27me3", "CTCF", "SUZ12")) %>%
        dplyr::left_join(meta_chip %>% dplyr::rename(Data = Factor),
                         by = c("ID_paper", "Material", "Data", "Replicate")) %>%
        distinct(Data, ID_paper, .keep_all = TRUE) %>%
        rowwise() %>% 
        mutate(bw = case_when(
            grepl("lustre", Path_bw) ~ Path_bw, # absolute path
            TRUE ~ here(Path_bw) # relative path
        ))
    
    config_atac <- config %>%
        filter(Data == "scATACseq") %>%
        dplyr::left_join(meta_chip, by = c("ID_paper", "Material")) %>%
        distinct(Data, ID_paper, Celltype, .keep_all = TRUE) %>%
        mutate(bw = here(paste0("data/scATACseq/pipeline_10X_ATAC/", ID_paper, "/celltype_bw/output/", Celltype, ".bw")),
               ID_paper = paste0(ID_paper, " ", Celltype))
    
    config_multiome_atac <- config %>% 
        filter(Data == "scMultiome_ATAC") %>%
        distinct(Data, ID_paper, Celltype, .keep_all = TRUE) %>%
        mutate(bw = here(paste0("R-4/data/scMultiome/pipeline_10X_Multiome/", ID_paper, "/celltype_bw/output/", Celltype, ".ATAC.bw")),
               ID_paper = paste0(ID_paper, " ", Celltype, " ATAC"))
    
    config_multiome_rna <- config %>% 
        filter(Data == "scMultiome_RNA") %>%
        distinct(Data, ID_paper, Celltype, .keep_all = TRUE) %>%
        mutate(bw = here(paste0("R-4/data/scMultiome/pipeline_10X_Multiome/", ID_paper, "/celltype_bw/output/", Celltype, ".RNA.bw")),
               ID_paper = paste0(ID_paper, " ", Celltype, " RNA"))
    
    if (!("Ymax" %in% colnames(config_chip))) config_chip$Ymax <- NA
    if (!("Ymax" %in% colnames(config_multiome_atac))) config_multiome_atac$Ymax <- NA
    
    config_all <- bind_rows(config_rna, config_chip, config_atac, config_multiome_atac, config_multiome_rna)
    
    if (!("Color" %in% colnames(config_all))) config_all$Color <- NA
    
    # check files exist
    config_all$Exists <- file.exists(config_all$bw)
    
    if (any(!config_all$Exists)) warning("Some .bw files not found: ", config_all$ID_paper[!config_all$Exists])
    
    # reorder to match config
    config_all <- config_all %>% arrange(Order)
    
    return(dplyr::select(config_all, Data, ID_paper, bw, Ymax, Color, Exists))
    
}


#' Automate the placement of one heatmap with legend
#' 
#' This function expects that bb_pageCreate() has already been called.
#' 
#' @param data data.frame as returned by bb_readHic()
#' @param y numeric, position of heatmap on y-axis
#' @param annotation character, string to use to label the heatmap in the top left
#' @param params object of class bb_params, containing shared params for the plot
#' including, at minimum, chrom, chromstart, chromend, and assembly
bb_placeHeatmapAndLegend <- function(counts, y, x = 2, annotation, params) {
     
    heatmap <- bb_plotHicTriangle(data = counts, params = params, x = x, y = y, width = 3, height = 1.5, just = "top",
                                  default.units = "inches")
    
    # annotate heatmap legend
    bb_annoHeatmapLegend(plot = heatmap, x = x + 1.5, y = y, width = 0.13, height = 1.2,
                         just = c("right", "top"))
    
    # plot label
    bb_plotText(annotation, x = x - 1.5, y = y, just = c("left", "top"))
    
}

bb_placeSignalAndLabel <- function(data, y, x, annotation = NULL, params, width = 3, height = 0.45, color = "navy", fontsize = 6, ...) {
    
    bb_plotSignal(data = data, params = params,
                  x = x, y = y, width = width, height = height,
                  just = c("left", "top"), default.units = "inches", linecolor = color, fill = color,
                  scale = TRUE, ...)
    
    if (!is.null(annotation)) bb_plotText(label = annotation, fonsize = fontsize, fontcolor = "black",
                                          x = x + 0.5*width, y = y, just = c("left", "top"), default.units = "inches")
    
}