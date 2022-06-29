library(RColorBrewer)
library(viridis)
library(magrittr)
library(ggplot2)

# General palettes -------------------------------------------------------------

# yellow -> red
ylrd    <- grDevices::colorRampPalette(brewer.pal(8, "OrRd"))(n = 100)

# blue -> red
rdbu    <- rev(grDevices::colorRampPalette(brewer.pal(8, "RdBu"))(n = 100))
rdbu2   <- colorRampPalette(c("navy", "white", "red"))(n = 100)
rdbu3   <- colorRampPalette(c("navy", "#2E2E97", "#4848A4", "white",
                              "#FF4E4E", "#FF3434", "red"))(n = 100)

# white -> purple
purples <- colorRampPalette(brewer.pal(9, "Purples"))(n = 100)

# gray -> red
gryrd   <- grDevices::colorRampPalette(c("gray83", "#E09797", "red"))(n = 200)

# modification on the viridis:: magma palette
# from Kinker et al (https://github.com/gabrielakinker/CCLE_heterogeneity/blob/master/custom_magma.R)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(viridis::magma(323, begin = 0.18)))
custom_viridis <- c(colorRampPalette(c("white", rev(viridis(323, begin = 0.15))[1]))(10), rev(viridis::viridis(323, begin = 0.18)))

# helper to show them
display_palette <- function(palette) {
    
    image(1:length(palette), 1, z = as.matrix(1:length(palette)), col = palette,
          xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
    axis(1, at = 1:length(palette), labels = names(palette))
    
}

# Custom palettes --------------------------------------------------------------

# borrowed from ArchR package
palette_region <- c("Distal"   = "#60BA64",
                    "Exonic"   = "#73C6FF",
                    "Intronic" = "#620FA3",
                    "Promoter" = "#FFC554")

# growth factor receptor alterations
palette_gfr <-c("ACVR1"          = "firebrick3",
                "ACVR1-G328E"    = "firebrick3",
                "ACVR1-G328V"    = "firebrick1",
                "ACVR1-G356D"    = "firebrick4",
                "ACVR1-R206H"    = "mediumvioletred",
                "PDGFRA"         = "darkorchid3",
                "PDGFRA mut"     = "darkorchid3",
                "PDGFRA mut&amp" = "darkorchid4",
                "PDGFRA-FIP1L1"  = "darkorchid1",
                "EGFR-L62R"      = "darksalmon",
                "BRAF"           = "deeppink4",
                "BRAF-V600E"     = "deeppink4",
                "None"           = "gray90")


palette_acvr1_simple <- c("WT" = 1, "Mutant" = 16)

palette_acvr1 <- c("ACVR1 mutant" = "orange", "ACVR1 KO" = "blue")

palette_tp53 <- c("TP53 mut/loss" = "navy", "NA" = "gray90")

# single-cell technology
palette_sc_tech <- c("10X Single Nuclei 3'" = "gray50",
                     "10X Single Cell 3'"   = "gray30",
                     "10X ATAC"             = "gray30",
                     "10X Multiome"         = "gray50",
                     "Smart-seq2"           = "red")

palette_sc_tech2 <- c("scRNA" = "#F8766D",
                      "scATAC" = "#00BA38",
                      "scMultiome_Singulator" = "#619CFF",
                      "scMultiome_OctoMACS" = "#1956bd",
                      "scMultiome_manual" = "#88a7db")


# molecular status
palette_molecular <- c("EZHIP"      = "darkorchid4",
                       "H3.1K27M"   = "orange",
                       "H3.2K27M"   = "lightpink",
                       "H3.1/2K27M" = "orange",
                       "H3.3K27M"   = "red3",
                       "H3.3G34R"   = "turquoise3",
                       "H3.3G34V"   = "turquoise4",
                       "H3.3G34R/V" = "turquoise4",
                       "H3WT"       = "gray50",
                       "WT"         = "gray50")


# location
palette_location <- c("PF"      = "#f7cebb",
                      "4V"      = "goldenrod3",
                      "Pons"    = "gold1",
                      "Thal."   = "darkolivegreen4",
                      "Spinal"  = "darkred",
                      "Cortex"  = "dodgerblue2")

# molecular groups
palette_groups <- c(
    "PFA-EZHIP-PF"           = "darkorchid4",
    "PFA-H3.1K27M-PF"        = "mediumorchid",
    "HGG-H3.1/2K27M-Pons"    = "orange",
    "HGG-H3.3K27M-4V"        = "#c96969",
    "HGG-H3.3K27M-Spinal"    = "#ff4747",
    "HGG-H3.3K27M-Pons"      = "red3",
    "HGG-H3.3K27M-Thal."     = "red4",
    "HGG-H3.3G34R/V-Cortex"  = "turquoise4",
    "HGG-H3WT-Cortex"        = "gray70",
    "Other"                  = "gray90")


palette_groups_simple <- c(
    "PFA-EZHIP"      = "darkorchid4",
    "PFA-H3.1K27M"   = "plum2",
    "HGG-H3.1/2K27M" = "orange",
    "HGG-H3.3K27M"   = "red3")

palette_genotype <- c("H3.1K27M"     = "orange",
                      "H3.3K27M"     = "red3",
                      "KO"           = "blue",
                      "H3.1K27M-KO"  = "blue",
                      "H3.3K27M-KO"  = "blue",
                      "EZHIP-PFA"    = "darkorchid4",
                      "H3.1K27M-PFA" = "mediumorchid",
                      "H3WT"         = "gray50")

palette_condition <- c("H3 WT + H3.1WT OE"         = "blue",
                       "H3.1WT OE"                 = "blue",
                       "H3 WT + H3.1K27M OE"       = "orange",
                       "H3 WT + H3.3K27M OE"       = "red3",
                       "H3 WT + H3.3K27R OE"       = "blue",
                       "H3.3K27R OE"               = "blue",
                       "H3.1K27M"                  = "orange",
                       "H3.1K27M OE"               = "orange",
                       "H3.1K27M KO + H3.3WT OE"   = "blue",
                       "H3.1K27M KO + H3.3K27M OE" = "red3",
                       "H3.3K27M"                  = "red3",
                       "H3.3K27M OE"               = "red3",
                       "H3.3K27M KO + H3.1WT OE"   = "blue",
                       "H3.3K27M KO + H3.1K27M OE" = "orange",
                       "ACVR1 mutant"              = "firebrick3",
                       "ACVR1 KO"                  = "gray50",
                       "H3.1K27M + ACVR1 mut"      = "firebrick3",
                       "H3.1K27M + ACVR1 KO"       = "gray50")

# shape palette
palette_cell_line <- c("DIPG21" = 21, "DIPG36" = 22, "DIPGIV" = 24, "BT245" = 21, "DIPGXIII" = 22, "HSJ019" = 24)

palette_cl_acvr1 <- c("ACVR1-R206H-KO" = "blue",
                      "ACVR1-R206H"    = "red3",
                      "ACVR1-G328E-KO" = "blue",
                      "ACVR1-G328E"    = "orange",
                      "ACVR1-G328V-KO" = "blue",
                      "ACVR1-G328V"    = "orange")
# Hox genes
# excuse the non-code-friendly spacing to match paralogs across clusters...
hoxa <- c("HOXA1", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7",          "HOXA9", "HOXA10", "HOXA11",           "HOXA13")
hoxb <- c("HOXB1", "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "HOXB7", "HOXB8", "HOXB9",                               "HOXB13")
hoxc <- c(                           "HOXC4", "HOXC5", "HOXC6",          "HOXC8", "HOXC9", "HOXC10", "HOXC11", "HOXC12", "HOXC13")
hoxd <- c("HOXD1",          "HOXD3", "HOXD4",                            "HOXD8", "HOXD9", "HOXD10", "HOXD11", "HOXD12", "HOXD13")

pal_hoxa <- c(rep("#15048E", 2), rep("#7B0FED", 3), rep("#DBAAF2", 5), rep("#CCCCCC", 1)) %>% set_names(hoxa)
pal_hoxb <- c(rep("#15048E", 2), rep("#7B0FED", 3), rep("#DBAAF2", 4), rep("#CCCCCC", 1)) %>% set_names(hoxb)
pal_hoxc <- c(rep("#7B0FED", 2), rep("#DBAAF2", 5), rep("#CCCCCC", 2)) %>% set_names(hoxc)
pal_hoxd <- c(rep("#CCCCCC", 1), rep("#7B0FED", 2), rep("#DBAAF2", 4), rep("#CCCCCC", 2)) %>% set_names(hoxd)
palette_hox <- c(pal_hoxa, pal_hoxb, pal_hoxc, pal_hoxd)

# track types
palette_tracks <- c(
    "RNAseq"          = "black",
    "Pol2"            = "black",
    "H3K27ac"         = "blue",
    "H3K4me3"         = "blue",
    "H3K27me3"        = "red",
    "EZH2"            = "red",
    "H2AK119Ub"       = "red",
    "CTCF"            = "gray20",
    "RING1B"          = "gray40",
    "SUZ12"           = "gray40",
    "scATACseq"       = "gray40",
    "scMultiome_ATAC" = "gray40",
    "scMultiome_RNA"  = "black")

# cell types
palette_type <- c("RGC"                  = "#ffcc00",
                  "Glial progenitors"    = "#d5d98b",
                  "OPC"                  = "#e0de53",
                  "Proliferating OPC"    = "#e6f957",
                  "Oligodendrocytes"     = "#b4e04e",
                  "Astrocytes"           = "#00a385",
                  "Ependymal"            = "#8ee5cf",
                  "Neuronal progenitors" = "#ffbda3",
                  "Neurons"              = "#135ca0",
                  "Immune"               = "gray50",
                  "Vascular & other"     = "gray70",
                  "Normal"               = "gray90")

if (grepl("oncohistones", here::here())) {
    
    load(here("data/scRNAseq/references/mouse_atlas_extended/palette_joint_mouse_extended.Rda"))
    load(here("data/scRNAseq/references/mouse_atlas_extended/palette_joint_mouse_extended_full.Rda"))    
}


# mouse timepoints
palette_timepoint <- c(rep(RColorBrewer::brewer.pal(7, "YlGnBu")[2:7], each = 2),
                       RColorBrewer::brewer.pal(4, "Reds")[2:4]) %>%
    set_names(c("E10.5", "E10", "E12.5", "E12", "E13.5", "E13", "E15.5", "E15",
                "E16.5", "E16", "E18.5", "E18", "P0", "P3", "P6"))

palette_timepoint <- c(palette_timepoint, "E12.5, E15.5" = palette_timepoint["E12.5"])

# human timepoints
palette_timepoint_human <- c(RColorBrewer::brewer.pal(9, "Blues")[3:9],
                             RColorBrewer::brewer.pal(7, "Reds")[2:7]) %>% 
    set_names(c("CS12", "CS13", "CS14", "CS15", "CS19", "CS20", "CS22", "GW14", "GW18", "GW19", "GW20", "GW22", "GW25"))


# Helper functions for aesthetics ----------------------------------------------

#' Apply a clean theme to a ggplot2 object
#'
#' @references https://github.com/sjessa/ggmin
theme_min <- function(base_size = 11, base_family = "",
                      border_colour = "grey90",
                      border_size = 1) {
    
    theme_light(base_size = 11, base_family = "") +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = border_colour, size = border_size),
            axis.ticks = element_line(colour = border_colour),
            strip.background = element_rect(fill = NA, colour = NA),
            strip.text.x = element_text(colour = "black", size = rel(1.2)),
            strip.text.y = element_text(colour = "black", size = rel(1.2)),
            title = element_text(size = rel(0.9)),
            axis.text = element_text(colour = "black", size = rel(0.8)),
            axis.title = element_text(colour = "black", size = rel(0.9)),
            legend.title = element_text(colour = "black", size = rel(0.9)),
            legend.key.size = unit(0.9, "lines"),
            legend.text = element_text(size = rel(0.7), colour = "black"),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.background = element_rect(colour = NA, fill = NA)
        )
}


#' A theme for when we just want to plot a color bar as a row
theme_row <- function() {
    
    theme(panel.border = element_blank(),
          axis.text.y  = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x  = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks   = element_blank())
    
}


#' Rotate the x axis labels in a ggplot
#'
#' @param angle Integer, value in degrees to rotate labels. Default: 90.
#'
#' @return A theme element to rotate labels
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(colour = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + rotateX()
rotate_x <- function(angle = 90) {
    
    theme(axis.text.x = element_text(angle = angle, hjust = 1))
    
}




#' Remove the legend in a ggplot
#'
#' @return A theme element to hide legend
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(colour = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + noLegend()
no_legend <- function() {
    
    theme(legend.position = "none")
    
}

bottom_legend <- function() {
    
    theme(legend.position = "bottom")
    
}

pal_qual1 <- c("#646BA8",
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

pal_qual2 <- c("red4",
               "darkslategray3",
               "dodgerblue1",
               "darkcyan",
               "gray79",
               "black", 
               "skyblue2",
               "dodgerblue4",
               "purple4",
               "maroon",
               "chocolate1",
               "bisque3",
               "bisque",
               "seagreen4",
               "lightgreen",
               "skyblue4",
               "mediumpurple3",
               "palevioletred1",
               "lightsalmon4",
               "darkgoldenrod1")

set_qual_pal <- function(n, set = 1) {
    
    if (set == 1) pal_qual <- pal_qual1
    else if (set == 2) pal_qual <- pal_qual2
    
    if (n <= length(pal_qual)) { pal_qual_ramped <- head(pal_qual, n)}
    else {
        
        pal_qual_ramped <- colorRampPalette(pal_qual)(n = n)
        
    }
    
    return(pal_qual_ramped)
    
}



#' Remove axis ticks and tick labels from a ggplot
#'
#' @return A theme element to remove ticks
no_ticks <- function() {
    
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
    
}
