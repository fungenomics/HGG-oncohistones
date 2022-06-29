
# plotgardener helpers -------------------------------------------------------------
# these are based on the helper functions at code/functions/BentoBox_helpers.R,
# updated following the package overhaul from BentoBox --> plotgardener


pg_placeSignalAndLabel <- function(data, y, x, chrom = NULL, annotation = NULL, params, width = 3, height = 0.45, color = "navy", fontsize = 6, ...) {
    
    plotSignal(data = data, chrom = chrom, params = params,
                  x = x, y = y, width = width, height = height,
                  just = c("left", "top"), default.units = "inches", linecolor = color, fill = color,
                  scale = TRUE, ...)
    
    if (!is.null(annotation)) plotText(label = annotation, fonsize = fontsize, fontcolor = "black",
                                          x = x + 0.5*width, y = y, just = c("left", "top"), default.units = "inches")
    
}

pg_placeSignalAndLabel2 <- function(data, y, x, annotation = NULL, params, width = 3, height = 0.45, color = "navy", fontsize = 6, ...) {
    
    plotSignal(data = data, params = params,
               x = x, y = y, width = width, height = height,
               just = c("left", "top"), default.units = "inches", linecolor = color, fill = color,
               scale = TRUE, ...)
    
    if (!is.null(annotation)) plotText(label = annotation, fonsize = fontsize, fontcolor = "black",
                                       x = x + 0.5*width, y = y, just = c("left", "top"), default.units = "inches")
    
}
