
# Collection of functions to work with the output of deepTools
# These functions are heavily based on code from Nicolas De Jay
# .../nicolas.dejay/from_hydra/2021/210928-stacked-heatmap-npz-manipulation-scripts

require(readr)
require(jsonlite)
require(data.table)
require(magrittr)
require(purrr)
require(R.utils)


#' Parse a .npz file output by deeptools computeMatrix
#' 
#' @param input_npz Character, path to .npz file from deepTools computeMatrix
parse_npz <- function(input_npz) {
    
    message("@ Loading npz header...")
    input_npz_header <- readr::read_lines(gzfile(input_npz), n_max = 1L)
    
    input_npz_header_parsed <- input_npz_header %>%
        # Remove trailing @
        sub("^@", '', .) %>%
        # Parse as JSON
        jsonlite::fromJSON()
    
    message("@ Loading npz body [slow step]...")
    
    # The body of the file contains abundances of the mark
    input_npz <- readr::read_tsv(gzfile(input_npz),
                                 comment = "@",
                                 col_names = FALSE)
    
    # Prevent feature starts and ends from being stored in
    # scientific notation.
    input_npz$X2 %<>% as.character()
    input_npz$X3 %<>% as.character()
    
    return(list("body" = input_npz,
                "header" = input_npz_header_parsed))
    
}


get_sample_indices <- function(npz) {
    
    message("@ Getting sample indices...")
    
    # Calculate the start and end column index for each sample.  Each interval
    # defines a fixed number of bins associated with each sample, over which the
    # abundance of a mark was averaged (presumably).
    sample_n <- length(npz$header$sample_boundaries) - 1
    
    sample_start_indices <-
        (npz$header$sample_boundaries + 7L) %>% head(n = -1L)
    
    sample_end_indices <-
        (npz$header$sample_boundaries + 6L) %>% tail(n = -1L)
    
    return(list("n"     = sample_n,
                "start" = sample_start_indices,
                "end"   = sample_end_indices))
    
}


split_npz_by_sample <- function(npz, keep_info_cols = TRUE) {
    
    message("@ Splitting data...")
    
    indices <- get_sample_indices(npz)
    
    # Subset to indices for each sample
    if (keep_info_cols) purrr::map(1:indices$n, ~ npz$body[, c(1:6, indices$start[[.x]]:indices$end[[.x]])])
    else purrr::map(1:indices$n, ~ npz$body[, indices$start[[.x]]:indices$end[[.x]]])
    
}


aggregate_npz <- function(npz) {
    
    message("@ Aggregating data...")
    
    indices <- get_sample_indices(npz)
    
    output_npz <- npz$body
    
    colSE <- function(x) {
        
        apply(x, 2, sd, na.rm = TRUE) / 
            nrow(x)
        
    }
    
    # Calculate average abundance across features for each bin per sample
    sample_averages <- purrr::map_dfc(1:indices$n,
                                      ~ npz$body[, indices$start[[.x]]:indices$end[[.x]]] %>%
                                          # Calculate the mean abundance over all bins
                                          colMeans(na.rm = TRUE)
    ) %>% 
        set_colnames(paste0("mean_", as.character(1:indices$n))) %>% 
        tibble::rowid_to_column(var = "bin")
    
    sample_se <- purrr::map_dfc(1:indices$n,
                                ~ npz$body[, indices$start[[.x]]:indices$end[[.x]]] %>%
                                    # Calculate the mean abundance over all bins
                                    colSE()
    ) %>% 
        set_colnames(paste0("SE_", as.character(1:indices$n)))
    
    return(bind_cols(sample_averages, sample_se))
    
}


#' Rescale a matrix to [0, 1]
#' 
#' @param x numeric matrix
scale_mat <- function(x){
    
    y <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    return(y)
}


#' Reformat columns of manipulated npz to match requirements for deeptools
#' 
#' @param x numeric vector
reformat_col <- function(x) {
    
    x %>% 
        # The original table only had 6 decimal points
        round(digits = 6) %>%
        as.character() %>%
        # Numpy (deeptools) requires "nan" instead of "NaN"
        dplyr::recode("NaN" = "nan")
    
}


#' Rescale npz values to between [0, 1] for each sample
#'
#' @param npz List as returned by \code{parse_npz()}
#'
#' @value Data frame, new \code{npz$body}
rescale_npz <- function(npz) {
    
    message("@ Preparing data for rescaling...")

    # split the npz sample values into one dataframe per sample
    npz_split <- split_npz_by_sample(npz, keep_info_cols = FALSE)
    
    message("@ Rescaling data to [0, 1] [slow step]...")
    
    # rescale values in each sample to [0, 1]
    npz_split_scaled <- map(npz_split, ~ .x %>% as.data.table %>% scale_mat)
    
    message("@ Reformatting data...")
    
    # reformat columns
    npz_split_scaled_reformatted <- map(npz_split_scaled,
                                        ~ apply(.x, MARGIN = 2, reformat_col))
    
    # generate output dataframe
    # (bind info cols with all of the reformatted data)
    output_npz <- dplyr::bind_cols(c(
        list(as.data.table(npz$body[, 1:6, drop = FALSE])),
        npz_split_scaled_reformatted
        ))
    
    output_npz <- as.data.frame(output_npz)
    
    return(output_npz)
    
}



write_npz <- function(npz_header, npz_body, output) {
    
    message("@ Preparing output file...")
    
    output_npz_header_parsed <- npz$header
    
    # R encodes scalars as one-element vectors, whereas numpy (deeptools) expects
    # bonafide scalars for most variables, except for the likes of "group_labels".
    # Similarly, it expects null instead of the empty hash {}.
    variables <- c("group_labels")
    for (i in variables)
        if (length(output_npz_header_parsed[[i]]) == 1)
            output_npz_header_parsed[[i]] %<>% list(.)
    
    output_npz_header <-
        output_npz_header_parsed %>%
        jsonlite::toJSON(auto_unbox = TRUE, null = "null") %>%
        paste0("@", .) %>%
        as.character() %>%
        data.frame()
    
    message("@ Saving data...")
    
    # write the header
    write.table(x            = output_npz_header,
                file         = paste0(output, ".raw"),
                append       = FALSE,
                quote        = FALSE,
                col.names    = FALSE,
                row.names    = FALSE)
    
    # append the values
    write.table(x            = npz_body,
                file         = paste0(output, ".raw"),
                append       = TRUE,
                quote        = FALSE,
                col.names    = FALSE,
                row.names    = FALSE,
                sep          = "\t")
    
    R.utils::gzip(paste0(output, ".raw"), destname = output, overwrite = TRUE, remove = TRUE)
    
}