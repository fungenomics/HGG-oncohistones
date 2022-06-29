# Helpers for internal testing, mainly using the ensurer package

require(ensurer)

# General ----------------------------------------------------------------------

find_missing_files <- function(paths) {
  
  exists <- file.exists(paths)
  
  if (all(exists)) message("All files exist")
  else {
    warning("Missing files!")
    return(paths[!exists])
  }
  
}

find_missing_files_df <- function(df, path_col) {
  
  exists <- file.exists(df[[path_col]])
  
  if (all(exists)) message("All files exist")
  else {
    warning("Missing files!")
    return(df[!exists,])
  }
  
}


# ensurer contracts ------------------------------------------------------------
# Uses the {ensurer} package: https://github.com/smbache/ensurer
# Define some contracts that can be reused

# ensure all inputs = FALSE
ensure_none <- ensures_that(all(map_lgl(., isFALSE)))

# ensure no inputs are NA
ensure_present <- ensures_that(all(. != "NA"))

# ensure all inputs (paths) exist
ensure_exists <- ensures_that(all(file.exists(.)))
