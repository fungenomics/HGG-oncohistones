# Following the data structures vignette at
# https://github.com/timoast/signac/blob/master/vignettes/data_structures.Rmd

library(motifmatchr)
library(TFBSTools)
library(JASPAR2020)

# get DNA sequence motif info
# it looks like this needs to be done only once, so could be saved?
pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species = "Homo sapiens", all_versions = FALSE)
)

save(pfm, file = here("R-4.0.0/data/scMultiome/references/pfm_JASPAR2020.Rda"))
