# 2022-03-01
# This script loads domain calls from H3K27me2, and
# 1) filters out segments with a score below a threshold
# 2) merges segments which are separated by less than X bp

library(here)
library(GenomicRanges)
library(rtracklayer)
library(purrr)

# config
threshold     <- 1.2     # segment score
collapse_dist <- 100     # base pairs

bed_files <- list.files(here("data/ChIPseq/domains_@mhulswit"), pattern = ".bed", full.names = TRUE)

clean_domain_calls <- function(path) {
    
    message("@@ ", basename(path))
    
    # load domain
    domains <- rtracklayer::import(path,  extraCols = c(segment = "character", value = "numeric"))
    
    # filter to segments above threshold
    domains_filt <- domains[mcols(domains)$value > threshold, ]
    
    # merge segments which are separated by less than <collapse_dist> bp
    domains_filt_merged <- GenomicRanges::reduce(domains_filt, min.gapwidth = collapse_dist)
    
    # write new BED file
    out <- gsub(".bed", ".filt.merged.bed", gsub("mhulswit/", "mhulswit/filtered/", path))
    rtracklayer::export(domains_filt_merged, out)
    
}

walk(bed_files, clean_domain_calls)

message("@ done.")
