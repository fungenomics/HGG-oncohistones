# Purpose: save TPM counts for bulk RNAseq cohort.
# Adapted from https://fungenomics.github.io/HGG-oncohistones/code/02-bulk_comparisons.html#prepare-tables

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(glue)
library(magrittr)

source(here("rr_helpers.R"))

out <- here("output/02")

# get sample order as Table 1
sample_order <- data.table::fread(here("data/metadata/bulk_info.samples.tsv"), data.table = FALSE) %>% .$Nickname

counts_RNA <- read.table(here("data/RNAseq/pipeline_l3/all_TPM/counts/Ensembl.ensGene.exon.tpm.tsv.gz"),
                         header = T, sep = "\t", check.names = FALSE) %>%
    tibble::rownames_to_column(var = "ID") %>%
    separate(ID, into = c("ENSID", "symbol"), sep = ":") %>%
    arrange(ENSID) %>%
    .[, c("ENSID", "symbol", sample_order)] %T>%
    rr_write_tsv(glue("{out}/TABLE_bulk_counts.TPM.tsv"),
                 "TPM counts for bulk RNAseq data")

# sanity checks
length(colnames(counts_RNA)) - 2

dim(counts_RNA)
