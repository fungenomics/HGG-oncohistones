#!/usr/bin/Rscript

# Adapted from Nisha Kabir
# Purpose:
# Collapse features from transcript level to gene level; keeping the longest interval

# libraries
library(tidyr)
library(dplyr)

transcripts <- read.delim("Ensembl.ensGene.whole.hg19.bed", header = FALSE, sep = "\t") %>% as.data.frame()

colnames(transcripts) <- c("chr", "start", "end", "attribute", "score", "strand")

transcripts <- transcripts %>%
	separate(attribute, into = c("enst", "ensg", "gene_symbol"), sep = ":") %>%
	# For each gene, defined by Ensembl gene ID get the longest interval
	group_by(ensg) %>%
	mutate(start = min(start)) %>%
	mutate(end = max(end)) %>%
	group_by(chr, start, end, ensg, gene_symbol, score, strand) %>%
	summarise() %>%
	unite(attribute, c("ensg", "gene_symbol"), sep = ":")

write.table(x = transcripts, file = "Ensembl.ensGene.hg19.collapsed.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

