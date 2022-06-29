#!/usr/bin/Rscript

# Adapted from Nisha Kabir
# Purpose:
# Collapse features from transcript level to gene level; keeping the longest interval

# libraries
library(tidyr)
library(dplyr)

transcripts <- read.delim("Ensembl.ensGene.mm10.whole.bed", header = FALSE, sep = "\t") %>% as.data.frame()

colnames(transcripts) <- c("chr", "start", "end", "attribute", "score", "strand")

transcripts <- transcripts %>%
	separate(attribute, into = c("id", "gene_id", "ensembl_transcript", "ensembl_gene", "gene_symbol"), sep = "; ") %>%
    separate(gene_symbol, into = c("drop", "gene_symbol"), sep = " ") %>% 
    separate(ensembl_gene, into = c("drop", "ensembl_gene"), sep = " ") %>% 
	# For each gene, defined by Ensembl gene ID get the longest interval
	group_by(ensembl_gene) %>%
	mutate(start = min(start)) %>%
	mutate(end = max(end)) %>%
	group_by(chr, start, end, ensembl_gene, gene_symbol, score, strand) %>%
	summarise() %>%
	unite(attribute, c("ensembl_gene", "gene_symbol"), sep = ":")

write.table(x = transcripts, file = "Ensembl.ensGene.mm10.collapsed.bed", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

