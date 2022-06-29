#!/bin/bash

# Adapted from Nisha Kabir
# Get the promoter coordinates for each gene, by taking the TSS +/- 2.5kbp
cat Ensembl.ensGene.hg19.collapsed.bed | \
	awk -F "\t" 'BEGIN{OFS="\t";} {if($6 == "+") {print $1,$2-2500,$2+2500,$4,$5,$6} else print $1,$3-2500,$3+2500,$4,$5,$6}' > \
	Ensembl.ensGene.hg19.collapsed.promoter5kb.bed

# Next, make a copy for autosomal chromosomes only (exclude random/M/X/Y)
# and ensure the regions to not exceed the chromsome bounds to avoid
# deepTools errors
CHR_BOUNDS="hg19_chr1_22.bed"

module load r/3.6.1
export R_LIBS_USER="../../renv/library/R-3.6/x86_64-pc-linux-gnu"

# 1. split promoters according to chromosome
# 2. for promoters on each chromosome, if the end
#    of the promoter exceeds the chromosome bound, set the
#    end of the promoter to be the end of the chromosome
# 3. put the promoters back together
R --vanilla -e "
  library(GenomicRanges)
  library(rtracklayer)
  promoters <- import('Ensembl.ensGene.hg19.collapsed.promoter5kb.bed')
  chr_bounds <- import('${CHR_BOUNDS}')
  promoters_split <- split(promoters, seqnames(promoters))
  promoters_split <- promoters_split[seqnames(chr_bounds)]
  chr_bounds <- split(chr_bounds, seqnames(chr_bounds))
  for (i in 1:length(chr_bounds)) print(i); end(promoters_split[[i]])[ end(promoters_split[[i]]) >= end(chr_bounds[[i]]) ] <- end(chr_bounds[[i]]) - 1
  promoters_fixed <- unlist(as(promoters_split, 'GRangesList'))
  export(object = promoters_fixed, con = 'Ensembl.ensGene.hg19.collapsed.promoter5kb.bounded.bed')
"
