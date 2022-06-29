# Author: Selin Jessa
# Date: 2021-05-21
#
# Purpose:
# 1. Modify the feature annotation used for RNAseq counts to quantify
# expression at each HOX transcript instead of each HOX gene.
#
# For this, we'll use a trick - since we want to collapse all transcripts by
# gene, *except* for the HOX transcripts, we'll replace the gene ID by the transcript
# ID for HOX transcripts only. Therefore each HOX transcript will be treated separately.
#
# 2. Extract GTF with all HOX transcripts for downstream analysis

library(here)
library(rtracklayer)

# load reference from the in-house bulk RNAseq pipeline
ensembl_gtf <- import("<pipeline>/code/refs/genomic_features/hg19/Ensembl.ensGene/Ensembl.ensGene.exon.gtf")

# get idx of HOX transcripts
hox_idx <- grepl(":HOX", ensembl_gtf$gene_id)

# for HOX transcripts, replace the gene id with the transcript ID so they are not collapsed
ensembl_gtf[hox_idx, ]$gene_id <- paste0(ensembl_gtf[hox_idx, ]$ensembl_transcript, ":", ensembl_gtf[hox_idx, ]$gene_symbol)

# sanity check it worked
ensembl_gtf[hox_idx, ]$gene_id

# save
rtracklayer::export(ensembl_gtf, here("data/RNAseq/references/Ensembl.ensGene.exon.plus_HOX_tx.gtf"))

# finally, on the command line, grep for HOX genes to generate a GTF specifically for the HOX transcripts
# $ grep ":HOX" data/RNAseq/references/Ensembl.ensGene.exon.plus_HOX_tx.gtf > \
#        final/data/RNAseq/references/Ensembl.ensGene.exon.HOX.gtf