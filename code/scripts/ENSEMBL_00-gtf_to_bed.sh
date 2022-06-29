#!/bin/bash

# Adapted from Nisha Kabir & Nicolas De Jay
# Purpose:
# Convert GTF used for gene annotaiton in the bulk pipeline (Ensemble.ensGene.whole.gtf) into a bed file
# use fields 1 (seqName (chr)), 4 (start), 5 (end), 7 (strand), and 9 (attribute)
# these fields become 1 (chr), 2 (start), 3 (end), 4 (attribute), 5 - score of 0, and 6 (strand)

GTF="<pipeline>/code/refs/genomic_features/hg19/Ensembl.ensGene/Ensembl.ensGene.whole.gtf"

cut -f1,4,5,7,9 $GTF | \
	sed 's/id.*ensembl_transcript "\(ENST[0-9]*\).*ensembl_gene "\(ENSG[0-9]*\).*gene_symbol "\(.*\)"/\1:\2:\3/' | 	
	awk -F "\t" 'BEGIN{OFS="\t";} {print $1,$2,$3,$5,"0",$4}' |
	sort > Ensembl.ensGene.hg19.whole.bed
