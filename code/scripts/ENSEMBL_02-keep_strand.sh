#!/bin/bash

# Adapted from Nisha Kabir and Steven Hebert
# Purpose:
# Get bedtools to keep strand info (as prep for homer which requires strand info)

awk '{if($6=="+"){ $4=$4":+"} else{$4=$4":-"} print $1"\t"$2"\t"$3"\t"$4}' Ensembl.ensGene.hg19.collapsed.bed > \
    Ensembl.ensGene.hg19.collapsed.strand.bed
