#!/usr/bin/env bash

#SBATCH --job-name="inferCNV"
#SBATCH --account=rrg-kleinman
#SBATCH --ntasks=1
#SBATCH --output=logs/infercnv-%j.out
#SBATCH --time=05:30:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4

module load nixpkgs/16.09
module load gcc/7.3.0
module load r/4.0.0
export R_LIBS_USER="../../../../../renv/library/R-4.0/x86_64-pc-linux-gnu"

# Run script to call CNVs and generate custom heatmap
Rscript ../../../../../../code/scripts/inferCNV.R
