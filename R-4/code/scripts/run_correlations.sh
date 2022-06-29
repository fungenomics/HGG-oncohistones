#!/usr/bin/env bash

#SBATCH --job-name="correlations"
#SBATCH --account=rrg-kleinman
#SBATCH --ntasks=1
#SBATCH --output=logs/correlations-%j.out
#SBATCH --time=00:20:00
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1

# -------------------------
# PREPARATION
# -------------------------

# Prep R environment
module load nixpkgs/16.09
module load gcc/7.3.0
module load r/4.0.0
export R_LIBS_USER="../../../../../renv/library/R-4.0/x86_64-pc-linux-gnu"

# -------------------------
# RUN
# -------------------------

Rscript ../../../../../../code/scripts/predict_celltype_cor.R

