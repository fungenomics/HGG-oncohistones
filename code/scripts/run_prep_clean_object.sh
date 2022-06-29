#!/usr/bin/env bash

#SBATCH --job-name="clean"
#SBATCH --account=rrg-kleinman
#SBATCH --ntasks=1
#SBATCH --output=.clean-%j.out
#SBATCH --time=00:15:00
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1


# -------------------------
# PREPARATION
# -------------------------

# Prep R environment
module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.1
export R_LIBS_USER="../../../../renv/library/R-3.6/x86_64-pc-linux-gnu"

# -------------------------
# RUN
# -------------------------

Rscript ../../../../code/scripts/prepare_clean_object_scRNA.R

