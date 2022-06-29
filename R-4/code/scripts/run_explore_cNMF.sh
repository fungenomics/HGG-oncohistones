#!/usr/bin/env bash

#SBATCH --job-name="explorecNMF"
#SBATCH --account=rrg-kleinman
#SBATCH --ntasks=1
#SBATCH --output=logs/explorecNMF-%j.out
#SBATCH --time=00:30:00
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

RMD="/lustre03/project/6004736/sjessa/from_beluga/HGG-oncohistones/final/R-4.0.0/code/scripts/explore_cNMF.Rmd"
WORKDIR=`pwd`
R --no-save -e \
  "rmarkdown::render(
    input = '$RMD',
    output_format = 'html_document',
    output_file = 'explore_cNMF.html',
    knit_root_dir = '$WORKDIR',
    output_dir = 'output_ngenes2000_niter100_malignant',
    params = list(output = 'output_ngenes2000_niter100_malignant', suffix = 'malignant'))"
