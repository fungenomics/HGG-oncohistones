#!/usr/bin/env bash

#SBATCH --job-name="clean"
#SBATCH --account=rrg-kleinman
#SBATCH --ntasks=1
#SBATCH --output=../../logs/clean-%j.out
#SBATCH --time=00:30:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

# for i in ../../data/scMultiome/pipeline_10X_Multiome/*; do id=`readlink -f $i`; echo "sbatch --export=ALL;DIR=$id run_prep_object_scMultiome.sh"; done


# -------------------------
# PREPARATION
# -------------------------

# Load requirements
module load StdEnv/2020
module load r/4.1.2
export R_LIBS_USER="../../../../renv/library/R-4.1/x86_64-pc-linux-gnu"

# -------------------------
# RUN
# -------------------------

Rscript ../../../../code/scripts/prepare_clean_object_scMultiome.R

