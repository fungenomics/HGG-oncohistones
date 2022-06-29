#!/usr/bin/env bash

#SBATCH --job-name="consensus"
#SBATCH --account=rrg-kleinman
#SBATCH --ntasks=1
#SBATCH --output=../../logs/consensus-%j.out
#SBATCH --time=00:20:00
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1

# for i in ../../data/scRNAseq/pipeline_10X/*; do id=`readlink -f $i`; echo "sbatch --export=ALL;DIR=$id run_add_consensus.sh"; done


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

cd $DIR

Rscript ../../../../code/scripts/add_consensus_projections.R

