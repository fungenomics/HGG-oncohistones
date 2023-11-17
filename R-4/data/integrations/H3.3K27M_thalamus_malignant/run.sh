#!/usr/bin/bash
#SBATCH --job-name="harmony"
#SBATCH --time=02:00:00
#SBATCH --output=logs/integrate_harmony-%j.out
#SBATCH --account=rrg-kleinman
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G

mkdir -p logs/

module load r/4.0.0

# use the project specific library
export R_LIBS_USER="../../../renv/library/R-4.0/x86_64-pc-linux-gnu"

# Rscript ../../../../code/scripts/integrate_harmony.R
Rscript ../../../../code/scripts/update_integration_metadata.R
