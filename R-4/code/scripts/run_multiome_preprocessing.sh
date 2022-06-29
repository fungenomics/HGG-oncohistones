#!/usr/bin/bash
#SBATCH --time=01:30:00
#SBATCH --output=logs/preprocessing.Rmd-%j.out
#SBATCH --account=rrg-kleinman
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G

# Load requirements
module load StdEnv/2020
module load r/4.1.2

# Load python & MACS2 for peak calling
export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use $MUGQIC_INSTALL_HOME/modulefiles
module load mugqic/python/3.7.3
module load mugqic/MACS2/2.2.7.1

# Set library
export R_LIBS_USER="../../../../../renv/library/R-4.1/x86_64-pc-linux-gnu"

# --------------------------
# Customizable parameters
# --------------------------

# Modify this to use a custom config file, which should be saved in this directory
# e.g. CONFIG="custom.preprocessing.config.tsv"
CONFIG="scMultiome_preprocessing.config.tsv"

# --------------------------
# Non-customizable parameters
# --------------------------

# R markdown file to render
RMD="../../../../../code/scripts/preprocessing_scMultiome.Rmd"

# Environment variable holding lab pipeline code/resources
SEURAT_V3_ASSETS="/lustre03/project/6004736/singlecell_pipeline/code/seurat_v3/main"

# Path to MACS2 for peak calling (corresponds to running `$ which macs2` after loading 
# the above modules)
MACS2="/cvmfs/soft.mugqic/CentOS6/software/MACS2/MACS2-2.2.7.1/bin/macs2"

# Path where $RMD should be evaluated
WORKDIR=`pwd`

# Run R Markdown file. This command:
# - constructs the render command in R
# - passes in the run, assets path, adn config file as parameters
# - runs & renders the  R Markdown
# - generates the output HTML prefixed by the run
R --no-save -e \
  "rmarkdown::render(
    input = '$RMD',
    output_format = 'html_document',
    output_file = 'preprocessing.html',
    params = list(config = '$CONFIG', macs2_path = '$MACS2', assets = '$SEURAT_V3_ASSETS'),
    knit_root_dir = '$WORKDIR',
    output_dir = '$WORKDIR')"
