#!/usr/bin/env bash

#SBATCH --job-name="cNMF"
#SBATCH --account=def-kleinman
#SBATCH --output=logs/cNMF-%j.out
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=6

# ------------------------------------------------------------------------------
# Set up
# ------------------------------------------------------------------------------

# Prep R environment
module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.1
export R_LIBS_USER="../../../../../../renv/library/R-3.6/x86_64-pc-linux-gnu"

# Generate raw counts TSV file with R (cell x gene)
R --vanilla -e "
  library(Seurat)
  load('../seurat.Rda')
  seurat <- subset(x = seurat, subset = Malignant_normal_consensus %in% c('Malignant', 'Likely malignant'))
  expression <- t(as.data.frame(GetAssayData(seurat, slot = 'counts')))
  if (!file.exists('expression_matrix.tsv.gz')) write.table(expression,
            file          = gzfile('expression_matrix.tsv.gz'),
            sep           = '\t',
            row.names     = TRUE,
            col.names     = TRUE,
            quote         = FALSE)
  "

# Prep cNMF environment
module load python/3.6.3
source ~/virtualenv/cNMF/bin/activate

# Path to cNMF code
CNMF="/lustre03/project/6004736/sjessa/software/cNMF"


# ------------------------------------------------------------------------------
# cNMF
# ------------------------------------------------------------------------------

# This part follows the overview at https://github.com/dylkot/cNMF

OUT="output_ngenes2000_niter100_malignant"
mkdir -p $OUT

# Step 1 - normalize the input matrix and prepare the run parameters
echo "@ Step 1"
python $CNMF/cnmf.py prepare --output-dir . --name $OUT -c expression_matrix.tsv.gz \
  -k 5 6 7 8 9 \
  --n-iter 100 \
  --total-workers 6 \
  --seed 100 \
  --numgenes 2000 \
  --beta-loss frobenius
 

# Step 2 - factorize the matrix
# Computationally intensive step -- run in parallel with two workers
echo "@ Step 2"
nohup parallel \
  python $CNMF/cnmf.py factorize --output-dir . --name $OUT \
  --worker-index {} ::: 0 1 2 3 4 5


# Step 3 - combine the results
echo "@ Step 3"
python $CNMF/cnmf.py combine --output-dir . --name $OUT

# remove intermediate files
rm $OUT/cnmf_tmp/$OUT.spectra.k_*.iter_*.df.npz


# Step 4 - select k
echo "@ Step 4"
python $CNMF/cnmf.py k_selection_plot --output-dir . --name $OUT


# Step 5 - get consensus estimates at chosen k
# - Manually choose K based on output of previous plot
# - Run once with threshold at 2.00 to see what distribution of 
#      average distances looks like
# - Run second time to filter outliers based on threshold selected from plot
#      output by previous command

echo "

@ done.
@ The final steps must be completed interactively.
@ Step 5: Issue the following commands and choose K:

salloc --time=00:10:00 --ntasks=1 --account=def-kleinman --mem-per-cpu=4G

module load python/3.6.3
source ~/virtualenv/cNMF/bin/activate

K=? # Choose the desired value

# Run first with a more permissive threshold
python ${CNMF}/cnmf.py consensus --output-dir . --name ${OUT} \
  --components \$K --local-density-threshold 2.00 --show-clustering

# Run a second time with a more stringent threshold, I'll uniformly use 0.02
python ${CNMF}/cnmf.py consensus --output-dir . --name ${OUT} \
  --components \$K --local-density-threshold 0.1 --show-clustering

"

