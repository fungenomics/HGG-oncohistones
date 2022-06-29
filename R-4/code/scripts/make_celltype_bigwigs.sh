#!/usr/bin/bash
#SBATCH --job-name="make_celltype_bw"
#SBATCH --account=rrg-kleinman
#SBATCH --cpus-per-task=6
#SBATCH --time=03:00:00
#SBATCH --mem=20G
#SBATCH --output=logs/make_celltype_bw-%j.out
#SBATCH --error=logs/make_celltype_bw-%j.out

# Purpose: This script extracts bigwig files from scATACseq data, with one
# bigwig per cluster, so they can be browsed in IGV
# Following the recommendation of the Signac developer (https://github.com/timoast/signac/issues/28)
# - Extract bams for each cluster (however, with the 10X tool subset-bam instead of sinto)
# - Create bigwig tracks from bams https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html
#
# Note: the main memory intensive part is the loading of the Seurat object; the rest
# of the script handling BAMs/bw typically requires <5G

# -----------------------------------------------
# 0. SETUP
# -----------------------------------------------

# make output dir
mkdir -p output

# -----------------------------------------------
# 1. Extract barcodes for each cluster from Seurat object
# -----------------------------------------------

echo "@ STEP 1: extract barcodes"

module load nixpkgs/16.09
module load gcc/7.3.0
module load r/4.0.0
export R_LIBS_USER="../../../../../renv/library/R-4.0/x86_64-pc-linux-gnu"

# Now we will load the Seurat object, splitting the 
# barcodes into groups based on their cluster assignment,
# which is stored in the @meta.data column "cluster_predicted.id"

R --vanilla -e "
  library(here)
  library(dplyr)
  library(Seurat)
  library(Signac)
  library(magrittr)
  source(here('code/functions/scRNAseq.R'))
  source(here('R-4.0.0/code/functions/scATACseq.R'))
  load('../seurat.Rda')
  seurat <- aggregate_cell_types_per_cluster(seurat)
  df <- data.frame('barcode' = colnames(seurat), 'cluster' = seurat\$Type_status_cluster, stringsAsFactors = F)
  df_list <- split(df, df\$cluster)
  names(df_list) <- gsub('/', '_', names(df_list)) # No '/' allowed in file names
  purrr::iwalk(df_list, ~ readr::write_tsv(.x %>% dplyr::select(barcode), glue:::glue('output/{.y}.barcodes.tsv'), col_names = FALSE))
  "

echo "@@ # of barcodes in each:"
wc -l output/*.barcodes.tsv

# -----------------------------------------------
# 2. Create bams for each cluster for each data type using subset-bam, and index
# -----------------------------------------------

# seff for this part:
# Job ID: 16360770
# Cluster: beluga
# User/Group: sjessa/rrg-kleinman
# State: FAILED (exit code 2)
# Nodes: 1
# Cores per node: 6
# CPU Utilized: 03:02:03
# CPU Efficiency: 29.34% of 10:20:30 core-walltime
# Job Wall-clock time: 01:43:25
# Memory Utilized: 1.18 GB
# Memory Efficiency: 11.78% of 10.00 GB

echo "@ STEP 2A: subset bams"

# Fill in the path to the BAM file from cellranger-atac_count
SAMPLE=`readlink -f ..`
SAMPLE=${SAMPLE##*/}
CR_DIR=`grep cellranger ../preprocessing/scMultiome_preprocessing.config.tsv | cut -f 2`
BAM_ATAC=$CR_DIR"/atac_possorted_bam.bam"
BAM_RNA=$CR_DIR"/gex_possorted_bam.bam"

if [[ ! -f $BAM_ATAC ]] ; then
    echo "@@ BAM ${BAM_ATAC} does not exist; exiting."
    exit
fi

if [[ ! -f $BAM_RNA ]] ; then
    echo "@@ BAM ${BAM_RNA} does not exist; exiting."
    exit
fi

echo "@@ Found ${BAM_ATAC}"
echo "@@ Found ${BAM_RNA}"

module load subset-bam/1.1
module load samtools

# Run the tools in a loop
for i in output/*.barcodes.tsv
do
  echo "@@ ${i}"
  cluster=${i##*/}
  cluster="${cluster%.barcodes.tsv}"
  subset-bam --bam $BAM_ATAC \
    --cell-barcodes $i \
    --cores 6 \
    --out-bam "output/${cluster}.ATAC.bam"
  subset-bam --bam $BAM_RNA \
    --cell-barcodes $i \
    --cores 6 \
    --out-bam "output/${cluster}.RNA.bam"
done

echo "@ STEP 2B: index bams"

for i in output/*.bam
do
  echo "@@ ${i}"
  samtools index ${i}
done

# -----------------------------------------------
# 3. Create bigwig tracks using deepTools
# -----------------------------------------------

# # seff for this part
# Job ID: 16365727
# Cluster: beluga
# User/Group: sjessa/rrg-kleinman
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 6
# CPU Utilized: 04:27:14
# CPU Efficiency: 86.65% of 05:08:24 core-walltime
# Job Wall-clock time: 00:51:24
# Memory Utilized: 977.20 MB
# Memory Efficiency: 23.86% of 4.00 GB

echo "@ STEP 3: make bw files"

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use $MUGQIC_INSTALL_HOME/modulefiles
module load mugqic/deepTools/3.5.0

for i in output/*.bam
do
  echo "@@ ${i}"
  cluster=${i##*/}
  cluster="${cluster%.bam}"
  bamCoverage -b $i \
    -o "output/${cluster}.bw" \
    -of bigwig \
    -p 6 \
    --binSize 1 \
    --ignoreDuplicates \
    --normalizeUsing RPKM \
    --ignoreForNormalization chrX chrY chrM    
done

echo "@ DONE."
echo "@ last git commit:"
git log --pretty=oneline | head -n 1
