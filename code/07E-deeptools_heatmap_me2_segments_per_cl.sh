#!/usr/bin/bash
#SBATCH --job-name=07E-heatmap
#SBATCH --account=rrg-kleinman
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem=40G
#SBATCH --output=../logs/07E-deeptools_heatmap_segments.%j.out
#SBATCH --error=../logs/07E-deeptools_heatmap_segments.%j.out

# ---------------------------------
# SETUP
# ---------------------------------

# Should be run after deeptools computeMatrix
# For each cell line, this script constructs one command per normalization type
# and runs them in parallel with GNU parallel

# Prep R environment
module load nixpkgs/16.09
module load gcc/7.3.0
module load r/3.6.1
export R_LIBS_USER="../renv/library/R-3.6/x86_64-pc-linux-gnu"


export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use $MUGQIC_INSTALL_HOME/modulefiles
module load mugqic/deepTools/3.5.0

# ---------------------------------
# RUN DEEPTOOLS
# ---------------------------------

QUANTILE="0.90"
PROJ_DIR=".."
IN_DIR="../output/07E/"
OUT_DIR="../output/07E/"
METADATA="${PROJ_DIR}/output/07A/ChIP_metadata_isogenic.deepTools_input.tsv"
key="_quantile${QUANTILE}_across_samples"
NORM="RxAdj"

CL_MASTER=("BT245" "DIPGIV" "DIPGXIII")

for CL in ${CL_MASTER[@]}; do

  echo $CL

  IN_DIR_CL="${IN_DIR}/${CL}/"
  OUT_DIR_CL="${OUT_DIR}/${CL}/"

  mkdir -p $OUT_DIR_CL
  
  input_prefix="${IN_DIR_CL}/${NORM}_segments_ext50000_me2"
  output_prefix="${OUT_DIR_CL}/${NORM}_segments_ext50000_me2"
  input_npz="${input_prefix}/matrix.npz"
  output_zmax="${output_prefix}/zMax${key}.txt"
  output_pdf="${output_prefix}/matrix.transformed${key}.pdf"
  output_sorted="${output_prefix}/matrix.transformed${key}.sorted.bed"
  
  mkdir -p $output_prefix

  if [ ! -f "$output_zmax" ] ; then
  # Get the X quantile for each sample
    R --vanilla -e "
      library(purrr)
      library(readr)
      library(data.table)
      source(file.path('$PROJ_DIR', 'code/functions/deepTools_helpers.R'))
      npz <- parse_npz('$input_npz')
      indices <- get_sample_indices(npz)
      npz_subset <- npz\$body[, indices\$start[[1]]:indices\$end[[indices\$n]]]
      zMax <- quantile(npz_subset, $QUANTILE, na.rm = TRUE)
      readr::write_lines(zMax, '$output_zmax')
    "
  else
    echo "@ already generated ${output_zmax}"
  fi
  
  # Read in zMaxes
  ZMAX=`cat $output_zmax`
  
  cat ${METADATA} | grep "$CL" |  grep H3K27me2 | cut -f 3 | sed 's/\s//g' > "${OUT_DIR_CL}/labels.txt"
  readarray LABELS < "${OUT_DIR_CL}/labels.txt"
  
  echo $output_pdf
  echo ${LABELS[@]}
  echo $ZMAX

  # Regenerate the heatmap
  # - We use "nearest" interpolation and a higher DPI to avoid too much
  #   data compression, see: https://www.biostars.org/p/322414/
  # - Set missing data to the lowest colour in the color scale to 
  #   avoid missing data plotted in black
  plotHeatmap -m $input_npz \
      -out $output_pdf \
      --outFileSortedRegions $output_sorted \
      --colorList "#4575b4,#91bfdb,#e0f3f8,#ffffbf,#fee090,#fc8d59,#d73027" \
      --startLabel "domain" \
      --endLabel "" \
      --sortRegions descend \
      --sortUsingSamples 1 2 \
      --sortUsing mean \
      --averageTypeSummaryPlot median \
      --missingDataColor "#4575b4" \
      --samplesLabel ${LABELS[@]} \
      --zMax $ZMAX \
      --verbose \
      --interpolationMethod nearest \
      --dpi 400 \
      --plotTitle $CL
      
done
    
end_time=$(date)
echo "@ Start time: $start_time"
echo "@ End time: $end_time"
echo "@ done."

