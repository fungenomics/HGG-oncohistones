#!/usr/bin/bash
#SBATCH --job-name=07D-heatmap
#SBATCH --account=rrg-kleinman
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem=40G
#SBATCH --output=../logs/07D-deeptools_heatmap_CGIs.%j.out
#SBATCH --error=../logs/07D-deeptools_heatmap_CGIs.%j.out

# ---------------------------------
# SETUP
# ---------------------------------

# Should be run after deeptools computeMatrix
# For each cell line, this script constructs one command per normalization type
# and runs them in parallel with GNU parallel

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use $MUGQIC_INSTALL_HOME/modulefiles
module load mugqic/deepTools/3.5.0

PROJ_DIR=".."
REFS_DIR="${PROJ_DIR}/data/ChIPseq/references"
DATA_DIR="${PROJ_DIR}/data/ChIPseq/"
OUT_DIR="${PROJ_DIR}/output/07D/"
METADATA="${PROJ_DIR}/output/07A/ChIP_metadata_isogenic.deepTools_input.tsv"
BLACKLIST="${REFS_DIR}/ENCODE.blacklist.wgEncodeHg19ConsensusSignalArtifactRegions.hg19.noOverlaps.bed"
ANCHOR="${REFS_DIR}/IGV.annotations.cpgIslands.hg19.bed"
readarray CHRS_EXCLUDE < "${REFS_DIR}/hg19_chrs_to_exclude.txt"

# ---------------------------------
# RUN DEEPTOOLS
# ---------------------------------

# List all cell lines to loop through
CL_MASTER=("BT245" "DIPGIV" "DIPGXIII")
NORM=("RxAdj")
EXT=20000

for CL in ${CL_MASTER[@]}; do

  OUT_DIR_CL="${OUT_DIR}/${CL}/"

  mkdir -p $OUT_DIR_CL
  
  cat ${METADATA} | grep "$CL" |  grep H3K27me3 | cut -f 3 | sed 's/\s//g' > "${OUT_DIR_CL}/labels.txt"
  readarray LABELS < "${OUT_DIR_CL}/labels.txt"

  echo $OUT_DIR_CL
  echo ${LABELS[@]}
  
  rm ${OUT_DIR_CL}/commands_heatmap.txt
  touch ${OUT_DIR_CL}/commands_heatmap.txt

  params="${NORM}_ext${EXT}_me3"
  input_prefix="${OUT_DIR_CL}/${params}"
  output_prefix="${OUT_DIR_CL}/${params}"
  input_npz="${input_prefix}/matrix.npz"
  output_sorted="${output_prefix}/matrix.sorted.bed"
  output_mat="${output_prefix}/matrix_heatmap.tab"
  readarray BW_FILES < "${input_prefix}/bw_files.txt"

  mkdir -p $output_prefix

  start_time=$(date)
  output_png="${output_prefix}/matrix_k1.pdf"
  echo "@ deeptools plotHeatmap: ${output_png}"
    
  # Create command to plot with blue-yellow-red color scheme
  echo "plotHeatmap -m ${input_npz} \
     -out ${output_png} \
     --outFileSortedRegions ${output_sorted} \
     --outFileNameMatrix ${output_mat} \
     --kmeans 1 \
     --colorList \"#4575b4,#91bfdb,#e0f3f8,#ffffbf,#fee090,#fc8d59,#d73027\" \
     --startLabel \"CGI\" \
     --endLabel \"\" \
     --sortRegions descend \
     --sortUsingSamples 1 2 \
     --sortUsing mean \
     --averageTypeSummaryPlot median \
     --samplesLabel "${LABELS[@]}" \
     --verbose" >> ${OUT_DIR_CL}/commands_heatmap.txt
       
  plotHeatmap -m ${input_npz} \
     -out ${output_png} \
     --outFileSortedRegions ${output_sorted} \
     --outFileNameMatrix ${output_mat} \
     --kmeans 1 \
     --colorList "#4575b4,#91bfdb,#e0f3f8,#ffffbf,#fee090,#fc8d59,#d73027" \
     --startLabel "CGI" \
     --endLabel "" \
     --sortRegions descend \
     --sortUsingSamples 1 2 \
     --sortUsing mean \
     --averageTypeSummaryPlot median \
     --samplesLabel "${LABELS[@]}" \
     --verbose

done

end_time=$(date)
echo "@ Start time: $start_time"
echo "@ End time: $end_time"
echo "@ done."



