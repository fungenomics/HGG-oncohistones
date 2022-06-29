#!/usr/bin/bash
#SBATCH --job-name=07D-deeptools
#SBATCH --account=rrg-kleinman
#SBATCH --cpus-per-task=20
#SBATCH --time=03:00:00
#SBATCH --mem=20G
#SBATCH --output=../logs/07D-deeptools_counts_CGIs.%j.out
#SBATCH --error=../logs/07D-deeptools_counts_CGIs.%j.out

# ---------------------------------
# SETUP
# ---------------------------------

# The deeptools scripts are based on code from Nicolas De Jay.
# This script sequentially runs computeMatrix for several samples/normalizations
# Each job is parallelized internally with deepTools

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use $MUGQIC_INSTALL_HOME/modulefiles
module load mugqic/deepTools/3.5.0

PROJ_DIR=".."
REFS_DIR="${PROJ_DIR}/data/ChIPseq/references"
DATA_DIR="${PROJ_DIR}/data/ChIPseq/"
OUT_DIR="../output/07D/"
METADATA="${PROJ_DIR}/output/07A/ChIP_metadata_isogenic.deepTools_input.tsv"
BLACKLIST="${REFS_DIR}/ENCODE.blacklist.wgEncodeHg19ConsensusSignalArtifactRegions.hg19.noOverlaps.bed"
ANCHOR="${REFS_DIR}/IGV.annotations.cpgIslands.hg19.bed"
readarray CHRS_EXCLUDE < "${REFS_DIR}/hg19_chrs_to_exclude.txt"

mkdir -p ${OUT_DIR}

# ---------------------------------
# RUN DEEPTOOLS
# ---------------------------------

# List all cell lines & window-sizes to loop through
CL_MASTER=("BT245" "DIPGXIII" "DIPGIV")
NORM_MASTER=("RxAdj")
EXT=20000

for CL in ${CL_MASTER[@]}; do

  OUT_DIR_CL="${OUT_DIR}/${CL}/"
  
  mkdir -p $OUT_DIR_CL

  # Get samples & paths to BW files
  # Log the files / samples used, and then read into an array
  cat ${METADATA} | grep "$CL" | grep H3K27me3 | cut -f 1,3 | awk -F '\t' 'NR>=1 {ID=$1"__"$2; print ID}' | sed 's/\s//g' > "${OUT_DIR_CL}/samples.txt"
  readarray SAMPLE_IDS < "${OUT_DIR_CL}/samples.txt"

  echo ${SAMPLE_IDS[@]}
  
  for NORM in ${NORM_MASTER[@]}; do
  
      params="${NORM}_ext${EXT}_me3"
    
      output_prefix="${OUT_DIR_CL}/${params}"
      output_npz="${output_prefix}/matrix.npz"
    
      mkdir -p $output_prefix
      
      echo $output_prefix
      
      if [ $NORM == Rx ] ; then
        BW_COL=5
      elif [ $NORM == RxAdj ]; then
        BW_COL=6
      elif [ $NORM == CPM_Rx ]; then
        BW_COL=7
      elif [ $NORM == CPM_RxAdj ]; then
        BW_COL=8
      fi

      cat ${METADATA} | grep "$CL" | grep H3K27me3 | cut -f $BW_COL  > "${output_prefix}/bw_files.txt"
      readarray BW_FILES < "${output_prefix}/bw_files.txt"
      
      echo ${BW_FILES[@]}
    
      echo "@ deeptools computeMatrix: ${output_npz}"
      start_time=$(date)
      
      if [ ! -f "${output_npz}.done" ] ; then
          echo "@ deeptools computeMatrix: generating ${output_npz}"
          echo "........................................................."
    
          echo "computeMatrix reference-point \
            -S ${BW_FILES[@]} \
            -R ${ANCHOR} \
            --referencePoint center \
            --beforeRegionStartLength ${EXT} \
            --afterRegionStartLength ${EXT} \
            --outFileName ${output_npz} \
            --blackListFileName ${BLACKLIST} \
            --skipZeros \
            -p 20" > ${output_prefix}/command.txt
        
          computeMatrix reference-point \
            -S ${BW_FILES[@]} \
            -R ${ANCHOR} \
            --referencePoint center \
            --beforeRegionStartLength ${EXT} \
            --afterRegionStartLength ${EXT} \
            --outFileName ${output_npz} \
            --blackListFileName ${BLACKLIST} \
            --skipZeros \
            -p 20
        
          echo "@ deeptools computeMatrix: completed ${output_npz}"
          touch "${output_npz}.done"
    else
        echo "@ deeptools computeMatrix: already generated ${output_npz}"
    fi
      
  done

done


end_time=$(date)
echo "@ Start time: $start_time"
echo "@ End time: $end_time"
echo "@ done."
