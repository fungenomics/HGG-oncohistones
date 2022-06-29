#!/usr/bin/bash
#SBATCH --job-name=07E-deeptools
#SBATCH --account=rrg-kleinman
#SBATCH --cpus-per-task=20
#SBATCH --time=3:00:00
#SBATCH --mem=10G
#SBATCH --output=../logs/07E-deeptools_counts_segments.%j.out
#SBATCH --error=../logs/07E-deeptools_counts_segments.%j.out

# ---------------------------------
# SETUP
# ---------------------------------

# This script sequentially runs computeMatrix for several samples and a given normalization
# Each job is parallelized internally with deepTools

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use $MUGQIC_INSTALL_HOME/modulefiles
module load mugqic/deepTools/3.5.0

PROJ_DIR=".."
REFS_DIR="${PROJ_DIR}/data/ChIPseq/references"
DATA_DIR="${PROJ_DIR}/data/ChIPseq/"
OUT_DIR="../output/07E/"
METADATA="${PROJ_DIR}/output/07A/ChIP_metadata_isogenic.deepTools_input.tsv"
BLACKLIST="${REFS_DIR}/ENCODE.blacklist.wgEncodeHg19ConsensusSignalArtifactRegions.hg19.noOverlaps.bed"
readarray CHRS_EXCLUDE < "${REFS_DIR}/hg19_chrs_to_exclude.txt"

mkdir -p ${OUT_DIR}

# ---------------------------------
# RUN DEEPTOOLS
# ---------------------------------

# List all cell lines to loop through
CL_MASTER=("BT245" "DIPGXIII" "DIPGIV")
NORM=("RxAdj")
EXT_MASTER=(50000)

for CL in ${CL_MASTER[@]}; do

  OUT_DIR_CL="${OUT_DIR}/${CL}/"
  
  mkdir -p $OUT_DIR_CL

  # Get samples & paths to BW files
  # Log the files / samples used, and then read into an array
  cat ${METADATA} | grep "$CL" | cut -f 1,3 | awk -F '\t' 'NR>=1 {ID=$1"__"$2; print ID}' | sed 's/\s//g' > "${OUT_DIR_CL}/samples.txt"
  readarray SAMPLE_IDS < "${OUT_DIR_CL}/samples.txt"

  echo ${SAMPLE_IDS[@]}
  
  for EXT in ${EXT_MASTER[@]}; do
    
    params="${NORM}_segments_ext${EXT}_me2"
    
    output_prefix="${OUT_DIR_CL}/${params}"
    output_npz="${output_prefix}/matrix.npz"
    
    mkdir -p $output_prefix
    echo $output_prefix

    if [ $CL == BT245 ] ; then
        ANCHOR="${DATA_DIR}/domains_@mhulswit/filtered/BT245_NKO_c1p5_H3K27me2_SE.filt.merged.bed"
    elif [ $CL == DIPGXIII ]; then
        ANCHOR="${DATA_DIR}/domains_@mhulswit/filtered/DIPGXIII_parental_C14_H3K27me2.filt.merged.bed"
    elif [ $CL == DIPGIV ]; then
        ANCHOR="${DATA_DIR}/domains_@mhulswit/filtered/DIPGIV-nko-c3p5-Rx_H3K27me2.filt.merged.bed"
    fi
    
    if [ $NORM == Rx ] ; then
        BW_COL=5
    elif [ $NORM == RxAdj ]; then
        BW_COL=6
    elif [ $NORM == CPM_Rx ]; then
        BW_COL=7
    elif [ $NORM == CPM_RxAdj ]; then
        BW_COL=8
    fi

    cat ${METADATA} | grep "$CL" | grep H3K27me2 | cut -f $BW_COL  > "${output_prefix}/bw_files.txt"
    readarray BW_FILES < "${output_prefix}/bw_files.txt"
      
    echo ${BW_FILES[@]}
    echo $ANCHOR
    
    echo "@ deeptools computeMatrix: ${output_npz}"
    start_time=$(date)
      
    if [ ! -f "${output_npz}.done" ] ; then
        echo "@ deeptools computeMatrix: generating ${output_npz}"
        echo "........................................................."
    
        echo "computeMatrix scale-regions \
            -S ${BW_FILES[@]} \
            -R ${ANCHOR} \
            --regionBodyLength 50000 \
            --beforeRegionStartLength ${EXT} \
            --afterRegionStartLength ${EXT} \
            --outFileName ${output_npz} \
            --blackListFileName ${BLACKLIST} \
            --skipZeros \
            --verbose \
            -p 20" > ${output_prefix}/command.txt
        
        computeMatrix scale-regions \
            -S ${BW_FILES[@]} \
            -R ${ANCHOR} \
            --regionBodyLength 50000 \
            --beforeRegionStartLength ${EXT} \
            --afterRegionStartLength ${EXT} \
            --outFileName ${output_npz} \
            --blackListFileName ${BLACKLIST} \
            --verbose \
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
