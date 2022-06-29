#!/usr/bin/bash
#SBATCH --job-name=07B-deeptools
#SBATCH --account=rrg-kleinman
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=10G
#SBATCH --output=../logs/07B-deeptools_counts.%j.out
#SBATCH --error=../logs/07B-deeptools_counts.%j.out

# ---------------------------------
# SETUP
# ---------------------------------

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use $MUGQIC_INSTALL_HOME/modulefiles
module load mugqic/deepTools/3.5.0

PROJ_DIR="/lustre06/project/6004736/sjessa/from_narval/HGG-oncohistones/stable/"
REFS_DIR="${PROJ_DIR}/data/ChIPseq/references"
DATA_DIR="${PROJ_DIR}/data/ChIPseq/"
OUT_DIR="../output/07B/"
METADATA="${PROJ_DIR}/output/07A/ChIP_metadata_parental.deeptools_input.tsv"
BLACKLIST="${REFS_DIR}/ENCODE.blacklist.wgEncodeHg19ConsensusSignalArtifactRegions.hg19.noOverlaps.bed"
readarray CHRS_EXCLUDE < "${REFS_DIR}/hg19_chrs_to_exclude.txt"

mkdir -p ${OUT_DIR}

# ---------------------------------
# PREP SAMPLES
# ---------------------------------

# Get samples & paths to BAM files
# Log the files / samples used, and then read into an array
cat ${METADATA} | cut -f 2,4,5,7,8,9,19 | awk -F '\t' 'NR>=2 {ID=$1"__"$6"__"$3; print ID}' > "${OUT_DIR}/samples.txt"
readarray SAMPLE_IDS < "${OUT_DIR}/samples.txt"

cat ${METADATA} | cut -f 2,4,5,7,8,9,19 | awk -F '\t' 'NR>=2 {BAM=$7; print BAM}' > "${OUT_DIR}/bam_files.txt"
readarray BAM_FILES < "${OUT_DIR}/bam_files.txt"


# ---------------------------------
# CONFIGURE SETTINGS
# ---------------------------------

# Based off code from Nicolas De Jay
mapq=0
bin_size=10000
normalization="blacklist_bs${bin_size}_raw_MAPQ${mapq}_noDup"

output_prefix="${OUT_DIR}/${normalization}"
output_npz="${output_prefix}/counts.npz"
output_tab="${output_prefix}/counts.tab"

mkdir -p $output_prefix

# ---------------------------------
# RUN DEEPTOOLS
# ---------------------------------

echo "@ deeptools multiBamSummary"
start_time=$(date)

echo "multiBamSummary bins \
        --binSize ${bin_size} \
        --outFileName ${output_npz} \
        --outRawCounts ${output_tab} \
        --blackListFileName $BLACKLIST \
        -b ${BAM_FILES[@]} \
        --labels ${SAMPLE_IDS[@]} \
        --ignoreDuplicates \
        --minMappingQuality $mapq \
        -p 4 -v"
        
multiBamSummary bins \
        --binSize ${bin_size} \
        --outFileName ${output_npz} \
        --outRawCounts ${output_tab} \
        --blackListFileName $BLACKLIST \
        -b ${BAM_FILES[@]} \
        --labels ${SAMPLE_IDS[@]} \
        --ignoreDuplicates \
        --minMappingQuality $mapq \
        -p 4 -v

end_time=$(date)
echo "@ Start time: $start_time"
echo "@ End time: $end_time"
echo "@ done."
