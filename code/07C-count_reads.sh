#!/usr/bin/bash
#SBATCH --job-name=07C-count
#SBATCH --account=rrg-kleinman
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=5G
#SBATCH --output=../logs/07C-count_reads.%j.out
#SBATCH --error=../logs/07C-count_reads.%j.out

# calculate total number of mapped reads from BAMs

# ---------------------------------
# SETUP
# ---------------------------------

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use $MUGQIC_INSTALL_HOME/modulefiles
module load samtools/1.10

PROJ_DIR=".."
OUT_DIR="../output/07C/"
METADATA="${PROJ_DIR}/output/07A/ChIP_metadata_parental.deeptools_input.tsv"

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
# COUNT READS
# ---------------------------------

echo "@ samtools view"
start_time=$(date)

echo -e "mapped_reads" > "${OUT_DIR}/n_mapped_reads.txt"

for i in "${BAM_FILES[@]}"
do
   : 
   echo $i
   samtools view -c -bh -F 260 $i >> "${OUT_DIR}/n_mapped_reads.txt"
done

end_time=$(date)
echo "@ Start time: $start_time"
echo "@ End time: $end_time"
echo "@ done."
