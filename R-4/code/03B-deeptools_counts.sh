#!/usr/bin/bash
#SBATCH --job-name=03B-deepTools
#SBATCH --account=rrg-kleinman
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --mem=5G
#SBATCH --output=../logs/03B-deeptools-%j.out
#SBATCH --error=../logs/03B-deeptools-%j.out

# ---------------------------------
# SETUP
# ---------------------------------

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
module use $MUGQIC_INSTALL_HOME/modulefiles
module load mugqic/deepTools/3.5.0

PROJ_DIR="../.."
REFS_DIR="${PROJ_DIR}/data/ChIPseq/references"
DATA_DIR="${PROJ_DIR}/data/Paired-Tag/Zhu_2021"
OUT_DIR="../output/05B/"
BW_FILES=(
  "${DATA_DIR}/bigwig/Paired-Tag_H3K27ac_OPC.bw"
  "${DATA_DIR}/bigwig/Paired-Tag_H3K27ac_Ependymal.bw"
  )

mkdir -p ${OUT_DIR}
echo ${BW_FILES[@]}

multiBigwigSummary BED-file \
	-b ${BW_FILES[@]} \
	-l OPC Ependymal \
	-o "${OUT_DIR}/MBS_Ensembl.ensGene.mm10.collapsed.promoter5kb.bed.npz" \
	--outRawCounts "${OUT_DIR}/MBS_Ensembl.ensGene.mm10.collapsed.promoter5kb.tab" \
	--BED "${REFS_DIR}/Ensembl.ensGene.mm10.collapsed.promoter5kb.bounded.bed" \
	-p 4 \
	--verbose

