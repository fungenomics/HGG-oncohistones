# Adapted from Steven Hebert
#
# parameters:
# - normalization: one of "Rx" or "CPM_Rx". If "Rx", then all the values in the bw will
#   be multiplied by <rx_factor>. If "CPM_Rx", the total signal across all bases
#   (a proxy for coverage) will be calculated and divided by 1,000,000,
#   and values in the bw will be divided by this scaling factor, and then multiplied
#   by <rx_factor>.
# - bw_file: name of bw file
# - rx_factor: numeric Rx value
#
# usage:
# $ generate_normalize_bw.sh <normalization> <bw_file> <rx_factor>

NORMALIZATION=$1
BW_FILE=$2
RX_FACTOR=$3

echo "@ --------------------"
echo "normalization: " $NORMALIZATION
echo "bw file: " $BW_FILE
echo "rx factor: " $RX_FACTOR

# Create logs directory
if [ ! -d "logs" ] ; then
  mkdir "logs"
fi

echo "#!/usr/bin/env bash

#SBATCH -J normalize_bw
#SBATCH -A rrg-kleinman
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -t 1:00:00
#SBATCH --output=logs/${BW_FILE%.bw}.%j.out
#SBATCH --error=logs/${BW_FILE%.bw}.%j.err
#SBATCH --mem=12G

cd \$SLURM_SUBMIT_DIR

module use /cvmfs/soft.mugqic/CentOS6/modulefiles
module add mugqic/ucsc/v387

CHR_SIZE_FILE=/lustre06/project/6004736/sjessa/from_narval/references/hg19.chrom.sizes

NORM=$NORMALIZATION

BW_NORM_FILE=normalized/${BW_FILE%.bw}.normalized.$NORMALIZATION.bw
BEDGRAPH_FILE=${BW_FILE%.bw}.$NORMALIZATION.bedGraph
BEDGRAPH_NORM_FILE=\${BEDGRAPH_FILE%.bedGraph}.normalized.bedGraph
BEDGRAPH_SORTED_FILE=\${BEDGRAPH_FILE%.bedGraph}.normalized.sorted.bedGraph

# Convert bw to bedGraph
bigWigToBedGraph $BW_FILE \$BEDGRAPH_FILE

# Calculate the sum of the coverage of the BW file
COVERAGE=\`awk -F'\t' '{sum+=\$4;} END{print sum;}' \$BEDGRAPH_FILE\`

# Convert coverage in scientific notation to decimal, then divide by 1 million to get
# per-million scaling factor 
# https://stackoverflow.com/a/56257197
SCALING_FACTOR=\`echo "\$COVERAGE" | awk -F"E" 'BEGIN{OFMT=\"%10.3f\"} {print \$1 * (10 ^ \$2) / 1000000}'\`

# Normalize file
if [ \$NORM == "Rx" ] || [ \$NORM == "RxAdj" ] ; then
  awk -v RX_FACTOR="$RX_FACTOR" -F'\t' 'BEGIN {OFS = FS} {if(\$4 = \$4*RX_FACTOR) print \$0}' \$BEDGRAPH_FILE > \$BEDGRAPH_NORM_FILE
elif [ \$NORM  == "CPM_Rx" ] || [ \$NORM  == "CPM_RxAdj" ] ; then
  awk -v SCALING_FACTOR=\"\$SCALING_FACTOR\" -v RX_FACTOR="$RX_FACTOR" -F'\t' 'BEGIN {OFS = FS} {if(\$4 = (\$4/SCALING_FACTOR)*RX_FACTOR) print \$0}' \$BEDGRAPH_FILE > \$BEDGRAPH_NORM_FILE
fi

# Sort bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n \$BEDGRAPH_NORM_FILE > \$BEDGRAPH_SORTED_FILE

# Convert bedGraph to bw
bedGraphToBigWig \$BEDGRAPH_SORTED_FILE \$CHR_SIZE_FILE \$BW_NORM_FILE

# Cleanup
rm -f \$BEDGRAPH_FILE \$BEDGRAPH_NORM_FILE \$BEDGRAPH_SORTED_FILE
" > ${BW_FILE%.bw}.normalize.$NORMALIZATION.sh
