#!/bin/bash

# Filtering chimeric sequences after OTU clustering using vsearch

# Variables
NCORES=12
PROJECT="Suthaus_2022"
MARKER="Full18S"
CELL="cellCombined"
SIM="sim95"
RAW_DATA="../../raw_data"
DENOISE_METHOD="RAD"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"
# Input directory
OTU_DIR="${RAW_DATA}/OTU_clust/${PROJECT}/${MARKER}/${CELL}/${SIM}/${DENOISE_METHOD}"
# Output directory
OTU_CHIM_FILT="${RAW_DATA}/OTU_nonchimeric/${PROJECT}/${MARKER}/${CELL}/${SIM}/${DENOISE_METHOD}"




################################################
## CHIMERA FILTERING USING INDIVIDUAL SAMPLES ##
################################################


SAMPLES=$(ls ${RAW_READS_DIR}/*.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }')



echo "Samples used:"
echo "$SAMPLES"


mkdir -p ${OTU_CHIM_FILT}/


for SAMPLE in ${SAMPLES}
do
  # Cleaning
  rm -f ${OTU_CHIM_FILT}/otu_${SAMPLE}.fasta \
  # Chimera removal
  echo "De novo chimera filtering: sample ${SAMPLE}"
  vsearch --uchime_denovo ${OTU_DIR}/${SAMPLE}_otu.fasta \
          --threads 0 \
          --nonchimeras ${OTU_CHIM_FILT}/${SAMPLE}_otu.fasta
done
