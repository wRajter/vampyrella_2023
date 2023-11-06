#!/bin/bash

# Filtering chimeric sequences after OTU clustering using vsearch

# Variables
NCORES=6
PROJECT="Suthaus_2022"
MARKER="Full18S"
SIM="sim_90"
RAW_DATA="../../raw_data"
DENOISE_METHOD="RAD"
# Input directory
CLUSTERED_DIR="${RAW_DATA}/clustered/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}"
# Output directory
CHIM_FILT_DIR="${RAW_DATA}/chimera_filtered/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}"




################################################
## CHIMERA FILTERING USING INDIVIDUAL SAMPLES ##
################################################


SAMPLES=$(ls ${CLUSTERED_DIR}/*.fasta | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }' | \
          awk -F '_' '{ print $1 "_" $2 }')


echo "Samples used:"
echo "$SAMPLES"


mkdir -p ${CHIM_FILT_DIR}/


for SAMPLE in ${SAMPLES}
do
  # Cleaning
  rm -f ${CHIM_FILT_DIR}/${SAMPLE}_otu.fasta
  # Chimera removal
  echo "De novo chimera filtering: sample ${SAMPLE}"
  vsearch --uchime_denovo ${CLUSTERED_DIR}/${SAMPLE}_otu.fasta \
          --threads ${NCORES} \
          --nonchimeras ${CHIM_FILT_DIR}/${SAMPLE}_otu.fasta
done
