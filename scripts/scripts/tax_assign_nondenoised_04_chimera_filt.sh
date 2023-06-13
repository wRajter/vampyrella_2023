#!/bin/bash

# Filtering chimeric sequences after OTU clustering using vsearch

# Variables
NCORES=12
PROJECT="Jamy_2022"
MARKER="rDNA"
CELL="cell"
SIM="sim97"
RAW_DATA="../../raw_data"
OTU_DIR="${RAW_DATA}/OTU_clust/${PROJECT}/${MARKER}/${CELL}/${SIM}"
PRECLUST_DIR="${RAW_DATA}/OTU_preclust/${PROJECT}/${MARKER}/${CELL}/${SIM}"
OTU_CHIM_FILT="${RAW_DATA}/OTU_nonchimeric/${PROJECT}/${MARKER}/${CELL}/${SIM}"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"
DENOISE_DIR="${RAW_DATA}/denoise/${PROJECT}/${MARKER}/${CELL}"



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
  vsearch --uchime_denovo ${OTU_DIR}/otu_${SAMPLE}.fasta \
          --threads 0 \
          --nonchimeras ${OTU_CHIM_FILT}/otu_${SAMPLE}.fasta
done
