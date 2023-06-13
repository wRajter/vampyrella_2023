#!/bin/bash

# Second filtering chimeric sequences after clustering sequences into OTUs using VSEARCH

# Variables
NCORES=12
PROJECT="Jamy_2022"
MARKER="rDNA"
CELL="cell"
SIM="sim99"
RAW_DATA="../../raw_data"
CLUST_DIR="${RAW_DATA}/OTU_clust/${PROJECT}/${MARKER}/${CELL}/${SIM}"
OTU_CHIM_FILT2="${RAW_DATA}/OTU_nonchimeric_second/${PROJECT}/${MARKER}/${CELL}/${SIM}"




#########################################
## USING INDIVIDUAL SAMPLES  (VSEARCH) ##
#########################################

# Get sample names for looping
SAMPLES=$(ls ${CLUST_DIR}/*.fasta | \
          awk -F '/' '{ print $NF }' | \
          awk -F '_' '{ print $3 }' |
          awk -F '.' '{ print $1 }')



mkdir -p ${OTU_CHIM_FILT2}/


for SAMPLE in ${SAMPLES}
do
  # Cleaning
  rm -f ${OTU_CHIM_FILT2}/otu_seqs_${SAMPLE}.fasta \
  # Chimera removal
  echo "De novo chimera filtering: sample ${SAMPLE}"
  vsearch --uchime_denovo ${CLUST_DIR}/otu_seqs_${SAMPLE}.fasta \
          --threads 0 \
          --nonchimeras ${OTU_CHIM_FILT2}/otu_seqs_${SAMPLE}.fasta
done
