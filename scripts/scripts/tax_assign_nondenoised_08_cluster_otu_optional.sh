#!/bin/bash

##########################################
## CLUSTERING ALL VAMP SEQUENCES (99%)" ##
##########################################

# De novo clustering of vamp-sequences into a 'new' OTUs using VSEARCH

# Variables
THREADS=12
PROJECT="Jamy_2022"
MARKER="rDNA"
CELL="cell"
THRESHOLD=0.99
SIM="sim99"
RAW_DATA="../../raw_data"
FILT_DIR="${RAW_DATA}/OTU_filtered/${PROJECT}/${MARKER}/${CELL}/${SIM}"
CLUST_DIR="${RAW_DATA}/OTU_clust/${PROJECT}/${MARKER}/${CELL}/${SIM}"
PRECLUST_DIR="${RAW_DATA}/OTU_preclust/${PROJECT}/${MARKER}/${CELL}/${SIM}"



###################
## USING VSEARCH ##
###################



mkdir -p ${CLUST_DIR}

# Get sample names for looping
SAMPLES=$(ls ${PRECLUST_DIR}/*.fasta | \
          awk -F '/' '{ print $NF }' | \
          awk -F '_' '{ print $3 }' |
          awk -F '.' '{ print $1 }')


# Loop over the samples
for SAMPLE in ${SAMPLES}
do

  # Cleaning
  rm -f ${CLUST_DIR}/otu_seqs_${SAMPLE}.fasta

  # Clustering
  vsearch --cluster_fast ${FILT_DIR}/nonrare_otu_${SAMPLE}.fasta \
          --id 0.99 \
          --threads "${THREADS}" \
          --consout ${CLUST_DIR}/otu_seqs_${SAMPLE}.fasta

done
