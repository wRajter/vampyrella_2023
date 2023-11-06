#!/bin/bash

##########################################
## CLUSTERING ALL VAMP SEQUENCES (99%)" ##
##########################################

# De novo clustering of vamp-sequences into a 'new' OTUs using VSEARCH

# Variables
THREADS=12
PROJECT="Suthaus_2022"
MARKER="rDNA"
CELL="cell"
SIM="sim90"
THRESHOLD=0.90
PRECLUST_SIM="sim99"
RAW_DATA="../../raw_data"
PRECLUST_DIR="${RAW_DATA}/OTU_preclust/${PROJECT}/${MARKER}/${CELL}/${PRECLUST_SIM}"
CLUST_DIR="${RAW_DATA}/OTU_clust/${PROJECT}/${MARKER}/${CELL}/${SIM}"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"



###################
## USING VSEARCH ##
###################

mkdir -p ${CLUST_DIR}

SAMPLES=$(ls ${RAW_READS_DIR}/*.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }')


echo "Samples used:"
echo "$SAMPLES"


# Loop over the samples
for SAMPLE in ${SAMPLES}
do

  # Cleaning
  rm -f ${CLUST_DIR}/otu_${SAMPLE}.fasta

  # Clustering
  vsearch --cluster_fast ${PRECLUST_DIR}/otu_${SAMPLE}.fasta \
          --id ${THRESHOLD} \
          --threads "${THREADS}" \
          --consout ${CLUST_DIR}/otu_${SAMPLE}.fasta

done
