#!/bin/bash


# Clustering of reads into OTUs using vsearch

# Variables
THREADS=12
PROJECT="Suthaus_2022"
MARKER="Full18S"
CELL="cellCombined"
THRESHOLD=0.95
SIM="sim95"
DENOISE_METHOD="RAD"
RAW_DATA="../../raw_data"
# input directory:
DENOISE_DIR="${RAW_DATA}/denoise/${PROJECT}/${MARKER}/${CELL}/${DENOISE_METHOD}"
# output directory:
CLUST_DIR="${RAW_DATA}/OTU_clust/${PROJECT}/${MARKER}/${CELL}/${SIM}/${DENOISE_METHOD}"



###############################
## OTU CLUSTERING PER SAMPLE ##
###############################


# per sample

mkdir -p ${CLUST_DIR}

SAMPLES=$(ls ${DENOISE_DIR}/*.fasta | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }' | \
          awk -F '_' '{ print $1 "_" $2 }')


echo "Samples used:"
echo "$SAMPLES"

for SAMPLE in ${SAMPLES}
do

  # Cleaning
  rm -f ${CLUST_DIR}/${SAMPLE}_otu.fasta

  # Clustering
  vsearch --cluster_fast ${DENOISE_DIR}/${SAMPLE}_asv.fasta \
          --id ${THRESHOLD} \
          --threads "${THREADS}" \
          --consout ${CLUST_DIR}/${SAMPLE}_otu.fasta

done
