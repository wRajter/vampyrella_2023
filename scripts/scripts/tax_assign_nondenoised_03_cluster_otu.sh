#!/bin/bash


# Clustering of reads into OTUs using vsearch

# Variables
THREADS=8
PROJECT="Suthaus_2022"
MARKER="rDNA"
CELL="cell"
THRESHOLD=0.97
SIM="sim97"
RAW_DATA="../../raw_data"
DENOISE_DIR="${RAW_DATA}/denoise/${PROJECT}/${MARKER}/${CELL}"
DEREP_DIR="${RAW_DATA}/dereplicate/${PROJECT}/${MARKER}/${CELL}"
CLUST_DIR="${RAW_DATA}/OTU_clust/${PROJECT}/${MARKER}/${CELL}/${SIM}"
PRECLUST_DIR="${RAW_DATA}/OTU_preclust/${PROJECT}/${MARKER}/${CELL}/${SIM}"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"




###############################
## OTU CLUSTERING PER SAMPLE ##
###############################


# per sample

mkdir -p ${CLUST_DIR}

SAMPLES=$(ls ${RAW_READS_DIR}/*.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }')


echo "Samples used:"
echo "$SAMPLES"

for SAMPLE in ${SAMPLES}
do

  # Cleaning
  rm -f ${CLUST_DIR}/otu_${SAMPLE}.fasta

  # Clustering
  vsearch --cluster_fast ${RAW_READS_DIR}/${SAMPLE}.hifi_reads.fastq.gz \
          --id ${THRESHOLD} \
          --threads "${THREADS}" \
          --consout ${CLUST_DIR}/otu_${SAMPLE}.fasta

done
