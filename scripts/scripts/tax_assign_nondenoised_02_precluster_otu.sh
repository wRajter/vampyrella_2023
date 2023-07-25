#!/bin/bash


# Clustering of reads into OTUs using vsearch

# Variables
THREADS=12
PROJECT="Jamy_2019"
MARKER="rDNA"
CELL="cell"
THRESHOLD=0.99
SIM="sim99"
RAW_DATA="../../raw_data"
PRECLUST_DIR="${RAW_DATA}/OTU_preclust/${PROJECT}/${MARKER}/${CELL}/${SIM}"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"



############################################
## ORGANIZING FILES AFTER DADA2 FILTERING ##
############################################

mkdir -p ${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/raw/

mv ${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/*.fastq.gz \
   ${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/raw/



###############################
## OTU CLUSTERING PER SAMPLE ##
###############################


# per sample

mkdir -p ${PRECLUST_DIR}

SAMPLES=$(ls ${RAW_READS_DIR}/*.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }')


echo "Samples used:"
echo "$SAMPLES"

for SAMPLE in ${SAMPLES}
do

  # Cleaning
  rm -f ${PRECLUST_DIR}/otu_${SAMPLE}.fasta

  # Clustering
  vsearch --cluster_fast ${RAW_READS_DIR}/${SAMPLE}.fastq.gz \
          --id ${THRESHOLD} \
          --threads "${THREADS}" \
          --consout ${PRECLUST_DIR}/otu_${SAMPLE}.fasta

done
