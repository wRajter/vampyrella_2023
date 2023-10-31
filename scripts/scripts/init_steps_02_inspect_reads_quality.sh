#!/bin/bash

# Inspecting reads quality and creating quality reports using fastqc and multiqc.
# Note: activate conda fastqc environment before running the script (conda activate fastqc).

# Variables
NCORES=12
PROJECT="Suthaus_2022"
MARKER="Full18S"
RAW_DATA="../../raw_data"
RESULTS="../../results"
RAW_READS="${RAW_DATA}/PacBio/cellCombined/${PROJECT}_${MARKER}/${CELL}/raw"
OUTPUT_DIR="${RESULTS}/multiqc/${PROJECT}/${MARKER}/"



echo "Creating reads quality report using multiqc"
mkdir -p ${OUTPUT_DIR}/
rm -f ${OUTPUT_DIR}/*

# creating report for all fastq files separatelly (optional):
# fastqc -t $NCORES ${RAW_READS}/*.fastq.gz -o ${OUTPUT_DIR}

#  aggregating the summary files into a single report
multiqc ${OUTPUT_DIR}/ -o ${OUTPUT_DIR}/
