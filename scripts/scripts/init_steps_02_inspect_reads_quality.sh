#!/bin/bash

# Inspecting reads quality and creating quality reports using multiqc.

# Variables
NCORES=6
PROJECT="Suthaus_2022"
MARKER="Full18S"
RAW_DATA="../../raw_data"
RESULTS_DIR="../../results"
RAW_READS="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/cellCombined"
OUTPUT_DIR="${RESULTS_DIR}/multiqc/${PROJECT}/${MARKER}"



echo "Creating reads quality report using multiqc"
# Cleanup of the output
mkdir -p ${OUTPUT_DIR}/
rm -f ${OUTPUT_DIR}/*

# Creating report for all fastq files separatelly:
fastqc -t $NCORES ${RAW_READS}/*.fastq.gz -o ${OUTPUT_DIR}

# Aggregating the summary files into a single report
multiqc ${OUTPUT_DIR}/ -o ${OUTPUT_DIR}/

# Cleanup: Delete individual FastQC reports, keeping only the MultiQC report:
rm -f ${OUTPUT_DIR}/*fastqc.html
rm -f ${OUTPUT_DIR}/*fastqc.zip
rm -rf ${OUTPUT_DIR}/multiqc_data/
