#!/bin/bash

# Inspecting reads quality and creating quality reports using fastqc and multiqc.
# Note: activate conda fastqc environment before running the script (conda activate fastqc).

# Variables
NCORES=12
PROJECT="Suthaus_2022"
CELL="cell"
MARKER="rDNA"
RAW_DATA="../../raw_data"
RAW_READS="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}"
OUTPUT_DIR="${RAW_DATA}/fastqc_out/${PROJECT}/${MARKER}/${CELL}"



echo "Creating reads quality report using multiqc"
mkdir -p ${OUTPUT_DIR}/
rm -f ${OUTPUT_DIR}/*


fastqc -t $NCORES ${RAW_READS}/*.fastq.gz -o ${OUTPUT_DIR} # creating report for all fastq files separatelly
multiqc ${OUTPUT_DIR}/ -o ${OUTPUT_DIR}/ #  aggregating the summary files into a single report
