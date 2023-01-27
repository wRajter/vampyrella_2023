#!/bin/bash

# Inspecting reads quality and creating quality reports using fastqc and multiqc
# activate conda fastqc environment before running the script

# variables
NCORES=12
CELL="cell1"
MARKER="Full18S"
PROJECT_FILE="/home/lubo/code/wRajter/vampyrella_2023"
RAW_READS="${PROJECT_FILE}/raw_data/PacBio/Suthaus${MARKER}/${CELL}/"
OUTPUT_DIR="/home/lubo/code/wRajter/vampyrella_2023/raw_data/fastqc_out/${MARKER}_${CELL}"



echo "Creating reads quality report using multiqc"
mkdir -p ${OUTPUT_DIR}
rm -f ${OUTPUT_DIR}/*


fastqc -t $NCORES ${RAW_READS}/*.fastq.gz -o ${OUTPUT_DIR} # creating report for all fastq files separatelly
multiqc ${OUTPUT_DIR}/ -o ${OUTPUT_DIR}/ #  aggregating the summary files into a single report
