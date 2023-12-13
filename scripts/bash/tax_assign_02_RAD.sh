#!/bin/bash

# Define project information
PROJECT="Jamy_2022"
MARKER="rDNA"
RAW_DATA="../../raw_data"
SUFFIX=".fastq.gz"
JL_SCRIPT="../julia/denoise_RAD.jl"



# Construct the reads directory path
READS_DIR="${RAW_DATA}/extracted_18S/${PROJECT}/${MARKER}"


# Decompress all .gz files in the directory
echo "Decompressing all .gz files in ${READS_DIR}"
for file in ${READS_DIR}/*.gz; do
  gunzip -k ${file}
done

# Call your Julia script here
echo "Running Julia script"
julia ${JL_SCRIPT}


# Remove the decompressed files
echo "Compressing all .fastq files in ${READS_DIR}"
for file in ${READS_DIR}/*.fastq; do
  rm -f $file
done

echo "Done"
