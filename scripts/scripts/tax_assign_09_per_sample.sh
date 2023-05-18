#!/bin/bash

# Filtering OTUs into files based on sample they belong to.
# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables:
PROJECT="Suthaus_2022"
CELL="cell"
MARKER="rDNA"
SIM="sim99"
RAW_DATA="../../raw_data"
PER_SAMPLE_DIR="${RAW_DATA}/per_sample_results/${PROJECT}/${MARKER}/${CELL}/${SIM}"
FILT_OTU_DIR="${RAW_DATA}/OTU_filtered/${PROJECT}/${MARKER}/${CELL}/${SIM}"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}"
# SAMPLES=$(ls ${RAW_READS_DIR}/*reads.fastq.gz | \
#           awk -F '/' '{ print $NF }' | \
#           awk -F '_' '{ print $1 }')

SAMPLES=$(ls ${RAW_READS_DIR}/*reads.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }')

SAMPLES=$(sed 's/_21R//g' <<<"$SAMPLES")


####################################
## FILTERING OTUs FOR EACH SAMPLE ##
####################################

echo "Working on ${PROJECT} ${MARKER} ${CELL} ${SIM}."

# Getting a fasta file with the OTU sequences for each sample

# Cleaning
mkdir -p ${PER_SAMPLE_DIR}/fasta/
rm -f ${PER_SAMPLE_DIR}/transposed_table.qz*

# Transposing the table (this will flip the sample and feature axes, necessary for the next step).
echo "Transposing OTU table"
qiime feature-table transpose \
  --i-table ${FILT_OTU_DIR}/otu_table_filtered.qza \
  --o-transposed-feature-table ${PER_SAMPLE_DIR}/transposed_table.qza


qiime metadata tabulate \
  --m-input-file ${PER_SAMPLE_DIR}/transposed_table.qza \
  --o-visualization ${PER_SAMPLE_DIR}/transposed_table.qzv

# Using the transposed table as feature metadata, and keep only the OTUs found in individual samples.
echo "Using transposed table for filtering OTUs into individual samples"


for SAMPLE in ${SAMPLES}
do
  # rm -f ${PER_SAMPLE_DIR}/otu_seqs_filtered_${SAMPLE}.qza \
  #       ${PER_SAMPLE_DIR}/otu_seqs_filtered_${SAMPLE}.qzv
  rm -f ${PER_SAMPLE_DIR}/fasta/otu_seqs_filtered_${SAMPLE}.fasta

  # filter asv sequence table
  qiime feature-table filter-seqs \
    --i-data ${FILT_OTU_DIR}/otu_seqs_filtered.qza \
    --m-metadata-file ${PER_SAMPLE_DIR}/transposed_table.qza \
    --p-where ${SAMPLE}_${MARKER}_${CELL} \
    --o-filtered-data ${PER_SAMPLE_DIR}/otu_seqs_filtered_${SAMPLE}.qza


  qiime feature-table tabulate-seqs \
    --i-data ${PER_SAMPLE_DIR}/otu_seqs_filtered_${SAMPLE}.qza \
    --o-visualization ${PER_SAMPLE_DIR}/otu_seqs_filtered_${SAMPLE}.qzv

  # 9 Convert ASV sequences from qza to fasta format
  echo "Converting the ${SAMPLE} qza file to fasta file"

  qiime tools export \
    --input-path ${PER_SAMPLE_DIR}/otu_seqs_filtered_${SAMPLE}.qza \
    --output-path ${PER_SAMPLE_DIR}/fasta/

  mv ${PER_SAMPLE_DIR}/fasta/dna-sequences.fasta \
     ${PER_SAMPLE_DIR}/fasta/otu_seqs_filtered_${SAMPLE}.fasta
done

# Cleaning
rm -f ${PER_SAMPLE_DIR}/transposed_table.qz*
