#!/bin/bash

# Filtering OTUs into files based on sample they belong to.
# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables:
CELL="cellCombined"
MARKER="Full18S"
SIM="sim97"
RAW_DATA="../../raw_data"
PER_SAMPLE_DIR="${RAW_DATA}/per_sample_results/${MARKER}/${CELL}/${SIM}"
OTU_CLUST_DIR="${RAW_DATA}/OTU_clust/${MARKER}/${CELL}/${SIM}"
SAMPLES="A3 \
         Mock \
         NH1 \
         NH4 \
         Sim17 \
         Sim22 \
         Th16 \
         Th38 \
         Th40 \
         X17007"


####################################
## FILTERING OTUs FOR EACH SAMPLE ##
####################################

echo "Working with ${SIM} OTUs from ${MARKER} marker and ${CELL} cell."

# Getting a fasta file with the OTU sequences for each sample

# Cleaning
mkdir -p ${PER_SAMPLE_DIR}/fasta/
rm -f ${PER_SAMPLE_DIR}/transposed_table.qz*

# Transposing the table (this will flip the sample and feature axes, necessary for the next step).
echo "Transposing OTU table"
qiime feature-table transpose \
  --i-table ${OTU_CLUST_DIR}/otu_table.qza \
  --o-transposed-feature-table ${PER_SAMPLE_DIR}/transposed_table.qza

qiime metadata tabulate \
  --m-input-file ${PER_SAMPLE_DIR}/transposed_table.qza \
  --o-visualization ${PER_SAMPLE_DIR}/transposed_table.qzv

# Using the transposed table as feature metadata, and keep only the OTUs found in individual samples.
echo "Using transposed table for filtering OTUs into individual samples"

for SAMPLE in ${SAMPLES}
do
  rm -f ${PER_SAMPLE_DIR}/otu_seqs_${SAMPLE}.qz
  rm -f ${PER_SAMPLE_DIR}/fasta/otu_seqs_${SAMPLE}.fasta

  # filter asv sequence table
  qiime feature-table filter-seqs \
    --i-data ${OTU_CLUST_DIR}/otu_seqs.qza \
    --m-metadata-file ${PER_SAMPLE_DIR}/transposed_table.qza \
    --p-where ${SAMPLE}_${MARKER}_${CELL} \
    --o-filtered-data ${PER_SAMPLE_DIR}/otu_seqs_${SAMPLE}.qza

  qiime feature-table tabulate-seqs \
    --i-data ${PER_SAMPLE_DIR}/otu_seqs_${SAMPLE}.qza \
    --o-visualization ${PER_SAMPLE_DIR}/otu_seqs_${SAMPLE}.qzv

  # 9 Convert ASV sequences from qza to fasta format
  echo "Converting the ${SAMPLE} qza file to fasta file"

  qiime tools export \
    --input-path ${PER_SAMPLE_DIR}/otu_seqs_${SAMPLE}.qza \
    --output-path ${PER_SAMPLE_DIR}/fasta/

  mv ${PER_SAMPLE_DIR}/fasta/dna-sequences.fasta \
     ${PER_SAMPLE_DIR}/fasta/otu_seqs_${SAMPLE}.fasta
done

# cleaning
rm -f ${PER_SAMPLE_DIR}/transposed_table.qz*
