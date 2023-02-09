#!/bin/bash

# De novo clustering of ASVs into OTUs using vsearch through Qiime2
# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
MARKER="Full18S"
CELL="cellCombined"
THRESHOLD=0.99
SIM="sim99"
RAW_DATA="../../raw_data"
DENOISE_DIR="${RAW_DATA}/denoise/${MARKER}/${CELL}"
CLUST_DIR="${RAW_DATA}/OTU_clust/${MARKER}/${CELL}/${SIM}"


################
## CLUSTERING ##
################

echo "Working on ${MARKER} marker and ${CELL} cell."

# Cleaning
mkdir -p ${CLUST_DIR}/
rm -f ${CLUST_DIR}/otu_table.qza \
      ${CLUST_DIR}/otu_seqs.qza


# Clustering
echo "Clustering ASVs into OTUs at ${THRESHOLD} similarity threshold..."

qiime vsearch cluster-features-de-novo \
  --i-table ${DENOISE_DIR}/asv_table.qza \
  --i-sequences ${DENOISE_DIR}/asv_seqs.qza \
  --p-perc-identity ${THRESHOLD} \
  --o-clustered-table ${CLUST_DIR}/otu_table.qza \
  --o-clustered-sequences ${CLUST_DIR}/otu_seqs.qza


#####################
## POST CLUSTERING ##
#####################

echo "Summarizing OTU table and sequences..."

# Cleaning
rm -f ${CLUST_DIR}/otu_table.qzv \
      ${CLUST_DIR}/otu_table_summarize.qzv \
      ${CLUST_DIR}/otu_seqs.qzv \
      ${CLUST_DIR}/otu_seqs.fasta


# visualize otu table
qiime metadata tabulate \
  --m-input-file ${CLUST_DIR}/otu_table.qza \
  --o-visualization ${CLUST_DIR}/otu_table.qzv

# summarize otu table
qiime feature-table summarize \
 --i-table ${CLUST_DIR}/otu_table.qza \
 --o-visualization ${CLUST_DIR}/otu_table_summarize.qzv

# visualizing otu sequences - converting qza to qzv
qiime feature-table tabulate-seqs \
  --i-data ${CLUST_DIR}/otu_seqs.qza \
  --o-visualization ${CLUST_DIR}/otu_seqs.qzv

# convert qiime2 otu sequences to fasta
qiime tools export \
  --input-path ${CLUST_DIR}/otu_seqs.qza \
  --output-path ${CLUST_DIR}/

mv ${CLUST_DIR}/dna-sequences.fasta \
   ${CLUST_DIR}/otu_seqs.fasta
