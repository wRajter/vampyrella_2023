#!/bin/bash

# de novo clustering of ASVs into OTUs using vsearch through qiime2

# variables
MARKER="Full18S"
CELL="cell2"
THRESHOLD=0.99
ID="sim99"
QIIME_OUT="/home/lubo/code/wRajter/vampyrella_2023/raw_data/qiime_output"
FEATURE_TABLE="${QIIME_OUT}/dada2_denoise/table_${MARKER}_${CELL}_ALLSAMPLES.qza"
REP_SEQS_DIR="${QIIME_OUT}/dada2_denoise/"
CLUST_DIR="${QIIME_OUT}/otu_clust/${MARKER}/${CELL}/${ID}"
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




# ################
# ## CLUSTERING ##
# ################

# cleaning
mkdir -p ${CLUST_DIR}
rm -f ${QIIME_OUT}/otu_clust/table_dn_ALLSAMPLES_${ID}.qza
rm -f ${QIIME_OUT}/otu_clust/otu_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qza

# clustering
qiime vsearch cluster-features-de-novo \
  --i-table ${FEATURE_TABLE} \
  --i-sequences ${REP_SEQS_DIR}/asv_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --p-perc-identity ${THRESHOLD} \
  --o-clustered-table ${CLUST_DIR}/table_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qza \
  --o-clustered-sequences ${CLUST_DIR}/otu_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qza


# #####################
# ## POST CLUSTETING ##
# #####################

# cleaning
mkdir -p ${CLUST_DIR}/fasta
rm -f ${CLUST_DIR}/otu_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qzv
rm -f ${CLUST_DIR}/table_summarize_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qzv
rm -f ${CLUST_DIR}/otu_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qzv
rm -f ${CLUST_DIR}/fasta/otu_${MARKER}_${CELL}_ALLSAMPLES_${ID}.fasta

# visualize otu table
qiime metadata tabulate \
  --m-input-file ${CLUST_DIR}/table_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qza \
  --o-visualization ${CLUST_DIR}/table_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qzv

# summarize otu table
qiime feature-table summarize \
 --i-table ${CLUST_DIR}/table_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qza \
 --o-visualization ${CLUST_DIR}/table_summarize_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qzv

# visualizing otu sequences - converting qza to qzv
qiime feature-table tabulate-seqs \
  --i-data ${CLUST_DIR}/otu_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qza \
  --o-visualization ${CLUST_DIR}/otu_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qzv

# convert qiime2 otu sequences to fasta
qiime tools export \
  --input-path ${CLUST_DIR}/otu_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qza \
  --output-path ${CLUST_DIR}/fasta/

mv ${CLUST_DIR}/fasta/dna-sequences.fasta \
   ${CLUST_DIR}/fasta/otu_${MARKER}_${CELL}_ALLSAMPLES_${ID}.fasta


####################################
## FILTERING OTUs FOR EACH SAMPLE ##
####################################

# Getting a fasta file with the OTU sequences for each sample

# Transposing the table (this will flip the sample and feature axes, necessary for the next step).
rm -f ${CLUST_DIR}/transposed_table_${MARKER}_${CELL}_ALLSAMPLES.qz*

echo "Transposing the ${MARKER} ${CELL} table"
qiime feature-table transpose \
  --i-table ${FEATURE_TABLE} \
  --o-transposed-feature-table ${CLUST_DIR}/transposed_table_${MARKER}_${CELL}_ALLSAMPLES.qza

qiime metadata tabulate \
  --m-input-file ${CLUST_DIR}/transposed_table_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --o-visualization ${CLUST_DIR}/transposed_table_${MARKER}_${CELL}_ALLSAMPLES.qzv

# Using the transposed table as feature metadata, and keep only the OTUs found in individual samples.

for SAMPLE in ${SAMPLES}
do
  rm -f ${CLUST_DIR}/otu_${MARKER}_${CELL}_${ID}_${SAMPLE}.qz*
  rm -f ${CLUST_DIR}/fasta/otu_${MARKER}_${CELL}_${ID}_${SAMPLE}.fasta

  # filter asv sequence table
  qiime feature-table filter-seqs \
    --i-data ${CLUST_DIR}/otu_${MARKER}_${CELL}_ALLSAMPLES_${ID}.qza \
    --m-metadata-file ${CLUST_DIR}/transposed_table_${MARKER}_${CELL}_ALLSAMPLES.qza \
    --p-where ${SAMPLE}_${MARKER}_${CELL} \
    --o-filtered-data ${CLUST_DIR}/otu_${MARKER}_${CELL}_${ID}_${SAMPLE}.qza

  qiime feature-table tabulate-seqs \
    --i-data ${CLUST_DIR}/otu_${MARKER}_${CELL}_${ID}_${SAMPLE}.qza \
    --o-visualization ${CLUST_DIR}/otu_${MARKER}_${ID}_${CELL}_${SAMPLE}.qzv

  # 9 Convert ASV sequences from qza to fasta format
  echo "Converting the ${MARKER} ${CELL} ${ID} ${SAMPLE} OTUs qza file to fasta file"

  qiime tools export \
    --input-path ${CLUST_DIR}/otu_${MARKER}_${CELL}_${ID}_${SAMPLE}.qza \
    --output-path ${CLUST_DIR}/fasta/

  mv ${CLUST_DIR}/fasta/dna-sequences.fasta \
     ${CLUST_DIR}/fasta/otu_${MARKER}_${CELL}_${ID}_${SAMPLE}.fasta
done
