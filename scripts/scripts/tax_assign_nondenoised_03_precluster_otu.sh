#!/bin/bash

########################################
## PRECLUSTERING ALL SEQUENCES (99%)" ##
########################################

# De novo clustering of ASVs into OTUs using vsearch through Qiime2
# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
THREADS=12
PROJECT="Jamy_2022"
MARKER="rDNA"
CELL="cell"
THRESHOLD=0.99
SIM="sim99"
RAW_DATA="../../raw_data"
DENOISE_DIR="${RAW_DATA}/denoise/${PROJECT}/${MARKER}/${CELL}"
CLUST_DIR="${RAW_DATA}/OTU_clust/${PROJECT}/${MARKER}/${CELL}/${SIM}"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"


#################
## USING QIIME ##
#################

# ######################################
# ## USING INDIVIDUAL SAMPLES (QIIME) ##
# ######################################

# SAMPLES=$(ls ${RAW_READS_DIR}/*.fastq.gz | \
#           awk -F '/' '{ print $NF }' | \
#           awk -F '_' '{ print $1 }' |
#           awk -F '.' '{ print $1 }')

# mkdir -p ${CLUST_DIR}/

# for SAMPLE in ${SAMPLES}
# do
#   echo "Working on ${PROJECT} project ${MARKER} marker, ${CELL}, and ${SAMPLE} sample."

#   # Cleaning
#   rm -f ${CLUST_DIR}/otu_table_${SAMPLE}.qza \
#         ${CLUST_DIR}/otu_seqs_${SAMPLE}.qza

#   # Clustering
#   qiime vsearch cluster-features-de-novo \
#     --i-table ${DENOISE_DIR}/asv_table_${SAMPLE}.qza \
#     --i-sequences ${DENOISE_DIR}/asv_seqs_${SAMPLE}.qza \
#     --p-perc-identity ${THRESHOLD} \
#     --o-clustered-table ${CLUST_DIR}/otu_table_${SAMPLE}.qza \
#     --o-clustered-sequences ${CLUST_DIR}/otu_seqs_${SAMPLE}.qza
# done





# ########################
# ## CLUSTERING (QIIME) ##
# ########################

# echo "Working on ${PROJECT} project ${MARKER} marker and ${CELL} cell."

# # Cleaning
# mkdir -p ${CLUST_DIR}/
# rm -f ${CLUST_DIR}/otu_table.qza \
#       ${CLUST_DIR}/otu_seqs.qza


# # Clustering
# echo "Clustering ASVs into OTUs at ${THRESHOLD} similarity threshold..."

# qiime vsearch cluster-features-de-novo \
#   --i-table ${DENOISE_DIR}/asv_table.qza \
#   --i-sequences ${DENOISE_DIR}/asv_seqs.qza \
#   --p-perc-identity ${THRESHOLD} \
#   --o-clustered-table ${CLUST_DIR}/otu_table.qza \
#   --o-clustered-sequences ${CLUST_DIR}/otu_seqs.qza


# #############################
# ## POST CLUSTERING (QIIME) ##
# #############################

# echo "Summarizing OTU table and sequences..."

# # Cleaning
# rm -f ${CLUST_DIR}/otu_table.qzv \
#       ${CLUST_DIR}/otu_table_summarize.qzv \
#       ${CLUST_DIR}/otu_seqs.qzv \
#       ${CLUST_DIR}/otu_seqs.fasta


# # visualize otu table
# qiime metadata tabulate \
#   --m-input-file ${CLUST_DIR}/otu_table.qza \
#   --o-visualization ${CLUST_DIR}/otu_table.qzv

# # summarize otu table
# qiime feature-table summarize \
#  --i-table ${CLUST_DIR}/otu_table.qza \
#  --o-visualization ${CLUST_DIR}/otu_table_summarize.qzv

# # visualizing otu sequences - converting qza to qzv
# qiime feature-table tabulate-seqs \
#   --i-data ${CLUST_DIR}/otu_seqs.qza \
#   --o-visualization ${CLUST_DIR}/otu_seqs.qzv

# # convert qiime2 otu sequences to fasta
# qiime tools export \
#   --input-path ${CLUST_DIR}/otu_seqs.qza \
#   --output-path ${CLUST_DIR}/

# mv ${CLUST_DIR}/dna-sequences.fasta \
#    ${CLUST_DIR}/otu_seqs.fasta



###################
## USING VSEARCH ##
###################

mkdir -p ${CLUST_DIR}

# SAMPLES=$(ls ${DENOISE_DIR}/*.fasta | \
#           awk -F '/' '{ print $NF }' | \
#           awk -F '_' '{ print $3 }' |
#           awk -F '.' '{ print $1 }')

# ERR6454478

vsearch --cluster_fast ${DENOISE_DIR}/asv_seqs_ERR6454478.fasta \
          --id 0.99 \
          --threads 0 \
          --consout ${CLUST_DIR}/otu_seqs_ERR6454478.fasta

# for SAMPLE in ${SAMPLES}
# do

#   # Cleaning
#   rm -f ${CLUST_DIR}/otu_seqs_${SAMPLE}.fasta

#   # Clustering
#   vsearch --cluster_fast ${DENOISE_DIR}/asv_seqs_${SAMPLE}.fasta \
#           --id 0.99 \
#           --threads "${THREADS}" \
#           --consout ${CLUST_DIR}/otu_seqs_${SAMPLE}.fasta

# done
