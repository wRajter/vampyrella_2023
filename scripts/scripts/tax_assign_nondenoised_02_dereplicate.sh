#!/bin/bash

# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
PROJECT="Suthaus_2022"
MARKER="rDNA"
CELL="cell"
RAW_DATA="../../raw_data"
READS_GZA_DIR="${RAW_DATA}/reads_qza/${PROJECT}/${MARKER}/${CELL}"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"
MANIFEST_FILE="${RAW_DATA}/manifest_files/${PROJECT}/PacBioCCSmanifest_${MARKER}_${CELL}.tsv"
MANIFEST_DIR="${RAW_DATA}/manifest_files/${PROJECT}"
DENOISE_DIR="${RAW_DATA}/denoise/${PROJECT}/${MARKER}/${CELL}"

# ##############################
# ## USING INDIVIDUAL SAMPLES ##
# ##############################

# SAMPLES=$(ls ${RAW_READS_DIR}/*.fastq.gz | \
#           awk -F '/' '{ print $NF }' | \
#           awk -F '_' '{ print $1 }' |
#           awk -F '.' '{ print $1 }')

# mkdir -p ${READS_GZA_DIR}/
# mkdir -p ${DENOISE_DIR}/

# for SAMPLE in ${SAMPLES}
# do
  # Cleaning
  # rm -f ${READS_GZA_DIR}/raw_reads_${SAMPLE}.qza \
  #       ${READS_GZA_DIR}/raw_reads_summary_${SAMPLE}.qzv \
  #       ${DENOISE_DIR}/asv_seqs_${SAMPLE}.fasta

  # # Converting FASTQs into Qiime2 (qza) artifact format
  # qiime tools import \
  #     --type SampleData[SequencesWithQuality] \
  #     --input-path ${MANIFEST_DIR}/PacBioCCSmanifest_${MARKER}_${CELL}_${SAMPLE}.tsv \
  #     --output-path ${READS_GZA_DIR}/raw_reads_${SAMPLE}.qza \
  #     --input-format SingleEndFastqManifestPhred33V2

  # # Creating reads summary file
  # qiime demux summarize \
  #   --i-data ${READS_GZA_DIR}/raw_reads_${SAMPLE}.qza \
  #   --o-visualization ${READS_GZA_DIR}/raw_reads_summary_${SAMPLE}.qzv

  # # Cleaning
  # rm -f ${DENOISE_DIR}/asv_table_${SAMPLE}.qza \
  #       ${DENOISE_DIR}/asv_seqs_${SAMPLE}.qza

  # qiime vsearch dereplicate-sequences \
  #   --i-sequences ${READS_GZA_DIR}/raw_reads_${SAMPLE}.qza \
  #   --o-dereplicated-table ${DENOISE_DIR}/asv_table_${SAMPLE}.qza \
  #   --o-dereplicated-sequences ${DENOISE_DIR}/asv_seqs_${SAMPLE}.qza

#   # Convert ASV sequences from qza to fasta format
#   echo "Converting the ASV sequence qza file to fasta file - sample ${SAMPLE}"

#   qiime tools export \
#     --input-path ${DENOISE_DIR}/asv_seqs_${SAMPLE}.qza \
#     --output-path ${DENOISE_DIR}/

#   mv ${DENOISE_DIR}/dna-sequences.fasta ${DENOISE_DIR}/asv_seqs_${SAMPLE}.fasta

# done


# #######################################
# ## PREPARING RAW READS FOR DENOISING ##
# #######################################

# echo "Working on ${MARKER} marker and ${CELL} cell."

# # Cleaning
# mkdir -p ${READS_GZA_DIR}/
# rm -f ${READS_GZA_DIR}/raw_reads.qza \
#       ${READS_GZA_DIR}/raw_reads_summary.qzv


# # Converting FASTQs into Qiime2 (qza) artifact format
# echo "Converting FASTQs into Qiime2 format and creating reads summary table..."

# qiime tools import \
#     --type SampleData[SequencesWithQuality] \
#     --input-path ${MANIFEST_FILE} \
#     --output-path ${READS_GZA_DIR}/raw_reads.qza \
#     --input-format SingleEndFastqManifestPhred33V2

# # Creating reads summary file
# qiime demux summarize \
#    --i-data ${READS_GZA_DIR}/raw_reads.qza \
#    --o-visualization ${READS_GZA_DIR}/raw_reads_summary.qzv



###################
## DEREPLICATING ##
###################

# if we are not denoising sequences, we can use this pipline:

# Cleaning
mkdir -p ${DENOISE_DIR}/
rm -f ${DENOISE_DIR}/asv_table.qza \
      ${DENOISE_DIR}/asv_seqs.qza \
      ${DENOISE_DIR}/asv_table.qzv \
      ${DENOISE_DIR}/asv_seqs.fasta \
      ${DENOISE_DIR}/asv_seqs.qzv \
      ${DENOISE_DIR}/params.log


qiime vsearch dereplicate-sequences \
  --i-sequences ${READS_GZA_DIR}/raw_reads.qza \
  --o-dereplicated-table ${DENOISE_DIR}/asv_table.qza \
  --o-dereplicated-sequences ${DENOISE_DIR}/asv_seqs.qza

qiime feature-table summarize \
 --i-table ${DENOISE_DIR}/asv_table.qza \
 --o-visualization ${DENOISE_DIR}/asv_table.qzv

# Creating ASV sequence qza file form qza file
echo "Creating ASV sequence qza file..."

qiime feature-table tabulate-seqs \
  --i-data ${DENOISE_DIR}/asv_seqs.qza \
  --o-visualization ${DENOISE_DIR}/asv_seqs.qzv


# Convert ASV sequences from qza to fasta format
echo "Converting the ASV sequence qza file to fasta file..."

qiime tools export \
  --input-path ${DENOISE_DIR}/asv_seqs.qza \
  --output-path ${DENOISE_DIR}/

mv ${DENOISE_DIR}/dna-sequences.fasta ${DENOISE_DIR}/asv_seqs.fasta
