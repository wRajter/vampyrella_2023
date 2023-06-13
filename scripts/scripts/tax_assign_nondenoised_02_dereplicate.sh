#!/bin/bash

# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
PROJECT="Jamy_2022"
MARKER="rDNA"
CELL="cell"
RAW_DATA="../../raw_data"
READS_GZA_DIR="${RAW_DATA}/reads_qza/${PROJECT}/${MARKER}/${CELL}"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"
MANIFEST_FILE="${RAW_DATA}/manifest_files/${PROJECT}/PacBioCCSmanifest_${MARKER}_${CELL}.tsv"
MANIFEST_DIR="${RAW_DATA}/manifest_files/${PROJECT}"
DEREP_DIR="${RAW_DATA}/dereplicate/${PROJECT}/${MARKER}/${CELL}"


# ##############################
# ## USING INDIVIDUAL SAMPLES ##
# ##############################

# SAMPLES=$(ls ${RAW_READS_DIR}/*.fastq.gz | \
#           awk -F '/' '{ print $NF }' | \
#           awk -F '_' '{ print $1 }' |
#           awk -F '.' '{ print $1 }')

# mkdir -p ${READS_GZA_DIR}/
# mkdir -p ${DEREP_DIR}/

# for SAMPLE in ${SAMPLES}
# do
  # Cleaning
  # rm -f ${READS_GZA_DIR}/raw_reads_${SAMPLE}.qza \
  #       ${READS_GZA_DIR}/raw_reads_summary_${SAMPLE}.qzv \
  #       ${DEREP_DIR}/asv_seqs_${SAMPLE}.fasta

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
  # rm -f ${DEREP_DIR}/asv_table_${SAMPLE}.qza \
  #       ${DEREP_DIR}/asv_seqs_${SAMPLE}.qza

  # qiime vsearch dereplicate-sequences \
  #   --i-sequences ${READS_GZA_DIR}/raw_reads_${SAMPLE}.qza \
  #   --o-dereplicated-table ${DEREP_DIR}/asv_table_${SAMPLE}.qza \
  #   --o-dereplicated-sequences ${DEREP_DIR}/asv_seqs_${SAMPLE}.qza

#   # Convert ASV sequences from qza to fasta format
#   echo "Converting the ASV sequence qza file to fasta file - sample ${SAMPLE}"

#   qiime tools export \
#     --input-path ${DEREP_DIR}/asv_seqs_${SAMPLE}.qza \
#     --output-path ${DEREP_DIR}/

#   mv ${DEREP_DIR}/dna-sequences.fasta ${DEREP_DIR}/asv_seqs_${SAMPLE}.fasta

# done


#######################################
## PREPARING RAW READS FOR DENOISING ##
#######################################

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
mkdir -p ${DEREP_DIR}/
rm -f ${DEREP_DIR}/asv_table.qza \
      ${DEREP_DIR}/asv_seqs.qza \
      ${DEREP_DIR}/asv_table.qzv \
      ${DEREP_DIR}/asv_seqs.fasta \
      ${DEREP_DIR}/asv_seqs.qzv \
      ${DEREP_DIR}/params.log


qiime vsearch dereplicate-sequences \
  --i-sequences ${READS_GZA_DIR}/raw_reads.qza \
  --o-dereplicated-table ${DEREP_DIR}/asv_table.qza \
  --o-dereplicated-sequences ${DEREP_DIR}/asv_seqs.qza

qiime feature-table summarize \
 --i-table ${DEREP_DIR}/asv_table.qza \
 --o-visualization ${DEREP_DIR}/asv_table.qzv

# Creating ASV sequence qza file form qza file
echo "Creating ASV sequence qza file..."

qiime feature-table tabulate-seqs \
  --i-data ${DEREP_DIR}/asv_seqs.qza \
  --o-visualization ${DEREP_DIR}/asv_seqs.qzv


# Convert ASV sequences from qza to fasta format
echo "Converting the ASV sequence qza file to fasta file..."

qiime tools export \
  --input-path ${DEREP_DIR}/asv_seqs.qza \
  --output-path ${DEREP_DIR}/

mv ${DEREP_DIR}/dna-sequences.fasta ${DEREP_DIR}/asv_seqs.fasta




###################
## USING VSEARCH ##
###################
