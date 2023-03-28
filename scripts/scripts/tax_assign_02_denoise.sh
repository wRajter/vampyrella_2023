#!/bin/bash

# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
NCORES=12
F_PRIMER="GGCAAGTCTGGTGCCAG"
R_PRIMER="GACGAGGCATTTGGCTACCTT"
LEARN=100000 # The number of reads to use when training the error model - recommended: 1000000
PROJECT="Jamy_2019"
MARKER="rDNA"
CELL="cell"
MIN_LENGTH=2500 # minimal length of the sequence, shorter sequences will be filtered out
MAX_LENGTH=6000 # maximal length of the sequence, longer sequences will be filtered out
RAW_DATA="../../raw_data"
PACBIO_READS="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}"
READS_GZA_DIR="${RAW_DATA}/reads_qza/${PROJECT}/${MARKER}/${CELL}"
MANIFEST_FILE="${RAW_DATA}/manifest_files/${PROJECT}/PacBioCCSmanifest_${MARKER}_${CELL}.tsv"
DENOISE_DIR="${RAW_DATA}/denoise/${PROJECT}/${MARKER}/${CELL}"



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




# ########################################################
# ## DENOISING THE READS INTO AMPLICON SEQUENCE VARIANT ##
# ########################################################

# # Cleaning
# mkdir -p ${DENOISE_DIR}/
# rm -f ${DENOISE_DIR}/asv_table.qza \
#       ${DENOISE_DIR}/asv_seqs.qza \
#       ${DENOISE_DIR}/asv_stats.qza

# # Running DADA2
# echo "Denoising raw reads and creating ASVs using DADA2..."
# qiime dada2 denoise-ccs \
#  --i-demultiplexed-seqs ${READS_GZA_DIR}/raw_reads.qza \
#  --p-n-reads-learn ${LEARN} \
#  --p-indels True \
#  --p-max-ee 2 \
#  --p-max-mismatch 2 \
#  --p-min-len ${MIN_LENGTH} \
#  --p-max-len ${MAX_LENGTH} \
#  --p-n-threads ${NCORES} \
#  --p-front ${F_PRIMER} \
#  --p-adapter ${R_PRIMER} \
#  --o-table ${DENOISE_DIR}/asv_table.qza \
#  --o-representative-sequences ${DENOISE_DIR}/asv_seqs.qza \
#  --o-denoising-stats ${DENOISE_DIR}/asv_stats.qza \
#  --verbose



############################
## NON-DENOISING APPROACH ##
############################

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

# Write the parameters used into a log file
echo "Forward primer: ${F_PRIMER}
      Reverse primer: ${R_PRIMER}
      The number of reads to use when training the error model: ${LEARN}" > ${DENOISE_DIR}/params.log


# ###################
# ## POST DENOISE  ##
# ###################

# # Cleaning
# rm -f ${DENOISE_DIR}/asv_stats.qzv \
#       ${DENOISE_DIR}/asv_stats.tsv \
#       ${DENOISE_DIR}/asv_table.qzv \
#       ${DENOISE_DIR}/asv_seqs.qzv \
#       ${DENOISE_DIR}/asv_seqs.fasta \
#       ${DENOISE_DIR}/params.log


# # Visualize the stats file from denoising
# echo "Creating stat qzv file..."

# qiime metadata tabulate \
#   --m-input-file ${DENOISE_DIR}/asv_stats.qza \
#   --o-visualization ${DENOISE_DIR}/asv_stats.qzv

# qiime tools export \
#   --input-path ${DENOISE_DIR}/asv_stats.qza \
#   --output-path ${DENOISE_DIR}/
# mv ${DENOISE_DIR}/stats.tsv ${DENOISE_DIR}/asv_stats.tsv

# # Summarizing DADA2 output
# echo "Creating ASVs summary..."

# qiime feature-table summarize \
#  --i-table ${DENOISE_DIR}/asv_table.qza \
#  --o-visualization ${DENOISE_DIR}/asv_table.qzv


# # Creating ASV sequence qza file form qza file
# echo "Creating ASV sequence qza file..."

# qiime feature-table tabulate-seqs \
#   --i-data ${DENOISE_DIR}/asv_seqs.qza \
#   --o-visualization ${DENOISE_DIR}/asv_seqs.qzv


# # Convert ASV sequences from qza to fasta format
# echo "Converting the ASV sequence qza file to fasta file..."

# qiime tools export \
#   --input-path ${DENOISE_DIR}/asv_seqs.qza \
#   --output-path ${DENOISE_DIR}/

# mv ${DENOISE_DIR}/dna-sequences.fasta ${DENOISE_DIR}/asv_seqs.fasta

# # Write the parameters used into a log file
# echo "Forward primer: ${F_PRIMER}
#       Reverse primer: ${R_PRIMER}
#       The number of reads to use when training the error model: ${LEARN}" > ${DENOISE_DIR}/params.log
