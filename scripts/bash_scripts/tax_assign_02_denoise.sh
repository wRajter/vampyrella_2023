#!/bin/bash

# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
NCORES=8
F_PRIMER="CTGGTTGATYCTGCCAGT"
R_PRIMER="TGATCCTTCTGCAGGTTCACCTAC"
LEARN=1000000 # The number of reads to use when training the error model - recommended: 1000000
CELL="cell2"
MARKER="Full18S"
RAW_DATA="../../raw_data"
PACBIO_READS="${RAW_DATA}/PacBio/Suthaus${MARKER}/${CELL}"
READS_GZA_DIR="${RAW_DATA}/reads_qza/${MARKER}/${CELL}"
MANIFEST_FILE="${RAW_DATA}/manifest_files/PacBioCCSmanifest_${MARKER}_${CELL}.tsv"
DENOISE_DIR="${RAW_DATA}/denoise/${MARKER}/${CELL}"


#######################################
## PREPARING RAW READS FOR DENOISING ##
#######################################

echo "Working on ${MARKER} marker and ${CELL} cell."

# Cleaning
mkdir -p ${READS_GZA_DIR}/
rm -f ${READS_GZA_DIR}/raw_reads.qza \
      ${READS_GZA_DIR}/raw_reads_summary.qzv


# Converting FASTQs into Qiime2 (qza) artifact format
echo "Converting FASTQs into Qiime2 format and creating reads summary table..."

qiime tools import \
    --type SampleData[SequencesWithQuality] \
    --input-path ${MANIFEST_FILE} \
    --output-path ${READS_GZA_DIR}/raw_reads.qza \
    --input-format SingleEndFastqManifestPhred33V2

# Creating reads summary file
qiime demux summarize \
   --i-data ${READS_GZA_DIR}/raw_reads.qza \
   --o-visualization ${READS_GZA_DIR}/raw_reads_summary.qzv



########################################################
## DENOISING THE READS INTO AMPLICON SEQUENCE VARIANT ##
########################################################

# Cleaning
mkdir -p ${DENOISE_DIR}/
rm -f ${DENOISE_DIR}/asv_table.qza \
      ${DENOISE_DIR}/asv_seqs.qza \
      ${DENOISE_DIR}/asv_stats.qza

# Running DADA2
echo "Denoising raw reads and creating ASVs using DADA2..."
qiime dada2 denoise-ccs --i-demultiplexed-seqs ${READS_GZA_DIR}/raw_reads.qza \
 --p-n-reads-learn ${LEARN} \
 --p-min-len 800 --p-max-len 3000 \
 --p-n-threads ${NCORES} \
 --p-front ${F_PRIMER} --p-adapter ${R_PRIMER} \
 --o-table ${DENOISE_DIR}/asv_table.qza \
 --o-representative-sequences ${DENOISE_DIR}/asv_seqs.qza \
 --o-denoising-stats ${DENOISE_DIR}/asv_stats.qza \
 --verbose



###################
## POST DENOISE  ##
###################

# Cleaning
rm -f ${DENOISE_DIR}/asv_stats.qzv \
      ${DENOISE_DIR}/asv_table.qzv \
      ${DENOISE_DIR}/asv_seqs.qzv \
      ${DENOISE_DIR}/asv_seqs.fasta


# Visualize the stats file from denoising
echo "Creating stat qzv file..."

qiime metadata tabulate \
  --m-input-file ${DENOISE_DIR}/asv_stats.qza \
  --o-visualization ${DENOISE_DIR}/asv_stats.qzv


# Summarizing DADA2 output
echo "Creating ASVs summary..."

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
