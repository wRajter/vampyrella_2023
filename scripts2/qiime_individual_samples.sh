#!/bin/bash

# Importing FASTQs as QIIME2 artifact

RAW_DATA="/home/lubo/code/wRajter/vampyrella_2023/raw_data/"
MANIFEST_DIR="/home/lubo/code/wRajter/vampyrella_2023/raw_data/qiime_input/individual_samples/"
MANIFEST=$(cd ${MANIFEST_DIR} && ls *.tsv)
OUT_DIR="/home/lubo/code/wRajter/vampyrella_2023/raw_data/qiime_output/individual_samples/reads_qza/"

# mkdir -p ${RAW_DATA}qiime_output/individual_samples/reads_qza/

# for FASTQ in ${MANIFEST}
# do
#   echo "working on ${FASTQ}"
#   qiime tools import \
#     --type SampleData[SequencesWithQuality] \
#     --input-path ${MANIFEST_DIR}${FASTQ} \
#     --output-path ${OUT_DIR}raw_reads_${FASTQ}.qza \
#     --input-format SingleEndFastqManifestPhred33V2
# done

# DENOISING THE READS INTO AMPLICON SEQUENCE VARIANT

# Variables

INPUT_DIR="/home/lubo/code/wRajter/vampyrella_2023/raw_data/qiime_output/individual_samples/reads_qza/"
GZA_READS=$(cd ${INPUT_DIR} && ls)
OUTPUT_DIR="${RAW_DATA}qiime_output/individual_samples/dada2_output/"

mkdir -p ${RAW_DATA}qiime_output/individual_samples/dada2_output


for GZA_READ in ${GZA_READS}
do
  echo "working on ${GZA_READ}"
  qiime dada2 denoise-ccs --i-demultiplexed-seqs ${INPUT_DIR}${GZA_READ} \
    --p-min-len 800 --p-max-len 3000 \
    --p-n-threads 12 \
    --p-front CTGGTTGATYCTGCCAGT --p-adapter TGATCCTTCTGCAGGTTCACCTAC \
    --o-table ${OUTPUT_DIR}table_${GZA_READ}.qza \
    --o-representative-sequences ${OUTPUT_DIR}representative_sequences_${GZA_READ}.qza \
    --o-denoising-stats ${OUTPUT_DIR}stats_${GZA_READ}.qza \
    --verbose
done



# mkdir reads_qza
