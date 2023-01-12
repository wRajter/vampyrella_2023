#!/bin/bash

## PREPARING DATA ##

# 1. Inspect read quality

# 2. Activating QIIME 2 conda environment

# conda activate qiime2-2022.2

# 3. Importing FASTQs as QIIME 2 artifact

# gzip -d raw_data/*.fastq.gz

mkdir reads_qza dada2_output

qiime tools import \
    --type SampleData[SequencesWithQuality] \
    --input-path ../raw_data/reads/PacBioCCSmanifest.tsv \
    --output-path reads_qza/raw_reads.qza \
    --input-format SingleEndFastqManifestPhred33V2

# 4. Summarizing trimmed FASTQs

qiime demux summarize \
   --i-data reads_qza/raw_reads.qza \
   --o-visualization reads_qza/raw_reads_summary.qzv




## DENOISING THE READS INTO AMPLICON SEQUENCE VARIANT

# 1. Running DADA2

qiime dada2 denoise-ccs --i-demultiplexed-seqs reads_qza/raw_reads.qza \
 --p-min-len 50 --p-max-len 5000 \
 --p-n-threads 4 \
 --p-front CTGGTTGATYCTGCCAGT --p-adapter TGATCCTTCTGCAGGTTCACCTAC \
 --o-table dada2_output/table.qza \
 --o-representative-sequences dada2_output/representative_sequences.qza \
 --o-denoising-stats dada2_output/stats.qza \
 --verbose



# 2. Summarizing DADA2 output

qiime feature-table summarize \
 --i-table dada2_output/table.qza \
 --o-visualization dada2_output/dada2_table_summary.qzv

# 3 convert ASV sequences from qza to fasta format

# qiime tools export --input-path dada2_output/representative_sequences.qza --output-path ASV/
