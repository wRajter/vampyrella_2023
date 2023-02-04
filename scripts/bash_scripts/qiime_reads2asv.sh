#!/bin/bash

# Activate conda qiime2-2022.11 environment before running the script

# Variables
NCORES=6
F_PRIMER="CTGGTTGATYCTGCCAGT"
R_PRIMER="TGATCCTTCTGCAGGTTCACCTAC"
LEARN=1000000 # The number of reads to use when training the error model - recommended: 1000000
CELL="cellCombined"
MARKER="Full18S"
PROJECT_FILE="/home/lubo/code/wRajter/vampyrella_2023"
RAW_READS="${PROJECT_FILE}/raw_data/PacBio/Suthaus${MARKER}/${CELL}/"
MANIFEST="${PROJECT_FILE}/raw_data/qiime_input/PacBioCCSmanifest_${MARKER}_${CELL}.tsv"
OUTPUT="${PROJECT_FILE}/raw_data/qiime_output"
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


##################################
## PREPARING READS AND METADATA ##
##################################



# 1. Importing FASTQs as QIIME2 artifact
mkdir -p ${OUTPUT}/reads_qza
rm -f ${OUTPUT}/reads_qza/raw_reads_${MARKER}_${CELL}.qza
rm -f ${OUTPUT}/reads_qza/raw_reads_summary_${MARKER}_${CELL}.qzv

echo "Importing ${MARKER} ${CELL} reads to QIIME2 artifact"
qiime tools import \
    --type SampleData[SequencesWithQuality] \
    --input-path ${MANIFEST} \
    --output-path ${OUTPUT}/reads_qza/raw_reads_${MARKER}_${CELL}.qza \
    --input-format SingleEndFastqManifestPhred33V2

# 2. Creating reads summary file
echo "Creating ${MARKER} ${CELL} reads summary file"
qiime demux summarize \
   --i-data ${OUTPUT}/reads_qza/raw_reads_${MARKER}_${CELL}.qza \
   --o-visualization ${OUTPUT}/reads_qza/raw_reads_summary_${MARKER}_${CELL}.qzv



# ########################################################
# ## DENOISING THE READS INTO AMPLICON SEQUENCE VARIANT ##
# ########################################################

# 3. Running DADA2
mkdir -p ${OUTPUT}/dada2_denoise
rm -f ${OUTPUT}/dada2_denoise/*_${MARKER}_${CELL}.gz*

echo "denoising ${MARKER} ${CELL} reads and creating ASVs using DADA2"
qiime dada2 denoise-ccs --i-demultiplexed-seqs ${OUTPUT}/reads_qza/raw_reads_${MARKER}_${CELL}.qza \
 --p-n-reads-learn ${LEARN} \
 --p-min-len 800 --p-max-len 3000 \
 --p-n-threads ${NCORES} \
 --p-front ${F_PRIMER} --p-adapter ${R_PRIMER} \
 --o-table ${OUTPUT}/dada2_denoise/table_${MARKER}_${CELL}_ALLSAMPLES.qza \
 --o-representative-sequences ${OUTPUT}/dada2_denoise/asv_${MARKER}_${CELL}_ALLSAMPLES.qza \
 --o-denoising-stats ${OUTPUT}/dada2_denoise/stats_${MARKER}_${CELL}_ALLSAMPLES.qza \
 --verbose


###################
## POST DENOISE  ##
###################

mkdir -p ${OUTPUT}/dada2_post_denoise
rm -f ${OUTPUT}/dada2_post_denoise/*_${MARKER}_${CELL}.gz*

# 4. Summarizing DADA2 output
echo "creating ${MARKER} ${CELL} ASVs summary"

qiime feature-table summarize \
 --i-table ${OUTPUT}/dada2_denoise/table_${MARKER}_${CELL}_ALLSAMPLES.qza \
 --o-visualization ${OUTPUT}/dada2_post_denoise/dada2_table_summary_${MARKER}_${CELL}_ALLSAMPLES.qzv

# 5. Visualization of the representative sequences
echo "creating ${MARKER} ${CELL} ASVs visualization"

qiime feature-table tabulate-seqs \
  --i-data ${OUTPUT}/dada2_denoise/asv_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --o-visualization ${OUTPUT}/dada2_post_denoise/asv_${MARKER}_${CELL}_ALLSAMPLES.qzv

# 6. Convert ASV sequences from qza to fasta format
echo "Converting the ${MARKER} ${CELL} ASVs qza file to fasta file"

mkdir -p ${OUTPUT}/dada2_denoise/fasta
rm -f ${OUTPUT}/dada2_post_denoise/fasta/asv_${MARKER}_${CELL}_ALLSAMPLES.fasta

qiime tools export \
  --input-path ${OUTPUT}/dada2_denoise/asv_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --output-path ${OUTPUT}/dada2_post_denoise/fasta/

mv ${OUTPUT}/dada2_post_denoise/fasta/dna-sequences.fasta ${OUTPUT}/dada2_post_denoise/fasta/asv_${MARKER}_${CELL}_ALLSAMPLES.fasta


########################################################
## FILTERING REPRESENTATIVE SEQUENCES FOR EACH SAMPLE ##
########################################################

# Getting a fasta file with the representative sequences for each sample


# 7. Transposing the table (this will flip the sample and feature axes, necessary for the next step).
rm -f ${OUTPUT}/dada2_post_denoise/transposed_table_${MARKER}_${CELL}_ALLSAMPLES.qz*

echo "Transposing the ${MARKER} ${CELL} table"
qiime feature-table transpose \
  --i-table ${OUTPUT}/dada2_denoise/table_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --o-transposed-feature-table ${OUTPUT}/dada2_post_denoise/transposed_table_${MARKER}_${CELL}_ALLSAMPLES.qza

qiime metadata tabulate \
  --m-input-file ${OUTPUT}/dada2_post_denoise/transposed_table_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --o-visualization ${OUTPUT}/dada2_post_denoise/transposed_table_${MARKER}_${CELL}_ALLSAMPLES.qzv


# 8. Using the transposed table as feature metadata, and keep only the ASVs found in individual samples.

for SAMPLE in ${SAMPLES}
do
  rm -f ${OUTPUT}/dada2_post_denoise/asv_${MARKER}_${CELL}_${SAMPLE}.qz*
  rm -f ${OUTPUT}/dada2_post_denoise/table_${MARKER}_${CELL}_${SAMPLE}.qza
  rm -f ${OUTPUT}/dada2_post_denoise/fasta/asv_${MARKER}_${CELL}_${SAMPLE}.fasta

  # filter asv sequence table
  qiime feature-table filter-seqs \
    --i-data ${OUTPUT}/dada2_denoise/asv_${MARKER}_${CELL}_ALLSAMPLES.qza \
    --m-metadata-file ${OUTPUT}/dada2_post_denoise/transposed_table_${MARKER}_${CELL}_ALLSAMPLES.qza \
    --p-where ${SAMPLE}_${MARKER}_${CELL} \
    --o-filtered-data ${OUTPUT}/dada2_post_denoise/asv_${MARKER}_${CELL}_${SAMPLE}.qza

  qiime feature-table tabulate-seqs \
    --i-data ${OUTPUT}/dada2_post_denoise/asv_${MARKER}_${CELL}_${SAMPLE}.qza \
    --o-visualization ${OUTPUT}/dada2_post_denoise/asv_${MARKER}_${CELL}_${SAMPLE}.qzv

  # 9 Convert ASV sequences from qza to fasta format
  echo "Converting the ${MARKER} ${CELL} ${SAMPLE} ASVs qza file to fasta file"

  qiime tools export \
    --input-path ${OUTPUT}/dada2_post_denoise/asv_${MARKER}_${CELL}_${SAMPLE}.qza \
    --output-path ${OUTPUT}/dada2_post_denoise/fasta/

  mv ${OUTPUT}/dada2_post_denoise/fasta/dna-sequences.fasta ${OUTPUT}/dada2_post_denoise/fasta/asv_${MARKER}_${CELL}_${SAMPLE}.fasta
done
