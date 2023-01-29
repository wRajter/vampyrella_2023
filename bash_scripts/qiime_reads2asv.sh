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
METADATA="${PROJECT_FILE}/raw_data/qiime_input/metadata_${MARKER}_${CELL}.tsv"
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


# 1. Transforming the metadata text file to the qiime artifact file

mkdir -p ${OUTPUT}/metadata_qzv
rm -f ${OUTPUT}/metadata_qzv/tabulated_metadata_${MARKER}_${CELL}.qzv

echo "transforming the metadata_${MARKER}_${CELL}.tsv file into the qiime QZV artifact file"
qiime metadata tabulate \
  --m-input-file ${METADATA} \
  --o-visualization ${OUTPUT}/metadata_qzv/tabulated_metadata_${MARKER}_${CELL}.qzv



# 2. Importing FASTQs as QIIME2 artifact
mkdir -p ${OUTPUT}/reads_qza
rm -f ${OUTPUT}/reads_qza/raw_reads_${MARKER}_${CELL}.qza
rm -f ${OUTPUT}/reads_qza/raw_reads_summary_${MARKER}_${CELL}.qzv

echo "Importing ${MARKER} ${CELL} reads to QIIME2 artifact"
qiime tools import \
    --type SampleData[SequencesWithQuality] \
    --input-path ${MANIFEST} \
    --output-path ${OUTPUT}/reads_qza/raw_reads_${MARKER}_${CELL}.qza \
    --input-format SingleEndFastqManifestPhred33V2

# 3. Creating reads summary file
echo "Creating ${MARKER} ${CELL} reads summary file"
qiime demux summarize \
   --i-data ${OUTPUT}/reads_qza/raw_reads_${MARKER}_${CELL}.qza \
   --o-visualization ${OUTPUT}/reads_qza/raw_reads_summary_${MARKER}_${CELL}.qzv



# ########################################################
# ## DENOISING THE READS INTO AMPLICON SEQUENCE VARIANT ##
# ########################################################

# 4. Running DADA2
mkdir -p ${OUTPUT}/dada2_output
rm -f ${OUTPUT}/dada2_output/*_${MARKER}_${CELL}.gz*

echo "denoising ${MARKER} ${CELL} reads and creating ASVs using DADA2"
qiime dada2 denoise-ccs --i-demultiplexed-seqs ${OUTPUT}/reads_qza/raw_reads_${MARKER}_${CELL}.qza \
 --p-n-reads-learn ${LEARN} \
 --p-min-len 800 --p-max-len 3000 \
 --p-n-threads ${NCORES} \
 --p-front ${F_PRIMER} --p-adapter ${R_PRIMER} \
 --o-table ${OUTPUT}/dada2_output/table_${MARKER}_${CELL}.qza \
 --o-representative-sequences ${OUTPUT}/dada2_output/representative_sequences_${MARKER}_${CELL}_ALLSAMPLES.qza \
 --o-denoising-stats ${OUTPUT}/dada2_output/stats_${MARKER}_${CELL}.qza \
 --verbose


# 5. Summarizing DADA2 output
echo "creating ${MARKER} ${CELL} ASVs summary"

qiime feature-table summarize \
 --i-table ${OUTPUT}/dada2_output/table_${MARKER}_${CELL}.qza \
 --o-visualization ${OUTPUT}/dada2_output/dada2_table_summary_${MARKER}_${CELL}.qzv

# 6. Visualization of the representative sequences
echo "creating ${MARKER} ${CELL} ASVs visualization"

qiime feature-table tabulate-seqs \
  --i-data ${OUTPUT}/dada2_output/representative_sequences_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --o-visualization ${OUTPUT}/dada2_output/representative_sequences_${MARKER}_${CELL}_ALLSAMPLES.qzv

# 7 Convert ASV sequences from qza to fasta format
echo "Converting the ${MARKER} ${CELL} ASVs qza file to fasta file"

mkdir -p ${OUTPUT}/ASV
rm -f ${OUTPUT}/ASV/representative_sequences_${MARKER}_${CELL}_ALLSAMPLES.fasta

qiime tools export \
  --input-path ${OUTPUT}/dada2_output/representative_sequences_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --output-path ${OUTPUT}/ASV/

mv ${OUTPUT}/ASV/dna-sequences.fasta ${OUTPUT}/ASV/representative_sequences_${MARKER}_${CELL}_ALLSAMPLES.fasta


########################################################
## FILTERING REPRESENTATIVE SEQUENCES FOR EACH SAMPLE ##
########################################################

# Getting a fasta file with the representative sequences for each sample


# 8. Transposing the table (this will flip the sample and feature axes, necessary for the next step).
rm -f ${OUTPUT}/dada2_output/transposed_table_${MARKER}_${CELL}.qz*

echo "Transposing the ${MARKER} ${CELL} table"
qiime feature-table transpose \
  --i-table ${OUTPUT}/dada2_output/table_${MARKER}_${CELL}.qza \
  --o-transposed-feature-table ${OUTPUT}/dada2_output/transposed_table_${MARKER}_${CELL}.qza

qiime metadata tabulate \
  --m-input-file ${OUTPUT}/dada2_output/transposed_table_${MARKER}_${CELL}.qza \
  --o-visualization ${OUTPUT}/dada2_output/transposed_table_${MARKER}_${CELL}.qzv


# 9. Using the transposed as feature metadata, and keep only the features found in samples 1L or 2L:

for SAMPLE in ${SAMPLES}
do
  rm -f representative_sequences_${MARKER}_${CELL}_${SAMPLE}.qzv
  mkdir -p ${OUTPUT}/ASV
  rm -f ${OUTPUT}/dada2_output/ASV/representative_sequences_${MARKER}_${CELL}_${SAMPLE}.fasta

  qiime feature-table filter-seqs \
    --i-data ${OUTPUT}/dada2_output/representative_sequences_${MARKER}_${CELL}_ALLSAMPLES.qza \
    --m-metadata-file ${OUTPUT}/dada2_output/transposed_table_${MARKER}_${CELL}.qza \
    --p-where ${SAMPLE}_${MARKER}_${CELL} \
    --o-filtered-data ${OUTPUT}/dada2_output/representative_sequences_${MARKER}_${CELL}_${SAMPLE}.qza

  qiime feature-table tabulate-seqs \
    --i-data ${OUTPUT}/dada2_output/representative_sequences_${MARKER}_${CELL}_${SAMPLE}.qza \
    --o-visualization ${OUTPUT}/dada2_output/representative_sequences_${MARKER}_${CELL}_${SAMPLE}.qzv

  # 10 Convert ASV sequences from qza to fasta format
  echo "Converting the ${MARKER} ${CELL} ${SAMPLE} ASVs qza file to fasta file"

  qiime tools export \
    --input-path ${OUTPUT}/dada2_output/representative_sequences_${MARKER}_${CELL}_${SAMPLE}.qza \
    --output-path ${OUTPUT}/ASV/

  mv ${OUTPUT}/ASV/dna-sequences.fasta ${OUTPUT}/ASV/representative_sequences_${MARKER}_${CELL}_${SAMPLE}.fasta
done
