#!/bin/bash
#SBATCH --mail-user=lrajter@uni-koeln.de
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=12
#SBATCH --mem=46g
#SBATCH --time=240:00:00
#SBATCH --account=ag-hess
#SBATCH --output=%x.o%A_%a

# Activate your Conda environment
source /home/lrajter/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2

# Variables
PROJECT="Jamy_2022"
MARKER="rDNA"
CELL="cell"
READS_GZA_DIR="reads_qza"
MANIFEST_FILE="PacBioCCSmanifest_rDNA_cell.tsv"
DENOISE_DIR="denoised"


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
