#!/bin/bash
#SBATCH --mail-user=lrajter@uni-koeln.de
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=12
#SBATCH --mem=46g
#SBATCH --time=240:00:00
#SBATCH --account=ag-hess
#SBATCH --output=%x.o%A_%a

source /home/lrajter/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2


for SAMPLE in "A3_18S_cell1" \
              "Mock_18S_cell1" \
              "NH1_18S_cell1" \
              "NH4_18S_cell1" \
              "Sim17_18S_cell1" \
              "Sim22_18S_cell1" \
              "Th16_18S_cell1" \
              "Th38_18S_cell1" \
              "Th40_18S_cell1" \
              "X17007_18S_cell1";
  do
    qiime feature-classifier classify-consensus-vsearch \
      --i-query representative_sequences_${SAMPLE}.qza \
      --i-reference-reads ref_sequences.qza \
      --i-reference-taxonomy ref_taxonomy.qza \
      --o-classification vsearch_taxonomy_${SAMPLE}.qza \
      --p-threads 12 \
      --verbose
