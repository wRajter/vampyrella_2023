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

qiime feature-classifier classify-consensus-vsearch \
  --i-query representative_sequences_NH1_18S_cell1.qza \
  --i-reference-reads ref_sequences.qza \
  --i-reference-taxonomy ref_taxonomy.qza \
  --o-classification vsearch_taxonomy.qza \
  --p-threads 12
  --output-dir assignment_results \
  --verbose
