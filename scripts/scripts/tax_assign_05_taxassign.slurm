#!/bin/bash

#SBATCH --mail-user=lrajter@uni-koeln.de
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=12
#SBATCH --mem=46g
#SBATCH --time=240:00:00
#SBATCH --account=ag-hess
#SBATCH --output=%x.o%A_%a


# Variables
MARKER="Full18S"
CELL="cell2"

source /home/lrajter/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2

qiime feature-classifier classify-consensus-vsearch \
  --i-query asv_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --i-reference-reads ref_sequences.qza \
  --i-reference-taxonomy ref_taxonomy.qza \
  --o-classification vsearch_taxonomy_${MARKER}_${CELL}_ALLSAMPLES.qza \
  --p-threads 12 \
  --verbose
