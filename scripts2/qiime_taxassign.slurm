#!/bin/bash

#SBATCH --mail-user=lrajter@uni-koeln.de
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=12
#SBATCH --mem=46g
#SBATCH --time=240:00:00
#SBATCH --account=ag-hess
#SBATCH --output=%x.o%A_%a


qiime feature-classifier classify-consensus-vsearch \
  --i-query representative_sequences.qza \
  --i-reference-reads ref_sequences.qza \
  --i-reference-taxonomy ref_taxonomy.qza \
  --o-classification vsearch_taxonomy_allsamples.qza \
  --p-threads 12 \
  --verbose
