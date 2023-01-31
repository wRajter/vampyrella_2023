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

source /home/lrajter/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2022.2


for SAMPLE in ${SAMPLES}
do
  qiime feature-classifier classify-consensus-vsearch \
    --i-query representative_sequences_${MARKER}_${CELL}_${SAMPLE}.qza \
    --i-reference-reads ref_sequences.qza \
    --i-reference-taxonomy ref_taxonomy.qza \
    --o-classification vsearch_taxonomy_${MARKER}_${CELL}_${SAMPLE}.qza \
    --p-threads 12 \
    --verbose
done
