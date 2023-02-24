#!/bin/bash

#SBATCH --mail-user=lrajter@uni-koeln.de
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=6
#SBATCH --mem=46g
#SBATCH --time=240:00:00
#SBATCH --account=ag-hess
#SBATCH --output=%x.o%A_%a

# Creating a reference tree based on the reference alignment using RAxML

# Variables
REF_ALIGNMENT="reference_data/${RAW_DATA}"
MODEL="GTR+G"
NCORES=6

# Activate conda environment
source /home/lrajter/miniconda3/etc/profile.d/conda.sh
conda activate phylo_placment

# Check alignment
raxml-ng --check \
         --msa ${REF_ALIGNMENT} \
         --model ${MODEL} \
         --prefix T1

# Compute tree
raxml-ng \
    --msa ${REF_ALIGNMENT} \
    --model ${MODEL} \
    --prefix T2 \
    --threads ${NCORES}
