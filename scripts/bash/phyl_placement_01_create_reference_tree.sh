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
TAXON="pr2"
RAW_DATA="../../raw_data"
REF_VERSION="gblocks2"
REF_ALIGNMENT="${RAW_DATA}/reference_alignments/vamp_phylo_placement/${TAXON}/reference_alignment_${REF_VERSION}/reference_alignment.fasta"
MODEL="GTR+G"
NCORES=12
TREE_VERSION="2023"
OTU_DIR="${RAW_DATA}/phyl_placement/reference_trees/${TAXON}/reference_tree_${REF_VERSION}"


mkdir -p ${OTU_DIR}

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


mv T1* T2* ${OTU_DIR}/
