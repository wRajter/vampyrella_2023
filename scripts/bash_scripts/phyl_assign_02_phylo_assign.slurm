#!/bin/bash

#SBATCH --mail-user=lrajter@uni-koeln.de
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=12
#SBATCH --mem=46g
#SBATCH --time=240:00:00
#SBATCH --account=ag-hess
#SBATCH --output=%x.o%A_%a

# Phylogenetic assignment of the environmental sequencies

# Variables
MARKER="Full18S"
CELL="cellCombined"
SIM="sim99"
REF_ALIGNMENT="reference_data/euk_ref_plus_vamp_18S_mafft_gblocks.phy"
REF_TREE="reference_data/euk_ref_plus_vamp_18S_mafft_gblocks.newick"
QUERY_SEQS="reference_data/otu_seqs_${MARKER}_${CELL}_${SIM}.fasta"

# Activate conda environment
# phylo_placment env contains these two packages:
    # raxml-ng=1.1.0
    # epa-ng=0.3.8
source /home/lrajter/miniconda3/etc/profile.d/conda.sh
conda activate phylo_placment


# Aligning query sequences based on the reference alignment and tree 
papara \
    -t ${REF_TREE} \
    -s ${REF_ALIGNMENT} \
    -q ${QUERY_SEQS} -r

# Splitting alignment
epa-ng --split \
        ${REF_ALIGNMENT} \
        papara_alignment.default

# Substitution model and its parameters
raxml-ng --evaluate \
         --msa reference.fasta \
         --tree ${REF_TREE} \
         --model GTR+G

# Phylogenetic placement
epa-ng \
    -t ${REF_TREE} \
    -s reference.fasta \
    -q query.fasta \
    --model reference.fasta.raxml.bestModel

