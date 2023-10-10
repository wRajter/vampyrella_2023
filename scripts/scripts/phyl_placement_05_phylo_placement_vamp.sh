#!/bin/bash


# Phylogenetic assignment of the environmental sequencies

# Variables
PROJECT="Suthaus_2022"
TAXON="vampyrellida"
RAW_DATA="../../raw_data"
PHYL_PLAC_DIR="${RAW_DATA}/phyl_placement/${PROJECT}"
REF_VERSION="2023"
REF_ALIGNMENT="${RAW_DATA}/reference_alignments/vamp_phylo_placement/${TAXON}/reference_alignment_${REF_VERSION}/reference_alignment.phy"
REF_TREE="${RAW_DATA}/phyl_placement/reference_trees/${TAXON}/reference_tree_${REF_VERSION}/T2.raxml.bestTree"
QUERY_SEQS="${PHYL_PLAC_DIR}/eukaryotes/downstream_analyses/extract_otus/allsamples_sequences/Vampyrellida.fasta"
OUT_DIR="${PHYL_PLAC_DIR}/${TAXON}/phyl_placement_analysis/allsamples"


# Activate conda phylo_placement environment that should contain these three packages:
    # raxml-ng=1.1.0
    # epa-ng=0.3.8
    # papara=2.5


mkdir -p ${OUT_DIR}/

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

mv epa_* \
   papara_* \
   query.fasta \
   reference.fasta \
   reference.fasta.raxml* \
   ${OUT_DIR}/
