#!/bin/bash


# Phylogenetic assignment of the environmental sequencies

# Variables
PROJECT="Suthaus_2022"
MARKER="Full18S"
CELL="cellCombined"
RAW_DATA="../../raw_data"
SIM="sim90"
DENOISE_METHOD="RAD"
TAXON="eukaryotes"
QUERY_DIR="${RAW_DATA}/OTU_nonchimeric/${PROJECT}/${MARKER}/${CELL}/${SIM}/${DENOISE_METHOD}"
REF_VERSION="2022"
REF_ALIGNMENT="${RAW_DATA}/reference_alignments/vamp_phylo_placement/${TAXON}/reference_alignment_${REF_VERSION}/reference_alignment.phy"
REF_TREE="${RAW_DATA}/phyl_placement/reference_trees/${TAXON}/reference_tree_${REF_VERSION}/T2.raxml.bestTree"
PLACEMENT_DIR="${RAW_DATA}/phyl_placement/${PROJECT}/${TAXON}/phyl_placement_analysis"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"

SAMPLES=$(ls ${QUERY_DIR} | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }' | \
          grep -v 'Mock_18S_otu') # If you want to remove mock community from the samples


echo "Samples used:"
echo "$SAMPLES"

# Activate conda phylo_placement environment that should contain these three packages:
    # raxml-ng=1.1.0
    # epa-ng=0.3.8
    # papara=2.5


for SAMPLE in ${SAMPLES}
do

  echo "Working on sample ${SAMPLE}."

  mkdir -p ${PLACEMENT_DIR}/${SAMPLE}/


  # Aligning query sequences based on the reference alignment and tree
  echo "Aligning..."
  papara \
    -t ${REF_TREE} \
    -s ${REF_ALIGNMENT} \
    -q ${QUERY_DIR}/${SAMPLE}.fasta -r

  # Splitting alignment
  epa-ng --split \
          ${REF_ALIGNMENT} \
          papara_alignment.default

  # Substitution model and its parameters
  echo "Computing substitution model..."
  raxml-ng --evaluate \
          --msa reference.fasta \
          --tree ${REF_TREE} \
          --model GTR+G

  # Phylogenetic placement
  echo "Phylogenetic placement..."
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
     ${PLACEMENT_DIR}/${SAMPLE}/

done
