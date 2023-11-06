#!/bin/bash

# Phylogenetic assignment of the environmental sequencies

# Variables
PROJECT="Suthaus_2022"
RAW_DATA="../../raw_data"
TAXON="vampyrellida"
MARKER="Full18S"
DENOISE_METHOD="RAD"
SIM="sim_90"
REF_VERSION="2023"
PHYL_PLAC_DIR="${RAW_DATA}/phyl_placement/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}"
REF_ALIGNMENT="${RAW_DATA}/reference_alignments/vamp_phylo_placement/${TAXON}/reference_alignment_${REF_VERSION}/reference_alignment.phy"
REF_TREE="${RAW_DATA}/phyl_placement/reference_trees/${TAXON}/reference_tree_${REF_VERSION}/T2.raxml.bestTree"
QUERY_DIR="${PHYL_PLAC_DIR}/eukaryotes/downstream_analyses/extract_otus"
PLACEMENT_DIR="${PHYL_PLAC_DIR}/${TAXON}/phyl_placement_analysis"


# samples
SAMPLES=$(ls ${QUERY_DIR} | \
          awk -F '_' '{ print $1 "_" $2 }' |
          sort -u)


echo "Samples used:"
echo "$SAMPLES"


for SAMPLE in ${SAMPLES}
do
  echo "Working on sample ${SAMPLE}."

  mkdir -p ${PLACEMENT_DIR}/${SAMPLE}/


  # Aligning query sequences based on the reference alignment and tree
  echo "Aligning..."
  papara \
    -t ${REF_TREE} \
    -s ${REF_ALIGNMENT} \
    -q ${QUERY_DIR}/${SAMPLE}_sequences/Vampyrellida.fasta -r

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
