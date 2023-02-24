#!/bin/bash


# Phylogenetic assignment of the environmental sequencies

# Variables
MARKER="Full18S"
CELL="cellCombined"
RAW_DATA="../../raw_data"
SIM="sim99"
TAXON="eukaryotes"
QUERY_DIR="${RAW_DATA}/per_sample_results/${MARKER}/${CELL}/${SIM}/fasta/"
REF_ALIGNMENT="${RAW_DATA}/reference_alignments/euk_ref/euk_ref_plus_vamp_18S_mafft_gblocks.phy"
REF_TREE="${RAW_DATA}/phyl_placement/${TAXON}/reference_tree/T2.raxml.bestTree"
RAW_READS_DIR="${RAW_DATA}/PacBio/Suthaus${MARKER}/${CELL}"
PLACEMENT_DIR="${RAW_DATA}/phyl_placement/${TAXON}/phyl_placement_analysis"
SAMPLES=$(ls ${RAW_READS_DIR}/*reads.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '_' '{ print $1 }')


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
    -q ${QUERY_DIR}/otu_seqs_filtered_${SAMPLE}.fasta -r

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
