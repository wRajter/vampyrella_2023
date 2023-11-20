#!/bin/bash

# Phylogenetic assignment of the environmental sequencies

# Variables
PROJECT="Suthaus_2022"
MARKER="Full18S"
RAW_DATA="../../raw_data"
RESULTS="../../results"
SIM="sim_90"
DENOISE_METHOD="RAD"
TAXON="eukaryotes"
# EUKARYOTES variables
QUERY_DIR="${RAW_DATA}/chimera_filtered/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}"
REF_ALIGNMENT_EUK="${RAW_DATA}/reference_alignments/vamp_phylo_placement/eukaryotes/reference_alignment_2022/reference_alignment.phy"
REF_TREE_EUK="${RAW_DATA}/phyl_placement/reference_trees/eukaryotes/reference_tree_2022/T2.raxml.bestTree"
LOG_FILES_EUK="${RAW_DATA}/phyl_placement/jplace_log_files/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}/eukaryotes"
JPLACE_DIR_EUK="${RESULTS}/phyl_placement/jplace/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}/eukaryotes"
TAXON_FILE_EUK="${RAW_DATA}/reference_alignments/vamp_phylo_placement/eukaryotes/reference_alignment_2022/taxon_vamp.tsv"
VAMP_FASTA_FILES="${RESULTS}/phyl_placement/vamp_specific/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}"
# VAMPYRELLIDS variables
REF_ALIGNMENT_VAMP="${RAW_DATA}/reference_alignments/vamp_phylo_placement/vampyrellida/reference_alignment_2023/reference_alignment.phy"
REF_TREE_VAMP="${RAW_DATA}/phyl_placement/reference_trees/vampyrellida/reference_tree_2023/T2.raxml.bestTree"
LOG_FILES_VAMP="${RAW_DATA}/phyl_placement/jplace_log_files/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}/vampyrellida"
JPLACE_DIR_VAMP="${RESULTS}/phyl_placement/jplace/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}/vampyrellida"
TAXON_FILE_VAMP="${RAW_DATA}/reference_alignments/vamp_phylo_placement/vampyrellida/reference_alignment_2023/taxon_file.tsv"

# Samples
SAMPLES=$(ls ${QUERY_DIR}/*.fasta | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }' | \
          awk -F '_' '{ print $1 "_" $2 }')


echo "Samples used:"
echo "$SAMPLES"



#########################################
# STEP 1: Phyl. placement of eukaryotes #
#########################################

for SAMPLE in ${SAMPLES}
do

  echo "Working on sample ${SAMPLE}."

  mkdir -p ${LOG_FILES_EUK}/
  mkdir -p ${JPLACE_DIR_EUK}/


  # Aligning query sequences based on the reference alignment and tree
  echo "Aligning..."
  papara -t ${REF_TREE_EUK} \
         -s ${REF_ALIGNMENT_EUK} \
         -q ${QUERY_DIR}/${SAMPLE}_otu.fasta -r

  # Splitting alignment
  epa-ng --split \
          ${REF_ALIGNMENT_EUK} \
          papara_alignment.default

  # Substitution model and its parameters
  echo "Computing substitution model..."
  raxml-ng --evaluate \
          --msa reference.fasta \
          --tree ${REF_TREE_EUK} \
          --model GTR+G

  # Phylogenetic placement
  echo "Phylogenetic placement..."
  epa-ng \
      -t ${REF_TREE_EUK} \
      -s reference.fasta \
      -q query.fasta \
      --model reference.fasta.raxml.bestModel

  # Move all the log files to the log directory
  mv epa_info.log ${LOG_FILES_EUK}/${SAMPLE}_epa_info.log
  mv papara_log.default ${LOG_FILES_EUK}/${SAMPLE}_papara_log.default
  mv reference.fasta.raxml.bestModel ${LOG_FILES_EUK}/${SAMPLE}_reference.fasta.raxml.bestModel
  mv reference.fasta.raxml.log ${LOG_FILES_EUK}/${SAMPLE}_reference.fasta.raxml.log

  # Move the jplace file to the result directory
  mv epa_result.jplace ${JPLACE_DIR_EUK}/epa_result_${SAMPLE}.jplace

  # Delete rest of the files
  rm -f epa_* \
        papara_* \
        reference.fasta \
        reference.fasta.raxml*


  # EXTRACTING VAMPYRELLIDA OTUs
  mkdir -p ${VAMP_FASTA_FILES}/

  gappa prepare extract \
    --jplace-path ${JPLACE_DIR_EUK}/epa_result_${SAMPLE}.jplace \
    --clade-list-file ${TAXON_FILE_EUK} \
    --fasta-path ${QUERY_DIR}/${SAMPLE}_otu.fasta \
    --allow-file-overwriting \
    --log-file ${LOG_FILES_EUK}/${SAMPLE}_extract_otus.log

  mv sequences/Vampyrellida.fasta ${VAMP_FASTA_FILES}/${SAMPLE}_otu.fasta

  # Remove query sequences
  rm query.fasta
  rm -r samples
  rm -r sequences

done

###########################################
# STEP 2: Phyl. placement of vampyrellids #
###########################################

echo "Samples used:"
echo "$SAMPLES"

for SAMPLE in ${SAMPLES}
do

  echo "Working on sample ${SAMPLE}."

  mkdir -p ${LOG_FILES_VAMP}/
  mkdir -p ${JPLACE_DIR_VAMP}/


  # Aligning query sequences based on the reference alignment and tree
  echo "Aligning..."
  papara -t ${REF_TREE_VAMP} \
         -s ${REF_ALIGNMENT_VAMP} \
         -q ${VAMP_FASTA_FILES}/${SAMPLE}_otu.fasta -r

  # Splitting alignment
  epa-ng --split \
          ${REF_ALIGNMENT_VAMP} \
          papara_alignment.default

  # Substitution model and its parameters
  echo "Computing substitution model..."
  raxml-ng --evaluate \
          --msa reference.fasta \
          --tree ${REF_TREE_VAMP} \
          --model GTR+G

  # Phylogenetic placement
  echo "Phylogenetic placement..."
  epa-ng \
      -t ${REF_TREE_VAMP} \
      -s reference.fasta \
      -q query.fasta \
      --model reference.fasta.raxml.bestModel

  # Move all the log files to the log directory
  mv epa_info.log ${LOG_FILES_VAMP}/${SAMPLE}_epa_info.log
  mv papara_log.default ${LOG_FILES_VAMP}/${SAMPLE}_papara_log.default
  mv reference.fasta.raxml.bestModel ${LOG_FILES_VAMP}/${SAMPLE}_reference.fasta.raxml.bestModel
  mv reference.fasta.raxml.log ${LOG_FILES_VAMP}/${SAMPLE}_reference.fasta.raxml.log

  # Move the jplace file to the result directory
  mv epa_result.jplace ${JPLACE_DIR_VAMP}/epa_result_${SAMPLE}.jplace

  # Delete rest of the files
  rm -f epa_* \
        papara_* \
        reference.fasta \
        reference.fasta.raxml* \
        query.fasta

done
