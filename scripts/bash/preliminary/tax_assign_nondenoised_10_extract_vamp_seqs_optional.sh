#!/bin/bash

# Extracting vampyrella-specific sequences after second chimera filtering


# Variables
NCORES=12
PROJECT="Jamy_2022"
MARKER="rDNA"
CELL="cell"
SIM="sim99"
RAW_DATA="../../raw_data"
CLUST_DIR="${RAW_DATA}/OTU_clust/${PROJECT}/${MARKER}/${CELL}/${SIM}"
OTU_CHIM_FILT2="${RAW_DATA}/OTU_nonchimeric_second/${PROJECT}/${MARKER}/${CELL}/${SIM}"
VAMP_SEQ_DIR="${RAW_DATA}/vamp_specific_seqs/${PROJECT}/${MARKER}/${CELL}/${SIM}"
ASSIGNMENT_DIR="${RAW_DATA}/tax_assign_results/${PROJECT}/${MARKER}/${CELL}/${SIM}"
SEQTK="${RAW_DATA}/packages/seqtk"




##########################
## PER SAMPLE (VSEARCH) ##
##########################


# Get sample names for looping
SAMPLES=$(ls ${CLUST_DIR}/*.fasta | \
          awk -F '/' '{ print $NF }' | \
          awk -F '_' '{ print $3 }' |
          awk -F '.' '{ print $1 }')


mkdir -p ${VAMP_SEQ_DIR}

# Loop through the samples
for SAMPLE in ${SAMPLES}
do
  # Cleaning
  rm -f ${VAMP_SEQ_DIR}/otu_seqs_vamp_nonchim_${SAMPLE}.fasta

  # Pulling the vampyrellids out
  grep 'Vampyrellida' ${ASSIGNMENT_DIR}/blast6_${SAMPLE}.tab | \
  awk '{print $1}' > \
  ${ASSIGNMENT_DIR}/vamp_ids_${SAMPLE}.txt

  ${SEQTK}/seqtk subseq \
    ${OTU_CHIM_FILT2}/otu_seqs_${SAMPLE}.fasta \
    ${ASSIGNMENT_DIR}/vamp_ids_${SAMPLE}.txt > \
    ${VAMP_SEQ_DIR}/otu_seqs_vamp_nonchim_${SAMPLE}.fasta

  rm -f ${ASSIGNMENT_DIR}/vamp_ids_${SAMPLE}.txt
done
