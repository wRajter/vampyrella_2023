#!/bin/bash

# Extracting the 18S rRNA region from the environmental sequences using Qiime2
# Extracted 18S will serve for taxonomic assignemnt
# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
PROJECT="Jamy_2019"
MARKER="rDNA"
CELL="cell"
SIM="sim99"
F_PRIMER="AACCTGGTTGATCCTGCCAG"
R_PRIMER="TGATCCTTCTGCAGGTTCACCTAC"
IDENT_THRESHOLD=0.7 # minimum combined primer match identity threshold
RAW_DATA="../../raw_data"
OTU_CHIM_FILT="${RAW_DATA}/OTU_nonchimeric/${PROJECT}/${MARKER}/${CELL}/${SIM}"
EXTRACTED_18S="${RAW_DATA}/extracted_18S/${PROJECT}/${MARKER}/${CELL}/${SIM}"

############################
## FILTERING AND TRIMMING ##
############################

echo "Working on ${PROJECT} project, ${MARKER} marker, ${CELL} cell, and ${SIM} similarity."

# Cleaning
mkdir -p ${EXTRACTED_18S}/
rm -f ${EXTRACTED_18S}/extracted_18S_seqs.qza \
      ${EXTRACTED_18S}/extracted_18S_seqs.qzv \
      ${EXTRACTED_18S}/extracted_18S_seqs.fasta \
      ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.fasta \
      ${EXTRACTED_18S}/trimming.log \
      ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.qza


# Trim the sequences based on the primers
cutadapt -a ${R_PRIMER} \
         -g ${F_PRIMER} \
         --minimum-length 1 \
         --trim-n \
         --discard-untrimmed \
         -o ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.fasta \
            ${OTU_CHIM_FILT}/otu_seqs_nonchimeric.fasta \
          > ${EXTRACTED_18S}/trimming.log


qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.fasta \
  --output-path ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.qza

# Visalize extracted sequences and convert them into Fasta files
qiime feature-table tabulate-seqs \
  --i-data ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.qza \
  --o-visualization ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.qzv
