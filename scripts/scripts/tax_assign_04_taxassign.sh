#!/bin/bash

# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
PROJECT="Jamy_2019"
MARKER="rDNA"
CELL="cell"
SIM="sim99"
RAW_DATA="../../raw_data"
QUERY_SEQ="${RAW_DATA}/OTU_clust/${PROJECT}/${MARKER}/${CELL}/${SIM}/otu_seqs.qza"
REF_SEQ="${RAW_DATA}/reference_alignments/pr2/ref_sequences.qza"
REF_TAX="${RAW_DATA}/reference_alignments/pr2/ref_taxonomy.qza"
ASSIGNMENT_DIR="${RAW_DATA}/tax_assign_results/${PROJECT}/${MARKER}/${CELL}/${SIM}"

# Taxonomic assignment analysis
qiime feature-classifier classify-consensus-vsearch \
  --i-query ${QUERY_SEQ} \
  --i-reference-reads ${REF_SEQ} \
  --i-reference-taxonomy ${REF_TAX} \
  --output-dir ${ASSIGNMENT_DIR} \
  --o-classification vsearch_taxonomy.qza \
  --p-threads 12 \
  --verbose
