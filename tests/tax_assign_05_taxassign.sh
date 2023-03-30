#!/bin/bash

# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
PROJECT="Jamy_2019"
MARKER="rDNA"
CELL="cell"
SIM="sim99"
QUERY_SEQS="../raw_data/extracted_18S/Jamy_2019/rDNA/cell/sim99/extracted_18S_seqs_trimmed.qza"


# Taxonomic assignment analysis
qiime feature-classifier classify-consensus-vsearch \
  --i-query sub_seqs.qza \
  --i-reference-reads sub_ref_alignment.qza \
  --i-reference-taxonomy sub_ref_taxonomy.qza \
  --o-classification vsearch_taxonomy.qza \
  --p-threads 12 \
  --o-search-results search_results_blast6 \
  --verbose
  # --p-perc-identity 0.6 \
  # --p-query-cov 0.6 \
