#!/bin/bash

# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
PROJECT="Jamy_2019"
MARKER="rDNA"
CELL="cell"
SIM="sim99"


# Taxonomic assignment analysis
qiime feature-classifier classify-consensus-vsearch \
  --i-query sub_seqs.qza \
  --i-reference-reads sub_ref_alignment.qza \
  --i-reference-taxonomy sub_ref_taxonomy.qza \
  --o-classification vsearch_taxonomy.qza \
  --p-perc-identity 0.6 \
  --p-query-cov 0.6 \
  --p-threads 12 \
  --o-search-results search_results_blast6 \
  --verbose
