#!/bin/bash

# Assign taxonomy to ASVs

# created based on:
# https://github.com/LangilleLab/microbiome_helper/wiki/PacBio-CCS-Amplicon-SOP-v2-(qiime2-2022.2)

# taxonomic classifier downloaded from:
# https://docs.qiime2.org/2022.2/data-resources/

# 1. Running taxonomic classification

qiime feature-classifier classify-sklearn \
   --i-reads dada2_output/representative_sequences.qza \
   --i-classifier ../raw_data/tax_assignment/silva-138-99-nb-weighted-classifier.qza \
   --p-n-jobs 1 \
   --output-dir taxa


# Vsearch taxonomy assignment through qiime

qiime feature-classifier classify-consensus-vsearch \
  --i-query /home/lubo/code/wRajter/vampyrella_2023/scripts2/dada2_output/representative_sequences.qza \
  --i-reference-reads ref_sequences.qza \
  --i-reference-taxonomy ref_taxonomy.qza \
  --o-classification vsearch_taxonomy.qza \
  --p-threads 12
  --output-dir assignment_results \
  --verbose

qiime metadata tabulate \
  --m-input-file vsearch_taxonomy.qza \
  --o-visualization assignment_results/tax_viz.qzv

# exporting the taxonomic assignment result to TSV table
qiime tools export \
   --input-path taxa/classification.qza --output-path taxa


# Filtering feature tables
qiime taxa filter-table \
  --i-table ../dada2_output/table.qza \
  --i-taxonomy vsearch_taxonomy.qza \
  --p-include "k__Eukaryota;d__Rhizaria;p__Cercozoa;c__Endomyxa;o__Vampyrellida" \
  --o-filtered-table vamp_table.qza


# Barplot
qiime taxa barplot \
  --i-table vamp_table.qza \
  --i-taxonomy vsearch_taxonomy.qza \
  --o-visualization eukaryota_barplot.qzv


# Sample filtering

# qiime feature-table filter-samples \
#   --i-table ../../dada2_output/table.qza \
#   --m-metadata-file A3_18S_cell1.tsv \
#   --o-filtered-table A3_18S_cell1_filtered_table.qza


# qiime taxa filter-table \
#   --i-table A3_18S_cell1_filtered_table.qza \
#   --i-taxonomy ../taxonomy.qza \
#   --p-include "k__Eukaryota" \
#   --o-filtered-table A3_18S_cell1_taxonomy.qza


# qiime metadata tabulate \
# --m-input-file A3_18S_cell1_taxonomy.qza \
# --o-visualization A3_18S_cell1_tax_viz.qzv
