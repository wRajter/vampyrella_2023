#!/bin/bash

# Creating qiime visualization and tsv table from the taxonomic assignemnt results

# Variables
PROJECT_DIR="/home/lubo/code/wRajter/vampyrella_2023"
ASSIGNMENT_DIR="${PROJECT_DIR}/raw_data/qiime_output/assignment_results"
CELL="cell1"
MARKER="Full18S"
SAMPLES="ALLSAMPLES \
         A3 \
         Mock \
         NH1 \
         NH4 \
         Sim17 \
         Sim22 \
         Th16 \
         Th38 \
         Th40 \
         X17007"


for SAMPLE in ${SAMPLES}
do
  # exporting the taxonomic assignment result to qiime visualization file
  rm -f ${ASSIGNMENT_DIR}/${MARKER}_${CELL}/tax_viz_${SAMPLE}.qzv
  qiime metadata tabulate \
    --m-input-file ${ASSIGNMENT_DIR}/${MARKER}_${CELL}/vsearch_taxonomy_${MARKER}_${CELL}_${SAMPLE}.qza \
    --o-visualization ${ASSIGNMENT_DIR}/${MARKER}_${CELL}/tax_viz_${SAMPLE}.qzv

  # exporting the taxonomic assignment result to TSV table
  qiime tools export \
    --input-path ${ASSIGNMENT_DIR}/${MARKER}_${CELL}/vsearch_taxonomy_${MARKER}_${CELL}_${SAMPLE}.qza \
    --output-path ${ASSIGNMENT_DIR}/${MARKER}_${CELL}/
  mv ${ASSIGNMENT_DIR}/${MARKER}_${CELL}/taxonomy.tsv ${ASSIGNMENT_DIR}/${MARKER}_${CELL}/taxonomy_${SAMPLE}.tsv
done




  # # Filtering feature tables
  # qiime taxa filter-table \
  #   --i-table ../dada2_output/table.qza \
  #   --i-taxonomy vsearch_taxonomy.qza \
  #   --p-include "k__Eukaryota;d__Rhizaria;p__Cercozoa;c__Endomyxa;o__Vampyrellida" \
  #   --o-filtered-table vamp_table.qza


  # # Barplot
  # qiime taxa barplot \
  #   --i-table vamp_table.qza \
  #   --i-taxonomy vsearch_taxonomy.qza \
  #   --o-visualization eukaryota_barplot.qzv


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
