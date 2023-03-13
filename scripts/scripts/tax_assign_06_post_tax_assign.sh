#!/bin/bash

# Creating qiime visualization and tsv table from the taxonomic assignemnt result

# Variables
PROJECT="Suthaus_2022"
MARKER="Full18S"
CELL="cellCombined"
SIM="sim97"
RAW_DATA="../../raw_data"
ASSIGNMENT_DIR="${RAW_DATA}/tax_assign_results/${PROJECT}/${MARKER}/${CELL}/${SIM}"


# Cleaning
rm -f ${ASSIGNMENT_DIR}/tax_viz.qzv \
      ${ASSIGNMENT_DIR}/taxonomy.tsv

# Exporting the taxonomic assignment result to qiime visualization file
qiime metadata tabulate \
  --m-input-file ${ASSIGNMENT_DIR}/vsearch_taxonomy.qza \
  --o-visualization ${ASSIGNMENT_DIR}/tax_viz.qzv

# Exporting the taxonomic assignment result to TSV table
qiime tools export \
  --input-path ${ASSIGNMENT_DIR}/vsearch_taxonomy.qza \
  --output-path ${ASSIGNMENT_DIR}/
