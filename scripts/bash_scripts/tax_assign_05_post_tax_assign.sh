#!/bin/bash

# Creating qiime visualization and tsv table from the taxonomic assignemnt result

# Variables
CELL="cellCombined"
MARKER="Full18S"
SIM="sim99"
RAW_DATA="../../raw_data"
ASSIGNMENT_DIR="${RAW_DATA}/tax_assign_results/${MARKER}/${CELL}/${SIM}"


# Exporting the taxonomic assignment result to qiime visualization file
rm -f ${ASSIGNMENT_DIR}/${MARKER}_${CELL}/tax_viz_${SAMPLE}.qzv
qiime metadata tabulate \
  --m-input-file ${ASSIGNMENT_DIR}/vsearch_taxonomy.qza \
  --o-visualization ${ASSIGNMENT_DIR}/tax_viz.qzv

# Exporting the taxonomic assignment result to TSV table
qiime tools export \
  --input-path ${ASSIGNMENT_DIR}/vsearch_taxonomy.qza \
  --output-path ${ASSIGNMENT_DIR}/
