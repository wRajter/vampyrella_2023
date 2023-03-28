#!/bin/bash

# Creating qiime visualization and tsv table from the taxonomic assignemnt result


# Exporting the taxonomic assignment result to qiime visualization file
qiime metadata tabulate \
  --m-input-file vsearch_taxonomy.qza \
  --o-visualization sub_tax_viz.qzv

# Exporting the taxonomic assignment result to TSV table
qiime tools export \
  --input-path vsearch_taxonomy.qza \
  --output-path .
