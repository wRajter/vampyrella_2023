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
