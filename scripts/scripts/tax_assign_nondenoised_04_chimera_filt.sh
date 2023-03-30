#!/bin/bash

# Filtering chimeric sequences after de novo clustering of ASVs into OTUs using vsearch through Qiime2
# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
PROJECT="Jamy_2019"
MARKER="rDNA"
CELL="cell"
SIM="sim99"
RAW_DATA="../../raw_data"
OTU_DIR="${RAW_DATA}/OTU_clust/${PROJECT}/${MARKER}/${CELL}/${SIM}"
OTU_CHIM_FILT="${RAW_DATA}/OTU_nonchimeric/${PROJECT}/${MARKER}/${CELL}/${SIM}"


###############################
## DE NOVO CHIMERA FILTERING ##
###############################

echo "Working on ${PROJECT} project, ${MARKER} marker, ${CELL} cell, and ${SIM} similarity."

# Cleaning
mkdir -p ${OTU_CHIM_FILT}/
rm -f ${OTU_CHIM_FILT}/chimeras.qza \
      ${OTU_CHIM_FILT}/nonchimeras.qza \
      ${OTU_CHIM_FILT}/stats.qza \
      ${OTU_CHIM_FILT}/stats.qzv


# Chimera removal
echo "De novo chimera filtering of clustered OTUs..."

qiime vsearch uchime-denovo \
  --i-table ${OTU_DIR}/otu_table.qza \
  --i-sequences ${OTU_DIR}/otu_seqs.qza \
  --o-chimeras ${OTU_CHIM_FILT}/chimeras.qza \
  --o-nonchimeras ${OTU_CHIM_FILT}/nonchimeras.qza \
  --o-stats ${OTU_CHIM_FILT}/stats.qza

# Visualize stats file
qiime metadata tabulate \
  --m-input-file ${OTU_CHIM_FILT}/stats.qza \
  --o-visualization ${OTU_CHIM_FILT}/stats.qzv


#######################################
## FILTERING OTU TABLE AND SEQUNECES ##
#######################################

# Cleaning
rm -f ${OTU_CHIM_FILT}/otu_table_nonchimeric.qza \
      ${OTU_CHIM_FILT}/otu_seqs_nonchimeric.qza

# Exclude chimeras in otu table
qiime feature-table filter-features \
  --i-table ${OTU_DIR}/otu_table.qza \
  --m-metadata-file ${OTU_CHIM_FILT}/nonchimeras.qza \
  --o-filtered-table ${OTU_CHIM_FILT}/otu_table_nonchimeric.qza
# Exclude chimeras in otu sequences
qiime feature-table filter-seqs \
  --i-data ${OTU_DIR}/otu_seqs.qza \
  --m-metadata-file ${OTU_CHIM_FILT}/nonchimeras.qza \
  --o-filtered-data ${OTU_CHIM_FILT}/otu_seqs_nonchimeric.qza


#############################################
## VISUALIZATION, STATS, & IMPORT TO FASTA ##
#############################################

# Cleaning
rm -f ${OTU_CHIM_FILT}/otu_table_nonchimeric.qzv \
      ${OTU_CHIM_FILT}/otu_table_nonchimeric_summarize.qzv \
      ${OTU_CHIM_FILT}/otu_seqs_nonchimeric.qzv \
      ${OTU_CHIM_FILT}/otu_seqs_nonchimeric.fasta

# visualize otu table
qiime metadata tabulate \
  --m-input-file ${OTU_CHIM_FILT}/otu_table_nonchimeric.qza \
  --o-visualization ${OTU_CHIM_FILT}/otu_table_nonchimeric.qzv

# summarize otu table
qiime feature-table summarize \
  --i-table ${OTU_CHIM_FILT}/otu_table_nonchimeric.qza \
  --o-visualization ${OTU_CHIM_FILT}/otu_table_nonchimeric_summarize.qzv

# visualizing otu sequences - converting qza to qzv
qiime feature-table tabulate-seqs \
  --i-data ${OTU_CHIM_FILT}/otu_seqs_nonchimeric.qza \
  --o-visualization ${OTU_CHIM_FILT}/otu_seqs_nonchimeric.qzv

# convert qiime2 otu sequences to fasta
qiime tools export \
  --input-path ${OTU_CHIM_FILT}/otu_seqs_nonchimeric.qza \
  --output-path ${OTU_CHIM_FILT}/

mv ${OTU_CHIM_FILT}/dna-sequences.fasta \
   ${OTU_CHIM_FILT}/otu_seqs_nonchimeric.fasta
