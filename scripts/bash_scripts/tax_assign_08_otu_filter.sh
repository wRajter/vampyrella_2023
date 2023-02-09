#!/bin/bash

# Filtering denoised ASVs
# Note: Activate conda qiime2-2022.11 environment before running the script

# Variables:
CELL="cellCombined"
MARKER="Full18S"
SIM="sim97"
RAW_DATA="../../raw_data"
OTU_CLUST_DIR="${RAW_DATA}/OTU_clust/${MARKER}/${CELL}/${SIM}"
TAX_ASSIGN_DIR="${RAW_DATA}/tax_assign_results/${MARKER}/${CELL}/${SIM}"
FILT_OTU_DIR="${RAW_DATA}/OTU_filtered/${MARKER}/${CELL}/${SIM}"


##########################
## FILTER OUT RARE OTUs ##
##########################

echo "Working on ${MARKER} ${CELL} ${SIM} sample."

# Filtering the denoised table to remove noise/singletons
echo "Filtering out rare OTUs..."

# Cleaning
mkdir -p ${FILT_OTU_DIR}/
rm -f ${FILT_OTU_DIR}/table_filt.qza

# Filtering out rare OTUs

read -p "Enter minimal OTU frequency threshold (more details in help.txt): " MIN_FREQ_THRESHOLD

qiime feature-table filter-features \
   --i-table ${OTU_CLUST_DIR}/otu_table.qza \
   --p-min-frequency ${MIN_FREQ_THRESHOLD} \
   --p-min-samples 1 \
   --o-filtered-table ${FILT_OTU_DIR}/table_filt.qza


##################################################
## FILTER OUT CONTAMINANT AND UNCLASSIFIED OTUs ##
##################################################

# Removing OTUs which are likely contaminants or noise based on the taxonomic labels
echo "Filtering out contaminant and unclassified OTUs..."

# Cleaning
rm -f ${FILT_OTU_DIR}/table_filt_contam.qza

# Filtering
# Note:
# In this step, we will remove any OTU,
# which contains word mitochondrial and chloroplast in its taxonomic label.
# Then, we will exclude any OTU that is unclassified at the phylum level
# as those sequences could be noise (e.g. possible chimeric sequences).

qiime taxa filter-table \
   --i-table ${FILT_OTU_DIR}/table_filt.qza \
   --i-taxonomy ${TAX_ASSIGN_DIR}/vsearch_taxonomy.qza \
   --p-include p__ \
   --p-exclude mitochondria,chloroplast \
   --o-filtered-table ${FILT_OTU_DIR}/table_filt_contam.qza


########################
## RAREFACTION CURVES ##
########################

# Visualizing rarefaction curves to determine the read depth plateaus and
# to see the richness of the samples

# Removing OTUs which are likely contaminants or noise based on the taxonomic labels
echo "Creating rarefaction curves..."

# Cleaning
rm -f ${FILT_OTU_DIR}/table_filt_contam_summary.qzv \
      ${FILT_OTU_DIR}/rarefaction_curves.qzv \
      ${FILT_OTU_DIR}/otu_table_filtered.qza

# Before creating rarefaction curves, we need to summarize the filtered table
# to determine the maximum depth across the samples.
qiime feature-table summarize \
   --i-table ${FILT_OTU_DIR}/table_filt_contam.qza \
   --o-visualization ${FILT_OTU_DIR}/table_filt_contam_summary.qzv

read -p "Enter maximum sequencing depth (more details in help.txt): " MAX_DEPTH

# Creating rarefaction curves
qiime diversity alpha-rarefaction \
   --i-table ${FILT_OTU_DIR}/table_filt_contam.qza \
   --p-max-depth ${MAX_DEPTH} \
   --p-steps 20 \
   --p-metrics 'observed_features' \
   --o-visualization ${FILT_OTU_DIR}/rarefaction_curves.qzv


# We can now examine the rarefication curves.
# If we want to exclude some sample(s) with poor sequncing depth,
# we can decide on a minimum depth cut-off.
# Otherwise we will retain all samples then
# we can simply make a copy of the QZA file with
# the final table filename (see loop below).

read -p "Do you want to exclude sample(s) based on sequencing depth (more details in help.txt) [Yes/No]: " FILT_SAMPLES

if [ "${FILT_SAMPLES}" = "Yes" ]
then
  read -p "Enter sequence depth cut-off (more details in help.txt): " SEQ_DEPTH_CUTOFF
  qiime feature-table filter-samples \
    --i-table ${FILT_OTU_DIR}/table_filt_contam.qza \
    --p-min-frequency ${SEQ_DEPTH_CUTOFF} \
    --o-filtered-table ${FILT_OTU_DIR}/otu_table_filtered.qza
elif [ "${FILT_SAMPLES}" = "No" ]
then
  cp ${FILT_OTU_DIR}/table_filt_contam.qza \
    ${FILT_OTU_DIR}/otu_table_filtered.qza
else
  echo "Please answer Yes/No"
fi


#########################################
## SUBSET AND SUMMARIZE FILTERED TABLE ##
#########################################

# Now that we have our final filtered table,
# we will need to subset the QZA of OTU sequences to the same set.

echo "Subsetting the OTU sequences based the filtered table..."

# Cleaning
rm -f ${FILT_OTU_DIR}/otu_seqs_filtered.qza \
      ${FILT_OTU_DIR}/otu_seqs_filtered.qzv \
      ${FILT_OTU_DIR}/otu_seqs_filtered.fasta \
      ${FILT_OTU_DIR}/table_filtered_summary.qzv


# Subsetting
qiime feature-table filter-seqs \
   --i-data ${OTU_CLUST_DIR}/otu_seqs.qza \
   --i-table ${FILT_OTU_DIR}/otu_table_filtered.qza \
   --o-filtered-data ${FILT_OTU_DIR}/otu_seqs_filtered.qza

qiime feature-table tabulate-seqs \
    --i-data ${FILT_OTU_DIR}/otu_seqs_filtered.qza \
    --o-visualization ${FILT_OTU_DIR}/otu_seqs_filtered.qzv

qiime tools export \
  --input-path ${FILT_OTU_DIR}/otu_seqs_filtered.qza \
  --output-path ${FILT_OTU_DIR}/

mv ${FILT_OTU_DIR}/dna-sequences.fasta \
   ${FILT_OTU_DIR}/otu_seqs_filtered.fasta

# Creating final summary table for filtered OTUs
qiime feature-table summarize \
   --i-table ${FILT_OTU_DIR}/otu_table_filtered.qza  \
   --o-visualization ${FILT_OTU_DIR}/table_filtered_summary.qzv

# Create log file with parameters used in the analyses
echo "Custome parameters used:
      Minimal OTU frequency threshold: ${MIN_FREQ_THRESHOLD}
      Maximal sequencing depth: ${MAX_DEPTH}" > ${FILT_OTU_DIR}/params.log

if [ -n "${SEQ_DEPTH_CUTOFF}" ]
then
  echo "Sequence depth cut-off: ${SEQ_DEPTH_CUTOFF}" >> ${FILT_OTU_DIR}/params.log
fi
