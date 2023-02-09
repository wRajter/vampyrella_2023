#!/bin/bash

# Filtering denoised ASVs
# Note: Activate conda qiime2-2022.11 environment before running the script

# Variables:
CELL="cellCombined"
MARKER="Full18S"
RAW_DATA="../../raw_data"
DADA2_OUT="${RAW_DATA}/qiime_output/dada2_output"
TAX_ASSIGN_DIR="${RAW_DATA}/qiime_output/assignment_results/${MARKER}_${CELL}"
FIL_ASSIGN_OUT="${RAW_DATA}/filt_assign_output/${MARKER}/${CELL}"
MIN_FREQ_THRESHOLD=13
MAX_DEPTH=31807
SEQ_DEPTH_CUTOFF=50


##########################
## FILTER OUT RARE ASVs ##
##########################

echo "Working on ${MARKER} marker and ${CELL} cell."

# Filtering the denoised table to remove noise/singletons
echo "Filtering out rare ASVs..."

# Cleaning
mkdir -p ${FIL_ASSIGN_OUT}
rm -f ${FIL_ASSIGN_OUT}/table_filt_${MARKER}_${CELL}.qza

# Filtering out rare ASVs
# Note:
# To calculate the cut-off for excluding rare ASVs,
# I use the mean sample depth.
# This number could be find in the dada2 summary table,
# which was created running the reads2asv script.
# Then, I multiply the mean sample depth by 0.001, and round to the nearest integer.
# This represents the min frequency threshold.

qiime feature-table filter-features \
   --i-table ${DADA2_OUT}/table_${MARKER}_${CELL}.qza \
   --p-min-frequency ${MIN_FREQ_THRESHOLD} \
   --p-min-samples 1 \
   --o-filtered-table ${FIL_ASSIGN_OUT}/table_filt_${MARKER}_${CELL}.qza


##################################################
## FILTER OUT CONTAMINANT AND UNCLASSIFIED ASVs ##
##################################################

# Removing ASVs which are likely contaminants or noise based on the taxonomic labels
echo "Filtering out contaminant and unclassified ASVs in the ${MARKER} ${CELL} sample..."

# Cleaning
rm -f ${FIL_ASSIGN_OUT}/table_filt_contam_${MARKER}_${CELL}.qza

# Filtering
# Note:
# In this step, we will remove any ASV,
# which contains word mitochondrial and chloroplast in its taxonomic label.
# Then, we will exclude any ASV that is unclassified at the phylum level
# as those sequences could be noise (e.g. possible chimeric sequences).

qiime taxa filter-table \
   --i-table ${FIL_ASSIGN_OUT}/table_filt_${MARKER}_${CELL}.qza \
   --i-taxonomy ${TAX_ASSIGN_DIR}/vsearch_taxonomy_${MARKER}_${CELL}_ALLSAMPLES.qza \
   --p-include p__ \
   --p-exclude mitochondria,chloroplast \
   --o-filtered-table ${FIL_ASSIGN_OUT}/table_filt_contam_${MARKER}_${CELL}.qza


########################
## RAREFACTION CURVES ##
########################

# Visualizing rarefaction curves to determine the read depth plateaus and
# to see the richness of the samples

# Removing ASVs which are likely contaminants or noise based on the taxonomic labels
echo "Creating rarefaction curves for the ${MARKER} ${CELL} sample..."

# Cleaning
rm -f ${FIL_ASSIGN_OUT}/table_filt_contam_summary_${MARKER}_${CELL}.qzv \
      ${FIL_ASSIGN_OUT}/rarefaction_curves_test_${MARKER}_${CELL}.qzv \
      ${FIL_ASSIGN_OUT}/table_final_${MARKER}_${CELL}.qza

# Before creating rarefaction curves, we need to summarize the filtered table.
# From this summary table, we can determine the maximum depth across your samples:
# the maximum sample total frequency.
# This number will be used as a threshold in the next step (--p-max-depth parameter).
qiime feature-table summarize \
   --i-table ${FIL_ASSIGN_OUT}/table_filt_contam_${MARKER}_${CELL}.qza \
   --o-visualization ${FIL_ASSIGN_OUT}/table_filt_contam_summary_${MARKER}_${CELL}.qzv

# Creating rarefaction curves
qiime diversity alpha-rarefaction \
   --i-table ${FIL_ASSIGN_OUT}/table_filt_contam_${MARKER}_${CELL}.qza \
   --p-max-depth ${MAX_DEPTH} \
   --p-steps 20 \
   --p-metrics 'observed_features' \
   --o-visualization ${FIL_ASSIGN_OUT}/rarefaction_curves_test_${MARKER}_${CELL}.qzv


# We can now examine the rarefication curves.
# If we want to exclude some sample(s) with poor sequncing depth,
# we can decide on a minimum depth cut-off.
# Once we decide on a hard cut-off we can exclude samples below
# this cut-off using the command below (where SEQ_DEPTH_CUTOFF is a placeholder
# for the minimum depth you select):

# qiime feature-table filter-samples \
#    --i-table ${FIL_ASSIGN_OUT}/table_filt_contam_${MARKER}_${CELL}.qza \
#    --p-min-frequency SEQ_DEPTH_CUTOFF \
#    --o-filtered-table ${FIL_ASSIGN_OUT}/table_final_${MARKER}_${CELL}.qza

# Alternatively, if we want to retain all samples then
# we can simply make a copy of the QZA file with
# the final table filename

cp ${FIL_ASSIGN_OUT}/table_filt_contam_${MARKER}_${CELL}.qza \
   ${FIL_ASSIGN_OUT}/table_final_${MARKER}_${CELL}.qza



#########################################
## SUBSET AND SUMMARIZE FILTERED TABLE ##
#########################################

# Now that we have our final filtered table,
# we will need to subset the QZA of ASV sequences to the same set.

echo "Subsetting the ASV sequences based the filtered table for ${MARKER} ${CELL} sample..."

# Cleaning
rm -f ${FIL_ASSIGN_OUT}/asv_seqs_final_${MARKER}_${CELL}.qza \
      ${FIL_ASSIGN_OUT}/table_final_summary_${MARKER}_${CELL}.qzv

# Subsetting
qiime feature-table filter-seqs \
   --i-data ${DADA2_OUT}/asv_${MARKER}_${CELL}_ALLSAMPLES.qza \
   --i-table ${FIL_ASSIGN_OUT}/table_final_${MARKER}_${CELL}.qza \
   --o-filtered-data ${FIL_ASSIGN_OUT}/asv_seqs_final_${MARKER}_${CELL}.qza

# Creating final summary table for filtered ASVs
qiime feature-table summarize \
   --i-table ${FIL_ASSIGN_OUT}/table_final_${MARKER}_${CELL}.qza  \
   --o-visualization ${FIL_ASSIGN_OUT}/table_final_summary_${MARKER}_${CELL}.qzv
