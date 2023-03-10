#!/bin/bash

# Downstream (phylogenetic placement) analyses in GAPPA

# Variables
PROJECT="Suthaus_2022"
TAXON="vampyrellida"
MARKER="Full18S"
CELL="cellCombined"
RAW_DATA="../../raw_data"
PHY_PLAC_DIR="${RAW_DATA}/phyl_placement/${TAXON}"
JPLACE_DIR="${PHY_PLAC_DIR}/phyl_placement_analysis"
OUT_DIR="${PHY_PLAC_DIR}/downstream_analyses"
TAXON_PATH="vamp_phylo_placement"
TAXON_FILE="${RAW_DATA}/reference_alignments/${TAXON_PATH}/taxon_file.tsv"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}"
SAMPLES=$(ls ${RAW_READS_DIR}/*reads.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '_' '{ print $1 }')

SAMPLES+=" allsamples" # add all samples combined into our sample names array




# ########################
# ## HEAT TREE ANALYSIS ##
# ########################

# mkdir -p ${OUT_DIR}/heat_tree/

# for SAMPLE in ${SAMPLES}
# do
#   echo "Working on sample ${SAMPLE}"
#   # log scale
#   gappa examine heat-tree \
#     --jplace-path ${JPLACE_DIR}/${SAMPLE}/epa_result.jplace \
#     --mass-norm absolute \
#     --log-scaling \
#     --out-dir ${OUT_DIR}/heat_tree/ \
#     --write-svg-tree \
#     --allow-file-overwriting \
#     --log-file ${OUT_DIR}/heat_tree/heat_tree_${SAMPLE}_log.log

#   mv ${OUT_DIR}/heat_tree/tree.svg \
#     ${OUT_DIR}/heat_tree/heat_tree_${SAMPLE}_log.svg

#   # normal scale
#   gappa examine heat-tree \
#     --jplace-path ${JPLACE_DIR}/${SAMPLE}/epa_result.jplace \
#     --mass-norm absolute \
#     --out-dir ${OUT_DIR}/heat_tree/ \
#     --write-svg-tree \
#     --allow-file-overwriting \
#     --log-file ${OUT_DIR}/heat_tree/heat_tree_${SAMPLE}_norm.log

#   mv ${OUT_DIR}/heat_tree/tree.svg \
#     ${OUT_DIR}/heat_tree/heat_tree_${SAMPLE}_norm.svg

# done



# ##########################
# ## TAXONOMIC ASSIGNMENT ##
# ##########################

# mkdir -p ${OUT_DIR}/tax_assignment/

# for SAMPLE in ${SAMPLES}
# do
#   echo "Working on sample ${SAMPLE}"
#   gappa examine assign \
#     --jplace-path ${JPLACE_DIR}/${SAMPLE}/epa_result.jplace \
#     --taxon-file ${TAXON_FILE} \
#     --out-dir ${OUT_DIR}/tax_assignment/ \
#     --per-query-results \
#     --allow-file-overwriting \
#     --log-file ${OUT_DIR}/tax_assignment/${SAMPLE}_tax_assignment.log \
#     --krona

#   /usr/local/bin/ktImportText ${OUT_DIR}/tax_assignment/krona.profile

#   # cleaning
#   rm ${OUT_DIR}/tax_assignment/labelled_tree.newick \
#      ${OUT_DIR}/tax_assignment/krona.profile
#   mv ${OUT_DIR}/tax_assignment/per_query.tsv \
#      ${OUT_DIR}/tax_assignment/${SAMPLE}_per_query.tsv
#   mv ${OUT_DIR}/tax_assignment/profile.tsv \
#      ${OUT_DIR}/tax_assignment/${SAMPLE}_profile.tsv
#   mv text.krona.html \
#      ${OUT_DIR}/tax_assignment/${SAMPLE}_text.krona.html

# done



###################
## LABELLED TREE ##
###################

# merge jplace files for all sampales except the mock community

# cleaning
mkdir -p ${OUT_DIR}/labelled_tree/
mkdir -p ${JPLACE_DIR}/allsamples_except_mock/
rm -f ${OUT_DIR}/labelled_tree/paths.txt

# # create a variable with the paths for all the samples except the mock community
# for SAMPLE in ${SAMPLES}
# do
#   if [ ${SAMPLE} != "Mock" ] && [ ${SAMPLE} != "allsamples" ] && [ ${SAMPLE} != "Th40" ]
#   then
#     ls ${JPLACE_DIR}/${SAMPLE}/*.jplace >> ${OUT_DIR}/labelled_tree/paths.txt
#   fi
# done

# merge jplace files
gappa edit merge \
  --jplace-path \
  ../../raw_data/phyl_placement/vampyrellida/phyl_placement_analysis/A3/epa_result.jplace \
  ../../raw_data/phyl_placement/vampyrellida/phyl_placement_analysis/NH1/epa_result.jplace \
  ../../raw_data/phyl_placement/vampyrellida/phyl_placement_analysis/NH4/epa_result.jplace \
  ../../raw_data/phyl_placement/vampyrellida/phyl_placement_analysis/Sim17/epa_result.jplace \
  ../../raw_data/phyl_placement/vampyrellida/phyl_placement_analysis/Sim22/epa_result.jplace \
  ../../raw_data/phyl_placement/vampyrellida/phyl_placement_analysis/Th16/epa_result.jplace \
  ../../raw_data/phyl_placement/vampyrellida/phyl_placement_analysis/Th38/epa_result.jplace \
  ../../raw_data/phyl_placement/vampyrellida/phyl_placement_analysis/X17007/epa_result.jplace \
  --out-dir ${JPLACE_DIR}/allsamples_except_mock/ \
  --log-file ${JPLACE_DIR}/allsamples_except_mock/allsamples_except_mock.log \
  --allow-file-overwriting

# creating the label tree
gappa examine graft \
  --jplace-path ${JPLACE_DIR}/allsamples_except_mock/merge.jplace \
  --fully-resolve \
  --out-dir ${OUT_DIR}/labelled_tree/ \
  --file-prefix allsamples_except_mock_ \
  --log-file ${OUT_DIR}/labelled_tree/allsamples_except_mock_labelled_tree.log \
  --allow-file-overwriting


# creating the label trees for each sample and all samples together

# for SAMPLE in ${SAMPLES}
# do
#   echo "Working on sample ${SAMPLE}"

#   gappa examine graft \
#     --jplace-path ${JPLACE_DIR}/${SAMPLE}/epa_result.jplace \
#     --fully-resolve \
#     --out-dir ${OUT_DIR}/labelled_tree/ \
#     --file-prefix ${SAMPLE}_ \
#     --log-file ${OUT_DIR}/labelled_tree/${SAMPLE}_labelled_tree.log \
#     --allow-file-overwriting
# done
