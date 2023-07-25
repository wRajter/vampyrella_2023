#!/bin/bash

# Downstream (phylogenetic placement) analyses in GAPPA

# Variables
PROJECT="all_seqs"
TAXON="eukaryotes"
MARKER="rDNA"
CELL="cell"
RAW_DATA="../../raw_data"
PHY_PLAC_DIR="${RAW_DATA}/phyl_placement/${PROJECT}/${TAXON}"
JPLACE_DIR="${PHY_PLAC_DIR}/phyl_placement_analysis"
OUT_DIR="${PHY_PLAC_DIR}/downstream_analyses"
REF_VERSION="2022"
TAXON_FILE_DIR="${RAW_DATA}/reference_alignments/vamp_phylo_placement/${TAXON}/reference_alignment_${REF_VERSION}"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"

# SAMPLES=$(ls ${RAW_READS_DIR}/*.fastq.gz | \
#           awk -F '/' '{ print $NF }' | \
#           awk -F '.' '{ print $1 }')


SAMPLES+=" allsamples" # add all samples combined into our sample names array


echo "Samples used:"
echo "$SAMPLES"

########################
## HEAT TREE ANALYSIS ##
########################

mkdir -p ${OUT_DIR}/heat_tree/

for SAMPLE in ${SAMPLES}
do
  echo "Working on sample ${SAMPLE}"
  # log scale
  gappa examine heat-tree \
    --jplace-path ${JPLACE_DIR}/${SAMPLE}/epa_result.jplace \
    --mass-norm absolute \
    --log-scaling \
    --out-dir ${OUT_DIR}/heat_tree/ \
    --write-svg-tree \
    --allow-file-overwriting \
    --log-file ${OUT_DIR}/heat_tree/heat_tree_${SAMPLE}_log.log

  mv ${OUT_DIR}/heat_tree/tree.svg \
    ${OUT_DIR}/heat_tree/heat_tree_${SAMPLE}_log.svg

  # normal scale
  gappa examine heat-tree \
    --jplace-path ${JPLACE_DIR}/${SAMPLE}/epa_result.jplace \
    --mass-norm absolute \
    --out-dir ${OUT_DIR}/heat_tree/ \
    --write-svg-tree \
    --allow-file-overwriting \
    --log-file ${OUT_DIR}/heat_tree/heat_tree_${SAMPLE}_norm.log

  mv ${OUT_DIR}/heat_tree/tree.svg \
    ${OUT_DIR}/heat_tree/heat_tree_${SAMPLE}_norm.svg

done



##########################
## TAXONOMIC ASSIGNMENT ##
##########################

mkdir -p ${OUT_DIR}/tax_assignment/

for SAMPLE in ${SAMPLES}
do
  echo "Working on sample ${SAMPLE}"
  gappa examine assign \
    --jplace-path ${JPLACE_DIR}/${SAMPLE}/epa_result.jplace \
    --taxon-file ${TAXON_FILE_DIR}/taxon_file.tsv \
    --out-dir ${OUT_DIR}/tax_assignment/ \
    --per-query-results \
    --allow-file-overwriting \
    --log-file ${OUT_DIR}/tax_assignment/${SAMPLE}_tax_assignment.log

  # cleaning
  mv ${OUT_DIR}/tax_assignment/per_query.tsv \
     ${OUT_DIR}/tax_assignment/${SAMPLE}_per_query.tsv
  mv ${OUT_DIR}/tax_assignment/profile.tsv \
     ${OUT_DIR}/tax_assignment/${SAMPLE}_profile.tsv

done



#####################
## EXTRACTING OTUs ##
#####################

mkdir -p ${OUT_DIR}/extract_otus


for SAMPLE in ${SAMPLES}
do
  echo "Working on sample ${SAMPLE}"
  gappa prepare extract \
    --jplace-path ${JPLACE_DIR}/${SAMPLE}/epa_result.jplace \
    --clade-list-file ${TAXON_FILE_DIR}/taxon_vamp.tsv \
    --fasta-path ${JPLACE_DIR}/${SAMPLE}/query.fasta \
    --allow-file-overwriting \
    --color-tree-file ${OUT_DIR}/extract_otus/${SAMPLE}_tree \
    --sequences-out-dir ${OUT_DIR}/extract_otus/${SAMPLE}_sequences \
    --samples-out-dir ${OUT_DIR}/extract_otus/${SAMPLE}_samples \
    --log-file ${OUT_DIR}/extract_otus/${SAMPLE}_extract_otus.log
done
