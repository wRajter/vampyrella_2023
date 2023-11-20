#!/bin/bash

# Downstream (phylogenetic placement) analyses in GAPPA

# Variables
PROJECT="Suthaus_2022"
TAXON="vampyrella"
MARKER="Full18S"
DENOISE_METHOD="RAD"
SIM="sim_90"
RAW_DATA="../../raw_data"
RESULTS="../../results"
JPLACE_DIR="${RESULTS}/phyl_placement/jplace/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}/vampyrellida"
TAXON_FILE="${RAW_DATA}/reference_alignments/vamp_phylo_placement/vampyrellida/reference_alignment_2023/taxon_file.tsv"
OUT_DIR="${RESULTS}/phyl_placement/tax_assignment/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}/${TAXON}"
LOG_FILES="${RAW_DATA}/phyl_placement/jplace_log_files/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}/vampyrellida"
QUERY_DIR="${RAW_DATA}/chimera_filtered/${PROJECT}/${MARKER}/${DENOISE_METHOD}/${SIM}"


# Samples
SAMPLES=$(ls ${QUERY_DIR}/*.fasta | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }' | \
          awk -F '_' '{ print $1 "_" $2 }')


echo "Samples used:"
echo "$SAMPLES"


##########################
## TAXONOMIC ASSIGNMENT ##
##########################

mkdir -p ${OUT_DIR}/

for SAMPLE in ${SAMPLES}
do
  echo "Working on sample ${SAMPLE}"
  gappa examine assign \
    --jplace-path ${JPLACE_DIR}/epa_result_${SAMPLE}.jplace \
    --taxon-file ${TAXON_FILE} \
    --out-dir ${OUT_DIR} \
    --per-query-results \
    --allow-file-overwriting \
    --log-file ${LOG_FILES}/${SAMPLE}_tax_assignment.log

  # cleaning
  mv ${OUT_DIR}/per_query.tsv \
     ${OUT_DIR}/${SAMPLE}_per_query.tsv
  mv ${OUT_DIR}/profile.tsv \
     ${OUT_DIR}/${SAMPLE}_profile.tsv
  rm -f ${OUT_DIR}/labelled_tree.newick

done
