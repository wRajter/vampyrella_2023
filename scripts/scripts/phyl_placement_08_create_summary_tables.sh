#!/bin/bash


# create a tsv table with OTUs and their stats:
  # OTU_ID
  # Sample_ID
  # LWR
  # Taxonomy

# Variables:
RAW_DATA="../../raw_data"
PROJECT="Suthaus_2022"
TAX_ASSIGN_DIR="${RAW_DATA}/phyl_placement/${PROJECT}/vampyrellida/downstream_analyses/tax_assignment"
SAMPLE_LABELS=$(ls ${TAX_ASSIGN_DIR}/*_per_query.tsv | \
                awk -F '/' '{ print $NF }' | \
                awk -F '_' '{ print $1 }' | \
                awk 'BEGIN { ORS = "\t" } { print $1 }')





###################
## PREPARE TABLE ##
###################

rm -f ${TAX_ASSIGN_DIR}/tax_assign_summary_table.tsv


for SAMPLE in ${SAMPLE_LABELS}
do
  echo "working on sample: ${SAMPLE}"
  tac ${TAX_ASSIGN_DIR}/${SAMPLE}_18S_otu_per_query.tsv | \
  sort -u -k1,1 | \
  head -n -1 | \
  awk -v sample=${SAMPLE} {'print $1"\t"sample"\t"$2"\t"$6'} >> ${TAX_ASSIGN_DIR}/tax_assign_summary_table.tsv
done

sed  -i "1i otu_id\tsample_id\tlwr\ttaxopath" ${TAX_ASSIGN_DIR}/tax_assign_summary_table.tsv
