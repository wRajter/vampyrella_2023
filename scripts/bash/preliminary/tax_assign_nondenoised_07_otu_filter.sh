#!/bin/bash

# Filtering rare OTUs


# Variables:
PROJECT="Jamy_2022"
CELL="cell"
MARKER="rDNA"
SIM="sim97"
RAW_DATA="../../raw_data"
TAX_ASSIGN_DIR="${RAW_DATA}/tax_assign_results/${PROJECT}/${MARKER}/${CELL}/${SIM}"
FILT_OTU_DIR="${RAW_DATA}/OTU_filtered/${PROJECT}/${MARKER}/${CELL}/${SIM}"
OTU_CHIM_FILT="${RAW_DATA}/OTU_nonchimeric/${PROJECT}/${MARKER}/${CELL}/${SIM}"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"



#########################
## FILTERING RARE OTUs ##
######################@##

# samples
SAMPLES=$(ls ${RAW_READS_DIR}/*.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }')


echo "Samples used:"
echo "$SAMPLES"


mkdir -p ${FILT_OTU_DIR}


for SAMPLE in ${SAMPLES}
do
  # cleaning
  rm ${FILT_OTU_DIR}/nonrare_otu_${SAMPLE}.fasta

  # calculate percentage of the singletons:
  TOTAL_SEQS=$(grep '>' ${OTU_CHIM_FILT}/otu_${SAMPLE}.fasta | wc -l)
  SINGLS=$(grep 'seqs=1' ${OTU_CHIM_FILT}/otu_${SAMPLE}.fasta | wc -l)
  PERC=$(awk " BEGIN { print $SINGLS / $TOTAL_SEQS * 100}")
  echo -e "\nSample: ${SAMPLE}\nTotal sequences: ${TOTAL_SEQS}\nSingletons: ${SINGLS}\nPercentage: ${PERC}%\n"

  TODO: filtrate out seqs that contains seq=1

done

















# head -2 ${OTU_PRECLUST_DIR}/otu_seqs_ERR6454461.fasta | tail -1
# sed '/seqs=1/,/\>/{/seqs\=1/!{/\>/!d}}' ${OTU_PRECLUST_DIR}/otu_seqs_ERR6454461.fasta | grep '>' | wc -l
# sed -e '/seqs=1/,/>/{//!d;};' ${OTU_PRECLUST_DIR}/otu_seqs_ERR6454461.fasta | grep '>' | wc -l

# see the nondenoised_filter.ipynb notebook for further filtering in a non-qiime way
