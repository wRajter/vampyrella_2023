#!/bin/bash

# script for the transforming the artifact files to fasta files
# from the qiime (vsearch) taxonomical assignment and
# extracting the vampyrellida-specific sequences into a new fasta files for each sample



# variables
PROJECT="Suthaus_2022"
TAXON="Vampyrellida"
MARKER="Full18S"
CELL="cellCombined"
SIM="sim99"
RAW_DATA="../../raw_data"
SUMMARY_TABLE="${RAW_DATA}/OTU_results/${PROJECT}/otu_summary_table_${MARKER}_${CELL}_${SIM}.tsv"
VAMP_SEQ_DIR="${RAW_DATA}/vamp_specific_seqs/${PROJECT}/${MARKER}/${CELL}/${SIM}"
SEQTK="${RAW_DATA}/packages/seqtk"
PER_SAMPLE_DIR="${RAW_DATA}/per_sample_results/${PROJECT}/${MARKER}/${CELL}/${SIM}"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}"
SAMPLES=$(ls ${RAW_READS_DIR}/*reads.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '_' '{ print $1 }')


mkdir -p ${VAMP_SEQ_DIR}/

for SAMPLE in ${SAMPLES}
do
  echo "Working on sample ${SAMPLE}"
  rm -f ${VAMP_SEQ_DIR}/otu_seqs_vamp_${SAMPLE}.fasta \
        ${VAMP_SEQ_DIR}/vamp_ids_${SAMPLE}.txt

  awk -vtaxon="Order" -vsample="$SAMPLE" \
      'NR==1\
      {for(i=1;i<=NF;i++)\
      {if($i==taxon)taxon_col=i; \
      if ($i==sample)sample_col=i;}} \
      NR>0 && $taxon_col=="Vampyrellida" && $sample_col>0 {print $1}' ${SUMMARY_TABLE} > ${VAMP_SEQ_DIR}/vamp_ids_${SAMPLE}.txt

  ${SEQTK}/seqtk subseq \
  ${PER_SAMPLE_DIR}/fasta/otu_seqs_filtered_${SAMPLE}.fasta \
  ${VAMP_SEQ_DIR}/vamp_ids_${SAMPLE}.txt > \
  ${VAMP_SEQ_DIR}/otu_seqs_vamp_${SAMPLE}.fasta

  rm -f ${VAMP_SEQ_DIR}/vamp_ids_${SAMPLE}.txt
done
