#!/bin/bash

# script for the transforming the artifact files to fasta files
# from the qiime (vsearch) taxonomical assignment and
# extracting the vampyrellida-specific sequences into a new fasta files for each sample



# variables
TAXON="Vampyrellida"
MARKER="Full18S"
CELL="cellCombined"
RAW_DATA="/home/lubo/code/wRajter/vampyrella_2023/raw_data"
# local machine at uni: "/home/lubomir/projects/vampyrella_2023/raw_data"
REP_SEQS="${RAW_DATA}/qiime_output/dada2_post_denoise"
OUTPUT="${REP_SEQS}/fasta"
TSV_FILES="${RAW_DATA}/qiime_output/assignment_results/${MARKER}_${CELL}"
SEQTK="${RAW_DATA}/packages/seqtk"
SAMPLES="ALLSAMPLES \
         A3 \
         Mock \
         NH1 \
         NH4 \
         Sim17 \
         Sim22 \
         Th16 \
         Th38 \
         Th40 \
         X17007"

for SAMPLE in ${SAMPLES}
do
    echo "Working on sample ${SAMPLE}"

    grep "${TAXON}" ${TSV_FILES}/taxonomy_${SAMPLE}.tsv | awk '{print $1}' > ${TSV_FILES}/vamp_ids_${MARKER}_${CELL}_${SAMPLE}.txt
    ${SEQTK}/seqtk subseq ${OUTPUT}/asv_${MARKER}_${CELL}_${SAMPLE}.fasta ${TSV_FILES}/vamp_ids_${MARKER}_${CELL}_${SAMPLE}.txt > ${OUTPUT}/${TAXON}_${MARKER}_${CELL}_${SAMPLE}.fasta
done
