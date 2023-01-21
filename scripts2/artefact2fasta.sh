#!/bin/bash

# script for the transforming the artifact files to fasta files
# from the qiime (vsearch) taxonomical assignment and
# extracting the vampyrellida-specific sequences into a new fasta files for each sample



# variables
RAW_DATA="/home/lubomir/projects/vampyrella_2023/raw_data"
REP_SEQS="${RAW_DATA}/qiime_output/individual_samples/dada2_output"
OUTPUT="${REP_SEQS}/ASVs"
TSV_FILES="${RAW_DATA}/qiime_output/individual_samples/assignment_results/tsv_files"
TAXON="Vampyrellida"
MARKER="18S"
CELL="cell1"


for SAMPLE in "A3_${MARKER}_${CELL}" \
              "Mock_${MARKER}_${CELL}" \
              "NH1_${MARKER}_${CELL}" \
              "NH4_${MARKER}_${CELL}" \
              "Sim17_${MARKER}_${CELL}" \
              "Sim22_${MARKER}_${CELL}" \
              "Th16_${MARKER}_${CELL}" \
              "Th38_${MARKER}_${CELL}" \
              "Th40_${MARKER}_${CELL}" \
              "X17007_${MARKER}_${CELL}";
    do
        echo "Working on sample ${SAMPLE}"
        qiime tools export \
            --input-path ${REP_SEQS}/representative_sequences_${SAMPLE}.qza \
            --output-path ${OUTPUT}/
        mv ${OUTPUT}/dna-sequences.fasta ${OUTPUT}/asv_${SAMPLE}.fasta

        grep "${TAXON}" ${TSV_FILES}/vsearch_taxonomy_${SAMPLE}.tsv | awk '{print $1}' > ${TSV_FILES}/vamp_ids_${SAMPLE}.txt
        seqtk subseq ${OUTPUT}/asv_${SAMPLE}.fasta ${TSV_FILES}/vamp_ids_${SAMPLE}.txt > ${OUTPUT}/${TAXON}_${SAMPLE}.fasta
    done
