#!/bin/bash


# Merging cell1 and cell2 reads fastq files into cellCombined file

# Variables
PROJECT_DIR="/home/lubo/code/wRajter/vampyrella_2023"
# local machine at uni: /home/lubomir/projects/vampyrella_2023
FASTQ_DIR="${PROJECT_DIR}/raw_data/PacBio/SuthausFull18S"
OUTPUT="${FASTQ_DIR}/cellCombined"
FILES="A3_18S.hifi_reads.fastq.gz \
       Mock_18S.hifi_reads.fastq.gz \
       NH1_18S.hifi_reads.fastq.gz \
       NH4_18S.hifi_reads.fastq.gz \
       Sim17_18S.hifi_reads.fastq.gz \
       Sim22_18S.hifi_reads.fastq.gz \
       Th16_18S.hifi_reads.fastq.gz \
       Th38_18S.hifi_reads.fastq.gz \
       Th40_18S.hifi_reads.fastq.gz \
       X17007_18S.hifi_reads.fastq.gz"

rm -rf ${OUTPUT}/
mkdir -p ${OUTPUT}

for FILE in ${FILES}
do
    cat ${FASTQ_DIR}/cell1/${FILE} ${FASTQ_DIR}/cell2/${FILE} > ${OUTPUT}/${FILE}
done
