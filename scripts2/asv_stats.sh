#!/bin/bash

# create a tsv table with ASVs and their stats:
  # ASV_ID
  # Sample_ID
  # Biome
  # Abundance
  # Num_Reads
  # Taxonomy
  # Percent similarity
  # Sequence


# Variables:
PROJECT_DIR="/home/lubo/code/wRajter/vampyrella_2023"
ASV_DIR="${PROJECT_DIR}/raw_data/qiime_output/ASV"
CELL="cell1"
MARKER="Full18S"
SAMPLES="A3 Mock"

for SAMPLE in ${SAMPLES}
do
  grep '>' ${ASV_DIR}/representative_sequences_${MARKER}_${CELL}_${SAMPLE}.fasta | \
  sed 's/>//' | \
  awk -v sample="${SAMPLE}" '{print $0"\t"sample}'
done



#TODO: copy asvs IDs as field 1
#TODO: copy asvs SampleID as field 2
#TODO: copy biome type as field 3
#TODO: copy abundance as field 4
#TODO: copy num_reads as field 5
#TODO: copy taxonomy as field 6
#TODO: copy Percent similarity as field 6
#TODO: copy sequences as field 7
