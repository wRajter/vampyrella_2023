#!/bin/bash


# Merging cell1 and cell2 reads fastq files into cellCombined file

# Variables
RAW_DATA="../../raw_data"
PROJECT="Suthaus_2022_Full18S"
FASTQ_DIR="${RAW_DATA}/PacBio/${PROJECT}"
OUTPUT="${FASTQ_DIR}/cellCombined"

# Save filenames to a variable
files=$(basename -a ${FASTQ_DIR}/cell1/*.fastq.gz)

# Cleaning
rm -r ${OUTPUT} 2>/dev/null
mkdir -p ${OUTPUT}


for file in ${files}
do
  file1="${FASTQ_DIR}/cell1/${file}"
  file2="${FASTQ_DIR}/cell2/${file}"

  if [ -e "$file1" ] && [ -e "$file2" ]; then
    cat $file1 $file2 > ${OUTPUT}/${file}
    echo "Successfully merged ${file} from cell1 and cell2 into cellCombined."
  else
    echo "Warning: Missing file in one of the cells for ${file}"
  fi
done
