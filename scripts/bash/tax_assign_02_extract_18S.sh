#!/bin/bash

# Extracting the 18S rRNA region from the environmental sequences using cutadapt
# Extracted 18S will serve for taxonomic assignemnt

# Variables

# Prompt the user for the PROJECT
echo -n "Enter the project name (e.g., Jamy_2019) and press [ENTER]: "
read PROJECT

# Prompt the user for the MARKER
echo -n "Enter the marker (e.g., rDNA) and press [ENTER]: "
read MARKER

F_PRIMER="AACCTGGTTGATCCTGCCAG"
R_PRIMER="TACAAAGGGCAGGGACGTAAT"
IDENT_THRESHOLD=0.8 # minimum combined primer match identity threshold
RAW_DATA="../../raw_data"
READS_DIR="${RAW_DATA}/dada2/${PROJECT}/${MARKER}/filtered"

# Outpup directory
EXTRACTED_18S="${RAW_DATA}/extracted_18S/${PROJECT}/${MARKER}"

# Reverse complement the reverse primer
R_PRIMER_RC=$(echo $R_PRIMER | rev | tr 'ATCGatcg' 'TAGCtagc')
echo "reverse complement of ${R_PRIMER} is ${R_PRIMER_RC}"

# Make output directory if it doesn't exist
mkdir -p ${EXTRACTED_18S}/


# Loop through all the fastq.gz files
for file in ${READS_DIR}/*.fastq.gz; do

  # Get sample name
  sample=$(echo $file | awk -F '/' '{ print $NF }' | awk -F '.' '{ print $1 }')
  echo "Sample: ${sample}"

  # Sanity check
  rm -f ${EXTRACTED_18S}/extracted_18S_${sample}.fasta \
        ${EXTRACTED_18S}/trimming_${sample}.log

  # Trim the sequences based on the primers
  cutadapt -a ${R_PRIMER_RC} \
           -g ${F_PRIMER} \
           -n 3 \
           -o ${EXTRACTED_18S}/extracted_18S_${sample}.fastq.gz \
              ${file} \
            > ${EXTRACTED_18S}/trimming_${sample}.log

done


#############################
## primers to test (5'-3') ##
#############################

# SSU_F04 GGCCAGCCGCGGTAAT
# SSU_R22 TGGTCCGTGTTTCAAGAC
# Sequence: GTCTTGAAACACGGACCA; Type: regular 3'; Length: 18; Trimmed: 923 times
# Sequence: TTAGCATGGAATAAT; Type: regular 5'; Length: 15; Trimmed: 748 times
# Sequence Count 952
# Min Length 1116
# Max Length	4401
# Mean Length 2346.57

# EukA AACCTGGTTGATCCTGCCAGT
# EukB TGCATGCCAAGCTCCTTCGCT
# Sequence: AGCGAAGGAGCTTGGCATGCA; Type: regular 3'; Length: 21; Trimmed: 0 times
# Sequence: AACCTGGTTGATCCTGCCAGT; Type: regular 5'; Length: 21; Trimmed: 0 times

# NS1 GTAGTCATATGCTTGTCTC
# NS8 AGCTGCGTTCTTCATCGA
# Sequence: TCGATGAAGAACGCAGCT; Type: regular 3'; Length: 18; Trimmed: 378 times
# Sequence: GTAGTCATATGCTTGTCTC; Type: regular 5'; Length: 19; Trimmed: 0 times

# Medlin:
# 18SF AACCTGGTTGATCCTGCCAG
# 18SR TAGTTCGCTGCCAGTCCCYRC
# Sequence: GRYGGGACTGGCAGCGAACTA; Type: regular 3'; Length: 21; Trimmed: 3 times
# Sequence: AACCTGGTTGATCCTGCCAG; Type: regular 5'; Length: 20; Trimmed: 927 times


# pr2:
# 18S-SSU-0817-5P TTAGCATGGAATAATRRAATAGGA
# 18S-SSU-1512-3P TACAAAGGGCAGGGACGTAAT
# Sequence: ATTACGTCCCTGCCCTTTGTA; Type: regular 3'; Length: 21; Trimmed: 919 times
# Sequence: TTAGCATGGAATAATRRAATAGGA; Type: regular 5'; Length: 24; Trimmed: 613 times
