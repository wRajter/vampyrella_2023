#!/bin/bash

# Extracting the 18S rRNA region from the environmental sequences using Qiime2
# Extracted 18S will serve for taxonomic assignemnt
# Note: Activate conda qiime2-2022.11 environment before running the script.

# Variables
PROJECT="Jamy_2019"
MARKER="rDNA"
CELL="cell"
SIM="sim99"
F_PRIMER="AACCTGGTTGATCCTGCCAG"
R_PRIMER="TACAAAGGGCAGGGACGTAAT"
IDENT_THRESHOLD=0.8 # minimum combined primer match identity threshold
RAW_DATA="../../raw_data"
OTU_CHIM_FILT="${RAW_DATA}/OTU_nonchimeric/${PROJECT}/${MARKER}/${CELL}/${SIM}"
EXTRACTED_18S="${RAW_DATA}/extracted_18S/${PROJECT}/${MARKER}/${CELL}/${SIM}"

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



############################
## FILTERING AND TRIMMING ##
############################

echo "Working on ${PROJECT} project, ${MARKER} marker, ${CELL} cell, and ${SIM} similarity."

# Cleaning
mkdir -p ${EXTRACTED_18S}/
rm -f ${EXTRACTED_18S}/extracted_18S_seqs.qza \
      ${EXTRACTED_18S}/extracted_18S_seqs.qzv \
      ${EXTRACTED_18S}/extracted_18S_seqs.fasta \
      ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.fasta \
      ${EXTRACTED_18S}/trimming.log \
      ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.qza


# Reverse complement the reverse primer
R_PRIMER_RC=$(echo $R_PRIMER | rev | tr 'ATCGatcg' 'TAGCtagc')
echo "reverse complement of ${R_PRIMER} is ${R_PRIMER_RC}"


# Trim the sequences based on the primers
cutadapt -a ${R_PRIMER_RC} \
         -g ${F_PRIMER} \
         -n 3 \
         -o ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.fasta \
            ${OTU_CHIM_FILT}/otu_seqs_nonchimeric.fasta \
          > ${EXTRACTED_18S}/trimming.log

        #  --minimum-length 1 \
        #  --trim-n \
        #  --discard-untrimmed \
        # --trimmed-only \


qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.fasta \
  --output-path ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.qza

# Visalize extracted sequences and convert them into Fasta files
qiime feature-table tabulate-seqs \
  --i-data ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.qza \
  --o-visualization ${EXTRACTED_18S}/extracted_18S_seqs_trimmed.qzv
