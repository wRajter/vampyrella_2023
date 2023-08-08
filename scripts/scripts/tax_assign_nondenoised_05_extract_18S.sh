#!/bin/bash

# Extracting the 18S rRNA region from the environmental sequences using cutadapt
# Extracted 18S will serve for taxonomic assignemnt

# Variables
PROJECT="Suthaus_2022"
MARKER="rDNA"
CELL="cell"
SIM="Sim99"
F_PRIMER="AACCTGGTTGATCCTGCCAG"
R_PRIMER="TACAAAGGGCAGGGACGTAAT"
IDENT_THRESHOLD=0.8 # minimum combined primer match identity threshold
RAW_DATA="../../raw_data"
OTU_CHIM_FILT="${RAW_DATA}/OTU_nonchimeric/${PROJECT}/${MARKER}/${CELL}/${SIM}"
EXTRACTED_18S="${RAW_DATA}/extracted_18S/${PROJECT}/${MARKER}/${CELL}/${SIM}"
RAW_READS_DIR="${RAW_DATA}/PacBio/${PROJECT}_${MARKER}/${CELL}/filtered"


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




################################
## EXTRACTING 18S PER SAMPLES ##
################################


mkdir -p ${EXTRACTED_18S}/


SAMPLES=$(ls ${RAW_READS_DIR}/*.fastq.gz | \
          awk -F '/' '{ print $NF }' | \
          awk -F '.' '{ print $1 }')

echo "Samples used:"
echo "$SAMPLES"


# Reverse complement the reverse primer
R_PRIMER_RC=$(echo $R_PRIMER | rev | tr 'ATCGatcg' 'TAGCtagc')
echo "reverse complement of ${R_PRIMER} is ${R_PRIMER_RC}"


for SAMPLE in ${SAMPLES}
do
  # Cleaning
  rm -f ${EXTRACTED_18S}/extracted_18S_${SAMPLE}.fasta \
        ${EXTRACTED_18S}/trimming_${SAMPLE}.log

  # Trim the sequences based on the primers
  cutadapt -a ${R_PRIMER_RC} \
           -g ${F_PRIMER} \
           -n 3 \
           -o ${EXTRACTED_18S}/extracted_18S_${SAMPLE}.fasta \
              ${OTU_CHIM_FILT}/otu_${SAMPLE}.fasta \
            > ${EXTRACTED_18S}/trimming_${SAMPLE}.log



  # Calculating how many 18S fragments were actually extracted from the long fragments
  NUM_SEQS=$(grep '>' ${OTU_CHIM_FILT}/otu_${SAMPLE}.fasta | wc -l)
  NUM_EXTRACT=$(grep '>' ${EXTRACTED_18S}/extracted_18S_${SAMPLE}.fasta | wc -l)
  PERCENTAGE=$(awk " BEGIN { print $NUM_SEQS / $NUM_EXTRACT * 100}")
  echo -e "\nSample: ${SAMPLE}\nNumber of sequences: ${NUM_SEQS}\nExtracted: ${NUM_EXTRACT}\nPercentage extracted: ${PERCENTAGE}%\n"

done
